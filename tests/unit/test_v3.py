"""Unit tests for the v3 core.

v2 had no unit tests — everything was validated by running the pipeline in CI.
The generator is the wrong place to repeat that: a wrong wildcard produces a
Snakefile that parses fine and builds the wrong DAG, which an end-to-end run
catches slowly and expensively if at all.

Run with: pixi run pytest tests/unit -q
"""

from __future__ import annotations

import html
import re
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from comparem2.catalogue import CATALOGUE  # noqa: E402
from comparem2 import cli as cli_mod  # noqa: E402
from comparem2.cli import default_databases, parse_overrides, slug  # noqa: E402
from comparem2.guidance import GUIDANCE, citations  # noqa: E402
from comparem2.report import (  # noqa: E402
    _PARTITION_VOCABULARY_MINIMUM,
    render_report,
)
from comparem2.snakefile import render  # noqa: E402
from comparem2.tools import Context, Registry, Scope, Tool  # noqa: E402

SAMPLES = ("A", "B", "C")


# --- sample names --------------------------------------------------

@pytest.mark.parametrize("stem,expected", [
    ("116_2 duplicate", "116_2_duplicate"),  # v2's own test data
    ("GCF_000009205.2", "GCF_000009205.2"),
    ("weird;name&here", "weird_name_here"),
    ("...", "sample"),
    ("", "sample"),
])
def test_slug(stem, expected):
    assert slug(stem) == expected


def test_slug_output_is_wildcard_safe():
    for bad in ["a b", "a;b", "a/b", "a*b", "a{b}"]:
        assert not (set(slug(bad)) & set(" ;/*{}"))


# --- the contract --------------------------------------------------

def test_catalogue_is_the_agreed_set():
    assert len(CATALOGUE) == 13
    assert "antismash" not in CATALOGUE  # dropped: breaks the single env
    assert "gapseq" not in CATALOGUE  # dropped in favour of carveme
    assert "carveme" in CATALOGUE  # the one capability no competitor has
    # Removed 2026-09-02: a read-based metagenome profiler cannot take
    # assemblies, so its command could never have produced output.
    assert "sylph" not in CATALOGUE


@pytest.mark.parametrize("setting,tool,expected", [
    # A single-dash flag has to work: skani's compression factor is `-c`, and
    # an earlier parser prepended "--" to whatever followed the first "--",
    # which rejected `skani-c=125` and turned `skani--c=125` into
    # `-c 70 --c 125` — neither replacing the default nor naming a real flag.
    ("skani-c=125", "skani", ("-c", "125")),
    ("treecluster--threshold=0.1", "treecluster", ("--threshold", "0.1")),
    # A tool name containing a dash must not be split on it.
    ("snp-dists--x=1", "snp-dists", ("--x", "1")),
    # An empty value is passed as a bare flag.
    ("bakta--force=", "bakta", ("--force", "")),
])
def test_overrides_keep_the_flag_spelling(setting, tool, expected):
    result = parse_overrides([setting])
    assert tool in result
    assert expected in result[tool]


def test_override_replaces_rather_than_appends():
    """Naming a flag replaces that flag and leaves the tool's others alone."""
    before = dict(CATALOGUE["mashtree"].params)
    after = dict(parse_overrides(["mashtree--kmerlength=17"])["mashtree"])
    assert after["--kmerlength"] == "17"
    assert before["--kmerlength"] != "17"
    # every other default survives, and no flag is duplicated
    assert set(after) == set(before)
    assert after["--genomesize"] == before["--genomesize"]


def test_override_rejects_unknown_tool_and_malformed_input():
    for bad in ["nope--x=1", "skani", "--threshold=1", ""]:
        with pytest.raises(SystemExit):
            parse_overrides([bad])


def test_closure_pulls_dependencies():
    assert [t.name for t in CATALOGUE.closure(["fasttree"])] == ["bakta", "panaroo", "fasttree"]
    assert [t.name for t in CATALOGUE.closure(["seqkit"])] == ["seqkit"]


def test_closure_rejects_unknown_tool():
    with pytest.raises(KeyError):
        CATALOGUE.closure(["not-a-tool"])


def test_every_tool_renders():
    for tool in CATALOGUE:
        ctx = Context(Path("out"), Path("db"), 4, SAMPLES,
                      sample="A" if tool.scope is Scope.GENOME else None)
        assert all(isinstance(a, str) for a in tool.command(ctx))
        assert list(tool.outputs(ctx)), f"{tool.name} declares no outputs"


def test_stdout_tools_declare_exactly_one_output():
    for tool in CATALOGUE:
        if tool.stdout_to_output:
            ctx = Context(Path("out"), Path("db"), 4, SAMPLES,
                          sample="A" if tool.scope is Scope.GENOME else None)
            assert len(list(tool.outputs(ctx))) == 1, tool.name


def test_install_size_excludes_unmeasured():
    reg = Registry([CATALOGUE["checkm2"]])
    assert reg.install_size() == 1_735_095_710
    assert reg.unmeasured() == []
    assert [d.name for d in CATALOGUE.unmeasured(["bakta"])] == ["bakta-light"]


def test_dependency_cycles_are_impossible_in_catalogue():
    # closure() terminates; a cycle would hang or repeat.
    for tool in CATALOGUE:
        names = [t.name for t in CATALOGUE.closure([tool.name])]
        assert len(names) == len(set(names))


# --- rule generation -----------------------------------------------

def test_genome_rules_use_wildcards_set_rules_do_not():
    text = render(CATALOGUE, ["seqkit", "mashtree"], Path("res"), Path("db"), SAMPLES)
    assert '"res/samples/{sample}/seqkit/contigs.tsv"' in text
    assert '"res/mashtree/mashtree.newick"' in text
    assert "{sample}" not in text.split("rule mashtree:")[1]


def test_threads_are_substituted_by_snakemake():
    text = render(CATALOGUE, ["mashtree"], Path("res"), Path("db"), SAMPLES)
    assert "--numcpus {threads}" in text
    assert "--numcpus 8" not in text


def test_shell_block_uses_wildcards_prefix():
    """Regression: bare {sample} in a shell block is a NameError at runtime.

    Snakemake resolves bare wildcards in input/output paths only; inside shell
    it exposes them as `wildcards.sample`. The generated Snakefile parsed fine
    and failed at execution, which is exactly the class of bug the end-to-end
    run is too slow to be the only guard against.
    """
    text = render(CATALOGUE, ["seqkit"], Path("res"), Path("db"), SAMPLES)
    shell = text.split("    shell:")[1]
    assert "{wildcards.sample}" in shell
    assert "{sample}" not in shell.replace("{wildcards.sample}", "")


def test_stdout_redirect_is_emitted():
    text = render(CATALOGUE, ["mashtree"], Path("res"), Path("db"), SAMPLES)
    assert "> res/mashtree/mashtree.newick" in text


def test_set_tool_depending_on_genome_tool_expands_all_samples():
    text = render(CATALOGUE, ["panaroo"], Path("res"), Path("db"), SAMPLES)
    for sample in SAMPLES:
        assert f'"res/samples/{sample}/bakta/{sample}.gff3"' in text
    # panaroo's own rule must not carry an unresolved wildcard
    assert "{sample}" not in text.split("rule panaroo:")[1]


def test_all_targets_cover_every_sample():
    text = render(CATALOGUE, ["seqkit"], Path("res"), Path("db"), SAMPLES)
    head = text.split("rule seqkit:")[0]
    for sample in SAMPLES:
        assert f'"res/samples/{sample}/seqkit/contigs.tsv"' in head


def test_only_isolated_tools_get_their_own_environment():
    """Shared environment is the default; isolation is the documented exception.

    checkm2 pins DIAMOND 2.1.x and cannot co-solve with a current Bakta, so it
    is the single isolated tool. v2 reached 25 environments by inverting this.
    """
    isolated = [t.name for t in CATALOGUE if t.isolated]
    assert isolated == ["checkm2"]

    text = render(CATALOGUE, None, Path("res"), Path("db"), SAMPLES)
    assert text.count("conda:") == len(isolated)
    assert "conda:" in text.split("rule checkm2:")[1].split("rule ")[0]
    assert "conda:" not in text.split("rule bakta:")[1].split("rule ")[0]


def test_per_rule_conda_forces_all():
    text = render(CATALOGUE, None, Path("res"), Path("db"), SAMPLES, per_rule_conda=True)
    assert text.count("conda:") == len(CATALOGUE)


def test_pinned_tools_carry_minimum_versions():
    """Unpinned bioconda solves silently select years-old broken builds."""
    pinned = {spec.split("::")[1] for t in CATALOGUE for spec in t.conda if ">=" in spec}
    assert any(p.startswith("bakta") for p in pinned)
    assert any(p.startswith("panaroo") for p in pinned)


# --- database downloads -------------------------------------------

def test_every_database_declares_how_to_fetch_itself():
    """A declared size with no fetch is the bug this replaced: cm2 announced
    143.2 GB and then downloaded nothing, failing later inside a tool."""
    for db in CATALOGUE.databases():
        assert db.fetch is not None, f"{db.name} declares no fetch"
        steps = db.fetch(Path("/db"))
        assert steps, db.name
        for step in steps:
            assert all(isinstance(a, str) for a in step), db.name


def test_database_ready_paths_are_distinct_and_under_the_root():
    paths = [db.ready_path(Path("/db")) for db in CATALOGUE.databases()]
    assert len(paths) == len(set(paths))
    for p in paths:
        assert str(p).startswith("/db/")


def test_download_rules_are_generated_and_wired_as_inputs():
    text = render(CATALOGUE, ["checkm2"], Path("/res"), Path("/db"), SAMPLES)
    assert "rule download_checkm2:" in text
    # the tool waits on the database, which is what puts fetching in the DAG
    checkm2 = text.split("rule checkm2:")[1]
    assert '"/db/checkm2/checkm2.dmnd",' in checkm2.split("output:")[0]


def test_download_rule_has_no_inputs_and_one_output():
    text = render(CATALOGUE, ["checkm2"], Path("/res"), Path("/db"), SAMPLES)
    block = text.split("rule download_checkm2:")[1].split("\nrule ")[0]
    assert "input:" not in block
    assert block.count('"/db/checkm2/checkm2.dmnd",') == 1


def test_only_selected_databases_get_download_rules():
    """--until seqkit needs nothing, so nothing should be fetched."""
    text = render(CATALOGUE, ["seqkit"], Path("/res"), Path("/db"), SAMPLES)
    assert "download_" not in text
    text = render(CATALOGUE, ["mashtree", "skani"], Path("/res"), Path("/db"), SAMPLES)
    assert "download_" not in text


def test_amrfinder_database_is_marked_out_of_tree():
    """`amrfinder -u` rejects -d, so its data cannot live under --databases.
    Recording that is the difference between a known limitation and a lie."""
    db = CATALOGUE["amrfinder"].database
    assert db.out_of_tree is True
    assert not any(db.out_of_tree for db in CATALOGUE.databases()
                   if db.name != "amrfinder")


def test_gtdbtk_gets_its_database_through_the_environment():
    """GTDB-Tk has no flag for its database — without the env var, --databases
    was silently ignored for the largest database in the pipeline."""
    text = render(CATALOGUE, ["gtdbtk"], Path("/res"), Path("/db"), SAMPLES)
    assert "export GTDBTK_DATA_PATH=/db/gtdb" in text
    # and no other tool exports anything it should not
    assert text.count("export ") == 1


def test_default_database_location_is_shared_not_cwd_relative(monkeypatch):
    """The default used to be `./databases`, so running the same command from
    two directories fetched a second copy of up to 143 GB."""
    monkeypatch.delenv("COMPAREM2_DATABASES", raising=False)
    default = default_databases()
    assert default == Path.home() / ".comparem2" / "databases"
    assert default.is_absolute()


def test_comparem2_databases_env_var_moves_the_default(monkeypatch):
    """A home quota makes home the wrong place for 143 GB on a cluster, and
    `-d` cannot be set once for every run."""
    monkeypatch.setenv("COMPAREM2_DATABASES", "/evo/postdoc/cm2-databases")
    assert default_databases() == Path("/evo/postdoc/cm2-databases")
    monkeypatch.setenv("COMPAREM2_DATABASES", "~/scratch/db")
    assert default_databases() == Path.home() / "scratch" / "db"


def _databases_main_used(monkeypatch, tmp_path, argv: list[str]) -> Path:
    """Run `main` far enough to see which database root it settled on.

    `--report-only` skips Snakemake entirely, so stubbing the two calls that
    receive the path is enough — and it catches the mismatch that was here
    before, where `render_report` got the unresolved `args.databases` while
    everything else got the resolved one.
    """
    seen: dict[str, Path] = {}

    def fake_prepare(registry, selected, workdir, databases, *a, **k):
        seen["prepare"] = databases
        return tmp_path / "Snakefile"

    def fake_report(registry, selected, workdir, databases, *a, **k):
        seen["report"] = databases
        return tmp_path / "report.html"

    monkeypatch.setattr(cli_mod, "prepare", fake_prepare)
    monkeypatch.setattr(cli_mod, "render_report", fake_report)
    (tmp_path / "a.fna").write_text(">c\nACGT\n")
    assert cli_mod.main([str(tmp_path / "a.fna"), "-o", str(tmp_path / "out"),
                         "--until", "seqkit", "--report-only", *argv]) == 0
    assert seen["prepare"] == seen["report"], "the two paths must agree"
    return seen["report"]


def test_explicit_flag_beats_the_env_var(monkeypatch, tmp_path):
    """Precedence is -d, then $COMPAREM2_DATABASES, then home."""
    monkeypatch.setenv("COMPAREM2_DATABASES", str(tmp_path / "from-env"))
    chosen = tmp_path / "from-flag"
    assert _databases_main_used(monkeypatch, tmp_path,
                                ["-d", str(chosen)]) == chosen


def test_env_var_is_used_when_the_flag_is_absent(monkeypatch, tmp_path):
    monkeypatch.setenv("COMPAREM2_DATABASES", str(tmp_path / "from-env"))
    assert _databases_main_used(monkeypatch, tmp_path, []) == tmp_path / "from-env"


def test_relative_database_path_is_resolved_absolute(monkeypatch, tmp_path):
    """Snakemake is handed `--directory <output>`, so a relative database path
    would resolve against the wrong directory once the run starts."""
    monkeypatch.chdir(tmp_path)
    used = _databases_main_used(monkeypatch, tmp_path, ["-d", "rel-db"])
    assert used.is_absolute()
    assert used == (tmp_path / "rel-db").resolve()


def test_tools_without_env_export_nothing():
    text = render(CATALOGUE, ["seqkit", "mashtree"], Path("/res"), Path("/db"), SAMPLES)
    assert "export " not in text


def test_generated_snakefile_parses_as_python_free_text():
    """Directives appear in order in every rule.

    Download rules have no `input:` — a fetch depends on nothing — so only the
    directives a rule actually declares are ordered.
    """
    text = render(CATALOGUE, None, Path("res"), Path("db"), SAMPLES)
    for block in text.split("\nrule ")[2:]:
        name = block.split(":")[0]
        present = [d for d in ("input:", "output:", "shell:") if d in block]
        if name.startswith("download_"):
            assert "input:" not in block, name
            assert present == ["output:", "shell:"], name
        else:
            assert present == ["input:", "output:", "shell:"], name
        assert [block.index(d) for d in present] == sorted(block.index(d) for d in present)


# --- report --------------------------------------------------------

def test_report_only_shows_tools_with_output(tmp_path):
    for sample in SAMPLES:
        d = tmp_path / "samples" / sample / "seqkit"
        d.mkdir(parents=True)
        (d / "contigs.tsv").write_text("c1\t1000\t50.0\nc2\t500\t40.0\n")
    path = render_report(CATALOGUE, ["seqkit", "mashtree"], tmp_path, Path("db"), SAMPLES)
    body = path.read_text()
    assert "seqkit" in body
    assert "mashtree" not in body.split("<footer>")[0]
    assert "1 of 2 tools produced output" in body


def test_report_computes_n50(tmp_path):
    d = tmp_path / "samples" / "A" / "seqkit"
    d.mkdir(parents=True)
    # 100 + 200 + 300 = 600; half is 300; sorted desc 300,200,100 -> N50 = 300
    (d / "contigs.tsv").write_text("c1\t100\t50.0\nc2\t300\t50.0\nc3\t200\t50.0\n")
    path = render_report(CATALOGUE, ["seqkit"], tmp_path, Path("db"), ("A",))
    body = path.read_text()
    assert "600" in body and "300" in body


def test_report_is_self_contained(tmp_path):
    """No fetched subresources, so the file works offline.

    Narrowed from "no href starting with http": citations link out to doi.org,
    and a hyperlink the reader may click is not an external asset. What must not
    appear is anything the browser fetches to render the page.
    """
    d = tmp_path / "samples" / "A" / "seqkit"
    d.mkdir(parents=True)
    (d / "contigs.tsv").write_text("c1\t100\t50.0\n")
    body = render_report(CATALOGUE, ["seqkit"], tmp_path, Path("db"), ("A",)).read_text()
    assert "<script" not in body
    assert "<link" not in body
    assert "<img" not in body
    assert "url(http" not in body  # no webfonts or background images via CSS
    assert "@import" not in body


def test_report_survives_no_results(tmp_path):
    body = render_report(CATALOGUE, ["seqkit"], tmp_path, Path("db"), ("A",)).read_text()
    assert "No results yet." in body


# --- guidance ------------------------------------------------------

def test_every_tool_has_guidance():
    """The invariant that lets guidance live outside catalogue.py.

    Without this, adding a tool silently ships a section nobody can interpret,
    which is the failure this whole module exists to prevent.
    """
    assert set(GUIDANCE) == {t.name for t in CATALOGUE}


def test_guidance_is_complete_and_attributed():
    for name, entry in GUIDANCE.items():
        assert entry.blurb and entry.method, name
        assert entry.reading, f"{name} explains no column"
        assert entry.caveats, f"{name} claims no limitations"
        for citation in (entry.citation, *entry.also):
            assert citation.doi and "/" in citation.doi, (name, citation.doi)
            assert citation.text, name


def test_guidance_carries_no_markup():
    """Guidance is escaped on render, so it must not contain pre-escaped HTML.

    A stray `&lt;` in the source would reach the page as a literal "&lt;".
    """
    for name, entry in GUIDANCE.items():
        blob = " ".join([entry.blurb, entry.method,
                         *(t for pair in entry.reading for t in pair),
                         *entry.caveats])
        for token in ("&lt;", "&gt;", "&amp;", "<p>", "<em>", "<br"):
            assert token not in blob, f"{name} contains {token!r}"


def test_report_includes_guidance_and_citations(tmp_path):
    d = tmp_path / "samples" / "A" / "seqkit"
    d.mkdir(parents=True)
    (d / "contigs.tsv").write_text("c1\t100\t50.0\n")
    body = render_report(CATALOGUE, ["seqkit"], tmp_path, Path("db"), ("A",)).read_text()
    assert "What this is, and how to read it" in body
    assert "Reading this section" in body
    assert "What this cannot tell you" in body
    # SeqKit's own DOI, and the methods list at the end.
    assert "10.1002/imt2.191" in body
    assert "Methods and citations" in body


def test_methods_lists_only_tools_that_ran(tmp_path):
    """The citation list is for pasting into a manuscript, so a tool that
    produced nothing must not appear in it."""
    d = tmp_path / "samples" / "A" / "seqkit"
    d.mkdir(parents=True)
    (d / "contigs.tsv").write_text("c1\t100\t50.0\n")
    body = render_report(CATALOGUE, ["seqkit", "carveme"], tmp_path,
                         Path("db"), ("A",)).read_text()
    refs = body.split('<ol class="refs">')[1]
    assert "10.1002/imt2.191" in refs  # seqkit ran
    assert "10.1093/nar/gky537" not in refs  # carveme did not


def test_citations_deduplicate_shared_methods():
    """DIAMOND is cited by checkm2, bakta and carveme; it should appear once."""
    papers = citations(["checkm2", "bakta", "carveme"])
    dois = [c.doi for c in papers]
    assert len(dois) == len(set(dois))
    assert dois.count("10.1038/s41592-021-01101-x") == 1


def _skani_report(tmp_path, with_af: bool):
    d = tmp_path / "skani"
    d.mkdir(parents=True)
    # Real shape, taken from a run on thylakoid: 116_2 and its duplicate are
    # identical, and the AF matrix is asymmetric.
    (d / "ani.tsv").write_text(
        "3\nA\t100.00\t100.00\t99.14\nB\t100.00\t100.00\t99.14\n"
        "C\t99.14\t99.14\t100.00\n")
    if with_af:
        (d / "ani.tsv.af").write_text(
            "3\nA\t100.00\t100.00\t89.86\nB\t100.00\t100.00\t89.86\n"
            "C\t74.39\t74.39\t100.00\n")
    return render_report(CATALOGUE, ["skani"], tmp_path, Path("db"),
                         ("A", "B", "C")).read_text()


def test_skani_renders_aligned_fraction_alongside_ani(tmp_path):
    """ANI alone is half the answer: skani emits an identity once alignment
    covers ~15% of a genome, so the coverage has to be on screen too."""
    body = _skani_report(tmp_path, with_af=True)
    assert "Average nucleotide identity" in body
    assert "Aligned fraction" in body
    assert "99.14" in body  # from the ANI matrix
    assert "74.39" in body  # only in the AF matrix
    assert "not symmetric" in body
    # each matrix shades against its own floor, not a shared one
    assert "shaded from 99.14% to 100%" in body.lower()
    assert "shaded from 74.39% to 100%" in body.lower()


def test_skani_section_survives_a_missing_af_file(tmp_path):
    """The .af sidecar is skani's, not ours — an older skani, or a partial run,
    must degrade to the ANI matrix rather than losing the section."""
    body = _skani_report(tmp_path, with_af=False)
    assert "Average nucleotide identity" in body
    assert "99.14" in body
    assert "Aligned fraction" not in body


def test_skani_declares_both_matrices_as_outputs():
    ctx = Context(Path("res"), Path("db"), 8, ("A", "B"))
    names = [p.name for p in CATALOGUE["skani"].outputs(ctx)]
    assert names == ["ani.tsv", "ani.tsv.af"]


def _panaroo_report(tmp_path, n, present):
    samples = tuple(f"g{i}" for i in range(n))
    d = tmp_path / "panaroo"
    d.mkdir(parents=True)
    pattern = ["1"] * present + ["0"] * (n - present)
    (d / "gene_presence_absence.Rtab").write_text(
        "Gene\t" + "\t".join(samples) + "\n" + "x\t" + "\t".join(pattern) + "\n")
    body = render_report(CATALOGUE, ["panaroo"], tmp_path, Path("db"), samples).read_text()
    # Scope to the partitions table: the report has three tables in this
    # section and a looser pattern picks up the per-genome one too.
    table = body.split("Pangenome partitions</h3>")[1].split("</table>")[0]
    rows = {html.unescape(k): v for k, v in
            re.findall(r"<td>([^<]+)</td><td[^>]*>([\d,]+)</td>", table)}
    return body, rows


@pytest.mark.parametrize("n,present,expected", [
    # Below 20 genomes the table reports exact counts, so every cluster lands in
    # a row that says precisely what it is rather than in a fraction-bin that
    # may be unreachable.
    (4, 4, "All 4 genomes (core)"),
    (4, 3, "3 of 4 genomes"),
    (4, 1, "1 genome only"),
    (7, 1, "1 genome only"),
    (19, 18, "18 of 19 genomes"),
])
def test_small_pangenomes_report_exact_counts(tmp_path, n, present, expected):
    body, rows = _panaroo_report(tmp_path, n, present)
    assert rows.get(expected) == "1", (expected, rows)
    assert [k for k, v in rows.items() if v == "1" and k != "Total"] == [expected]
    # The unreachable bins must not appear as table rows. They are still
    # explained in the guidance prose, which is why this checks the rows and
    # not the whole document.
    assert not [k for k in rows if "Soft core" in k or "Cloud" in k]


@pytest.mark.parametrize("present,expected", [
    (20, "Core (≥99%)"),
    (19, "Soft core (95–99%)"),  # 0.95 — reachable exactly at n = 20
    (10, "Shell (15–95%)"),
    (2, "Cloud (<15%)"),  # 0.10
])
def test_large_pangenomes_use_the_conventional_bins(tmp_path, present, expected):
    body, rows = _panaroo_report(tmp_path, 20, present)
    assert rows.get(expected) == "1", (expected, rows)
    assert [k for k, v in rows.items() if v == "1" and k != "Total"] == [expected]


def test_partition_switch_over_is_where_all_four_bins_are_reachable():
    """The threshold is not arbitrary: it is the smallest N at which every one
    of the four conventional bins can hold a cluster."""
    n = _PARTITION_VOCABULARY_MINIMUM
    assert (n - 1) / n >= 0.95  # Soft core reachable
    assert 1 / n < 0.15  # Cloud reachable
    assert (n - 2) / n < 0.95  # and Shell still has room below Soft core
    # one genome fewer and Soft core collapses
    assert (n - 2) / (n - 1) < 0.95


def test_guidance_escapes_angle_brackets_into_entities(tmp_path):
    """CheckM2's thresholds are written with bare '<', which must survive as an
    entity rather than being swallowed as a tag."""
    d = tmp_path / "checkm2"
    d.mkdir(parents=True)
    (d / "quality_report.tsv").write_text(
        "Name\tCompleteness\tContamination\nA\t99.1\t0.4\n")
    body = render_report(CATALOGUE, ["checkm2"], tmp_path, Path("db"), ("A",)).read_text()
    assert "&lt;5%" in body
    assert "<5%" not in body


# --- TUI -----------------------------------------------------------

@pytest.mark.asyncio
async def test_tui_lists_tools_and_shows_install_cost():
    pytest.importorskip("textual")
    from comparem2.tui import ComparemTUI
    from textual.widgets import DataTable

    app = ComparemTUI([], Path("results"), Path("databases"), SAMPLES, 4)
    async with app.run_test() as pilot:
        assert app.query_one(DataTable).row_count == len(CATALOGUE)
        # The whole point of the cost line: the 141 GB is visible before running.
        assert "143.2 GB" in app.cost_text
        assert "unknown size" in app.cost_text

        await pilot.press("n")
        assert "0 tools selected" in app.cost_text
        await pilot.press("a")
        assert f"{len(CATALOGUE)} tools selected" in app.cost_text


@pytest.mark.asyncio
async def test_tui_toggle_and_dependency_marking():
    pytest.importorskip("textual")
    from comparem2.tui import ComparemTUI

    app = ComparemTUI([], Path("results"), Path("databases"), SAMPLES, 4)
    async with app.run_test() as pilot:
        await pilot.press("n")
        assert app.selected == set()
        await pilot.press("a")
        await pilot.press("space")  # toggle the row under the cursor off
        assert len(app.selected) == len(CATALOGUE) - 1


def test_tui_does_not_override_textual_reserved_names():
    """Regression: `on_event` is Textual's internal dispatch.

    Overriding it swallowed every framework message and the app died with
    ScreenStackError on the first keypress.
    """
    pytest.importorskip("textual")
    from textual.app import App

    from comparem2.tui import ComparemTUI

    for name in ("on_event", "on_mount", "compose"):
        overridden = ComparemTUI.__dict__.get(name)
        if overridden is not None:
            assert name != "on_event", "on_event is reserved by Textual"
        assert hasattr(App, "on_event")


def test_tui_selection_marks_are_not_rich_markup():
    """Regression: the marks were `[x]`, `[+]` and `[ ]`, and drew as nothing.

    A DataTable cell given a `str` is parsed as Rich markup, and `[x]` is a tag
    rather than text. The selection column was blank — in an interface whose
    whole job is showing what is selected.
    """
    from comparem2.tui import MARK_DEP, MARK_OFF, MARK_ON

    marks = (MARK_ON, MARK_DEP, MARK_OFF)
    for mark in marks:
        assert "[" not in mark and "]" not in mark, f"{mark!r} parses as markup"
    assert len(set(marks)) == 3, "the three states must be distinguishable"


@pytest.mark.asyncio
async def test_tui_rows_reflect_the_seeded_selection_at_startup():
    """`on_mount` refreshed the cost line but not the rows.

    With the selection seeded from `--until`, eleven rows claimed to be
    `pending` when they were not selected at all.
    """
    pytest.importorskip("textual")
    from textual.widgets import DataTable

    from comparem2.tui import MARK_OFF, MARK_ON, ComparemTUI

    app = ComparemTUI([], Path("results"), Path("databases"), SAMPLES, 4,
                      selected=["mashtree"])
    async with app.run_test():
        table = app.query_one(DataTable)
        assert table.get_cell("mashtree", "sel") == MARK_ON
        assert table.get_cell("gtdbtk", "sel") == MARK_OFF
        assert table.get_cell("gtdbtk", "status") == "not selected"


def test_prepare_writes_the_snakefile_and_every_env_file(tmp_path):
    """One entry point for both the CLI and the TUI.

    Writing the Snakefile without the env files is a silent failure: the
    Snakefile parses, the DAG builds, the job runs, and Snakemake dies
    afterwards inside its metadata bookkeeping.
    """
    from comparem2.snakefile import prepare, render_envs

    snakefile = prepare(CATALOGUE, ["checkm2"], tmp_path, tmp_path / "db", SAMPLES)

    assert snakefile == tmp_path / ".comparem2" / "Snakefile"
    assert "rule checkm2:" in snakefile.read_text()
    written = {p.name for p in (tmp_path / ".comparem2" / "envs").iterdir()}
    assert written == set(render_envs(CATALOGUE, ["checkm2"]))
    # checkm2 is the one isolated tool, so its rule carries a `conda:`
    # directive and this file is not optional.
    assert "checkm2.yaml" in written


@pytest.mark.asyncio
async def test_tui_writes_envs_and_withholds_the_report_when_nothing_ran(
        tmp_path, monkeypatch):
    """Regression, found on thylakoid: two defects in one run.

    The TUI built the Snakefile itself instead of calling prepare(), so no env
    file was written and checkm2's `conda:` directive pointed at nothing —
    Snakemake aborted the whole workflow, not the job. It then printed a green
    "Report:" line for a run that had produced no output at all.
    """
    pytest.importorskip("textual")
    from comparem2 import tui as tui_mod

    reports: list[object] = []
    monkeypatch.setattr(tui_mod, "run", lambda *a, **k: iter(
        [tui_mod.Event("job_error", rule="checkm2", message="boom")]))
    monkeypatch.setattr(tui_mod, "render_report",
                        lambda *a, **k: reports.append(a) or Path("report.html"))

    app = tui_mod.ComparemTUI([], tmp_path, tmp_path / "db", SAMPLES, 4,
                              selected=["checkm2"])
    async with app.run_test() as pilot:
        await pilot.press("r")
        await app.workers.wait_for_complete()
        await pilot.pause()
        assert app.state["checkm2"] == tui_mod.FAILED

    assert (tmp_path / ".comparem2" / "envs" / "checkm2.yaml").exists()
    assert reports == [], "a run that produced nothing must not claim a report"


@pytest.mark.asyncio
async def test_tui_settles_every_tool_to_a_terminal_state(tmp_path, monkeypatch):
    """A tool that never started is `not run`, not `pending` for ever.

    Twelve rows reading `pending` after the workflow had already aborted made
    total failure look like a run still in progress.
    """
    pytest.importorskip("textual")
    from comparem2 import tui as tui_mod

    monkeypatch.setattr(tui_mod, "run", lambda *a, **k: iter(
        [tui_mod.Event("job_finished", rule="mashtree")]))
    monkeypatch.setattr(tui_mod, "render_report", lambda *a, **k: Path("report.html"))

    app = tui_mod.ComparemTUI([], tmp_path, tmp_path / "db", SAMPLES, 4,
                              selected=["mashtree", "treecluster"])
    async with app.run_test() as pilot:
        await pilot.press("r")
        await app.workers.wait_for_complete()
        await pilot.pause()
        assert app.state["mashtree"] == tui_mod.DONE
        assert app.state["treecluster"] == tui_mod.NOT_RUN


@pytest.mark.asyncio
async def test_tui_seeds_its_selection_from_until():
    """`--tui --until mashtree` must not open with gtdbtk's 141.4 GB armed.

    The TUI ignored `--until` and selected the whole catalogue, which put a
    141.4 GB download one keypress away for a user who had asked for one tool.
    """
    pytest.importorskip("textual")
    from comparem2.tui import ComparemTUI

    app = ComparemTUI([], Path("results"), Path("databases"), SAMPLES, 4,
                      selected=["mashtree"])
    async with app.run_test():
        assert app.selected == {"mashtree"}
        assert "1 tools selected" in app.cost_text
        assert "no databases" in app.cost_text


def test_runner_names_the_rule_that_finished(tmp_path):
    """Regression: `job_finished` carries `job_id`, not `jobid`, and no rule.

    `job_info` is the only record carrying `rule_name`, so the runner has to
    remember the mapping. Without it every finish event arrived unattributable,
    the TUI dropped them all, and a workflow that completed 6 of 6 steps on
    thylakoid reported that nothing had run.

    This one needs real Snakemake: the defect was in the field names, which a
    fake event stream would simply have agreed with.
    """
    pytest.importorskip("snakemake")
    from comparem2.runner import run

    snakefile = tmp_path / "Snakefile"
    snakefile.write_text(
        'rule all:\n    input: "a.txt"\n\n'
        'rule a:\n    output: "a.txt"\n    shell: "touch {output}"\n')

    events = list(run(snakefile, 1, workdir=tmp_path))

    finished = {e.rule for e in events if e.kind == "job_finished"}
    assert "a" in finished, f"unnamed job_finished among {[e.kind for e in events]}"
    assert (tmp_path / "a.txt").exists()


def test_runner_dry_run_selects_the_dryrun_executor(tmp_path):
    """`ExecutionSettings(dryrun=True)` raises TypeError — it is not a field.

    The broad `except Exception` in run() turned that into an `error` event, so
    a dry run reported itself as a failed workflow rather than as a bug.
    """
    pytest.importorskip("snakemake")
    from comparem2.runner import run

    target = tmp_path / "done.txt"
    snakefile = tmp_path / "Snakefile"
    snakefile.write_text(
        f'rule all:\n    output: "{target}"\n    shell: "touch {{output}}"\n')

    kinds = [e.kind for e in run(snakefile, 1, workdir=tmp_path, dry_run=True)]

    assert "done" in kinds and "error" not in kinds
    assert not target.exists(), "a dry run must not execute the rule"


# --- contig figure -------------------------------------------------

def test_viridis_endpoints_and_monotonic():
    from comparem2.report import viridis, VIRIDIS
    assert viridis(0.0) == VIRIDIS[0]
    assert viridis(1.0) == VIRIDIS[-1]
    assert viridis(-5) == VIRIDIS[0] and viridis(5) == VIRIDIS[-1]  # clamped
    assert all(viridis(i / 20).startswith("#") and len(viridis(i / 20)) == 7
               for i in range(21))


def test_axis_ticks_are_round_numbers():
    """Quartering the max produced labels like '819.114 kb'."""
    from comparem2.report import _ticks, _bp
    assert _ticks(3_276_455) == [0, 1_000_000, 2_000_000, 3_000_000]
    assert [_bp(t) for t in _ticks(3_276_455)] == ["0", "1 Mb", "2 Mb", "3 Mb"]
    assert _ticks(0) == [0.0]
    assert all(t <= 4200 for t in _ticks(4200))


def test_contig_figure_one_rect_per_record():
    from comparem2.report import draw_contigs
    svg = draw_contigs([
        ("A", 40.0, [(1000, 40.0), (500, 38.0)]),
        ("B", 35.0, [(1200, 35.0)]),
    ])
    # 3 contigs + legend swatch
    assert svg.count("<rect") == 4
    assert svg.count("<title>") == 3  # hover text per contig
    assert "A (40.0%)" in svg and "B (35.0%)" in svg


def test_contig_figure_widths_are_proportional():
    from comparem2.report import draw_contigs
    import re as _re
    svg = draw_contigs([("A", 40.0, [(1000, 40.0), (500, 40.0)])])
    widths = [float(w) for w in _re.findall(r'<rect[^>]*width="([0-9.]+)"', svg)]
    assert abs(widths[0] - 2 * widths[1]) < 0.01


def test_contig_figure_tiny_contig_stays_visible():
    from comparem2.report import draw_contigs
    import re as _re
    svg = draw_contigs([("A", 40.0, [(10_000_000, 40.0), (1, 40.0)])])
    widths = [float(w) for w in _re.findall(r'<rect[^>]*width="([0-9.]+)"', svg)]
    assert min(widths) >= 0.4


def test_seqkit_section_has_figure_and_table(tmp_path):
    for sample, recs in (("A", "c1\t1000\t40.0\nc2\t500\t38.0\n"),
                         ("B", "c1\t1200\t35.0\n")):
        d = tmp_path / "samples" / sample / "seqkit"
        d.mkdir(parents=True)
        (d / "contigs.tsv").write_text(recs)
    body = render_report(CATALOGUE, ["seqkit"], tmp_path, Path("db"), ("A", "B")).read_text()
    assert "contig sizes coloured by GC content" in body
    assert "<h3>Assembly statistics</h3>" in body
    # length-weighted mean GC for A: (1000*40 + 500*38)/1500 = 39.33
    assert "39.3" in body


# --- pangenome matrix ----------------------------------------------

def _write_rtab(tmp_path, genomes, rows):
    d = tmp_path / "panaroo"
    d.mkdir(parents=True, exist_ok=True)
    lines = ["Gene\t" + "\t".join(genomes)]
    for i, pattern in enumerate(rows):
        lines.append(f"gene{i}\t" + "\t".join("1" if p else "0" for p in pattern))
    (d / "gene_presence_absence.Rtab").write_text("\n".join(lines) + "\n")


def test_pangenome_matrix_compresses_identical_patterns(tmp_path):
    """Genes with the same presence pattern collapse into one block.

    With N genomes there are at most 2**N - 1 patterns, so this is lossless and
    keeps the SVG small where a per-gene column would emit tens of thousands.
    """
    from comparem2.report import draw_pangenome
    genomes = ["A", "B", "C"]
    # 100 core, 10 shared by A+B, 1 unique to C
    patterns = [((True, True, True), 100), ((True, True, False), 10),
                ((False, False, True), 1)]
    svg = draw_pangenome(genomes, patterns)
    # 3 + 2 + 1 filled cells
    assert svg.count("<rect") == 6
    assert svg.count("<text") == len(genomes) + 2  # labels + two axis captions


def test_pangenome_rare_patterns_stay_visible(tmp_path):
    from comparem2.report import draw_pangenome
    svg = draw_pangenome(["A", "B"], [((True, True), 100000), ((False, True), 1)])
    widths = [float(w) for w in __import__("re").findall(r'width="([0-9.]+)"', svg)]
    assert min(widths) >= 0.6, "a one-gene pattern must not collapse to zero width"


@pytest.mark.parametrize("value,expected", [
    ("2,091", True),    # was left-aligned: float("2,091") raises
    ("116_2", False),   # was right-aligned: float("116_2") == 1162.0
    ("116_2_duplicate", False),
    ("55.3%", True),
    ("584", True),
    ("0", True),
    ("-3", True),
    ("E8202", False),
    ("", False),
])
def test_numeric_alignment_predicate(value, expected):
    """Regression: `float()` misclassified in both directions.

    It rejects thousands separators, so formatted numbers were left-aligned
    while unformatted ones in the same column were right-aligned. And Python
    accepts underscores in numeric literals, so the sample name `116_2` was
    right-aligned as though it were 1162.
    """
    from comparem2.report import _numeric
    assert _numeric(value) is expected


def test_raw_table_escapes_nothing_but_headers():
    from comparem2.report import _raw_table
    out = _raw_table([['<span>x</span>', "1"]], header=["A", "B"], numeric_columns={1})
    assert "<span>x</span>" in out          # markup preserved
    assert '<td class="n">1</td>' in out    # numeric column right-aligned


def test_panaroo_section_has_labelled_subsections(tmp_path):
    """Three unlabelled tables in a row read as one table. They need headings."""
    _write_rtab(tmp_path, ["A", "B"], [(True, True), (True, False)])
    body = render_report(CATALOGUE, ["panaroo"], tmp_path, Path("db"), ("A", "B")).read_text()
    for heading in ("Shared gene content", "Pangenome partitions", "Genes per genome"):
        assert f"<h3>{heading}</h3>" in body


def test_panaroo_overlap_uses_dot_patterns(tmp_path):
    _write_rtab(tmp_path, ["A", "B", "C"], [(True, True, True), (True, False, True)])
    body = render_report(CATALOGUE, ["panaroo"], tmp_path, Path("db"),
                         ("A", "B", "C")).read_text()
    assert "●●●" in body   # present in all three
    assert "●·●" in body   # present in A and C only
    assert "Reading order:" in body


def test_panaroo_section_reports_overlaps(tmp_path):
    genomes = ["A", "B", "C"]
    rows = ([(True, True, True)] * 5 + [(True, True, False)] * 3
            + [(False, False, True)] * 2)
    _write_rtab(tmp_path, genomes, rows)
    body = render_report(CATALOGUE, ["panaroo"], tmp_path, Path("db"),
                         tuple(genomes)).read_text()
    assert "10 gene clusters across 3 genomes" in body
    assert "3 distinct presence patterns" in body
    assert "all genomes" in body      # the 5 core genes
    assert "A, B" in body             # the 3 shared by A and B
    assert "pangenome presence matrix" in body


def test_panaroo_section_survives_single_genome(tmp_path):
    _write_rtab(tmp_path, ["A"], [(True,), (True,)])
    body = render_report(CATALOGUE, ["panaroo"], tmp_path, Path("db"), ("A",)).read_text()
    assert "2 gene clusters across 1 genomes" in body


# --- newick / tree drawing -----------------------------------------

def test_newick_parses_topology_and_lengths():
    """Regression: a delimiter set containing ':' truncated every label.

    The parser stopped at the branch-length separator, so `pos` was left on
    ':', the sibling loop broke after one child, and a four-leaf tree silently
    became a one-leaf chain. It drew without error — only the leaf count showed it.
    """
    from comparem2.report import parse_newick, _leaves
    root = parse_newick("((A:0.1,B:0.2):0.3,C:0.4,D:0.5);")
    assert len(root.children) == 3
    assert [leaf.name for leaf in _leaves(root)] == ["A", "B", "C", "D"]
    assert [leaf.length for leaf in _leaves(root)] == [0.1, 0.2, 0.4, 0.5]
    assert root.children[0].length == 0.3


def test_newick_handles_quotes_and_missing_lengths():
    from comparem2.report import parse_newick, _leaves
    assert [n.name for n in _leaves(parse_newick("('A b',C);"))] == ["A b", "C"]
    assert _leaves(parse_newick("(A,B);"))[0].length == 0.0


def test_draw_tree_labels_every_leaf():
    from comparem2.report import parse_newick, draw_tree, _leaves
    root = parse_newick("((A:0.1,B:0.2):0.3,C:0.4);")
    svg = draw_tree(root)
    assert svg.count("<text") == len(_leaves(root)) == 3
    assert svg.startswith("<svg")


def test_draw_tree_colours_by_cluster():
    from comparem2.report import parse_newick, draw_tree
    svg = draw_tree(parse_newick("(A:0.1,B:0.2,C:0.3);"),
                    {"A": "1", "B": "1", "C": "2"})
    import re
    assert len(set(re.findall(r'fill="(#[0-9a-f]{6})"', svg))) == 2


def test_malformed_newick_degrades_to_text(tmp_path):
    d = tmp_path / "mashtree"
    d.mkdir(parents=True)
    (d / "mashtree.newick").write_text("((((not a tree")
    body = render_report(CATALOGUE, ["mashtree"], tmp_path, Path("db"), ("A",)).read_text()
    assert "mashtree" in body  # section survives


def test_every_tool_has_a_renderer_or_fallback(tmp_path):
    # The fallback must cover anything without a specific section, so that no
    # tool can run and display nothing — v2's failure with 5 tools.
    from comparem2.report import SECTIONS, _fallback
    for tool in CATALOGUE:
        assert SECTIONS.get(tool.name, _fallback) is not None
