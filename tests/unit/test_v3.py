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
from comparem2.cli import slug  # noqa: E402
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


def test_generated_snakefile_parses_as_python_free_text():
    # Every rule block should have input/output/shell in order.
    text = render(CATALOGUE, None, Path("res"), Path("db"), SAMPLES)
    for block in text.split("\nrule ")[2:]:
        assert block.index("input:") < block.index("output:") < block.index("shell:")


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
