# CLAUDE.md

Guidance for Claude Code working in this repository.

**This branch (`v3`) is a rewrite.** v2 has been removed from it — no
`workflow/`, no `dynamic_report/`, no R, no 25 conda environments, no
`./comparem2` launcher, no `config/config.yaml`. If you are looking for any of
those, they are on `master`. Do not reintroduce v2 patterns here.

Three files carry the context, and they are deliberately separate:

- **`DESIGN.md`** — what v3 is and why it is shaped this way. No dates. Read
  this before making a design decision, especially the *Rules that must not be
  quietly undone* section.
- **`DECISIONS.md`** — the dated log of how it got here, including decisions
  reversed within a day and a *What went wrong* section. Read this before
  undoing one of those rules, or to check whether an idea has already been
  tried.
- **`STATUS.md`** — what is currently true of a real run: which tools have
  actually been executed, measured database sizes, where things live on
  thylakoid, and what is known broken.

## What this is

CompareM2 v3: a Snakemake-driven pipeline that runs 13 analysis tools over a set
of microbial assemblies and produces one self-contained HTML report.
**Linux-only** — the tools are `linux-64`, so `pixi install` will not work on
macOS. Unit tests are pure Python and run anywhere.

**Publication (describes v2):** Kobel et al. "CompareM2 is a genomes-to-report
pipeline for comparing microbial genomes." *Bioinformatics* 41(9), btaf517
(2025). doi:10.1093/bioinformatics/btaf517.

## Architecture

Declarative specs, generated workflow. There is no hand-written Snakefile.

```
src/comparem2/
  tools.py      the contract: Tool, Database, Context, Registry, Scope
  catalogue.py  the 13 Tool specs — command lines and outputs. THE source of truth.
  guidance.py   what each tool does and how to read its output, for the report
  snakefile.py  generates a Snakefile (and envs/*.yaml) from the specs
  cli.py        argument parsing, input canonicalisation, hands off to Snakemake
  runner.py     drives Snakemake via its API, emits structured events
  tui.py        Textual interface over those events
  report.py     renders the HTML report
```

The flow: `cli.py` canonicalises every input to
`<workdir>/samples/<sample>/<sample>.fna`, `snakefile.py` renders a Snakefile
from `CATALOGUE`, Snakemake executes it, `report.py` renders the result.

**Adding or changing a tool means editing `catalogue.py` and `guidance.py`.**
Nothing else should need to know a tool exists. A test enforces that the two
stay in step.

### Two things that follow from the design

- **Wildcard rendering.** A `Context` built with `sample="{sample}"` and
  `threads="{threads}"` makes every path and command come out templated, which
  is how the same code renders both a real command and a Snakemake rule. Inside
  a `shell:` block, wildcards are `{wildcards.sample}`, not `{sample}` — the
  bare form parses fine and fails at runtime.
- **Declared outputs.** A tool declares its outputs rather than writing wherever
  it likes. That is what buys resumability, progress reporting, and the report's
  knowledge of what to render.

## Development

```bash
# Unit tests — pure Python, no pixi needed, runs on macOS
python -m pytest tests/unit -q          # or: pixi run pytest

# On linux, with the tool environment
pixi install
pixi run cm2 --help
pixi run test-fast                      # 4 genomes, no databases required
pixi run cm2 <assemblies>... --dry-run
```

Isolated tools need an absolute launcher, because Snakemake's shell does not
inherit an interactive PATH:

```bash
pixi run cm2 <assemblies>... \
  --isolated-launcher "/home/thylakoid/.pixi/bin/pixi run -e {tool}"
```

Moving a pixi project invalidates its environments — conda bakes the absolute
prefix into shebangs and RPATHs — so `rm -rf .pixi && pixi install` after any
move.

### Testing

`tests/unit/test_v3.py`, 75 tests, ~0.7 s. This is the primary instrument: the
codebase is a generator, and a wrong wildcard produces a Snakefile that parses
cleanly and builds the wrong DAG, which an end-to-end run catches slowly if at
all. CI (`.github/workflows/unit.yaml`) runs it on 3.11–3.13 without pixi.

Test genomes are shipped zipped under `tests/`; `pixi run unpack` extracts them.
Unpacked `.fna` files are gitignored on purpose.

`tests/E._faecium/` contains `116_2.fna` and `116_2 duplicate.fna` — the same
genome twice, and a filename with a space. Any tool that treats them differently
is wrong, and it is the standing cross-check: 0.00000 mash distance, 100.00%
ANI, 0 SNPs, identical CDS counts.

## Conventions that matter

- **Pin a minimum version for every tool.** An unconstrained bioconda spec lets
  the solver reach back years to satisfy some other package's constraint. Both
  bakta and panaroo resolved to builds that installed cleanly and crashed on
  first use. **`pixi install` succeeding says nothing about whether the pipeline
  works.**
- **Commands are argument lists, never shell strings.** A tool that writes to
  stdout declares `stdout_to_output=True` and the redirect is added by whatever
  runs it.
- **`isolated=True` is an exception that must carry its reason in the spec.**
  Exactly one tool has it (checkm2, which pins DIAMOND 2.1.x against bakta's
  2.2.x). v2 reached 25 environments by making this the default.
- **Databases declare a measured size**, taken from `content-length`. `None`
  means unmeasured — never guess, and say "unmeasured" when totalling.
- **Passthrough parameters**: `--set tool--flag=value` on the CLI, `params` on
  the spec. v2 spelled this `set_tool--flag: value` in `config.yaml`.
- **Report guidance is quoted from papers and checked.** Numbers in
  `guidance.py` were copied from the tool's own paper and verified against the
  PDF text. If you add a number there, it needs the same treatment; if a claim
  is methodological caution rather than a paper's finding, say so in the
  sentence.

## Verification status

`STATUS.md` carries the table of which tool command lines have actually been
executed. **It tracks execution, never installation** — a clean `pixi install`
says nothing about whether a tool runs. Commands drafted against documented
interfaces and never executed are the standing risk here; 12 of 13 are now
verified, and GTDB-Tk is the exception.
