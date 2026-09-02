"""Small steps a rule needs around a tool's own command.

Most tools are one command. GTDB-Tk is not: it needs a batchfile written before
it runs, and it writes its bacterial and archaeal results to separate files that
have to become the one table the pipeline declared. The first is a file with
content, the second is text handling — neither is expressible as an argument
list, which is the discipline every command in `catalogue.py` follows.

A file with content is `Tool.files`, rendered from the Context and written by
`prepare()`. The reshaping afterwards is `Tool.post`, and it comes here rather
than into an `awk` invocation because this is testable and that is not.

**Invoked through an absolute `sys.executable`**, never a bare `python`. Under a
conda deployment the rule's environment holds the tool and its interpreter, not
CompareM2, so `python -m comparem2.steps` would find the wrong Python. The
absolute path is the one that generated the workflow, which is the one that has
this module.
"""

from __future__ import annotations

import argparse
import glob
from pathlib import Path


def merge_tsv(out: Path, patterns: list[str]) -> None:
    """Concatenate TSVs that share a header into one, header kept once.

    Globs rather than named paths, deliberately. GTDB-Tk's documentation does
    not say whether `classify_wf` leaves its summaries in `--out_dir` or in a
    `classify/` subdirectory below it, and the answer has moved between major
    versions, so the caller passes both candidates and whichever exists wins.

    A domain with no genomes in it is a header-only file, or absent entirely —
    an all-bacterial set has no `ar53` summary — and both contribute nothing.
    Matching *nothing at all* is an error: silently writing an empty table
    would let the rule report success for a tool that classified no genomes.
    """
    paths = sorted({p for pattern in patterns for p in glob.glob(pattern)})
    if not paths:
        raise SystemExit(f"merge-tsv: no input matched {' '.join(patterns)}")

    header: str | None = None
    rows: list[str] = []
    for path in paths:
        lines = Path(path).read_text().splitlines()
        if not lines:
            continue
        if header is None:
            header = lines[0]
        elif lines[0] != header:
            raise SystemExit(
                f"merge-tsv: {path} has a different header from {paths[0]}")
        rows += [line for line in lines[1:] if line.strip()]

    if header is None:
        raise SystemExit(f"merge-tsv: every input was empty: {' '.join(paths)}")

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join([header, *rows]) + "\n")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="python -m comparem2.steps",
        description="Steps a generated rule runs around a tool's command.")
    sub = p.add_subparsers(dest="step", required=True)

    merge = sub.add_parser("merge-tsv", help="concatenate TSVs sharing a header")
    merge.add_argument("--out", type=Path, required=True)
    merge.add_argument("patterns", nargs="+",
                       help="paths or globs; at least one must match")

    args = p.parse_args(argv)
    if args.step == "merge-tsv":
        merge_tsv(args.out, args.patterns)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
