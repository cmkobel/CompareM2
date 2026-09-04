"""Six *Enterococcus faecium* plasmids, shipped so that `--demo` needs nothing.

Plasmids rather than chromosomes because of what they cost: 452 KB compressed
against 3.4 MB for the four-chromosome set, on a package that is otherwise
109 KB. That buys a demo every user can run without downloading a genome, and
the four analyses `--demo` runs — contig statistics, an alignment-free tree,
clusters and all-against-all ANI — are as meaningful on plasmids as on
chromosomes.

**They are plasmids, so the whole catalogue is not meaningful on them.** CheckM2
would report a completeness near zero, correctly and uselessly, because it is
looking for a bacterial chromosome's marker genes. `--demo` therefore fixes its
own tool list rather than running everything.

The seventh input is the sixth one again under another name. It costs no bytes,
it exercises the space in `116_2 duplicate.fna` through sample-name
canonicalisation, and it gives the report something worth looking at: two
inputs that must come out at 0.00000 mash distance and 100.00% ANI. That pair is
the project's standing cross-check — see `tests/E._faecium`.
"""

from __future__ import annotations

import shutil
import zipfile
from importlib import resources
from pathlib import Path

ARCHIVE = "plasmids.zip"

# Named rather than counted, so a change to the archive is a test failure and
# not a quietly different demo.
PLASMIDS = ("116_2.fna", "Dallas_55.fna", "E8202.fna", "EF_VRE.fna",
            "ISMMS_VRE_1.fna", "VB3240.fna")

# The duplicated input: which file, and the name the copy takes. The space is
# deliberate.
DUPLICATE_OF = "116_2.fna"
DUPLICATE_AS = "116_2 duplicate.fna"

# What `--demo` runs. Every one of these needs no database and no network, which
# is the whole point: the demo must not download 62.5 GB to say hello.
TOOLS = ("seqkit", "mashtree", "treecluster", "skani")


def extract(into: Path) -> list[Path]:
    """Write the demo assemblies into `into` and return their paths.

    Overwrites rather than skipping. A demo that silently reused a directory
    someone had edited would be a demo of that edit.
    """
    into.mkdir(parents=True, exist_ok=True)
    archive = resources.files(__package__) / ARCHIVE
    with resources.as_file(archive) as path:
        with zipfile.ZipFile(path) as zf:
            names = set(zf.namelist())
            missing = [n for n in PLASMIDS if n not in names]
            if missing:
                raise SystemExit(
                    f"the bundled demo archive is missing {', '.join(missing)}; "
                    "this is a packaging fault, not something you did")
            for name in PLASMIDS:
                (into / name).write_bytes(zf.read(name))

    shutil.copyfile(into / DUPLICATE_OF, into / DUPLICATE_AS)
    return [into / n for n in (*PLASMIDS, DUPLICATE_AS)]
