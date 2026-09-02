# The bioconda recipe

`meta.yaml` here is a **draft**. The authoritative copy lives in
[bioconda-recipes](https://github.com/bioconda/bioconda-recipes) at
`recipes/comparem2/meta.yaml`, which already exists — publishing v3 is a
version bump to a recipe `cmkobel` already maintains, not a new-recipe review.

**The reasoning lives in this file, not in the recipe.** Recipes in
bioconda-recipes are terse: of forty `noarch: python` recipes sampled from a
clone, most carry no comments at all and the wordiest had nine lines in
fifty-three. Anything that needs explaining is explained here.

## Why the autobump bot cannot do this

Bioconda's auto-updater changes the version, the source URL, the checksum and
the build number. Its own documentation says it "does not yet know how to
monitor changed dependencies", and it does not touch build scripts or test
sections. This bump changes all three: `noarch: generic` becomes
`noarch: python` with a pip install and entry points, the run dependencies are
replaced wholesale, and the test section is new.

So it needs a hand-written PR — and worse, **the bot is a hazard here**. Tagging
`v3.0.0` may prompt it to open a routine-looking version bump that keeps v2's
`snakemake-minimal <8`, `mamba <2` and `noarch: generic`, which would build a
package that installs nothing usable. Get the real PR in first, or close the
bot's.

The current published recipe is **2.16.2**, build 0, `noarch: generic`, with run
dependencies `snakemake-minimal <8`, `pulp <2.8`, `python <3.12`, `mamba <2`,
`pandas`.

## What changes from v2

| | v2.16.2 | v3.0.0 |
| --- | --- | --- |
| build | `noarch: generic` — shipped the `./comparem2` bash launcher | `noarch: python`, `pip install`, two entry points |
| snakemake | `snakemake-minimal <8` | `snakemake-minimal >=9,<10` |
| conda frontend | `mamba <2` | `conda` — Snakemake 9 deprecated alternative frontends |
| python | `<3.12` | `>=3.11` |
| dropped | `pulp`, `pandas` | neither is imported by v3 |
| added | | `textual`, the two executor plugins |
| kept | `run_exports` pinning at major version | the PR template asks for it; semantic versioning, so `max_pin="x"` |

The tool set is *not* in either recipe. v2 deployed 25 environments at run
time; v3 deploys 14. That is the same model, and it is why a recipe for a
13-tool pipeline is this short.

## Release steps

1. **Version.** Set `__version__` in `src/comparem2/__init__.py` and
   `version` in `pixi.toml` to `3.0.0` — a unit test fails if they disagree.
   Update `citation.cff` (`version: 2` today) at the same time.
2. **Tag and push.** `git tag v3.0.0 && git push origin v3.0.0`. GitHub
   generates the source tarball the recipe fetches.
3. **Hash it.**
   ```bash
   curl -sL https://github.com/cmkobel/CompareM2/archive/refs/tags/v3.0.0.tar.gz \
     | shasum -a 256
   ```
4. **Fork and edit.** In a fork of bioconda-recipes, replace
   `recipes/comparem2/meta.yaml` with this file, version and sha256 filled in.
5. **Lint.** `bioconda-utils lint --packages comparem2` (or let CI do it).
6. **PR.** CI builds the package and the BioContainer. Then comment
   `@BiocondaBot please add label` to request the merge label.

## About the container

Bioconda builds a quay.io BioContainer for every package automatically, so one
appears whether or not we want it. It contains **the pipeline and no analysis
tools**: it can render a report and run `--dry-run`, and a real run inside it
would deploy the tool environments on first use, needing network and a writable
`--conda-prefix`.

**A hand-built image with all thirteen tools is not planned** (decided
2026-09-02). pixi and conda both install this well enough that the image would
be a third thing to keep in step, and it could not be one environment anyway —
CheckM2 pins DIAMOND 2.1.x against Bakta's 2.2.x, so it would be two
environments inside one image. If it is ever needed, that is the shape.
