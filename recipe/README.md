# The bioconda recipe

`meta.yaml` here is a **draft**. The authoritative copy lives in
[bioconda-recipes](https://github.com/bioconda/bioconda-recipes) at
`recipes/comparem2/meta.yaml`, which already exists — publishing v3 is a
version bump to a recipe `cmkobel` already maintains, not a new-recipe review.

**The reasoning lives in this file, not in the recipe.** Recipes in
bioconda-recipes are terse: of forty `noarch: python` recipes sampled from a
clone, most carry no comments at all and the wordiest had nine lines in
fifty-three. Anything that needs explaining is explained here.

## The autobump bot: a hazard for v2→v3, and the right tool afterwards

Bioconda's auto-updater changes the version, the source URL, the checksum and
the build number. Its own documentation says it "does not yet know how to
monitor changed dependencies", and it does not touch build scripts or test
sections.

**The v2→v3 bump changed all three**, which is why it needed a hand-written PR:
`noarch: generic` became `noarch: python` with a pip install and entry points,
the run dependencies were replaced wholesale, and the test section was new. The
bot was an active hazard then — it would have opened a routine-looking version
bump keeping v2's `snakemake-minimal <8`, `mamba <2` and `noarch: generic`,
building a package that installs nothing usable. It did exactly that:
[PR #68821](https://github.com/bioconda/bioconda-recipes/pull/68821), opened
minutes after the `v3.0.0` tag and now closed.

**From v3.0.0 onward the bot is what does the release.** The published recipe is
already the v3 shape, so a version-and-checksum bump is the whole change, which
is precisely what the bot is for — and it is enabled. **Pushing the tag is the
release**; the PR arrives on its own. Only a change to the dependencies, the
build script or the test section needs a hand-written PR again, and this file is
where to draft it.

### It arrives on its own, and then waits three days

Learned 2026-09-07, when v3.1.0's PR
[#68842](https://github.com/bioconda/bioconda-recipes/pull/68842) had sat open
and green since 2026-09-04T20:15:57Z. Nothing was wrong with it. Bioconda's
[`.mergify.yml`](https://github.com/bioconda/bioconda-recipes/blob/master/.mergify.yml)
is one 19-line rule that approves and squash-merges an autobump PR when every
commit is BiocondaBot's, the three build checks pass, and:

```yaml
- created-at<3 days ago   # i.e. the PR must be older than three days
```

So **the release lands three days after the PR is opened**, and no label or
comment is needed — the 2.x bumps that merged in 36 minutes predate this rule.
Two consequences worth knowing before a tag:

- **The clock runs on the PR, not on the commit.** The bot pushes each new
  version onto the same `bump/comparem2` branch and rewrites the PR title —
  verified on a live example, pybiolib #68885, created 08:14 for 1.4.474 and
  carrying 1.4.481 by 12:11 the same day, one `created-at`. A tag pushed while
  an unmerged bump PR is open therefore inherits that PR's age instead of
  starting a new three days. The bot picked up the v3.1.0 tag about nine
  minutes after it was pushed, and runs roughly hourly.
- **Do not push to that branch by hand.** `commits[*].author=BiocondaBot` is a
  condition, so one commit of ours — or a maintainer's merge of master, as on
  dotmatch #68860 — disqualifies the PR from auto-merge and it needs a human
  after that.

Before v3.0.0 the published recipe was **2.16.2**, build 0, `noarch: generic`,
with run dependencies `snakemake-minimal <8`, `pulp <2.8`, `python <3.12`,
`mamba <2`, `pandas`.

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
time; v3 deploys 6, grouped by dependency ecosystem rather than one per tool.
That is the same model, more cheaply, and it is why a recipe for a
fourteen-tool pipeline is this short.

## Release steps

Steps 1–3 are the release. Steps 4–6 are only for a release that also changes
the recipe's shape; otherwise autobump does them.

1. **Version.** Set `__version__` in `src/comparem2/__init__.py` and `version`
   in `pixi.toml` to the new number — a unit test fails if they disagree.
   Update `citation.cff` (`version` and `date-released`) at the same time.
2. **Tag and push.** `git tag vX.Y.Z && git push origin vX.Y.Z`. GitHub
   generates the source tarball the recipe fetches, and autobump opens the
   bump PR from the tag — then holds it three days, so check whether an older
   bump PR is already open and can carry this version instead.
3. **Hash it**, to check the bot against and to update the draft here.
   ```bash
   curl -sL https://github.com/cmkobel/CompareM2/archive/refs/tags/vX.Y.Z.tar.gz \
     | shasum -a 256
   ```
4. **Fork and edit**, *if the recipe itself changed*. In a fork of
   bioconda-recipes, replace `recipes/comparem2/meta.yaml` with this file,
   version and sha256 filled in.
5. **Lint.** `bioconda-utils lint --packages comparem2` (or let CI do it).
6. **PR.** CI builds the package and the BioContainer. Then comment
   `@BiocondaBot please add label` to request the merge label.

> **Do not move a tag once a recipe is building against it.** `v3.0.0` was
> force-moved inside the two-minute window in which bioconda was building, and
> whether the build fetched the old bytes or a cached archive is not
> determinable after the fact. The check is one line: read the PR's
> `.merged_at` first.

## About the container

Bioconda builds a quay.io BioContainer for every package automatically, so one
appears whether or not we want it. It contains **the pipeline and no analysis
tools**: it can render a report and run `--dry-run`, and a real run inside it
would deploy the tool environments on first use, needing network and a writable
`--conda-prefix`.

**A hand-built image with all fourteen tools is not planned** (decided
2026-09-02). pixi and conda both install this well enough that the image would
be a third thing to keep in step, and it could not be one environment anyway —
CheckM2 pins DIAMOND 2.1.x against Bakta's 2.2.x, so it would be two
environments inside one image. If it is ever needed, that is the shape.
