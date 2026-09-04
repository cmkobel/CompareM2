# intbitset 4.1.2 + osx-arm64 builds

**Sent 2026-09-04:
[conda-forge/intbitset-feedstock#21](https://github.com/conda-forge/intbitset-feedstock/pull/21)**,
from `cmkobel:osx-arm64-and-4.1.2`. Two commits: the recipe change, and the
rerender separately so the human change stays reviewable on its own. The local
branch is `~/postdoc/cm2-macos/upstream/intbitset-feedstock`.

Everything below is the working note the PR was written from; the PR body
itself is much shorter. Two changes were made after the first push, both from
the conda-forge linter, which failed the PR on the first of them:

- **`{{ stdlib('c') }}` in `build:`** — now required of any recipe using a
  compiler ([conda-forge.github.io#2102](https://github.com/conda-forge/conda-forge.github.io/issues/2102)).
  This changes the variant config, so it needed a second rerender; py3.13 was
  rebuilt locally afterwards and still passes. The build string moves from
  `py313h2f2c7d1_0` to `py313hd55261d_0` because `c_stdlib` joins the hash.
- **`files.pythonhosted.org` instead of `pypi.io`** — the linter's advisory,
  not a failure.

The linter is green after both: *"found it was in an excellent condition."*

Measured 2026-09-03 on macOS 26.2 / Apple silicon. This is the first of the two
changes CompareM2 needs for panaroo on macOS; the other is
[bioconda-panaroo-pr.md](bioconda-panaroo-pr.md). Neither is sufficient alone.

## Why

`intbitset` has no `osx-arm64` build. conda-forge ships 2.4.0, 2.4.1 and 3.0.2
for `linux-64`, `osx-64` and `win-64` only. Two consequences, and only the
second is about macOS:

- **panaroo cannot be installed on Apple silicon at all**, because `intbitset`
  is in its run dependencies. That costs three tools rather than one: panaroo's
  core gene alignment is the input to both `snp-dists` and `fasttree`.
- **`python 3.12` + panaroo is unsatisfiable on `linux-64` today**, because
  intbitset has no py3.12 build. This is what holds CompareM2's default
  environment at python 3.11.16, on Linux, right now.

The feedstock is stalled rather than hostile: last human merge **2023-09-19**,
**seven** open autotick-bot PRs (#14–#20), sole maintainer `gtonkinhill` — who
also wrote panaroo.

**But the maintainer is not absent, only this feedstock is.** Checked
2026-09-04: his public activity runs to **2026-09-03**, he merged two panaroo
PRs on 2026-07-02, and he pushed to `bioconda/bioconda-recipes` on 2026-08-23.
So the three-year feedstock silence is not evidence that a ping will go
unanswered, and it argues against escalating to `@conda-forge/core` early.
`.github/CODEOWNERS` is `* @gtonkinhill`, so opening the PR requested his
review automatically — he has already been notified once, without anyone
@-mentioning him.

## The cause is one missing line, and no rerender will fix it

`conda-forge.yml` carries no `provider` entry. Per conda-forge's current
[conda-forge.yml reference](https://conda-forge.org/docs/maintainer/conda_forge_yml/),
x86_64 platforms are enabled by default and the others are not — an empty
`provider` is equivalent to `osx_arm64: None`. So conda-smithy emits no arm
configs however often it runs, which is why bot PR #20 rerendered on
2026-04-28 and produced zero.

Adding

```yaml
provider:
  osx_arm64: default
```

and rerendering emits `osx_arm64_python3.{11,12,13,14}.yaml`. **Verified
2026-09-03** against a local clone with conda-smithy **2026.9.1** and
conda-forge-pinning **2026.09.03.11.37.26**.

## Why the PR has to carry the version bump too

`3.0.2` does not compile on py3.12. A rerender at 3.0.2 would emit configs for
3.11–3.14 and three of the four would fail, so the provider change alone lands
red.

`4.1.2` builds. From the sdist on `osx-arm64`, all four Pythons produce
`Mach-O 64-bit bundle arm64`, and the upstream test suite passes **13,930
tests** — 4.81 s on py3.13, 6.08 s on py3.14. Upstream's changelog from 3.0.2
to 4.1.2 is Python-version support, dropped EOL Pythons and regenerated C only;
there is no API change, and panaroo uses the constructor, `add`, `discard`,
`remove`, `copy`, `intersection`, `|=`, `&`, `len`, `in` and iteration.

This also supersedes five stale bot version PRs (#15 3.1.0, #17 4.0.0, #19
4.1.0, #20 4.1.2) and the two rebuild PRs (#14 py3.12, #16 py3.13, #18 py3.14).

## The other recipe changes, and how sure each one is

| change | why | confidence |
| ------ | --- | ---------- |
| drop `six` from `run:` | 4.1.2 declares no runtime dependencies. PyPI `requires_dist` lists only the `tests` extra, and the only `six` left in the tree is a comment at `tests/test_intbitset.py:33` | **certain** — leaving it in is simply wrong |
| add `setuptools` to `host:`, add `--no-build-isolation` | the current conda-forge pattern for a setup.py-only project; without it pip's isolated build environment has to reach PyPI for setuptools, which CI does not allow, and py3.12+ no longer bundles it | **plausible, not verified** — see below |
| `about: home:` → `inveniosoftware-contrib` | upstream moved orgs; the recipe still points at `inveniosoftware` | incidental, drive-by |

**On the setuptools change.** Bot PR #20 is red on `cp313` across `linux_64`,
`osx_64` and `win_64` while 3.10/3.11/3.12 are green. Its logs have expired
(the Actions API returns HTTP 410), so the actual failure is *not known*. The
setuptools/`--no-build-isolation` pair is the most likely explanation and is
correct regardless. **This draft must not claim it as the diagnosed cause** —
what is verified is that the recipe as written builds, not why #20 does not.

## The recipe builds, on all four arm configs

Run locally 2026-09-03 with conda-build 26.7.1 against the four `.ci_support`
configs the rerender emits. Every one exited 0 and produced a real `.conda`,
and the `imports: intbitset` test ran:

```
intbitset-4.1.2-py311hcc37be2_0.conda
intbitset-4.1.2-py312h72574e9_0.conda
intbitset-4.1.2-py313h2f2c7d1_0.conda   <- the config #20 fails on
intbitset-4.1.2-py314h1686f52_0.conda
```

Peak memory 228.1 MB. This is a native arm build on Apple silicon, which is
what conda-forge's runners now do for `osx_arm64: default`.

## Together with the panaroo change, this is sufficient

Not "should be" — checked. With those four artifacts in a local channel and
bioconda's panaroo recipe carrying only the one-line change in
[bioconda-panaroo-pr.md](bioconda-panaroo-pr.md), CompareM2's own specs
`panaroo>=1.5`, `snp-dists>=1.2.0` and `fasttree>=2.2.0` solve on `osx-arm64`
from conda alone — **278 packages**, no pip, no `--no-deps`, no local pin
beyond the intbitset build itself. panaroo 1.8.0 `py_1`, intbitset 4.1.2
`py313h2f2c7d1_0`, **python 3.13.15**, and zero mentions of prokka.

That python version is the second half of the point: 4.1.2 lifts the cap that
holds CompareM2's Linux environment at 3.11.16.

## What this draft still needs before it is sent

- A decision on the merge path. conda-forge's documented route for a stalled
  feedstock is: open the PR, ping the maintainer, wait about a week, then ask
  `@conda-forge/core` to merge. `@conda-forge-admin, please add user @cmkobel`
  is the other route, but conda-forge's own FAQ says that is *not* the
  recommended entry point — lead with a PR.
- `@conda-forge-admin, please rerender` as the first PR comment, rather than
  committing the local rerender. The bot rerenders against current pinnings,
  which keeps the diff to two files and avoids shipping a rerender that is
  stale by the time anyone reviews it.

## Rebuilding the numbers

```bash
# arm builds from sdist, all four Pythons, plus the upstream suite
for v in 3.11 3.12 3.13 3.14; do
  uv venv --python $v v$v
  uv pip install --python v$v/bin/python --no-binary :all: --no-cache intbitset==4.1.2
  v$v/bin/python -c "from intbitset import intbitset; print(list(intbitset([9,2,5])))"
done
curl -sL https://files.pythonhosted.org/packages/source/i/intbitset/intbitset-4.1.2.tar.gz | tar -xz
cd intbitset-4.1.2 && ../v3.13/bin/python -m pytest -q

# the rerender claim
conda-smithy rerender --no-check-uptodate && ls .ci_support/ | grep osx_arm64
```
