# Draft: add linux-aarch64 to intbitset-feedstock

**Sent 2026-10-05 as [intbitset-feedstock#25](https://github.com/conda-forge/intbitset-feedstock/pull/25)**, from `cmkobel/intbitset-feedstock` branch
`linux-aarch64`: two commits on conda-forge `main` `17f282b`: `29f8860`
*Enable linux-aarch64 builds* and `8542061`, the rerender with conda-smithy
2026.9.23. PR text: `~/postdoc/cm2-macos/linux-aarch64/pr-body.md`. Destination: a PR against
[conda-forge/intbitset-feedstock](https://github.com/conda-forge/intbitset-feedstock),
the same feedstock as [intbitset-feedstock-pr.md](intbitset-feedstock-pr.md).
Working material, outside git: `~/postdoc/cm2-macos/linux-aarch64/` — the
feedstock clone with the change committed and rerendered, the built `.conda`,
the solve probes and the run outputs.

## What it buys

panaroo on Linux on ARM. intbitset is the **only** thing missing there:
conda-forge builds it for linux-64, osx-64, osx-arm64 and win-64, and has never
built linux-aarch64. Unlike osx-arm64, prokka is not in the way — so
[bioconda-panaroo-pr.md](bioconda-panaroo-pr.md) (#69901) is not needed for
this platform, and this change alone is sufficient.

Solve-only, `pixi lock` for linux-aarch64 against conda-forge and bioconda:
panaroo's recipe `run:` list minus `intbitset`, plus `snp-dists>=1.2.0` and
`fasttree>=2.2.0`, solves **with** prokka (453 packages) and without it (147),
python 3.13.15 both ways. `intbitset` on its own: `No candidates were found`.

Windows on ARM would reach this through WSL2, which runs Linux of the host's
architecture — general knowledge, not tested here.

## The change

One line in `conda-forge.yml`, next to #21's:

```diff
 provider:
+  linux_aarch64: default
   osx_arm64: default
```

`conda-smithy rerender` then adds four configs, `linux_aarch64` py3.11–3.14,
built natively on GitHub's `ubuntu-24.04-arm` runners in
`quay.io/condaforge/linux-anvil-aarch64:alma10`. No recipe change; upstream
already ships manylinux2014 aarch64 wheels of 4.1.2 for cp310–cp314.

## Evidence that the result works

Built on the laptop (Apple silicon, Docker 29.3.1, native aarch64 containers)
in the image CI uses, by the feedstock's own `.scripts/build_steps.sh`. All four
configs from commit `8542061` exit 0 and pass the recipe's import test, 49–62 s
each: `py311h5a00dd8_0`, `py312h84d9e0c_0`, `py313hc935821_0`,
`py314h35c850b_0` — the same build strings as a first round from the
conda-smithy 2026.9.1 rerender, whose aarch64 configs were identical.

Then, in `condaforge/miniforge3` on arm64: bioconda's **published** panaroo
1.8.0 — prokka line and all — with `snp-dists` and `fasttree`, the local
intbitset being the only unpublished package. 456 packages, created in 59.3 s.
Same command and prodigal GFFs as the macOS runs in `../STATUS.md`:

| | |
| --- | --- |
| panaroo wall clock | 77.6 s |
| gene clusters | 3,569 |
| duplicate pair, clusters differing | **0 of 3,569** |
| `snp-dists` 1.2.0, duplicate pair | **0** SNPs (4,111 to E8202) |
| `FastTree -nt -gtr` duplicate pair | branch length 0.0 |
| `gene_presence_absence.Rtab` | sha256 `0d0bebcc…`, **byte-identical** to the 09-03 and 09-28 osx-arm64 runs |

`build-locally.py` itself fails on a Mac, and not because of the recipe: it
bind-mounts the feedstock from case-insensitive APFS, and ncurses' terminfo
holds names that differ only by case (`[Errno 22] … terminfo/32/2621A`), so the
build environment never finishes installing. Running `build_steps.sh` with the
feedstock in a Docker volume avoids it.

Aside: prokka installs on linux-aarch64, but `tbl2asn-forever`'s binary there
is an x86-64 ELF (`file`), so prokka's tbl2asn step cannot run natively. Not
executed. panaroo never calls it.

## Windows is not reachable this way

`pixi lock` fails on win-64 (first at `python-edlib`) and win-arm64 (first at
`numba`). mafft, prank, mash, cd-hit and python-edlib have no Windows build on
either architecture — bioconda does not build for Windows — and win-arm64 lacks
numba and intbitset as well.

## What is left

- Gerry Tonkin-Hill is the feedstock's sole maintainer, as for #21, and has
  not been told it is open.
- Build number left at 0 — Carl's call, 2026-10-05. conda-forge's PR template
  asks for "Bumped the build number (if the version is unchanged)"; the PR text
  leaves that box unchecked and offers to bump. Not checked what CI's upload
  does with the existing 0-builds of the other twelve configs.
