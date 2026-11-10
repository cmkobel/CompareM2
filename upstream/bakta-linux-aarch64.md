# Draft: bakta on linux-aarch64 installs, then crashes on the first genome

**Unsent, and the destination is not settled.** Measured 2026-11-10 on the
laptop, an M1 running Asahi Linux (linux-aarch64, 16 KiB pages). Working
material, outside git: `~/cm2-aarch64/` — the failed pipeline run (`full13/`),
the bakta 1.11.4 runs with pyhmmer 0.11.4 (`bakta1114-pyhmmer011/`), and the
patched-prefix run (`full13b/`).

## What happens

`catalogue.py`'s `ANNOTATION_ENV` (`bakta>=1.10`, `ncbi-amrfinderplus>=4.2.7`)
solves on linux-aarch64 to **bakta 1.11.4 with pyhmmer 0.12.3**. On linux-64
and osx-arm64 the same specs give bakta 1.12.2. The pipeline's bakta rule then
fails on every genome, after gene calling:

```
  File ".../bakta/features/orf.py", line 25, in detect_spurious
    orf = orf_by_aa_digest[hit.name.decode()]
AttributeError: 'str' object has no attribute 'decode'. Did you mean: 'encode'?
```

## Two causes, one in each of two recipes

1. **bakta 1.12 needs `blast >=2.17`, and bioconda has no linux-aarch64 build
   of blast 2.17.0.** The recipe has the platform commented out:
   `# CircleCI arm.large runner times out as the build takes too long`.
   bioconda-recipes#62276 ("re-add ARM builds", merged 2026-01-30) brought
   back osx-arm64 only, and its thread ends with the aarch64 build still
   timing out. The infrastructure issue is bioconda-recipes#58217, open. The
   newest linux-aarch64 blast is 2.16.0 build 5. So the solver steps back to
   the newest bakta that accepts 2.16, which is 1.11.4.
2. **bioconda's bakta 1.11.4 allows any `pyhmmer >=0.10.15`**, with no
   ceiling. pyhmmer 0.12.0 (2026-01-21) made `name` and the other text
   attributes `str` instead of `bytes`, a breaking change listed as such in
   its changelog. bakta followed in 1.12.0 ("Update PyHMMER to v0.12.0",
   bakta#421): `orf.py` line 25 is `hit.name.decode()` at v1.11.4 and
   `hit.name` at v1.12.0. So 1.11.4 installs cleanly beside the pyhmmer that
   breaks it.

## bakta 1.11.4 itself is fine

With `pyhmmer <0.12`, the solver picks 0.11.4, and bakta 1.11.4 on linux-aarch64
runs the catalogue's command against the db v6.0 light database the pipeline
fetches (`software-min` 1.11). Four genomes, `tests/E._faecium`, exit 0 at
93–110 s each, 8 threads. Against thylakoid's bakta 1.12.1 (linux-64,
`results_two_envs`, 2026-09-03):

- **every count in the `.txt` summary identical** for all four genomes, from
  tRNAs to oriTs. 116_2: 68 tRNA, 1 tmRNA, 18 rRNA, 13 ncRNA, 34 ncRNA
  regions, 2,587 CDS + 1 sORF, 305 hypotheticals, 2 oriC. E8202 3,180 CDS,
  SRR24 2,815.
- **the `.faa` byte-identical** for all four genomes, and the duplicate pair
  is identical to each other.

## Routes, cheapest first

- **A ceiling on bakta <=1.11's pyhmmer, as a repodata patch.**
  `pyhmmer <0.12` on bakta 1.10 and 1.11, which all call `hit.name.decode()`
  (read at v1.10.0, v1.10.4, v1.11.0 and v1.11.4). With that,
  linux-aarch64 resolves a bakta that works, with no rebuild. Where it goes is
  the open question: `bioconda/bioconda-repodata-patches` exists, but its last
  push was 2021-03-15. Find out from bioconda core whether that is still how
  they patch, before writing anything.
- **blast 2.17 for linux-aarch64.** This is the real fix, and it gives
  linux-aarch64 the same bakta as everyone else. It is blocked on bioconda's
  CI, not on the recipe, and nothing here can speed up NCBI's build. Watch
  #58217. Do not open a duplicate.
- **Nothing in `catalogue.py`.** A `pyhmmer <0.12` there would push linux-64
  and osx-arm64 down to bakta 1.11.4 too, because bioconda's bakta 1.12.2
  requires `pyhmmer >=0.12.0` (and `blast >=2.17.0`). A floor of
  `bakta>=1.12` would make the environment fail to build on linux-aarch64
  rather than fail at run time. That is more honest, but it also takes
  AMRFinder with it, which shares the environment and works. Leave the floor
  where it is, so that either route above fixes linux-aarch64 with no change
  here.

## What is left

- Read the other 1.10.x and 1.11.x tags, if the patch is to name builds
  rather than a version range.
- Ask bioconda where repodata patches go now.
- Nothing has been sent.
