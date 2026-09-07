# Known-good environment locks

**These are evidence, not yet part of the build.** Nothing in the pipeline
reads them. They are here because on 2026-09-07 the floors-only specs in
`catalogue.py` stopped solving, and this is the only surviving record of an
environment that actually ran all thirteen tools.

Exported with `conda list --explicit` from the environments on thylakoid that
produced the verified runs in `STATUS.md`:

| File | Packages | Source prefix on thylakoid |
| --- | ---: | --- |
| `main.linux-64.lock` | 394 | `/evo/postdoc/cm2-envs-two/f35bbb1ff167437785dcb4a2729c2beb_` |
| `checkm2.linux-64.lock` | 130 | `/evo/postdoc/cm2-envs-two/d2bf6f91528e69362e223a294d337592_` |

`@EXPLICIT` format: full package URLs with hashes, no solving. `conda create
--file <lock>` reproduces the environment exactly, on `linux-64` only — which
is the only platform the tools exist for anyway.

Measured 2026-09-07 on GenomeDK, where the floors-only spec fails after
4 min 07 s: `conda create --dry-run --file main.linux-64.lock` returns **exit 0
in 6.9 s**.

Regenerate against a prefix that is known to work:

```bash
conda list --explicit -p <prefix> > envs/locks/main.linux-64.lock
```

See `STATUS.md`, *Known broken: the main environment no longer solves*, for what
went wrong and what the options are.
