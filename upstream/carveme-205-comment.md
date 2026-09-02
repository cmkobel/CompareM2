# Draft: comment on cdanielmachado/carveme#205

**Destination:** a comment on [#205 "Carveme not finishing with SCIP solver"](https://github.com/cdanielmachado/carveme/issues/205)
(open, unanswered since 2024-06-21). Not a new issue — that one is this bug.
[#220](https://github.com/cdanielmachado/carveme/issues/220) has a 2026-04-23
comment ("1.6.6 takes forever", singularity) that is probably also this; worth
linking, not worth claiming, since #220's own report is about CPLEX.

**Unsent.** Before sending: decide whether to open the DIAMOND-output-path
paragraph as its own issue instead, and offer the reproducer files (see
[README.md](README.md)).

---

Diagnosed on 1.6.6 — this reproduces for me, and the cause is the SCIP *build*
rather than CarveMe's use of it.

conda-forge's `scip` links PaPILO; the PyPI `pyscipopt` wheel's SCIP does not
(`Model().printExternalCodeVersions()` shows the difference). Since bioconda's
`carveme` depends on `scip`, every conda, bioconda and container install has
PaPILO — and on the carving MILP it runs to the hardcoded 600 s limit in
`reconstruction/carving.py` and returns whatever it has. One `.faa` of 2588
proteins, 24 idle cores, identical DIAMOND hits and universe model in both rows:

| | MILP | reported | model | annotated reactions dropped |
| --- | ---: | --- | --- | ---: |
| as installed | 601 s | time limit | 1193 rx / 750 met | 253 of 1069 |
| `presolving/milp/maxrounds=0` | 43 s | gap limit | 1742 rx / 1175 met | 44 |

A second genome (3185 proteins) repeats it: 908 s and 261 annotated reactions
dropped, against 17 s and 45. So this is not only slow — users are getting
models a quarter smaller in gene-supported reactions than the run intends, with
nothing in the output to say so.

Three lines in `reconstruction/carving.py` fix it, beside the limits already set
there:

```python
if solver.__class__.__name__ == 'SCIPSolver':
    solver.problem.setParam('limits/time', 600)
    solver.problem.setParam('limits/gap', 0.001)
    solver.problem.setParam('presolving/milp/maxrounds', 0)   # PaPILO; see below
```

`reconstruction/gapfilling.py` sets the same two limits and would want the same
line — and since each gap-fill medium is its own MILP with its own 600 s, a
`--gapfill A,B,C` run multiplies the cost, which may be part of what you are
seeing here.

Until then, a wrapper works. Patch at `solve`, because `minmax_reduction` sets
its parameters immediately before calling it, so anything set earlier is
overwritten:

```python
from reframed.solvers.scip_solver import SCIPSolver
_solve = SCIPSolver.solve
def solve(self, *a, **k):
    self.problem.setParam("presolving/milp/maxrounds", 0)
    return _solve(self, *a, **k)
SCIPSolver.solve = solve

import sys; sys.argv[0] = "carve"
from carveme.cli.carve import main
raise SystemExit(main())
```

One thing worth flagging separately, because it is the deeper issue. On this
problem SCIP reports three different objective values depending on
configuration — 913.5 at the time limit as installed, 943.5 "optimal" when
handed a starting solution, 947.5 "optimal" with the presolver off — and all of
those points are feasible only to tolerance: max constraint violation ~5e-11,
but integrality violation exactly 1e-6 against the default `feastol` of 1e-6,
and rounding their binaries to exact integers leaves an infeasible LP.
`minmax_reduction` couples each flux to its indicator with `bigM=1e3` against
`eps=1e-3`, six orders of magnitude apart, so which model a user gets depends on
their solver build. Scaling those apart, or tightening `numerics/feastol`, may
matter more than the presolver does.

Also, unrelated and small: `run_blast` writes to `os.path.splitext(inputfile)[0]
+ '.tsv'`, so `carve path/to/proteins.faa` overwrites `path/to/proteins.tsv` if
one exists. In a pipeline that annotates with Bakta and carves the resulting
`.faa` in place, that silently replaces Bakta's feature table — which is exactly
what happened to us, and went unnoticed because nothing downstream reads it.
Happy to open it as its own issue.

I can send the written-out MILP and the two solution files if they are useful.
