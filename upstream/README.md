# Upstream reports

Drafts of things this project found in software it depends on, kept here so they
do not evaporate with the session that wrote them. **Nothing here has been
sent.** Each file says what it is, where it goes, and what it still needs.

| draft | destination | status |
| ----- | ----------- | ------ |
| [carveme-205-comment.md](carveme-205-comment.md) | comment on [cdanielmachado/carveme#205](https://github.com/cdanielmachado/carveme/issues/205) | unsent |
| [scip-question.md](scip-question.md) | issue or discussion on [scipopt/scip](https://github.com/scipopt/scip) | unsent |

The finding behind both is in [../STATUS.md](../STATUS.md) (*CarveMe was nine
minutes for the wrong reason*), [../DECISIONS.md](../DECISIONS.md) and
`../src/comparem2/carve_scip.py`.

Deliberately **not** filed: conda-forge's `scip` feedstock and bioconda's
`carveme` recipe. PaPILO is a legitimate part of the feedstock's build, there is
no PaPILO-free conda SCIP to pin to, and the fix belongs in CarveMe.

## The reproducer

Both drafts offer three attachments: the written-out MILP and two solutions.
They are generated files, not kept in git. To rebuild them from a protein FASTA
in an environment that has carveme:

```python
# write the MILP exactly as `carve` would build it, then stop
import pandas as pd
from carveme import config, project_dir
from carveme.reconstruction.carving import minmax_reduction
from carveme.reconstruction.diamond import load_diamond_results
from carveme.reconstruction.scoring import reaction_scoring
from reframed import load_cbmodel
from reframed.core.transformation import apply_bounds
from reframed.solvers.scip_solver import SCIPSolver

universe = load_cbmodel(project_dir + config.get('generated', 'default_universe'),
                        flavor='bigg')
universe.id = "probe"
apply_bounds(universe)
gprs = pd.read_csv(project_dir + config.get('generated', 'bigg_gprs'))
gprs = gprs[gprs.reaction.isin(list(universe.reactions))]
scores, _ = reaction_scoring(load_diamond_results("hits.tsv"), gprs)

def solve(self, *a, **k):
    self.problem.writeProblem("carveme.lp")
    raise SystemExit(0)

SCIPSolver.solve = solve
minmax_reduction(universe, dict(scores[['reaction', 'normalized_score']].values),
                 default_score=-1.0, uptake_score=0.0, soft_score=1.0)
```

`hits.tsv` is what `carve` leaves beside its input: `diamond blastp -d
<carveme>/data/generated/bigg_proteins.dmnd -q proteins.faa -o hits.tsv
--more-sensitive --top 10 --quiet`.

The two solutions come from the LP itself, no CarveMe needed:

```bash
scip -c "read carveme.lp set presolving milp maxrounds 0 set limits gap 0.001 \
         optimize write solution carveme_947.sol quit"     # the 947.4997 point
scip -c "read carveme.lp set limits gap 0.001 optimize \
         write solution carveme_943.sol quit"              # with the wheel's SCIP
```

The copies used for the numbers in these drafts were produced on thylakoid from
`tests/E._faecium/116_2.fna` on 2026-09-02 and staged, gzipped, in
`~/Downloads/carveme-scip-repro/` on the laptop — which is not a durable
location, hence the recipe above.
