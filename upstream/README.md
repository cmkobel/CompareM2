# Upstream reports

Drafts of things this project found in software it depends on, kept here so they
do not evaporate with the session that wrote them. Each file says what it is,
where it goes, and what it still needs. **One has been sent** — the status
column is the record of which.

| draft | destination | status |
| ----- | ----------- | ------ |
| [carveme-205-comment.md](carveme-205-comment.md) | comment on [cdanielmachado/carveme#205](https://github.com/cdanielmachado/carveme/issues/205) | unsent |
| [scip-question.md](scip-question.md) | issue or discussion on [scipopt/scip](https://github.com/scipopt/scip) | unsent |
| [panaroo-intbitset-note.md](panaroo-intbitset-note.md) | issue or discussion on [gtonkinhill/panaroo](https://github.com/gtonkinhill/panaroo) | unsent, and incomplete on purpose — see its last section |
| [intbitset-feedstock-pr.md](intbitset-feedstock-pr.md) | PR against [conda-forge/intbitset-feedstock](https://github.com/conda-forge/intbitset-feedstock) | **sent 2026-09-04** — [PR #21](https://github.com/conda-forge/intbitset-feedstock/pull/21) |
| [bioconda-panaroo-pr.md](bioconda-panaroo-pr.md) | PR against [bioconda/bioconda-recipes](https://github.com/bioconda/bioconda-recipes) `recipes/panaroo` | unsent; patch prepared and verified to apply |

The first two findings are in [../STATUS.md](../STATUS.md) (*CarveMe was nine
minutes for the wrong reason*), [../DECISIONS.md](../DECISIONS.md) and
`../src/comparem2/carve_scip.py`.

## The two panaroo-on-macOS PRs are a pair

`intbitset-feedstock-pr.md` and `bioconda-panaroo-pr.md` are what it takes to
make panaroo installable on Apple silicon, and **neither works alone**. Verified
together on 2026-09-03: with both applied, `panaroo>=1.5`, `snp-dists>=1.2.0`
and `fasttree>=2.2.0` solve on `osx-arm64` from conda alone. Sequence them in
either order, but do not send one and stop.

`panaroo-intbitset-note.md` is a different question — whether panaroo should
carry the dependency at all — and its answer was *probably yes, leave it*. It is
not a prerequisite for either PR.

Working material for both, including the prepared feedstock branch, the applied
patch and the local arm builds, is in `~/postdoc/cm2-macos/` — outside git,
because it holds a conda-bld tree.

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
