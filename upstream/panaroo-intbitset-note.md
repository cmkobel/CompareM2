# Draft: what replacing intbitset would cost panaroo

**Unsent.** Destination: an issue or discussion on
[gtonkinhill/panaroo](https://github.com/gtonkinhill/panaroo). Measured
2026-09-03 on macOS 26.2 / Apple silicon, panaroo 1.8.0.

Written because the question came up here and the answer is not obvious in
either direction. It is **not** a request to make the change — the benchmark
below argues against the easy version of it.

## Why the dependency is worth a second look

`intbitset` is the only compiled dependency panaroo carries, and its
conda-forge feedstock is stalled: last human merge 2023-09-19, seven open bot
PRs, and builds only for python 3.8–3.11. Two consequences fall on panaroo
rather than on intbitset:

- `python 3.12` + panaroo is unsatisfiable on `linux-64`. intbitset is what
  caps panaroo's python on conda, and it is why this project's environment sits
  at python 3.11.16.
- The feedstock never opted into `osx_arm64`, so panaroo cannot be installed on
  Apple silicon at all. (Fixable in the feedstock — two lines in
  `conda-forge.yml` plus a version bump — see `../STATUS.md`.)

Both are packaging problems with packaging fixes. Dropping the dependency would
also solve them, which is why it is worth knowing what it would cost.

## The surface is small

Panaroo 1.8.0 uses eight operations: construction from an iterable of ints (13
call sites, no other form), `add`, `discard`, `remove`, `copy`, `intersection`,
`|=`, `&`, plus `len`, `in` and iteration. Every one has identical `set`
semantics. The C extension iterates in sorted order and panaroo serialises
members by joining them in iteration order, so a replacement has to sort.

```python
class intbitset(set):
    def __init__(self, iterable=(), *a, **k): super().__init__(iterable)
    def __iter__(self): return iter(sorted(set.__iter__(self)))
    def copy(self): return intbitset(set.__iter__(self))
```

Shadowing the C extension with that (`PYTHONPATH` ahead of site-packages) and
running `panaroo --clean-mode strict -a core -t 4 --remove-invalid-genes` over
three *E. faecium* assemblies: `gene_presence_absence.Rtab`,
`gene_presence_absence.csv`, `gene_presence_absence_roary.csv`,
`summary_statistics.txt`, `pan_genome_reference.fa` and
`struct_presence_absence.Rtab` are **byte-identical** to the C-extension run.

`core_gene_alignment.aln` and `final_graph.gml` differ in record order and
`seqIDs` order, with identical sequences — but **two runs of the unmodified
C-extension panaroo differ the same way**, so that is existing run-to-run
nondeterminism, not the replacement. Wall clock 3:16.96 against 3:25.33; at
three genomes the data structure is irrelevant.

## And that is exactly why the easy version fails

Synthetic benchmark mirroring panaroo's usage — 10,000 pangenome nodes (30 %
core, present in every genome; the rest in up to 10), 200,000 `members |=
intbitset([mem])` as `merge_nodes` does, 200,000 `len(a & b)` as `cdhit` does,
then one iteration pass for serialisation. Peak RSS is for holding the 10,000
membership objects.

| genomes | intbitset | `set` | int bitmask | intersect: intbitset → `set` → int |
| ------- | --------- | ----- | ----------- | ---------------------------------- |
| 100 | 1.3 MB | 29.2 MB | 0.5 MB | 0.10 s → 0.27 s → 0.08 s |
| 1,000 | 2.3 MB | 170.0 MB | 2.1 MB | 0.12 s → 0.66 s → 0.09 s |
| 5,000 | 7.4 MB | **1,863.6 MB** | 13.5 MB | 0.14 s → 2.02 s → 0.14 s |

At 5,000 genomes the `set` version needs **250× the memory** for node
membership alone, before edges — which are more numerous. Panaroo is built for
that scale, so `set` is not a serious proposal.

A membership bitmask in a python `int` is a different matter: memory stays
within 2× of the C extension (13.5 MB against 7.4 MB at 5,000 genomes) and
union and intersection match it (0.15 s, 0.14 s), because CPython's big-int
`|` and `&` are the same word-at-a-time loop intbitset runs. The cost lands on
iteration, which the naive decode does bit by bit: 4.74 s against 0.33 s at
5,000 genomes, a 14× regression on every serialisation path. A byte-lookup
decode should close most of that, unmeasured.

So the honest summary is: the dependency is replaceable in principle, `set` is
not the way, and an int-bitmask replacement is a real piece of work — writing
and maintaining a bitset, touching every site that iterates or serialises
members — rather than deleting an import.

## What this draft still needs before it is sent

- A byte-lookup iteration decode, measured. Without it the int-bitmask column
  is an argument with a hole in it.
- A panaroo run at a genome count where any of this matters. Everything above
  at three genomes says only that the replacement is *correct*, not that it is
  affordable; the table is synthetic.
- The recipe half stated plainly: bioconda's panaroo lists `intbitset` in
  `requirements: run:`, so patching the source changes nothing for
  `conda install panaroo` until that line goes too.

## Rebuilding the numbers

The benchmark, to run against each variant in an environment that has panaroo:

```python
# bench.py — usage: python bench.py {c,set,int} <n_genomes>
import random, resource, sys, time

variant, n_genomes = sys.argv[1], int(sys.argv[2])
N_NODES, N_OPS = 10_000, 200_000
random.seed(0)

if variant == "c":
    from intbitset import intbitset as make
elif variant == "set":
    make = set
elif variant == "int":
    make = lambda xs: sum(1 << x for x in set(xs))

rss0 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
nodes = []
for i in range(N_NODES):
    members = range(n_genomes) if i % 10 < 3 else random.sample(
        range(n_genomes), min(10, n_genomes))
    nodes.append(make(members))
rss_built = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

t = time.perf_counter()                                  # merge_nodes.py
for _ in range(N_OPS):
    i, mem = random.randrange(N_NODES), random.randrange(n_genomes)
    nodes[i] |= (1 << mem) if variant == "int" else make([mem])
t_union = time.perf_counter() - t

t = time.perf_counter()                                  # cdhit.py
for _ in range(N_OPS):
    a, b = nodes[random.randrange(N_NODES)], nodes[random.randrange(N_NODES)]
    n = (a & b).bit_count() if variant == "int" else len(a & b)
t_inter = time.perf_counter() - t

t = time.perf_counter()                                  # __main__.py join
for n in nodes:
    if variant == "int":
        sum(1 for i in range(n.bit_length()) if n >> i & 1)
    else:
        sum(1 for _ in n)
t_iter = time.perf_counter() - t

print(f"{variant} G={n_genomes} union {t_union:.2f}s intersect {t_inter:.2f}s "
      f"iterate {t_iter:.2f}s RSS {(rss_built - rss0) / 2**20:.1f} MB")
```

The equivalence check is the `class intbitset(set)` above, saved as
`intbitset.py` in a directory put on `PYTHONPATH` ahead of site-packages, then
the same panaroo command run twice — once with it, once without — and the
outputs compared. Compare the C extension against *itself* too, or the ordering
differences will look like a regression.
