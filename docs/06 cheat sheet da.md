# Snydeark til CompareM2

Kun Linux. `cm2` er et alias for `comparem2` overalt. I et git-checkout eller et
pixi-workspace skal hver kommando nedenfor læses som `pixi run comparem2 …`.

## 1. Peg på genomerne

```bash
cd /sti/til/genomer        # relative stier løses herfra, også under pixi run
ls *.fna                   # skallen udfolder mønsteret; det er præcis det, der sendes videre
```

Nævner du slet ingen assemblies, bruger CompareM2 de `*.fna`-, `*.fa`-,
`*.fasta`- og `*.fas`-filer, der ligger i denne mappe, og skriver hvilken
mappe det var.

## 2. Kør

```bash
comparem2                          # de *.fna-, *.fa-, *.fasta- og *.fas-filer, der ligger her
comparem2 *.fna                    # alle 14 værktøjer -> ./results_comparem2/report.html
comparem2 *.fna -o run1            # outputmappe
comparem2 --demo                   # 6 medfølgende plasmider; ingen genomer, ingen databaser
comparem2 *.fna --dry-run          # hvad der ville blive kørt
comparem2 *.fna --tui              # tastaturflade; vælg med mellemrum/a, kør med r
```

Før der hentes noget som helst, skriver programmet, hvad det kommer til at
koste:

```
4 assemblies, 14 tools
to download: checkm2, gtdb, bakta-light, amrfinder (62.5 GB + 2 of unknown size) -> ~/.comparem2/databases
tool environments: 6 in ~/.comparem2/envs (none built yet)
```

Den første kørsel er stille i omkring et minut, mens de seks conda-miljøer
bliver løst. Det er forventet.

## 3. Kør et udvalg: den store besparelse

`--until` navngiver værktøjer og trækker selv deres forudsætninger med.

```bash
comparem2 *.fna --until fasttree               # kører bakta, panaroo, fasttree
comparem2 *.fna --until seqkit mashtree treecluster skani   # hurtigt, slet ingen databaser

# alt undtagen GTDB-downloaden på 60,8 GB (~3 GB i alt i stedet for 62,5 GB)
comparem2 *.fna --until seqkit checkm2 bakta amrfinder mlst mashtree treecluster \
                       skani panaroo snp-dists fasttree carveme biosynthesis
```

Værktøjsnavne: `seqkit checkm2 gtdbtk bakta amrfinder mlst mashtree treecluster
skani panaroo snp-dists fasttree carveme biosynthesis`.

Afhængigheder: `amrfinder`, `panaroo`, `carveme` ← `bakta`; `snp-dists`,
`fasttree` ← `panaroo`; `treecluster` ← `mashtree`; `biosynthesis` ← `carveme`.

Afhængigheder vælges automatisk.

## 4. Output

```
results_comparem2/
├── report.html                    selvstændig fil; overlever at blive mailet hjem fra en cluster
├── <værktøj>/…                    resultater for hele sættet
├── samples/<navn>/<værktøj>/…     resultater per genom
├── logs/… og samples/<navn>/logs/…
└── .comparem2/Snakefile           den genererede workflow, hvis du skal fejlsøge
```

Prøvenavne dannes af filnavnets stamme, hvor alt uden for `[A-Za-z0-9._-]`
erstattes af `_`, og omdøbningen bliver skrevet ud.

```bash
comparem2 *.fna --report-only      # gentegn rapporten ud fra det, der allerede ligger på disken
comparem2 *.fna --keep-going       # stop ikke de øvrige værktøjer, fordi ét fejler
```

## 5. `--set`: videregiv argumenter til et værktøj

`--set <værktøj><flag>=<værdi>`, kan gentages. Stav flaget præcis som værktøjet
selv gør, bindestreger og det hele: TreeClusters lange flag har to, skanis `-c`
har én. At navngive ét flag erstatter kun det flag; værktøjets øvrige
standardværdier bliver stående.

```bash
# de to, der er værd at kende
comparem2 *.fna --set skani-c=125                  # komplette isolater; 70 passer til fragmenterede MAG'er
comparem2 *.fna --set treecluster--threshold=0.1   # standard 0.05; tærsklen afgør klyngedannelsen

# flere på én gang
comparem2 *.fna \
  --set treecluster--method=avg_clade \
  --set treecluster--threshold=0.1 \
  --set skani-c=125 \
  --set mashtree--genomesize=3000000 \
  --set bakta--gram=+

# et flag uden værdi sendes bart videre
comparem2 *.fna --set bakta--force=
```

Standardværdierne, du overskriver: `bakta --force`; `mashtree --genomesize
5000000 --mindepth 5 --kmerlength 21 --sketch-size 10000`; `treecluster --method
max_clade --threshold 0.05`; `skani -c 70`. Hvert værktøjs fulde kommandolinje
står under [hvilke analyser laver den](30 what analyses does it do.md),
genereret ud fra specifikationerne, så den ikke kan komme ud af trit med det,
der faktisk køres.

---

*Oversat fra `docs/06 cheat sheet.md`, CompareM2 v3.4.0, 22. september 2026.
Kommandoer, flag, stier og programmets egne udskrifter står uoversat, fordi det
er dem, du skal skrive og læse på skærmen. Som PDF til udlevering:
[`CHEATSHEET.da.pdf`](assets/CHEATSHEET.da.pdf), to A4-sider, gentegnet med
opskriften øverst i `docs/assets/cheatsheet-print.css`. Dokumentationen på
[comparem2.readthedocs.io](https://comparem2.readthedocs.io) er kun på engelsk.*
