#!/usr/bin/env bash
# Builds the trimmed GTDB-Tk reference package the weekly `gtdbtk` job runs
# against, because the real one is 60.8 GB to download and 92 GB extracted, and
# a free runner has about 14 GB of disk.
#
# What it keeps, and why that is enough: the four E. faecium test genomes are
# all classified by GTDB-Tk's ANI screen, and when every genome is, classify_wf
# skips identify and align ("Identify and Align steps will be skipped.") and
# never reaches pplacer. So `msa/`, `pplacer/` and `split/` (14 GB) are left
# out, and the 80 GB skani database is replaced by one sketched from five
# GTDB r232 species representatives. Everything else is copied from the real
# package unchanged, ~157 MB, mostly the marker HMMs.
#
# What that means for the job: it covers installing GTDB-Tk, the reference
# version check, the skani screen, the summaries, and the pipeline's own rule
# around them. It does not cover prodigal, hmmer or pplacer. A genome that
# skani cannot place would make GTDB-Tk look for the files left out.
#
# Built on thylakoid on 2026-10-27, from the r232 package fetched 2026-09-02,
# with the skani in the pipeline's gtdbtk environment (0.3.2; GTDB built r232's
# with 0.3.1). Against it, GTDB-Tk 2.7.2 gave a summary identical in all 20
# columns to the real package's, for all four genomes.
#
#   tests/environments/gtdb_mock.sh <real r232 package dir> <skani> <out.tar.gz>

set -euo pipefail

real=$1 skani=$2 out=$3

# The species representatives, as GTDB r232 names them. E. faecium is the one
# the test genomes match; E. lactis is what the real package reports as their
# nearest other species (94.35% ANI), so `other_related_references` comes out
# the same. The other three give the screen something to reject.
REPRESENTATIVES=(
    GCF_029023785.1   # s__Enterococcus_B faecium
    GCF_015751045.1   # s__Enterococcus_B lactis
    GCF_029024925.1   # s__Enterococcus faecalis
    GCF_003697165.2   # s__Escherichia coli
    GCF_001027105.1   # s__Staphylococcus aureus
)

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
pkg="$work/release232"   # the real tarball's wrapper, so the same
mkdir -p "$pkg/skani" "$work/refs"   # `--strip-components=1` unpacks both

cp -r "$real"/{markers,masks,metadata,mrca_red,radii,taxonomy} "$pkg/"

# GTDB-Tk names a reference by its file name less `_genomic.fna.gz`, so the
# files keep NCBI's names. The assembly directory carries a suffix the
# accession does not, so it is read from the listing.
for acc in "${REPRESENTATIVES[@]}"; do
    n=${acc#GCF_}; n=${n%.*}
    dir="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${n:0:3}/${n:3:3}/${n:6:3}/"
    asm=$(curl -fsSL "$dir" | grep -o "${acc}_[^/\"]*" | head -1)
    curl -fsSL -o "$work/refs/${acc}_genomic.fna.gz" "$dir$asm/${asm}_genomic.fna.gz"
done

"$skani" sketch "$work"/refs/*_genomic.fna.gz -o "$pkg/skani/database"

cat > "$pkg/README" <<EOF
A trimmed GTDB-Tk r232 reference package, for CompareM2's continuous
integration only. It cannot classify anything but the five species sketched
under skani/database, and has no msa/, pplacer/ or split/: do not use it for
real work. Built by tests/environments/gtdb_mock.sh in
https://github.com/cmkobel/CompareM2 on $(date -u +%Y-%m-%d), with
$("$skani" --version).

Derived from the GTDB-Tk r232 reference package, https://gtdb.ecogenomic.org
(Parks et al.), and like it licensed under the Creative Commons
Attribution-ShareAlike 4.0 International License,
https://creativecommons.org/licenses/by-sa/4.0/. Changes: msa/, pplacer/ and
split/ removed, and skani/database rebuilt from five NCBI RefSeq assemblies,
GTDB's representatives for the species named in gtdb_mock.sh.
EOF

tar -czf "$out" -C "$work" release232
md5sum "$out"
