#!/bin/bash
set +euo pipefail
# This script downloads the newest version from github and tests relevant output files for correct hashes.



echo ">>> Test started <<<"



tmpdir=$(echo -n _; echo ${RANDOM} | md5sum | head -c 8) # the underscore in front makes it easy to glob all test directories.
mkdir -p $tmpdir

echo ">>> Cloning <<<"
git clone https://github.com/cmkobel/assemblycomparator2.git $tmpdir

echo ">>> Entering tmpdir $tmpdir <<<"
cd $tmpdir


ASSCOM2_BASE=$(realpath .)

echo ">>> ASSCOM2_BASE is $ASSCOM2_BASE <<<"


echo ">>> Creating starter environment with prefix starter${tmpdir} <<<"

cd $ASSCOM2_BASE && mamba env create -f environment.yaml -p ${ASSCOM2_BASE}/starter${tmpdir}

echo ">>> These are the starter environment software versions <<<"
conda run \
    --live-stream \
    --prefix ${ASSCOM2_BASE}/starter${tmpdir} \
        python --version; mamba --version; conda --version



echo ">>> Entering test genome directory <<<"
cd tests/E._faecium
pwd

echo ">>> Running <<<"

set -euo
conda run \
    --live-stream \
    --prefix ${ASSCOM2_BASE}/starter${tmpdir} \
        snakemake --snakefile ${ASSCOM2_BASE}/snakefile \
        --profile ${ASSCOM2_BASE}/profiles/local/ \
        --configfile ${ASSCOM2_BASE}/config.yaml
set +euo

echo ">>> Run has ended with exitcode $exitcode <<<"



echo ">>> Done. Please consider deleting the test ASSCOM2_BASE ${ASSCOM2_BASE} <<<"