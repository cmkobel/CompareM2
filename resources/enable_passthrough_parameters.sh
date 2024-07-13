#!/bin/bash


# In order to enable the passthrough parameters, hyphens must be allowed. With this oneliner, the regex pattern is modified from [a-zA-Z_]\w*$ to [a-zA-Z_][\w-]*\w$ in snakemake's parse_config().

sed -i 's|valid = re.compile(r"\[a-zA-Z_\]\\w\*\$")|valid=re.compile(r"\[a-zA-Z_\]\[\\w-\]\*\\w\$")|g' $CONDA_PREFIX/lib/python3.11/site-packages/snakemake/__init__.py

