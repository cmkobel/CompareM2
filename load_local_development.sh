# This file is used for developing new features into asscom2.
# Use this file by executing this script in the current shell.
# $ . load_local_development.sh 
# Then you can change the snakefile, the conda-yamls, the report and finally update the Dockerfile before committing the new feature.

# Add cwd to the path to we can run this asscom2 from any dir.
export PATH=$(realpath .):$PATH

# Set cwd as the asscom2 base directory.
export ASSCOM2_BASE=$(realpath .)

# You should use the conda profile so you can change the upstream dockerfile. Local development is always better because you get a faster development cycle.
export ASSCOM2_PROFILE=$(realpath profiles/conda/local/)

