## Site setup template, sourced by TTpipe.sh and unlock.sh.
## Make python3 (>= 3.11, for Snakemake 8) and conda (>= 24.7.1) available here, e.g. `module load ...`
## or your own Miniforge (set TTPIPE_CONDA_ROOT).
source "${TTPIPE_CONDA_ROOT:-$HOME/miniforge3}/etc/profile.d/conda.sh"
conda activate base
