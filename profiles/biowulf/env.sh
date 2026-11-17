## Site setup for NIH Biowulf, sourced by TTpipe.sh and unlock.sh.
## Needs python3 (>= 3.11, for Snakemake 8) and conda (>= 24.7.1).
## Biowulf's `conda` module may be older than Snakemake 8.27 needs; if the first
## run complains, install your own Miniforge and set TTPIPE_CONDA_ROOT:
##   wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
##   bash Miniforge3-Linux-x86_64.sh -b -p /data/$USER/miniforge3
if [ -f "${TTPIPE_CONDA_ROOT:-/nonexistent}/etc/profile.d/conda.sh" ]; then
    source "${TTPIPE_CONDA_ROOT}/etc/profile.d/conda.sh"
else
    module load conda
    source "$(conda info --base)/etc/profile.d/conda.sh"
fi
conda activate base
