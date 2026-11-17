## Site setup for Case Western Pioneer, sourced by TTpipe.sh and unlock.sh.
## Pioneer's Miniconda3 modules ship conda 23.x, too old for Snakemake 8.27, so use your own
## Miniforge (one-time install, ~0.5 GB); set TTPIPE_CONDA_ROOT if it is elsewhere:
##   wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
##   bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3
source "${TTPIPE_CONDA_ROOT:-$HOME/miniforge3}/etc/profile.d/conda.sh"
conda activate base
