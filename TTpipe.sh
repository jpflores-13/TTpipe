#! /bin/bash -login
#SBATCH -J TTpipe
#SBATCH -t 10-00:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=16gb    # building the conda envs on a first run needs > 4 GB
#SBATCH -o "%x-%j.out"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80

## TTpipe launcher (MicroCPipe / ChromSimPipe pattern).
##
## Usage:
##   sbatch TTpipe.sh                         # SLURM cluster
##   bash   TTpipe.sh                         # profiles/local (no scheduler)
##   bash   TTpipe.sh -n                      # dry run: list jobs only
##   sbatch TTpipe.sh --configfile other.yaml # use a different config
##   sbatch TTpipe.sh --profile profiles/<site>   # override hpcProfile
##
## The compute site is a profile folder: its env.sh loads software, its
## config.yaml tells Snakemake how to run jobs. It is chosen by, in order:
## --profile, $TTPIPE_PROFILE, then hpcProfile in the config.
## Per-cluster settings that stay out of git (e.g. in ~/.bashrc):
##   export TTPIPE_PROFILE=profiles/biowulf        # this cluster's profile
##   export TTPIPE_OUTPUT_ROOT=/data/$USER/ttpipe  # where a relative outputDir goes
## Email: `export TTPIPE_EMAIL=you@university.edu` and
##   - the workflow mails you when it starts, finishes or fails (`notify:` in the config);
##   - SLURM also mails about the launcher job itself if you submit with
##     `sbatch --mail-user=$TTPIPE_EMAIL TTpipe.sh` (the --mail-type below is set here).
## Any other arguments are passed straight to snakemake.

## Exit if any command fails
set -e

## Config file and profile (overridable with --configfile / --profile)
CONFIG="config/TTpipeConfig.yaml"
PROFILE=""
ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --configfile) CONFIG="$2";  shift 2 ;;
        --profile)    PROFILE="$2"; shift 2 ;;
        *)            ARGS+=("$1"); shift ;;
    esac
done

## Read a top-level scalar from the config, e.g. cfg outputDir
cfg() { sed -nE "s/^$1:[[:space:]]*['\"]?([^'\"#[:space:]]+).*/\1/p" "${CONFIG}" | head -1; }
PROFILE="${PROFILE:-${TTPIPE_PROFILE:-$(cfg hpcProfile)}}"
PROFILE="${PROFILE%/}"

## Output directory, resolved the same way as workflows/rules/common.smk
OUTDIR="$(cfg outputDir)"
OUTDIR="${OUTDIR:-output}"
OUTDIR="${OUTDIR/#\~/$HOME}"
[[ "${OUTDIR}" = /* ]] || OUTDIR="${TTPIPE_OUTPUT_ROOT:-$PWD}/${OUTDIR}"

## Load this site's software (python3, conda)
[[ -f "${PROFILE}/env.sh" ]] || { echo "No ${PROFILE}/env.sh: set --profile, \$TTPIPE_PROFILE or hpcProfile in ${CONFIG}" >&2; exit 1; }
source "${PROFILE}/env.sh"

## Create and activate virtual environment with requirements
[ -d env ] || python3 -m venv env
source env/bin/activate
python3 -m pip install -q -r config/requirements.txt

## SLURM job logs (only for profiles that use the SLURM executor)
EXTRA=()
if grep -qE "^executor:[[:space:]]*slurm" "${PROFILE}/config.yaml"; then
    SLURM_LOGS="${OUTDIR}/logs_slurm"
    mkdir -p "${SLURM_LOGS}"
    EXTRA+=(--slurm-logdir "${SLURM_LOGS}" --max-jobs-per-second 5 --max-status-checks-per-second 0.5)
fi

## Execute TTpipe snakemake workflow
snakemake -s workflows/TTpipe.snakefile --profile "${PROFILE}" --configfile "${CONFIG}" --rerun-incomplete --keep-going -p "${EXTRA[@]}" "${ARGS[@]}"

## Success message
echo "Workflow completed successfully!"
