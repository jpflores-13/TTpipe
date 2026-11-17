#! /bin/bash -login

## Release the Snakemake lock after a failed or interrupted run, then
## resubmit with TTpipe.sh.
##
## Usage:
##   bash unlock.sh                                  # default config
##   bash unlock.sh --configfile other.yaml --profile profiles/<site>
## The profile is chosen as in TTpipe.sh: --profile, $TTPIPE_PROFILE,
## then hpcProfile in the config.

## Exit if any command fails
set -e

CONFIG="config/TTpipeConfig.yaml"
PROFILE=""
SNAKEFILE="workflows/TTpipe.snakefile"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --configfile) CONFIG="$2";  shift 2 ;;
        --profile)    PROFILE="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done

cfg() { sed -nE "s/^$1:[[:space:]]*['\"]?([^'\"#[:space:]]+).*/\1/p" "${CONFIG}" | head -1; }
PROFILE="${PROFILE:-${TTPIPE_PROFILE:-$(cfg hpcProfile)}}"

## Load this site's software
[[ -f "${PROFILE%/}/env.sh" ]] || { echo "No ${PROFILE%/}/env.sh: set --profile, \$TTPIPE_PROFILE or hpcProfile in ${CONFIG}" >&2; exit 1; }
source "${PROFILE%/}/env.sh"

## Activate virtual environment
source env/bin/activate

## Unlock snakemake workflow
snakemake --unlock -s "${SNAKEFILE}" --cores 1 --configfile "${CONFIG}"

## Deactivate virtual environment
deactivate

## Success message
echo -e "\033[0;32mDirectory unlocked, ready to rerun."
