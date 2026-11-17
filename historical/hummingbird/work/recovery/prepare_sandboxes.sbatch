#!/bin/bash
#SBATCH --job-name=ipr_prepare_containers
#SBATCH --partition=128x24
#SBATCH --account=128x24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --no-requeue
set -euo pipefail
umask 027
run='/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261003T043208Z_b2a74f5f'
recovery="$run/recovery_20261004T163308Z"
source "$run/scripts/environment.sh"
tmp_root=$(mktemp -d "/tmp/ipr-${USER}-${SLURM_JOB_ID}.XXXXXX")
cleanup() {
    rc=$?
    trap - EXIT
    printf 'Extraction test exit code: %s\n' "$rc"
    case "$tmp_root" in /tmp/ipr-"${USER}"-"${SLURM_JOB_ID}".*) rm -rf -- "$tmp_root" ;; esac
    exit "$rc"
}
trap cleanup EXIT
export SINGULARITY_TMPDIR="$tmp_root/singularity"
export NXF_TEMP="$tmp_root/nextflow"
export TMPDIR="$tmp_root/tmp"
mkdir -p "$SINGULARITY_TMPDIR" "$NXF_TEMP" "$TMPDIR"
hostname
python3 "$recovery/prepare_sandboxes.py"
