#!/bin/bash
#SBATCH --job-name=ipr_container_test
#SBATCH --partition=128x24
#SBATCH --account=128x24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --no-requeue
set -euo pipefail
umask 027
run='/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261003T043208Z_b2a74f5f'
source "$run/scripts/environment.sh"
test -n "${SLURM_JOB_ID:-}"
tmp_root=$(mktemp -d "/tmp/ipr-${USER}-${SLURM_JOB_ID}.XXXXXX")
cleanup() {
    rc=$?
    trap - EXIT
    printf 'Container test exit code: %s\n' "$rc"
    case "$tmp_root" in /tmp/ipr-"${USER}"-"${SLURM_JOB_ID}".*) rm -rf -- "$tmp_root" ;; esac
    exit "$rc"
}
trap cleanup EXIT
export SINGULARITY_TMPDIR="$tmp_root/singularity"
export NXF_TEMP="$tmp_root/nextflow"
export TMPDIR="$tmp_root/tmp"
mkdir -p "$SINGULARITY_TMPDIR" "$NXF_TEMP" "$TMPDIR"
hostname
df -PT "$tmp_root"
stat -f -c 'Local temporary filesystem: %T' "$tmp_root"
count=0
while IFS= read -r container; do
    printf 'TEST %s\n' "$container"
    timeout 60 singularity exec --no-home -B /mnt/beegfs -B /hb "$container" /bin/sh -c 'test -r "$1" && test -r "$2" && test -w "$3"' sh "$run/input/all_significant_proteins.fasta" '/hb/groups/cornejo_training/StrepA_project_revamped/tools/interproscan6/data/interpro/109.0/databases.json' "$run/work"
    printf 'PASS %s\n' "$container"
    count=$((count + 1))
done < "$run/provenance/container_paths.txt"
printf 'ALL %s CONTAINERS PASSED\n' "$count"
