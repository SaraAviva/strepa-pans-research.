#!/bin/bash
#SBATCH --job-name=ipr_extract_test
#SBATCH --partition=128x24
#SBATCH --account=128x24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
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
image=$(head -1 "$run/provenance/container_paths.txt")
sandbox="$recovery/sandboxes/$(basename "$image" .img)"
test ! -e "$sandbox"
singularity sif list "$image"
offset=$(singularity sif list "$image" | awk -F '|' '/FS \(Squashfs\/\*System\/amd64\)/ {gsub(/ /,"",$4); split($4,a,"-"); print a[1]}')
[[ "$offset" =~ ^[0-9]+$ ]]
LD_LIBRARY_PATH="$recovery/extraction_tools/lib" "$recovery/extraction_tools/unsquashfs" -no-progress -no-xattrs -processors 1 -data-queue 32 -frag-queue 32 -offset "$offset" -dest "$sandbox" "$image"
timeout 90 singularity exec --no-home -B /mnt/beegfs -B /hb "$sandbox" /bin/sh -c 'test -r "$1" && test -r "$2" && test -w "$3" && hmmsearch -h | head -3' sh "$run/input/all_significant_proteins.fasta" '/hb/groups/cornejo_training/StrepA_project_revamped/tools/interproscan6/data/interpro/109.0/databases.json' "$run/work"
printf 'EXTRACTED CONTAINER EXECUTION PASSED\n'
sha256sum "$image"
