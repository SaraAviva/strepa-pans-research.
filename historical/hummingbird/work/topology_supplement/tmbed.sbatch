#!/bin/bash
#SBATCH --job-name=topology_tmbed_66
#SBATCH --partition=128x24
#SBATCH --account=128x24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --chdir=/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261006T234134Z_topology_c655cc0c/launch
#SBATCH --output=/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261006T234134Z_topology_c655cc0c/logs/tmbed-%j.out
#SBATCH --error=/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261006T234134Z_topology_c655cc0c/logs/tmbed-%j.err
#SBATCH --mail-user=researcher@example.org
#SBATCH --mail-type=END,FAIL
#SBATCH --no-requeue
set -euo pipefail
umask 027
set +u
module load miniconda3 singularity-ce/singularity-ce.4.1.4
set -u
task_tmp=$(mktemp -d "/tmp/topology-${USER}-${SLURM_JOB_ID}.XXXXXX")
export SINGULARITY_TMPDIR="$task_tmp/singularity"
export SINGULARITY_CACHEDIR="$task_tmp/cache"
export TMPDIR="$task_tmp/tmp"
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
mkdir -p "$SINGULARITY_TMPDIR" "$SINGULARITY_CACHEDIR" "$TMPDIR"
cleanup() {
 rc=$?
 trap - EXIT
 printf 'Component finished with exit %s at ' "$rc"
 date -u
 case "$task_tmp" in /tmp/topology-"${USER}"-"${SLURM_JOB_ID}".*) rm -rf -- "$task_tmp" ;; esac
 exit "$rc"
}
trap cleanup EXIT
hostname
python3 /hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261006T234134Z_topology_c655cc0c/scripts/run_components.py /hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261006T234134Z_topology_c655cc0c tmbed
