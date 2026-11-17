#!/bin/bash
#SBATCH --job-name=topology_verify66
#SBATCH --partition=128x24
#SBATCH --account=128x24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --dependency=afterany:796204
#SBATCH --chdir=/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261006T234134Z_topology_c655cc0c/launch
#SBATCH --output=/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261006T234134Z_topology_c655cc0c/logs/verify-%j.out
#SBATCH --error=/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261006T234134Z_topology_c655cc0c/logs/verify-%j.err
#SBATCH --mail-user=researcher@example.org
#SBATCH --mail-type=END,FAIL
#SBATCH --no-requeue
set -euo pipefail
umask 027
set +u
module load miniconda3
set -u
python3 -B /hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261006T234134Z_topology_c655cc0c/scripts/verify_supplement.py /hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261006T234134Z_topology_c655cc0c
