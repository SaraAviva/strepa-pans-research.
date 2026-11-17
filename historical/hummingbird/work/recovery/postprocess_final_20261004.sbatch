#!/bin/bash
#SBATCH --job-name=ipr_verify66_final_b2a74f5f
#SBATCH --partition=128x24
#SBATCH --account=128x24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --no-requeue
#SBATCH --chdir=/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261003T043208Z_b2a74f5f
#SBATCH --output=/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261003T043208Z_b2a74f5f/logs/postprocess-%j.out
#SBATCH --error=/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261003T043208Z_b2a74f5f/logs/postprocess-%j.err
#SBATCH --mail-user=researcher@example.org
#SBATCH --mail-type=END,FAIL
set -euo pipefail
umask 027
test -n "${SLURM_JOB_ID:-}"
/usr/bin/python3 - <<'PY'
import subprocess
rows=subprocess.check_output(['sacct','-X','-n','-P','-j','792566','--format=JobID,State,ExitCode'],universal_newlines=True).splitlines()
assert '792566|COMPLETED|0:0' in [r.rstrip('|') for r in rows],rows
PY
exec /usr/bin/python3 '/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261003T043208Z_b2a74f5f/scripts/postprocess_final_20261004.py'
