#!/bin/bash
#SBATCH --job-name=ipr_all66_b2a74f5f
#SBATCH --partition=128x24
#SBATCH --account=128x24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --no-requeue
#SBATCH --chdir=/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261003T043208Z_b2a74f5f/launch
#SBATCH --output=/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261003T043208Z_b2a74f5f/logs/analysis-%j.out
#SBATCH --error=/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261003T043208Z_b2a74f5f/logs/analysis-%j.err
set -euo pipefail
umask 027
run='/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261003T043208Z_b2a74f5f'
test -n "${SLURM_JOB_ID:-}"
test "${SLURM_CPUS_PER_TASK:-0}" -eq 8
source "$run/scripts/environment.sh"
exec 9>"$run/provenance/analysis.lock"
flock -n 9 || { echo 'Another analysis process holds the run lock' >&2; exit 75; }
test ! -e "$run/provenance/workflow_completion.json"
{
 date -u
 hostname
 module list 2>&1
 nextflow -version
 singularity --version
 java -version 2>&1
 ldd /hb/software/utils/singularity-ce/gnu13-4.1.4/libexec/singularity/bin/squashfuse_ll
} > "$run/provenance/runtime-${SLURM_JOB_ID}.txt" 2>&1
/usr/bin/python3 "$run/scripts/preflight.py"
while IFS= read -r container; do
 printf 'Checking container %s
' "$container"
 singularity exec --no-home -B /mnt/beegfs -B /hb "$container" /bin/sh -c 'test -r "$1" && test -r "$2"' sh "$run/input/all_significant_proteins.fasta" '/hb/groups/cornejo_training/StrepA_project_revamped/tools/interproscan6/data/interpro/109.0/databases.json'
done < "$run/provenance/container_paths.txt"
while IFS= read -r container; do sha256sum "$container"; done < "$run/provenance/container_paths.txt" > "$run/provenance/container_hashes.sha256"
pipeline="$run/runtime/pipeline"
cmd=(nextflow -log "$run/logs/nextflow.log" -C "$pipeline/nextflow.config,$run/scripts/run.config" run "$pipeline" -profile singularity -ansi-log false -w "$run/work" -with-trace "$run/logs/trace.tsv" -with-report "$run/logs/report.html" -with-timeline "$run/logs/timeline.html" --input "$run/input/all_significant_proteins.fasta" --datadir '/hb/groups/cornejo_training/StrepA_project_revamped/tools/interproscan6/data' --interpro 109.0 --applications 'AntiFam,CATH-Gene3D,CATH-FunFam,CDD,COILS,HAMAP,MobiDB-lite,NCBIFAM,PANTHER,Pfam,PIRSF,PIRSR,PRINTS,PROSITE-patterns,PROSITE-profiles,SFLD,SMART,SUPERFAMILY' --formats tsv,json,xml --outdir "$run/results" --outprefix all_significant_proteins --goterms --pathways --no-matches-api --skip-interpro-version-check --cpus 1 --max-workers 6)
printf '%q ' "${cmd[@]}" > "$run/provenance/nextflow_command.txt"
printf '
' >> "$run/provenance/nextflow_command.txt"
set +e
"${cmd[@]}"
workflow_rc=$?
set -e
/usr/bin/python3 - "$run" "$workflow_rc" "$SLURM_JOB_ID" <<'PY'
import datetime,json,pathlib,sys
p=pathlib.Path(sys.argv[1])/'provenance/workflow_completion.json'
with p.open('x') as f:
 json.dump(dict(exit_code=int(sys.argv[2]),slurm_job_id=sys.argv[3],finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f,indent=2)
 f.write('\n')
PY
exit "$workflow_rc"
