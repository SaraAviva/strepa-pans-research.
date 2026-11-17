import json
import pathlib
import shutil

here=pathlib.Path(__file__).resolve().parent
run=json.loads((here/'setup_result.json').read_text())['run_dir']
package=here/'package'
(package/'scripts').mkdir(parents=True,exist_ok=True)
(package/'runtime/bin').mkdir(parents=True,exist_ok=True)
(package/'provenance').mkdir(parents=True,exist_ok=True)
cache='/hb/groups/cornejo_training/StrepA_project_revamped/tools/interproscan6/singularity_cache'
db='/hb/groups/cornejo_training/StrepA_project_revamped/tools/interproscan6/data'
groups={
 'interpro-hmmer-3.3.img':'SEARCH_ANTIFAM|SEARCH_NCBIFAM|SEARCH_PFAM|SEARCH_PIRSR|SEARCH_PIRSF|SEARCH_SUPERFAMILY',
 'interpro-cath-0.16.10.img':'SEARCH_GENE3D|SEARCH_FUNFAM',
 'interpro-cdd-1.0.img':'SEARCH_CDD',
 'interpro-ncoils-2.2.1.img':'RUN_COILS',
 'interpro-pftools-3.2.12.img':'RUN_PFSEARCH|RUN_PSSCAN',
 'interpro-idrpred-1.0.3.img':'RUN_MOBIDBLITE',
 'interpro-panther-1.0.img':'SEARCH_PANTHER|RUN_TREEGRAFTER',
 'interpro-fingerprintscan-3.597.ebiftp.img':'RUN_PRINTS',
 'interpro-sfld-ebiftp.img':'SEARCH_SFLD',
 'interpro-smart-1.0.img':'PREFILTER_SMART|SEARCH_SMART',
 'interpro-groovy-4.0.27-1.img':'ADD_XREFS|REPRESENTATIVE_LOCATIONS',
 'interpro-download-1.0.img':'DOWNLOAD_INTERPRO|DOWNLOAD_DATABASE'}
def put(name,text,executable=False):
    p=package/name
    with p.open('x') as f: f.write(text.replace('@RUN@',run).replace('@DB@',db))
    if executable:p.chmod(0o750)

put('provenance/container_paths.txt','\n'.join(cache+'/'+x for x in groups)+'\n')
config='''// Reuse installed 6.0.1 configuration with bounded local execution in SLURM.
executor {
    name = 'local'
    cpus = 6
    memory = '24 GB'
    queueSize = 6
}
process {
    executor = 'local'
    cpus = 1
    errorStrategy = 'terminate'
    maxRetries = 0
    withLabel: dynamic { cpus = 1 }
    withName: 'DOWNLOAD_INTERPRO|DOWNLOAD_DATABASE' {
        beforeScript = 'echo "Database downloads are forbidden for this run" >&2; exit 74'
        errorStrategy = 'terminate'
        maxRetries = 0
    }
'''
for img,processes in groups.items():
    config+="    withName: '"+processes+"' { container = '"+cache+'/'+img+"' }\n"
config+='''}
singularity {
    enabled = true
    autoMounts = true
    cacheDir = '@RUN@/runtime/singularity_cache'
    runOptions = '-B /mnt/beegfs -B /hb'
}
trace {
    fields = 'task_id,hash,native_id,name,status,exit,submit,duration,realtime,cpus,memory,peak_rss,container,workdir'
    raw = true
}
'''
put('scripts/run.config',config)
put('runtime/bin/singularity','''#!/bin/bash
# Existing user-installed FUSE support; no software installation or security changes.
export LD_LIBRARY_PATH="/home/RESEARCH_USER/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PATH="/home/RESEARCH_USER/bin:$PATH"
exec /hb/software/utils/singularity-ce/gnu13-4.1.4/bin/singularity "$@"
''',True)
put('scripts/environment.sh','''#!/bin/bash
set +u
module load miniconda3 nextflow/25.10.4 singularity-ce/singularity-ce.4.1.4
export PATH="@RUN@/runtime/bin:$PATH"
export NXF_HOME="@RUN@/runtime/nextflow_home"
export NXF_TEMP="@RUN@/runtime/tmp"
export NXF_OFFLINE=true
export NXF_VER=25.10.4
export NXF_OPTS='-Xms256m -Xmx4g -XX:ActiveProcessorCount=8'
export SINGULARITY_CACHEDIR="@RUN@/runtime/singularity_cache"
export SINGULARITY_TMPDIR="@RUN@/runtime/tmp"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
set -u
''',True)
put('scripts/analysis.sbatch','''#!/bin/bash
#SBATCH --job-name=ipr_all66_b2a74f5f
#SBATCH --partition=128x24
#SBATCH --account=128x24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --no-requeue
#SBATCH --chdir=@RUN@/launch
#SBATCH --output=@RUN@/logs/analysis-%j.out
#SBATCH --error=@RUN@/logs/analysis-%j.err
set -euo pipefail
umask 027
run='@RUN@'
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
 printf 'Checking container %s\n' "$container"
 singularity exec --no-home -B /mnt/beegfs -B /hb "$container" /bin/sh -c 'test -r "$1" && test -r "$2"' sh "$run/input/all_significant_proteins.fasta" '@DB@/interpro/109.0/databases.json'
done < "$run/provenance/container_paths.txt"
while IFS= read -r container; do sha256sum "$container"; done < "$run/provenance/container_paths.txt" > "$run/provenance/container_hashes.sha256"
pipeline="$run/runtime/pipeline"
cmd=(nextflow -log "$run/logs/nextflow.log" -C "$pipeline/nextflow.config,$run/scripts/run.config" run "$pipeline" -profile singularity -ansi-log false -w "$run/work" -with-trace "$run/logs/trace.tsv" -with-report "$run/logs/report.html" -with-timeline "$run/logs/timeline.html" --input "$run/input/all_significant_proteins.fasta" --datadir '@DB@' --interpro 109.0 --applications 'AntiFam,CATH-Gene3D,CATH-FunFam,CDD,COILS,HAMAP,MobiDB-lite,NCBIFAM,PANTHER,Pfam,PIRSF,PIRSR,PRINTS,PROSITE-patterns,PROSITE-profiles,SFLD,SMART,SUPERFAMILY' --formats tsv,json,xml --outdir "$run/results" --outprefix all_significant_proteins --goterms --pathways --no-matches-api --skip-interpro-version-check --cpus 1 --max-workers 6)
printf '%q ' "${cmd[@]}" > "$run/provenance/nextflow_command.txt"
printf '\n' >> "$run/provenance/nextflow_command.txt"
set +e
"${cmd[@]}"
workflow_rc=$?
set -e
/usr/bin/python3 - "$run" "$workflow_rc" "$SLURM_JOB_ID" <<'PY'
import datetime,json,pathlib,sys
p=pathlib.Path(sys.argv[1])/'provenance/workflow_completion.json'
with p.open('x') as f:
 json.dump(dict(exit_code=int(sys.argv[2]),slurm_job_id=sys.argv[3],finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f,indent=2)
 f.write('\\n')
PY
exit "$workflow_rc"
''',True)
put('scripts/postprocess.sbatch','''#!/bin/bash
#SBATCH --job-name=ipr_verify66_b2a74f5f
#SBATCH --partition=128x24
#SBATCH --account=128x24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --no-requeue
#SBATCH --chdir=@RUN@
#SBATCH --output=@RUN@/logs/postprocess-%j.out
#SBATCH --error=@RUN@/logs/postprocess-%j.err
set -euo pipefail
umask 027
test -n "${SLURM_JOB_ID:-}"
exec /usr/bin/python3 '@RUN@/scripts/postprocess.py'
''',True)
put('scripts/status.sh','''#!/bin/bash
set -euo pipefail
run='@RUN@'
printf 'Run directory: %s\\n' "$run"
ids=()
for role in analysis postprocess; do
 if [ -f "$run/provenance/${role}_job_id.txt" ]; then
  id=$(cat "$run/provenance/${role}_job_id.txt")
  printf '%s job: %s\\n' "$role" "$id"
  ids+=("$id")
 fi
done
if [ "${#ids[@]}" -gt 0 ]; then
 jobs=$(IFS=,; echo "${ids[*]}")
 squeue -j "$jobs" -o '%.18i %.12P %.28j %.10T %.10M %R' || true
 sacct -j "$jobs" --format=JobID,JobName%28,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS,NodeList
fi
if [ -f "$run/summary/RUN_SUMMARY.txt" ]; then cat "$run/summary/RUN_SUMMARY.txt"; fi
if [ -f "$run/logs/nextflow.log" ]; then tail -18 "$run/logs/nextflow.log"; fi
''',True)
put('scripts/submit.py','''#!/usr/bin/env python3
import datetime,json,pathlib,re,shlex,subprocess
run=pathlib.Path(__file__).resolve().parent.parent
def submit(role,args):
 command=['sbatch','--parsable']+args+[str(run/'scripts'/(role+'.sbatch'))]
 marker=run/'provenance'/(role+'_submission_command.txt')
 with marker.open('x') as f:
  f.write(shlex.join(command)+'\\n'); f.flush()
 result=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
 with (run/'provenance'/(role+'_submission_receipt.json')).open('x') as f:
  json.dump(dict(command=command,returncode=result.returncode,stdout=result.stdout,stderr=result.stderr,time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f,indent=2)
 if result.returncode!=0:
  raise RuntimeError(role+' submission failed: '+result.stderr)
 jobid=result.stdout.strip().split(';')[0]
 if not re.fullmatch(r'[0-9]+',jobid):
  raise RuntimeError('Unclear submission acknowledgement; inspect squeue/sacct before retrying: '+result.stdout)
 with (run/'provenance'/(role+'_job_id.txt')).open('x') as f:
  f.write(jobid+'\\n'); f.flush()
 print(role+' job ID: '+jobid,flush=True)
 return jobid
analysis=submit('analysis',[])
submit('postprocess',['--dependency=afterok:'+analysis,'--kill-on-invalid-dep=yes'])
''',True)
for name in ['preflight.py','postprocess.py']:
    shutil.copy2(here/name,package/'scripts'/name)
shutil.copy2(here/'setup_run.py',package/'provenance/setup_run.py')
print(run)
