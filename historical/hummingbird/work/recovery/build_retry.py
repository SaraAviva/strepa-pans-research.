import ast
import pathlib
import re
import subprocess

local = pathlib.Path(__file__).resolve().parent
original = local/'original/scripts'
out = local/'retry_scripts'
out.mkdir(exist_ok=True)
run = '/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261003T043208Z_b2a74f5f'
recovery = run+'/recovery_20261004T163308Z'

analysis = (original/'analysis.sbatch').read_text()
analysis = analysis.replace('ipr_all66_b2a74f5f','ipr_all66_retry_b2a74f5f')
analysis = analysis.replace('set -euo pipefail', '#SBATCH --mail-user=researcher@example.org\n#SBATCH --mail-type=FAIL\nset -euo pipefail', 1)
insertion = r'''
# Compute-node temporary space and an extracted copy of each unchanged image.
# This avoids the unavailable FUSE mount/extractor combination on compute nodes.
job_tmp=$(mktemp -d "/tmp/ipr-${USER}-${SLURM_JOB_ID}.XXXXXX")
stage=environment
finish_attempt() {
 rc=$?
 trap - EXIT
 /usr/bin/python3 - "$run" "$SLURM_JOB_ID" "$stage" "$rc" <<'PY'
import datetime,json,pathlib,sys
p=pathlib.Path(sys.argv[1])/'provenance'/('attempt-'+sys.argv[2]+'.json')
with p.open('x') as f:
 json.dump(dict(slurm_job_id=sys.argv[2],last_stage=sys.argv[3],exit_code=int(sys.argv[4]),finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f,indent=2)
PY
 case "$job_tmp" in /tmp/ipr-"${USER}"-"${SLURM_JOB_ID}".*) rm -rf -- "$job_tmp" ;; esac
 exit "$rc"
}
trap finish_attempt EXIT
export SINGULARITY_TMPDIR="$job_tmp/singularity"
export NXF_TEMP="$job_tmp/nextflow"
export TMPDIR="$job_tmp/tmp"
mkdir -p "$SINGULARITY_TMPDIR" "$NXF_TEMP" "$TMPDIR"
'''
analysis = analysis.replace('source "$run/scripts/environment.sh"\n', 'source "$run/scripts/environment.sh"\n'+insertion)
analysis = analysis.replace(' ldd /hb/software/utils/singularity-ce/gnu13-4.1.4/libexec/singularity/bin/squashfuse_ll', ' df -PT "$job_tmp"')
analysis = analysis.replace('/usr/bin/python3 "$run/scripts/preflight.py"', 'stage=input_and_container_preflight\n/usr/bin/python3 "$run/scripts/preflight.py"')
analysis = analysis.replace('done < "$run/provenance/container_paths.txt"', 'done < "'+recovery+'/sandbox_paths.txt"', 1)
analysis = analysis.replace('$run/scripts/run.config','$run/scripts/run_retry_20261004.config')
analysis = analysis.replace('set +e\n"${cmd[@]}"', 'stage=interproscan\nset +e\n"${cmd[@]}"')
(out/'analysis_retry_20261004.sbatch').write_text(analysis)

config = (original/'run.config').read_text()
source_prefix = '/hb/groups/cornejo_training/StrepA_project_revamped/tools/interproscan6/singularity_cache/'
config = re.sub(re.escape(source_prefix)+r'([^\s\x27]+)\.img', lambda m:recovery+'/sandboxes/'+m.group(1), config)
(out/'run_retry_20261004.config').write_text(config)

post = (original/'postprocess.py').read_text()
post = post.replace("provenance/run_metadata.json", "provenance/retry_20261004_run_metadata.json")
post = post.replace("'provenance/'+role+'_job_id.txt'", "'provenance/retry_20261004_'+role+'_job_id.txt'")
post = post.replace('provenance/analysis_submission_command.txt', 'provenance/retry_20261004_analysis_submission_command.txt')
post = post.replace('provenance/postprocess_submission_command.txt', 'provenance/retry_20261004_postprocess_submission_command.txt')
ast.parse(post)
(out/'postprocess_retry_20261004.py').write_text(post)
postbatch = (original/'postprocess.sbatch').read_text().replace('ipr_verify66_b2a74f5f','ipr_verify66_retry_b2a74f5f')
postbatch = postbatch.replace('set -euo pipefail', '#SBATCH --mail-user=researcher@example.org\n#SBATCH --mail-type=END,FAIL\nset -euo pipefail', 1)
postbatch = postbatch.replace('/scripts/postprocess.py', '/scripts/postprocess_retry_20261004.py')
(out/'postprocess_retry_20261004.sbatch').write_text(postbatch)

status = pathlib.Path('work/package/scripts/status.sh').read_text()
status = status.replace('${role}_job_id.txt','retry_20261004_${role}_job_id.txt')
(out/'status_retry_20261004.sh').write_text(status)
for path in out.iterdir():
    if path.suffix in ('.sbatch','.sh'):
        subprocess.run(['bash','-n',str(path)],check=True)
print('Prepared retry scripts; shell and Python syntax checks passed.')
