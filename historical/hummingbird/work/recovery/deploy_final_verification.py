import hashlib,json,pathlib,subprocess

base=pathlib.Path(__file__).resolve().parent
run='/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261003T043208Z_b2a74f5f'
sbatch=(base/'retry_scripts/postprocess_retry_20261004.sbatch').read_text().replace('ipr_verify66_retry_b2a74f5f','ipr_verify66_final_b2a74f5f').replace('postprocess_retry_20261004.py','postprocess_final_20261004.py')
gate="""/usr/bin/python3 - <<'PY'
import subprocess
rows=subprocess.check_output(['sacct','-X','-n','-P','-j','792566','--format=JobID,State,ExitCode'],universal_newlines=True).splitlines()
assert '792566|COMPLETED|0:0' in [r.rstrip('|') for r in rows],rows
PY
"""
sbatch=sbatch.replace('exec /usr/bin/python3',gate+'exec /usr/bin/python3')
(base/'postprocess_final_20261004.sbatch').write_text(sbatch)
note=dict(analysis_job='792566',failed_verification_job='792567',
 causes=['MD5 hexadecimal text was compared case-sensitively.','The success phrase is in analysis stdout, not nextflow.log.','CDD emitted an explicit empty map for group_2090, which REPORT_NO_MATCHES excludes but COMBINE_MATCHES retains.'],
 corrections=['Normalize only validated 32-character hexadecimal MD5 strings; retain exact sequence SHA-256 and length checks.','Require the correct analysis stdout success phrase and Nextflow all-processes-finished evidence.','Accept explicit empty maps only from successful REPORT_NO_MATCHES or COMBINE_MATCHES tasks in this run trace; final JSON must also explicitly contain an empty matches array.'],
 validation={'regression_cases':json.loads((base/'final_verifier_test_results.json').read_text()),'actual_records_verified_locally':66,'actual_with_matches':51,'actual_no_hits':15,'successful_workflow_tasks':51},
 scheduling='Analysis 792566 is already COMPLETED 0:0 in SLURM accounting and has left the active controller, so a new afterok dependency cannot refer to it. This verification-only recheck explicitly gates on accounting COMPLETED 0:0 both before submission and inside its allocation. The original verification 792567 used afterok. No analysis resubmission.',
 preserved='The failed summary and overwritten provenance snapshots are archived in recovery_20261004T163308Z/verification_792567_failed; original scripts, logs, and raw analysis results remain intact.')
files={'scripts/postprocess_final_20261004.py':(base/'postprocess_final_20261004.py').read_text(),
       'scripts/postprocess_final_20261004.sbatch':sbatch,
       'provenance/final_verifier_correction.json':json.dumps(note,indent=2)+'\n'}
remote="""import pathlib,json,subprocess,shutil,hashlib,datetime,shlex,re
run=pathlib.Path(RUN)
files=FILES
def state(job):
 out=subprocess.check_output(['sacct','-X','-n','-P','-j',job,'--format=JobID,State,ExitCode'],universal_newlines=True)
 rows=[x.split('|') for x in out.splitlines() if x.startswith(job+'|')]
 assert len(rows)==1,out
 return rows[0][1:3]
assert state('792566')==['COMPLETED','0:0']
assert state('792567')==['FAILED','2:0']
queue=subprocess.check_output(['squeue','-h','-u','research_user','-o','%i %j %T'],universal_newlines=True)
assert not queue.strip(),queue
assert not (run/'provenance/final_postprocess_job_id.txt').exists()
assert all(not (run/p).exists() for p in files)
old=json.loads((run/'summary/run_summary.json').read_text())
assert old['verified_complete'] is False and old['input_records']==66
archive=run/'recovery_20261004T163308Z/verification_792567_failed'
archive.mkdir()
shutil.move(str(run/'summary'),str(archive/'summary'))
(run/'summary').mkdir()
for name in ['analysis_accounting_at_postprocess.txt','results_and_deliverables.sha256']:
 if (run/'provenance'/name).exists(): shutil.copy2(str(run/'provenance'/name),str(archive/name))
for rel,contents in files.items():
 with (run/rel).open('x') as f: f.write(contents)
for p in [run/'scripts/postprocess_final_20261004.py',run/'scripts/postprocess_final_20261004.sbatch']:
 p.chmod(0o750)
command=['sbatch','--hold','--parsable',str(run/'scripts/postprocess_final_20261004.sbatch')]
with (run/'provenance/final_postprocess_submission_command.txt').open('x') as f: f.write(shlex.join(command)+'\\n')
proc=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
receipt=dict(command=command,returncode=proc.returncode,stdout=proc.stdout,stderr=proc.stderr,submitted_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
with (run/'provenance/final_postprocess_submission_receipt.json').open('x') as f: json.dump(receipt,f,indent=2)
assert proc.returncode==0,receipt
job=proc.stdout.strip().split(';')[0]
assert re.fullmatch(r'[0-9]+',job),receipt
with (run/'provenance/final_postprocess_job_id.txt').open('x') as f:f.write(job+'\\n')
release=subprocess.run(['scontrol','release',job],stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
with (run/'provenance/final_postprocess_release_receipt.json').open('x') as f:json.dump(dict(job_id=job,returncode=release.returncode,stdout=release.stdout,stderr=release.stderr),f,indent=2)
print(json.dumps(dict(job_id=job,submission=receipt,release_returncode=release.returncode)),flush=True)
assert release.returncode==0,release.stderr
""".replace('RUN',repr(run),1).replace('FILES',repr(files),1)
proc=subprocess.run(['ssh','-S','/path/to/ssh-session/control','-o','BatchMode=yes','researcher@example.org','python3 -'],input=remote,capture_output=True,text=True)
(base/'final_verification_submission.json').write_text(json.dumps({'returncode':proc.returncode,'stdout':proc.stdout,'stderr':proc.stderr},indent=2)+'\n')
print(proc.stdout,proc.stderr)
assert proc.returncode==0
