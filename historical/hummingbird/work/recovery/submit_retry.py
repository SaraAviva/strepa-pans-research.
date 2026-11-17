#!/usr/bin/env python3
import datetime
import hashlib
import json
import pathlib
import re
import shlex
import subprocess

run = pathlib.Path('/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261003T043208Z_b2a74f5f')
recovery = run/'recovery_20261004T163308Z'
prefix = 'retry_20261004_'
def state(job):
    text = subprocess.check_output(['sacct','-X','-n','-P','-j',str(job),'--format=JobID,State,ExitCode'],universal_newlines=True)
    rows = [line.split('|') for line in text.splitlines() if line.startswith(str(job)+'|')]
    assert len(rows) == 1, text
    return rows[0][1:3]

assert state(785450)[0] == 'FAILED'
assert state(785451)[0].startswith('CANCELLED')
assert state(792553) == ['COMPLETED','0:0']
assert state(792556) == ['COMPLETED','0:0']
assert not (run/'provenance/workflow_completion.json').exists()
assert not any((run/'results').iterdir())
assert not any((run/'summary').iterdir())
queue = subprocess.check_output(['squeue','-h','-u','research_user','-o','%i %j %T'],universal_newlines=True)
assert not queue.strip(), 'Active jobs require inspection before resubmission: '+queue
proof = json.loads((recovery/'container_preparation.json').read_text())
images = (run/'provenance/container_paths.txt').read_text().splitlines()
assert len(proof['containers']) == len(images)
assert {p['image'] for p in proof['containers']} == set(images)
assert all(p['execution_passed'] and pathlib.Path(p['sandbox']).is_dir() for p in proof['containers'])
with (recovery/'sandbox_paths.txt').open('x') as f:
    f.write('\n'.join(p['sandbox'] for p in proof['containers'])+'\n')
meta = json.loads((run/'provenance/retry_20261004_run_metadata.json').read_text())
assert hashlib.sha256((run/'input/all_significant_proteins.fasta').read_bytes()).hexdigest() == meta['input_sha256']
assert hashlib.sha256(pathlib.Path(meta['source']).read_bytes()).hexdigest() == meta['input_sha256']

def submit(role, extra):
    command = ['sbatch','--parsable']+extra+[str(run/'scripts'/(role+'_retry_20261004.sbatch'))]
    with (run/'provenance'/(prefix+role+'_submission_command.txt')).open('x') as f:
        f.write(shlex.join(command)+'\n')
    result = subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    receipt = dict(command=command,returncode=result.returncode,stdout=result.stdout,stderr=result.stderr,
                   submitted_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    with (run/'provenance'/(prefix+role+'_submission_receipt.json')).open('x') as f:
        json.dump(receipt,f,indent=2); f.write('\n')
    assert result.returncode == 0, receipt
    job = result.stdout.strip().split(';')[0]
    assert re.fullmatch(r'[0-9]+',job), 'Ambiguous acknowledgement: inspect queue before retrying'
    with (run/'provenance'/(prefix+role+'_job_id.txt')).open('x') as f:
        f.write(job+'\n')
    print(json.dumps(dict(role=role,job_id=job,receipt=receipt)),flush=True)
    return job

analysis = submit('analysis',[])
submit('postprocess',['--dependency=afterok:'+analysis,'--kill-on-invalid-dep=yes'])
