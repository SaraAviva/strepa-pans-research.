#!/usr/bin/env python3
"""Extract unchanged existing images and verify them on an allocated compute node."""
import datetime
import hashlib
import json
import os
import pathlib
import re
import subprocess

run = pathlib.Path('/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261003T043208Z_b2a74f5f')
recovery = run/'recovery_20261004T163308Z'
assert os.environ.get('SLURM_JOB_ID'), 'Must run in a SLURM allocation'
assert not (recovery/'container_preparation.json').exists()
images = [pathlib.Path(p) for p in (run/'provenance/container_paths.txt').read_text().splitlines()]
records = []
for index, source in enumerate(images):
    sandbox = recovery/'sandboxes'/source.stem
    print('PREPARING '+str(source), flush=True)
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if sandbox.exists():
        assert index == 0, 'Unexpected existing sandbox: '+str(sandbox)
        first_log = (recovery/'extraction-smoke-792553.out').read_text()
        assert 'EXTRACTED CONTAINER EXECUTION PASSED' in first_log
        assert digest in first_log
    else:
        listing = subprocess.check_output(['singularity','sif','list',str(source)],universal_newlines=True)
        parts = [line.split('|') for line in listing.splitlines() if 'FS (Squashfs/*System/amd64)' in line]
        assert len(parts) == 1, listing
        offset = int(parts[0][3].strip().split('-')[0])
        env = dict(os.environ, LD_LIBRARY_PATH=str(recovery/'extraction_tools/lib'))
        command = [str(recovery/'extraction_tools/unsquashfs'),'-no-progress','-no-xattrs',
                   '-processors','1','-data-queue','32','-frag-queue','32',
                   '-offset',str(offset),'-dest',str(sandbox),str(source)]
        subprocess.run(command, env=env, check=True, timeout=300)
    probe = ['singularity','exec','--no-home','-B','/mnt/beegfs','-B','/hb',str(sandbox),
             '/bin/sh','-c','test -r "$1" && test -r "$2" && test -w "$3"', 'sh',
             str(run/'input/all_significant_proteins.fasta'),
             '/hb/groups/cornejo_training/StrepA_project_revamped/tools/interproscan6/data/interpro/109.0/databases.json',
             str(run/'work')]
    subprocess.run(probe, check=True, timeout=90)
    assert hashlib.sha256(source.read_bytes()).hexdigest() == digest
    records.append(dict(image=str(source),sha256=digest,sandbox=str(sandbox),execution_passed=True))
    print('PASS '+str(sandbox), flush=True)
report = dict(slurm_job_id=os.environ['SLURM_JOB_ID'],hostname=subprocess.check_output(['hostname'],universal_newlines=True).strip(),
              completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
              extraction='Existing login-node unsquashfs and dependency libraries copied byte-identically into this run; unchanged SIF system partitions extracted without privileged xattrs.',
              containers=records)
with (recovery/'container_preparation.json').open('x') as f:
    json.dump(report,f,indent=2); f.write('\n')
print('ALL '+str(len(records))+' CONTAINERS PASSED', flush=True)
