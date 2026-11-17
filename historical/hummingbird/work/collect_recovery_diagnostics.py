#!/usr/bin/env python3
"""Read-only remote diagnostics for the failed all-protein InterProScan run."""
import datetime
import json
import pathlib
import subprocess

run = pathlib.Path('/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261003T043208Z_b2a74f5f')
report = {'collected_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
          'run': str(run), 'commands': [], 'files': {}, 'directories': {}}

for command in [
    ['hostname'], ['id'],
    ['squeue', '-u', 'research_user', '-o', '%.18i %.12P %.32j %.12T %.10M %R'],
    ['sacct', '-S', '2026-10-02', '-j', '785450,785451', '--parsable2',
     '--format=JobID,JobName%40,State,ExitCode,Elapsed,Start,End,AllocCPUS,ReqMem,MaxRSS,NodeList,Reason'],
    ['sinfo', '-o', '%P %a %l %D %C %m'],
]:
    try:
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True, timeout=30)
        report['commands'].append({'argv': command, 'returncode': result.returncode,
                                   'stdout': result.stdout, 'stderr': result.stderr})
    except Exception as exc:
        report['commands'].append({'argv': command, 'error': str(exc)})

paths = [run/'runtime/bin/singularity', run/'launch/.nextflow.log',
         run/'provenance/previous_runtime_failure.nextflow.log']
for directory in ['logs', 'scripts', 'results', 'summary', 'provenance', 'launch']:
    parent = run/directory
    try:
        report['directories'][directory] = [
            {'name': p.name, 'bytes': p.stat().st_size, 'directory': p.is_dir()}
            for p in sorted(parent.iterdir())]
    except Exception as exc:
        report['directories'][directory] = {'error': str(exc)}
    if directory in ('logs', 'scripts'):
        paths.extend(p for p in parent.iterdir() if p.is_file() and p.suffix not in ('.html', '.png'))
    if directory == 'provenance':
        for pattern in ['runtime-*.txt', 'database_preflight_*.json', 'workflow_completion.json',
                        'run_metadata.json', '*job_id.txt', 'nextflow_command.txt']:
            paths.extend(parent.glob(pattern))

for path in sorted(set(paths)):
    if path.is_file():
        with path.open('rb') as stream:
            size = path.stat().st_size
            if size > 100000:
                stream.seek(size - 100000)
            content = stream.read().decode('utf-8', errors='replace')
        report['files'][str(path.relative_to(run))] = {
            'total_bytes': size, 'tail_only': size > 100000, 'content': content}

print(json.dumps(report, indent=2))
