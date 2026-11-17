#!/usr/bin/env python3
"""Verify publication-copy integrity and syntax without executing research code."""
from pathlib import Path
import ast
import hashlib
import json
import shutil
import subprocess
import sys

root = Path(__file__).resolve().parents[1]
manifest = json.loads((root / 'provenance/source-manifest.json').read_text())
failures = []
counts = {'hashes': 0, 'python_syntax': 0, 'shell_syntax': 0, 'input_hashes': 0, 'source_pdf_hashes': 0}
for entry in manifest['files']:
    path = root / entry['path']
    if not path.is_file():
        failures.append({'path': entry['path'], 'problem': 'Missing file'})
        continue
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != entry['published_copy_sha256']:
        failures.append({'path': entry['path'], 'problem': 'Publication-copy hash mismatch'})
    counts['hashes'] += 1
    if path.suffix == '.py':
        try:
            ast.parse(raw.decode(), filename=entry['path'])
            counts['python_syntax'] += 1
        except SyntaxError as exc:
            failures.append({'path': entry['path'], 'problem': str(exc)})
    if path.suffix in {'.sh', '.sbatch'} or path.name == 'singularity':
        bash = shutil.which('bash')
        if bash:
            result = subprocess.run([bash, '-n', str(path)], capture_output=True, text=True)
            if result.returncode:
                failures.append({'path': entry['path'], 'problem': result.stderr})
            else:
                counts['shell_syntax'] += 1
        else:
            failures.append({'path': entry['path'], 'problem': 'Bash unavailable; syntax check skipped'})
for entry in json.loads((root/'provenance/located-inputs.json').read_text()):
    path = root/entry['path']
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != entry['sha256']:
        failures.append({'path': entry['path'], 'problem': 'Input missing or hash mismatch'})
    counts['input_hashes'] += 1
pdf = root/'provenance/source-documents/original_2015_hp_analysis_code.pdf'
for entry in manifest['files']:
    if 'source_document_sha256' in entry:
        if not pdf.is_file() or hashlib.sha256(pdf.read_bytes()).hexdigest() != entry['source_document_sha256']:
            failures.append({'path': pdf.relative_to(root).as_posix(), 'problem': 'Source PDF missing or hash mismatch'})
        counts['source_pdf_hashes'] += 1
print(json.dumps({'checks': counts, 'failures': failures}, indent=2))
sys.exit(bool(failures))
