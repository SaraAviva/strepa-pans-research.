#!/usr/bin/env python3
"""Run the recovered synthetic cases in isolation, without a real-data rerun."""
from pathlib import Path
import json
import shutil
import subprocess
import sys
import tempfile

root = Path(__file__).resolve().parents[1]
source = root / 'historical/hummingbird/work'
marker = '# Recheck the actual 66-record result with the same script.'
original = (source / 'recovery/test_final_verifier.py').read_text()
if original.count(marker) != 1:
    raise SystemExit('The original test boundary changed; inspect it before adapting the runner.')
synthetic_cases = original.split(marker, 1)[0]

with tempfile.TemporaryDirectory(prefix='strepa-verifier-tests-') as temp:
    work = Path(temp) / 'work'
    recovery = work / 'recovery'
    recovery.mkdir(parents=True)
    shutil.copy2(source / 'test_postprocess.py', work / 'test_postprocess.py')
    shutil.copy2(source / 'recovery/postprocess_final_20261004.py',
                 recovery / 'postprocess_final_20261004.py')
    test = recovery / 'test_final_verifier.py'
    test.write_text(synthetic_cases)
    result = subprocess.run([sys.executable, str(test)], capture_output=True, text=True)
    if result.returncode:
        sys.stderr.write(result.stdout + result.stderr)
        raise SystemExit(result.returncode)
    checks = json.loads((recovery / 'final_verifier_test_results.json').read_text())
    if len(checks) != 10 or not all(c['passed'] for c in checks):
        raise SystemExit('Unexpected test count or unsuccessful original test case.')
    print(json.dumps({'synthetic_cases_passed': len(checks),
                      'cases': [c['case'] for c in checks],
                      'real_data_recheck': 'Not run',
                      'cluster_or_network_jobs': 'None'}, indent=2))
