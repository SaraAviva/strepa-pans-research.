import json,pathlib,subprocess
base=pathlib.Path(__file__).resolve().parent
remote=r'''import pathlib,json,subprocess,hashlib,tarfile
run=pathlib.Path('/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261003T043208Z_b2a74f5f')
summary=json.loads((run/'summary/run_summary.json').read_text())
assert summary['verified_complete'] and summary['processed_records']==66 and not summary['issues']
accounting=subprocess.check_output(['sacct','-n','-P','-j','792566,792576','--format=JobID,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS,NodeList'],universal_newlines=True)
states={row.split('|')[0]:row.split('|')[1:3] for row in accounting.splitlines()}
assert states['792566']==states['792576']==['COMPLETED','0:0'],states
with (run/'provenance/final_slurm_accounting.txt').open('x') as f:f.write(accounting)
status="""#!/bin/bash
set -euo pipefail
run='/hb/groups/cornejo_training/StrepA_project_revamped/interproscan_all_significant_20261003T043208Z_b2a74f5f'
sacct -j 792566,792576 --format=JobID,JobName%28,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS,NodeList
cat "$run/summary/RUN_SUMMARY.txt"
"""
with (run/'scripts/status_final.sh').open('x') as f:f.write(status)
(run/'scripts/status_final.sh').chmod(0o750)
readme="""InterProScan results — VERIFIED COMPLETE

All 66 original proteins were processed and reconciled with JSON and XML results.
51 proteins have matches; 15 completed with no hits. No failed/unverified records.
44 have family/domain/functional signatures; 7 have only disorder/coiled-coil predictions.
43 have integrated InterPro entries; 30 have GO terms.
Matches are annotation evidence, not proof of virulence or novelty.

Start with summary/per_protein_summary.csv (one row per original protein).
The equivalent tab-separated table is summary/per_protein_summary.tsv.
summary/RUN_SUMMARY.txt gives versions, analyses, commands, checksums and paths.
results/ contains unchanged raw TSV, JSON and XML, including full coordinates.
input/ contains the unchanged supplied FASTA; provenance/ contains audit records.
The original failures and verifier corrections are preserved in logs/ and recovery_20261004T163308Z/.
Selected explicit no-match intermediate evidence is retained under work/.

Successful analysis: SLURM 792566, 13 minutes 30 seconds, COMPLETED 0:0.
Successful final verification: SLURM 792576, COMPLETED 0:0.
All 51 workflow tasks completed. No analysis rerun was needed to repair verification.

Hummingbird run directory (Evan can access the results here):
"""+str(run)+"\n\nPer-protein summary on Hummingbird:\n"+str(run/'summary/per_protein_summary.csv')+"\n\nRaw outputs on Hummingbird:\n"+str(run/'results')+"\n\nStatus command:\nbash "+str(run/'scripts/status_final.sh')+"\n\nNo requested analysis or verification remains outstanding. No messages were sent to Evan.\n"
with (run/'README_RESULTS.txt').open('x') as f:f.write(readme)
paths=[]
for folder in ['input','results','summary','provenance','scripts','logs']:
 paths.extend(p for p in (run/folder).rglob('*') if p.is_file())
archive=run/'recovery_20261004T163308Z/verification_792567_failed'
paths.extend(p for p in archive.rglob('*') if p.is_file())
paths.append(run/'recovery_20261004T163308Z/container_preparation.json')
paths.append(run/'README_RESULTS.txt')
for evidence in summary['no_match_evidence'].values():
 for item in evidence:
  p=run/item['path']
  assert hashlib.sha256(p.read_bytes()).hexdigest()==item['sha256']
  paths.append(p)
paths=sorted(set(paths))
manifest=run/'FINAL_DELIVERY_SHA256.txt'
with manifest.open('x') as f:
 for p in paths:f.write(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(run))+'\n')
paths.append(manifest)
bundle=run/'final_results_20261004.tar.gz'
with tarfile.open(str(bundle),'x:gz',dereference=True) as tar:
 for p in paths:tar.add(str(p),arcname=str(p.relative_to(run)),recursive=False)
print(json.dumps(dict(archive=str(bundle),bytes=bundle.stat().st_size,files=len(paths),sha256=hashlib.sha256(bundle.read_bytes()).hexdigest(),jobs=states)))
'''
proc=subprocess.run(['ssh','-S','/path/to/ssh-session/control','-o','BatchMode=yes','researcher@example.org','python3 -'],input=remote,capture_output=True,text=True)
(base/'final_bundle_receipt.json').write_text(json.dumps({'returncode':proc.returncode,'stdout':proc.stdout,'stderr':proc.stderr},indent=2)+'\n')
print(proc.stdout,proc.stderr)
assert proc.returncode==0
