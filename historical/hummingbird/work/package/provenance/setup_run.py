import csv
import hashlib
import json
import os
import pathlib
import shutil
import uuid
from collections import Counter, defaultdict
from datetime import datetime, timezone

os.umask(0o027)
parent = pathlib.Path('/hb/groups/cornejo_training/StrepA_project_revamped')
source = parent / 'significant_proteins/all_significant_proteins.fasta'
expected = '3285311fbbe7d7088ff63a24ea515e680cabd1e0d3519803d47b24deae04aedb'
data = source.read_bytes()
assert data and hashlib.sha256(data).hexdigest() == expected, 'Input changed since inspection'
run = parent / ('interproscan_all_significant_' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '_' + uuid.uuid4().hex[:8])
run.mkdir(mode=0o750)
for name in ['scripts', 'input', 'logs', 'provenance', 'results', 'summary', 'launch', 'work', 'runtime', 'runtime/bin', 'runtime/tmp', 'runtime/nextflow_home', 'runtime/singularity_cache']:
    (run / name).mkdir(mode=0o750)
snapshot = run / 'input/all_significant_proteins.fasta'
with snapshot.open('xb') as f:
    f.write(data)
assert hashlib.sha256(snapshot.read_bytes()).hexdigest() == expected
assert hashlib.sha256(source.read_bytes()).hexdigest() == expected
snapshot.chmod(0o440)
records = []
header = None
sequence = []
for line in data.decode('utf-8').splitlines():
    if line.startswith('>'):
        if header is not None:
            records.append((header, ''.join(sequence)))
        header, sequence = line[1:], []
    else:
        assert header is not None or not line.strip(), 'Content before FASTA header'
        sequence.append(line)
if header is not None:
    records.append((header, ''.join(sequence)))
ids = Counter(h.split()[0] if h.split() else '' for h, s in records)
groups = defaultdict(list)
allowed = set('ACDEFGHIKLMNPQRSTVWYUOXBZJ')
audit_rows = []
for n, (h, s) in enumerate(records, 1):
    ident = h.split()[0] if h.split() else ''
    groups[s].append(ident)
    audit_rows.append(dict(record_number=n, original_id=ident, original_header=h, length=len(s), sequence_md5=hashlib.md5(s.encode()).hexdigest(), sequence_sha256=hashlib.sha256(s.encode()).hexdigest(), gap_count=sum(s.count(c) for c in '-.'), invalid_symbols=''.join(sorted(set(s)-allowed)), x_count=s.count('X')))
audit = dict(source=str(source), snapshot=str(snapshot), sha256=expected, bytes=len(data), records=len(records), unique_ids=len(ids), total_residues=sum(len(s) for h,s in records), minimum_length=min(len(s) for h,s in records), maximum_length=max(len(s) for h,s in records), empty_records=[r['original_id'] for r in audit_rows if not r['length']], duplicate_ids={k:v for k,v in ids.items() if v>1}, duplicate_sequence_groups=[v for v in groups.values() if len(v)>1], gap_characters=sum(r['gap_count'] for r in audit_rows), invalid_records=[r for r in audit_rows if r['invalid_symbols']], symbol_counts=dict(sorted(Counter(''.join(s for h,s in records)).items())), alphabet_source='runtime/pipeline/lib/uk/ac/ebi/interpro/FastaFile.groovy:PROTEIN_ALPHABET', input_modifications='None. Snapshot is byte-for-byte identical to supplied FASTA.')
(run/'provenance/input_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
with (run/'provenance/input_records.tsv').open('x',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(audit_rows[0]),delimiter='\t'); writer.writeheader(); writer.writerows(audit_rows)
(run/'provenance/input.sha256').write_text(expected+'  input/all_significant_proteins.fasta\n')
assert len(records)==66 and len(ids)==66 and '' not in ids
assert not audit['empty_records'] and not audit['invalid_records'] and not audit['gap_characters']
installed = pathlib.Path('/home/RESEARCH_USER/.nextflow/assets/ebi-pf-team/interproscan6')
pipeline = run/'runtime/pipeline'
pipeline.mkdir()
for name in ['conf','lib','modules','subworkflows','bin','workflows']:
    shutil.copytree(installed/name,pipeline/name)
for name in ['nextflow.config','main.nf','README.md','LICENSE']:
    shutil.copy2(installed/name,pipeline/name)
manifest=[]
for f in sorted(pipeline.rglob('*')):
    if f.is_file():
        original=installed/f.relative_to(pipeline)
        digest=hashlib.sha256(f.read_bytes()).hexdigest()
        assert digest==hashlib.sha256(original.read_bytes()).hexdigest()
        manifest.append(digest+'  '+str(f.relative_to(run)))
(run/'provenance/pipeline.sha256').write_text('\n'.join(manifest)+'\n')
shutil.copy2('/home/RESEARCH_USER/interproscan6_run/.nextflow.log',run/'provenance/previous_success.nextflow.log')
shutil.copy2('/home/RESEARCH_USER/pans_interproscan_gt90/.nextflow.log',run/'provenance/previous_runtime_failure.nextflow.log')
dbroot=parent/'tools/interproscan6/data'
shutil.copy2(dbroot/'interpro/109.0/databases.json',run/'provenance/database_versions.json')
meta=dict(run_dir=str(run),created_utc=datetime.now(timezone.utc).isoformat(),source=str(source),input_sha256=expected,software_version='6.0.1',pipeline_original=str(installed),pipeline_git_commit='9fe11bce400d185cfc194d583c180e2423f41c37',nextflow_version='25.10.4',singularity_version='4.1.4',interpro_version='109.0',database_root=str(dbroot),execution_model='Nextflow local tasks inside one SLURM allocation; no nested SLURM jobs',resources=dict(cpus=8,memory_gb=32,time_limit='02:00:00',local_task_cpus=6,local_task_memory_gb=24,local_queue_size=6),enabled_analyses=['AntiFam','CATH-Gene3D','CATH-FunFam','CDD','COILS','HAMAP','MobiDB-lite','NCBIFAM','PANTHER','Pfam','PIRSF','PIRSR','PRINTS','PROSITE patterns','PROSITE profiles','SFLD','SMART','SUPERFAMILY'],skipped_analyses={'Phobius':'Licensed application directory unconfigured in installed setup','DeepTMHMM':'ML disabled; licensed application unconfigured','InterPro-N':'ML disabled; licensed application unconfigured','SignalP-Euk':'ML disabled; licensed application unconfigured','SignalP-Prok':'ML disabled; licensed application unconfigured','TMbed':'ML disabled in established standard workflow'},decisions=['All 66 supplied records included; no selection or deduplication by this run.','Force local computation with Matches API disabled.','Reuse installed release 109.0; disable online version check using prior successful exact-version run as compatibility evidence.','Pin every container to an existing local image; fail download processes if unexpectedly invoked.','Use existing personal FUSE runtime support via a run-local wrapper, addressing the documented prior missing-library failure.'])
(run/'provenance/run_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
(run/'provenance/change_log.txt').write_text('Initial setup: original FASTA snapshotted and verified; installed pipeline copied without source modifications; previous logs retained.\n')
print(json.dumps(dict(run_dir=str(run),input_sha256=expected,records=len(records),snapshot_verified=True)))
