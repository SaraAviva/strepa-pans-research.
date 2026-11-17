#!/usr/bin/env python3
"""Reconcile original records against structured outputs and workflow evidence."""
import csv
import hashlib
import json
import os
import pathlib
import re
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from urllib.parse import quote

run = pathlib.Path(__file__).resolve().parent.parent
meta = json.loads((run/'provenance/run_metadata.json').read_text())
audit = json.loads((run/'provenance/input_audit.json').read_text())
records = list(csv.DictReader((run/'provenance/input_records.tsv').open(), delimiter='\t'))
issues = []

def read_json(path):
    try:
        return json.loads(path.read_text())
    except Exception as exc:
        issues.append(str(path)+': '+str(exc))
        return {}

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()

def dump(path, value):
    with path.open('x') as f:
        json.dump(value,f,indent=2); f.write('\n')

if (run/'summary/run_summary.json').exists():
    sys.exit('Summary already exists; preserve it and write any re-verification to a new directory.')
if sha(run/'input/all_significant_proteins.fasta') != audit['sha256']:
    issues.append('Input snapshot checksum mismatch')
try:
    original_still_matches = sha(pathlib.Path(meta['source'])) == audit['sha256']
except OSError:
    original_still_matches = None

completion=read_json(run/'provenance/workflow_completion.json')
if completion.get('exit_code') != 0:
    issues.append('Successful workflow exit is not established')
logpath=run/'logs/nextflow.log'
logtext=logpath.read_text(errors='replace') if logpath.exists() else ''
if 'InterProScan completed successfully' not in logtext or 'all processes finished' not in logtext:
    issues.append('Nextflow success and all-processes-finished log evidence missing')
tracepath=run/'logs/trace.tsv'
trace=list(csv.DictReader(tracepath.open(),delimiter='\t')) if tracepath.exists() else []
if not trace:
    issues.append('Workflow trace missing or empty')
bad_tasks=[r for r in trace if r.get('status') not in ('COMPLETED','CACHED') or r.get('exit') not in ('0','-',None,'')]
if bad_tasks:
    issues.append('One or more workflow tasks failed or remain unverified')
required=['ANTIFAM:PARSE_ANTIFAM','CATH:PARSE_CATHGENE3D','CATH:PARSE_FUNFAM','CDD:PARSE_CDD','COILS:PARSE_COILS','HAMAP:PARSE_PFSEARCH','MOBIDBLITE:PARSE_MOBIDBLITE','NCBIFAM:PARSE_NCBIFAM','PANTHER:PARSE_PANTHER','PFAM:PARSE_PFAM','PIRSF:PARSE_PIRSF','PIRSR:PARSE_PIRSR','PRINTS:PARSE_PRINTS','PROSITE_PATTERNS:PARSE_PSSCAN','PROSITE_PROFILES:PARSE_PFSEARCH','SFLD:PARSE_SFLD','SMART:PARSE_SMART','SUPERFAMILY:PARSE_SUPERFAMILY']
names=[r.get('name','') for r in trace]
missing_analyses=[name for name in required if not any(':'+name in n for n in names)]
if missing_analyses:
    issues.append('Missing analysis completion evidence: '+', '.join(missing_analyses))
if any('DOWNLOAD_' in n for n in names):
    issues.append('Unexpected database download task appeared')

doc=read_json(run/'results/all_significant_proteins.json')
if doc.get('interproscan-version') != meta['software_version'] or doc.get('interpro-version') != meta['interpro_version']:
    issues.append('Structured-output versions do not match pinned versions')
by_id=defaultdict(list)
for result in doc.get('results',[]):
    for xref in result.get('xref',[]):
        by_id[xref.get('id')].append(result)
expected_ids={r['original_id'] for r in records}
extra_ids=set(by_id)-expected_ids
if extra_ids:
    issues.append('Unexpected structured-result IDs: '+repr(sorted(extra_ids)))
nohit_md5=set()
for p in (run/'work').glob('*/*/no_matches.json'):
    nohit_md5.update(read_json(p))

tsvpath=run/'results/all_significant_proteins.tsv'
tsv_ids=set()
expected={r['original_id']:r for r in records}
if not tsvpath.exists():
    issues.append('Raw TSV missing')
else:
    with tsvpath.open() as f:
        for lineno,line in enumerate(f,1):
            cols=line.rstrip('\n').split('\t')
            if len(cols)<13:
                issues.append('Malformed TSV row '+str(lineno)); continue
            ident=cols[0]; tsv_ids.add(ident)
            r=expected.get(ident)
            if r is None or cols[1]!=r['sequence_md5'] or cols[2]!=r['length']:
                issues.append('TSV ID/sequence mismatch on row '+str(lineno))
            try:
                if not 1<=int(cols[6])<=int(cols[7])<=int(cols[2]):
                    issues.append('TSV coordinates out of range on row '+str(lineno))
            except ValueError:
                issues.append('Non-numeric TSV coordinates on row '+str(lineno))

rows=[]
functional_exclusions={'mobidb-lite','coils','phobius','deeptmhmm','signalp','tmbed','antifam'}
member_endpoints={'pfam':'pfam','panther':'panther','cdd':'cdd','cath-gene3d':'cathgene3d','cath-funfam':'cathfunfam','hamap':'hamap','ncbifam':'ncbifam','pirsf':'pirsf','prints':'prints','prosite patterns':'prosite','prosite profiles':'prosite','sfld':'sfld','smart':'smart','superfamily':'ssf'}
for original in records:
    ident=original['original_id']; evidence=[]; local_issues=[]
    candidates=by_id.get(ident,[])
    result=candidates[0] if len(candidates)==1 else None
    if result is None:
        local_issues.append('Missing or ambiguous structured result')
    else:
        seq=result.get('sequence','')
        if len(seq)!=int(original['length']) or hashlib.sha256(seq.encode()).hexdigest()!=original['sequence_sha256'] or result.get('md5')!=original['sequence_md5']:
            local_issues.append('Structured sequence does not match original')
    matches=result.get('matches',[]) if result is not None else []
    if result is not None and not isinstance(matches,list):
        local_issues.append('Invalid matches array'); matches=[]
    if matches and ident not in tsv_ids:
        local_issues.append('Annotated structured result missing from TSV')
    if not matches and ident in tsv_ids:
        local_issues.append('TSV hits disagree with structured result')
    if result is not None and not matches and original['sequence_md5'] not in nohit_md5:
        local_issues.append('Explicit no-match intermediate evidence missing')
    libs=set(); signatures=set(); entries=set(); gos=set(); links=set()
    for match in matches:
        sig=match.get('signature') or {}; release=sig.get('signatureLibraryRelease') or {}
        lib=release.get('library','unknown'); libs.add(lib)
        accession=sig.get('accession',''); desc=sig.get('description') or sig.get('name') or ''
        signatures.add(lib+'|'+accession+'|'+desc)
        entry=sig.get('entry') or {}
        if entry.get('accession'):
            entries.add(entry['accession']+'|'+(entry.get('name') or '')+'|'+(entry.get('description') or ''))
            links.add('https://www.ebi.ac.uk/interpro/entry/InterPro/'+quote(entry['accession'],safe='')+'/')
        if accession and lib.lower() in member_endpoints:
            links.add('https://www.ebi.ac.uk/interpro/entry/'+member_endpoints[lib.lower()]+'/'+quote(accession,safe='')+'/')
        gos.update(re.findall(r'GO:\d{7}',json.dumps(match)))
        for loc in match.get('locations',[]):
            if not isinstance(loc.get('start'),int) or not isinstance(loc.get('end'),int) or not 1<=loc['start']<=loc['end']<=int(original['length']):
                local_issues.append('Structured annotation coordinates out of range')
    links.update('https://www.ebi.ac.uk/QuickGO/term/'+g for g in gos)
    lowerlibs={s.lower() for s in libs}
    if not matches:
        category='no_matches'
    elif lowerlibs <= {'mobidb-lite'}:
        category='disorder_only'
    elif lowerlibs <= {'mobidb-lite','coils'}:
        category='disorder_or_coiled_coil_only'
    elif lowerlibs-functional_exclusions:
        category='family_domain_or_functional_signature'
    elif 'antifam' in lowerlibs:
        category='includes_possible_spurious_ORF_signature'
    else:
        category='other_prediction'
    status='failed_or_unverified' if issues or local_issues else ('completed_with_annotations' if matches else 'completed_no_hits')
    rows.append(dict(record_number=original['record_number'],original_id=ident,length=original['length'],processing_status=status,evidence_category=category,has_functional_or_domain_signature=bool(lowerlibs-functional_exclusions),has_antifam_match='antifam' in lowerlibs,match_count=len(matches),member_databases='; '.join(sorted(libs)),member_signatures_descriptions='; '.join(sorted(signatures)),integrated_interpro_entries='; '.join(sorted(entries)),go_terms='; '.join(sorted(gos)),identifier_links='; '.join(sorted(links)),sequence_md5=original['sequence_md5'],verification_notes='; '.join(local_issues)))
    if local_issues:
        issues.extend(ident+': '+x for x in local_issues)

# Any global integrity discrepancy prevents a successful overall declaration.
counts=Counter(r['processing_status'] for r in rows)
for extension,delimiter in [('tsv','\t'),('csv',',')]:
    with (run/('summary/per_protein_summary.'+extension)).open('x',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter=delimiter)
        writer.writeheader(); writer.writerows(rows)
job_ids={}
for role in ['analysis','postprocess']:
    p=run/('provenance/'+role+'_job_id.txt')
    if p.exists(): job_ids[role]=p.read_text().strip()
accounting={}
if job_ids.get('analysis'):
    proc=subprocess.run(['sacct','-j',job_ids['analysis'],'-X','-n','-P','--format=JobID,State,ExitCode,Elapsed,AllocCPUS,ReqMem,NodeList'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    accounting=dict(returncode=proc.returncode,stdout=proc.stdout,stderr=proc.stderr)
    (run/'provenance/analysis_accounting_at_postprocess.txt').write_text(proc.stdout+proc.stderr)
    main_rows=[s.split('|') for s in proc.stdout.splitlines() if s.split('|')[0]==job_ids['analysis']]
    if main_rows and (main_rows[0][1]!='COMPLETED' or main_rows[0][2]!='0:0'):
        issues.append('SLURM accounting does not report analysis COMPLETED 0:0')
hashlines=[]
for folder in ['input','results','summary','scripts']:
    for p in sorted((run/folder).rglob('*')):
        if p.is_file(): hashlines.append(sha(p)+'  '+str(p.relative_to(run)))
(run/'provenance/results_and_deliverables.sha256').write_text('\n'.join(hashlines)+'\n')
summary=dict(verified_complete=not issues,input_records=len(records),processed_records=counts['completed_with_annotations']+counts['completed_no_hits'],failed_or_unverified=counts['failed_or_unverified'],completed_with_annotations=counts['completed_with_annotations'],completed_no_hits=counts['completed_no_hits'],proteins_with_any_match=sum(r['match_count']>0 for r in rows),proteins_with_integrated_interpro=sum(bool(r['integrated_interpro_entries']) for r in rows),proteins_with_go_terms=sum(bool(r['go_terms']) for r in rows),evidence_categories=dict(Counter(r['evidence_category'] for r in rows)),input_sha256=audit['sha256'],original_source_still_matches_snapshot=original_still_matches,workflow_tasks=len(trace),failed_workflow_tasks=len(bad_tasks),missing_analysis_completion=missing_analyses,job_ids=job_ids,software=meta,verified_at_utc=datetime.now(timezone.utc).isoformat(),issues=issues,results_dir=str(run/'results'),per_protein_tsv=str(run/'summary/per_protein_summary.tsv'),per_protein_csv=str(run/'summary/per_protein_summary.csv'),interpretation='Matches are annotation evidence, not proof of virulence. Disorder/coiled-coil-only results are distinguished. No new findings relative to prior runs are claimed.')
dump(run/'summary/run_summary.json',summary)
text=[('VERIFIED COMPLETE' if not issues else 'INCOMPLETE / VERIFICATION ISSUES'),'',f"Run: {run}",f"Input: {meta['source']}",f"Input SHA-256: {audit['sha256']}",f"Input records: {len(records)}",f"Processed: {summary['processed_records']}; failed/unverified: {summary['failed_or_unverified']}",f"With annotations: {summary['completed_with_annotations']}; completed no hits: {summary['completed_no_hits']}",f"With integrated InterPro entries: {summary['proteins_with_integrated_interpro']}; with GO terms: {summary['proteins_with_go_terms']}",f"Evidence categories: {summary['evidence_categories']}",f"Software: InterProScan {meta['software_version']}; InterPro {meta['interpro_version']}; Nextflow {meta['nextflow_version']}; SingularityCE {meta['singularity_version']}",f"Enabled: {', '.join(meta['enabled_analyses'])}",'Skipped: '+json.dumps(meta['skipped_analyses']),f"Jobs: {job_ids}",f"Workflow tasks: {len(trace)}; unsuccessful: {len(bad_tasks)}",'Commands: provenance/analysis_submission_command.txt, provenance/postprocess_submission_command.txt, provenance/nextflow_command.txt','Database versions: provenance/database_versions.json','Raw outputs: results/all_significant_proteins.tsv, results/all_significant_proteins.json, results/all_significant_proteins.xml','Per-protein summaries: summary/per_protein_summary.tsv and .csv','Checksums: provenance/input.sha256, pipeline.sha256, container_hashes.sha256, results_and_deliverables.sha256','Full provenance, scripts, logs, and isolated Nextflow work are retained.',summary['interpretation'],'Original source still matches snapshot: '+str(original_still_matches),'Unresolved issues: '+('; '.join(issues) if issues else 'None')]
(run/'summary/RUN_SUMMARY.txt').write_text('\n'.join(text)+'\n')
print('\n'.join(text))
sys.exit(0 if not issues else 2)
