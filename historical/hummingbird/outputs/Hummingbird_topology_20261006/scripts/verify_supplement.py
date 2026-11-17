#!/usr/bin/env python3
"""Independently reconcile supplemental raw results, sequence identity and job evidence."""
import argparse, collections, csv, datetime, hashlib, json, pathlib, subprocess, sys, traceback
import xml.etree.ElementTree as ET
from run_components_v4 import EXPECTED_SHA, fasta, sha, verify_phobius, verify_tmbed

def table(path):
    rows = list(csv.DictReader(path.open(), delimiter='\t'))
    result = {r['protein_id']:r for r in rows}
    assert len(result) == len(rows), 'Duplicate summary IDs: '+str(path)
    return result

def scheduler(job):
    text = subprocess.check_output(['sacct','-nX','-j',str(job),'-o','JobIDRaw,State,ExitCode','-P'], text=True)
    matches = [line.split('|') for line in text.splitlines() if line.split('|')[0] == str(job)]
    assert len(matches) == 1, 'Missing accounting for '+str(job)
    assert matches[0][1:3] == ['COMPLETED','0:0'], 'Job did not complete successfully: '+str(matches)
    return dict(zip(['job_id','state','exit_code'], matches[0][:3]))

def verify_ebi(run, seqs):
    out = run/'results/ebi'
    receipt = json.loads((run/'provenance/ebi_receipt.json').read_text())
    assert receipt['last_status'] == 'FINISHED' and receipt['input_sha256'] == EXPECTED_SHA
    assert set(receipt['applications']) == {'SignalP-Prok','InterPro-N'}
    for name, record in receipt['retrieved_files'].items():
        assert sha(out/name) == record['sha256'], 'EBI file hash mismatch: '+name
    log = (out/'log.txt').read_text()
    assert 'InterProScan completed successfully' in log
    assert 'RUN_SIGNALP_CPU_PROK' in log and 'RUN_INTERPRO_N_CPU' in log
    assert fasta(out/'sequence.txt') == seqs
    doc = json.loads((out/'json.txt').read_text())
    jseen = {}; canonical = {}
    for p in doc['results']:
        assert len(p['xref']) == 1
        ident = p['xref'][0]['id']; assert ident in seqs and ident not in jseen
        assert p['sequence'] == seqs[ident]
        assert p['md5'].lower() == hashlib.md5(seqs[ident].encode()).hexdigest()
        values = []
        for m in p['matches']:
            assert m['source'] in ['SignalP','InterPro-N']
            for loc in m['locations']:
                assert 1 <= loc['start'] <= loc['end'] <= len(seqs[ident])
                values.append((m['source'], m['signature']['accession'], loc['start'], loc['end']))
        canonical[ident] = sorted(values); jseen[ident] = p
    assert set(jseen) == set(seqs) and len(jseen) == 66
    xseen = set()
    for p in ET.parse(out/'xml.txt').getroot().findall('protein'):
        xrefs = p.findall('xref'); assert len(xrefs) == 1
        ident = xrefs[0].get('id'); assert ident in seqs and ident not in xseen
        xseen.add(ident); sequence = p.find('sequence')
        assert sequence.text == seqs[ident]
        assert sequence.get('md5').lower() == hashlib.md5(seqs[ident].encode()).hexdigest()
        values = []
        for m in p.findall('./matches/match'):
            for loc in m.findall('./locations/location'):
                values.append((m.get('source'),m.find('signature').get('ac'),int(loc.get('start')),int(loc.get('end'))))
        assert sorted(values) == canonical[ident], 'JSON/XML coordinate disagreement: '+ident
    assert xseen == set(seqs)
    tsv = collections.defaultdict(list)
    for row in csv.reader((out/'tsv.txt').open(), delimiter='\t'):
        assert len(row) == 15 and row[0] in seqs
        assert row[1].lower() == hashlib.md5(seqs[row[0]].encode()).hexdigest()
        assert int(row[2]) == len(seqs[row[0]])
        tsv[row[0]].append((row[9],row[4],int(row[6]),int(row[7])))
    for ident in seqs: assert sorted(tsv[ident]) == canonical[ident], 'TSV/JSON mismatch: '+ident
    rows = table(out/'per_protein.tsv'); assert set(rows) == set(seqs)
    for ident, row in rows.items():
        counts = collections.Counter(m['source'] for m in jseen[ident]['matches'])
        assert int(row['signalp_hits']) == counts['SignalP']
        assert int(row['interpro_n_matches']) == counts['InterPro-N']
        assert int(row['length']) == len(seqs[ident])
        assert row['sequence_sha256'] == hashlib.sha256(seqs[ident].encode()).hexdigest()
    return rows, {'job_id':receipt['job_id'],'records':len(rows),'interproscan_version':doc['interproscan-version'],
                  'interpro_version':doc['interpro-version'],'json_xml_tsv_agree':True,
                  'signalp_positive_proteins':sum(int(r['signalp_hits'])>0 for r in rows.values()),
                  'interpro_n_matches':sum(int(r['interpro_n_matches']) for r in rows.values())}

def verify_component(run, name, seqs, check_scheduler):
    out = run/'results'/name
    verification = json.loads((out/'verification.json').read_text())
    assert verification['verified_complete'] and verification['input_sha256'] == EXPECTED_SHA
    assert verification['processed_records'] == verification['input_records'] == 66
    for filename, digest in verification['output_sha256'].items(): assert sha(out/filename) == digest
    if name == 'phobius': features = verify_phobius(out/'phobius.out', seqs)
    else: features, _ = verify_tmbed(out/'tmbed.pred', seqs)
    rows = table(out/'per_protein.tsv'); assert set(rows) == set(seqs)
    for ident, row in rows.items():
        assert int(row['length']) == len(seqs[ident]) and row['status'] == 'verified_complete'
        assert row['sequence_sha256'] == hashlib.sha256(seqs[ident].encode()).hexdigest()
        tm = [x for x in features[ident] if x[0] == 'TRANSMEM'] if name == 'phobius' else [x for x in features[ident] if x[0] in 'HhBb']
        sig = [x for x in features[ident] if x[0] == ('SIGNAL' if name == 'phobius' else 'S')]
        assert int(row['tm_segments']) == len(tm)
        assert row['tm_coordinates'] == '; '.join('%d-%d'%(x[1],x[2]) for x in tm)
        assert row['signal_peptide'] == str(bool(sig))
        assert row['signal_coordinates'] == '; '.join('%d-%d'%(x[1],x[2]) for x in sig)
    if check_scheduler: verification['accounting'] = scheduler(verification['slurm_job_id'])
    return rows, verification

def verify_run(run, available_only=False):
    src = run/'input/all_significant_proteins.fasta'
    assert sha(src) == EXPECTED_SHA
    seqs = fasta(src); assert len(seqs) == 66
    ph, pv = verify_component(run,'phobius',seqs,not available_only)
    eb, ev = verify_ebi(run,seqs)
    tm, tv = ({}, {'status':'NOT_CHECKED'}) if available_only else verify_component(run,'tmbed',seqs,True)
    rows = []
    for ident, seq in seqs.items():
        row = {'protein_id':ident,'length':len(seq),'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest()}
        for prefix, values in [('phobius',ph),('tmbed',tm)]:
            if values:
                for field in ['tm_segments','tm_coordinates','signal_peptide','signal_coordinates']:
                    row[prefix+'_'+field] = values[ident][field]
        for field in ['signalp_hits','signalp_types','signalp_regions','signalp_scores','interpro_n_matches']:
            row[field] = eb[ident][field]
        rows.append(row)
    summary = {'verified_complete':not available_only,'status':'AVAILABLE_COMPONENTS_VERIFIED' if available_only else 'FOUR_ANALYSES_VERIFIED',
               'input_records':66,'accounted_records':len(rows),'input_sha256':EXPECTED_SHA,'run_directory':str(run),
               'verified_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
               'phobius':pv,'ebi':ev,'tmbed':tv,
               'unavailable':{'DeepTMHMM':'Separate licensed standalone package unavailable; not run'},
               'not_applicable':{'SignalP-Euk':'Bacterial dataset; SignalP-Prok used'},
               'interpretation':'Sequence predictions do not establish experimental localization, secretion or virulence.'}
    return summary, rows

if __name__ == '__main__':
    ap=argparse.ArgumentParser();ap.add_argument('run',type=pathlib.Path);ap.add_argument('--available-only',action='store_true');args=ap.parse_args()
    try:
        summary, rows=verify_run(args.run,args.available_only)
        if not args.available_only:
            out=args.run/'summary';out.mkdir(exist_ok=True)
            assert not (out/'run_summary.json').exists(), 'Existing final summary; review before overwriting'
            for extension, delimiter in [('tsv','\t'),('csv',',')]:
                with (out/('per_protein_summary.'+extension)).open('x',newline='') as f:
                    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter=delimiter);w.writeheader();w.writerows(rows)
            summary['output_sha256']={str(p.relative_to(args.run)):sha(p) for folder in ['results/phobius','results/tmbed','results/ebi','summary'] for p in (args.run/folder).iterdir() if p.is_file()}
            with (out/'run_summary.json').open('x') as f:json.dump(summary,f,indent=2)
        print(json.dumps(summary,indent=2))
    except Exception as error:
        out=args.run/'provenance'/('verification_failure_'+str(__import__('os').environ.get('SLURM_JOB_ID','local'))+'.json')
        out.write_text(json.dumps({'verified_complete':False,'error':str(error),'traceback':traceback.format_exc(),'run':str(args.run)},indent=2))
        traceback.print_exc();sys.exit(2)
