import ast,csv,hashlib,json,pathlib,shutil,subprocess,sys,xml.etree.ElementTree as ET
here=pathlib.Path(__file__).resolve().parent
script=here/'postprocess.py'
tree=ast.parse(script.read_text())
required=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='required' for t in n.targets))
base=here/'verification-tests'
base.mkdir(exist_ok=True)
checks=[]
for case in ['complete','missing_structured_record','failed_task','missing_nohit_evidence']:
    run=base/case
    run.mkdir()
    for folder in ['scripts','provenance','input','summary','results','logs','work/aa/bb']:
        (run/folder).mkdir(parents=True,exist_ok=True)
    shutil.copy2(script,run/'scripts/postprocess.py')
    proteins=[('protein_hit','MKTAA'),('protein_nohit','MKVGG')]
    fasta=run/'input/all_significant_proteins.fasta'
    fasta.write_text(''.join('>'+i+'\n'+s+'\n' for i,s in proteins))
    digest=hashlib.sha256(fasta.read_bytes()).hexdigest()
    original=[]
    for n,(ident,seq) in enumerate(proteins,1):
        original.append(dict(record_number=n,original_id=ident,length=len(seq),sequence_md5=hashlib.md5(seq.encode()).hexdigest(),sequence_sha256=hashlib.sha256(seq.encode()).hexdigest()))
    with (run/'provenance/input_records.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(original[0]),delimiter='\t');w.writeheader();w.writerows(original)
    metadata=dict(source=str(fasta),input_sha256=digest,software_version='6.0.1',interpro_version='109.0',nextflow_version='25.10.4',singularity_version='4.1.4',enabled_analyses=['Pfam'],skipped_analyses={})
    (run/'provenance/run_metadata.json').write_text(json.dumps(metadata))
    (run/'provenance/input_audit.json').write_text(json.dumps(dict(sha256=digest)))
    (run/'provenance/workflow_completion.json').write_text(json.dumps(dict(exit_code=0)))
    (run/'logs/nextflow.log').write_text('InterProScan completed successfully\nall processes finished\n')
    with (run/'logs/trace.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['name','status','exit'],delimiter='\t');w.writeheader()
        for i,name in enumerate(required):
            fail=case=='failed_task' and i==0
            w.writerow(dict(name='INTERPROSCAN:SCAN_SEQUENCES:'+name+' (1)',status='FAILED' if fail else 'COMPLETED',exit='1' if fail else '0'))
    outputs=[]
    for idx,(ident,seq) in enumerate(proteins):
        matches=[] if idx else [dict(signature=dict(accession='PF00001',name='Test domain',description='Fixture only',signatureLibraryRelease=dict(library='Pfam',version='38.2'),entry=dict(accession='IPR000001',name='Fixture entry',description='Test only',goXRefs=[dict(id='GO:0008150')])),locations=[dict(start=1,end=4)])]
        outputs.append(dict(sequence=seq,md5=original[idx]['sequence_md5'],xref=[dict(id=ident,name=ident)],matches=matches))
    doc={'interproscan-version':'6.0.1','interpro-version':'109.0','results':outputs[:1] if case=='missing_structured_record' else outputs}
    (run/'results/all_significant_proteins.json').write_text(json.dumps(doc))
    xml=ET.Element('results',{'interproscan-version':'6.0.1','interpro-version':'109.0'})
    for idx,(ident,seq) in enumerate(proteins):
        p=ET.SubElement(xml,'protein');ET.SubElement(p,'sequence',{'md5':original[idx]['sequence_md5']}).text=seq;ET.SubElement(p,'xref',{'id':ident})
    ET.ElementTree(xml).write(run/'results/all_significant_proteins.xml')
    (run/'results/all_significant_proteins.tsv').write_text('\t'.join(['protein_hit',original[0]['sequence_md5'],'5','Pfam','PF00001','Test domain','1','4','1e-5','T','03-10-2026','IPR000001','Test entry'])+'\n')
    if case!='missing_nohit_evidence':
        (run/'work/aa/bb/no_matches.json').write_text(json.dumps({original[1]['sequence_md5']:{}}))
    proc=subprocess.run([sys.executable,str(run/'scripts/postprocess.py')],stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    if not (run/'summary/run_summary.json').exists():
        raise AssertionError(proc.stderr)
    summary=json.loads((run/'summary/run_summary.json').read_text())
    rows=list(csv.DictReader((run/'summary/per_protein_summary.tsv').open(),delimiter='\t'))
    if case=='complete':
        assert proc.returncode==0 and summary['verified_complete'],proc.stdout+proc.stderr
        assert [r['processing_status'] for r in rows]==['completed_with_annotations','completed_no_hits']
        assert summary['processed_records']==2 and summary['proteins_with_integrated_interpro']==1
    else:
        assert proc.returncode==2 and not summary['verified_complete']
        assert rows[1]['processing_status']=='failed_or_unverified'
    checks.append(dict(case=case,passed=True,returncode=proc.returncode,issues=summary['issues']))
path=here/'postprocessor_test_results.json'
path.write_text(json.dumps(checks,indent=2)+'\n')
print(path.read_text())
