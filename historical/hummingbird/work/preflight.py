#!/usr/bin/env python3
import hashlib
import json
import os
import pathlib
import sys

run=pathlib.Path(__file__).resolve().parent.parent
meta=json.loads((run/'provenance/run_metadata.json').read_text())
db=pathlib.Path(meta['database_root'])
snapshot=run/'input/all_significant_proteins.fasta'
assert hashlib.sha256(snapshot.read_bytes()).hexdigest()==meta['input_sha256']
required={
 'interpro/109.0':['databases.json','entries.json','goterms.json','goterms.ipr.json','pathways.json','pathways.ipr.json'],
 'antifam/8.0':['AntiFam.hmm'],
 'cath/4.3.0/gene3d':['gene3d_main.hmm','model_to_family_map.tsv','discontinuous_regs.pkl.py3'],
 'cath/4.3.0/funfam':[''],
 'cdd/3.21':['db','data'],
 'hamap/2026_01':['profiles'],
 'ncbifam/19.0':['ncbifam.hmm'],
 'panther/19.0':['famhmm/binHmm','Tree_MSF','PAINT_Annotations'],
 'pfam/38.2':['pfam_a.hmm','pfam_a.dat'],
 'pirsf/3.10':['sf_hmm_all','pirsf.dat'],
 'pirsr/2025_05':['sr_hmm_all','sr_uru.json'],
 'prints/42.0':['prints.pval','FingerPRINTShierarchy.db'],
 'prosite/2026_01':['prosite_patterns.dat','evaluator.dat','prosite_profiles','skip_flagged_profiles.txt'],
 'sfld/4':['sfld.hmm','sfld_hierarchy.tsv','sfld_sites.annot'],
 'smart/9.0':['hmmer3/smart.HMMs','hmmer2'],
 'superfamily/1.75':['hmmlib_1.75','self_hits.tab','dir.cla.scop.txt_1.75','model.tab','pdbj95d']}
inventory=[]
for directory,children in required.items():
    for child in children:
        p=db/directory/child
        assert p.exists() and os.access(str(p),os.R_OK), 'Missing/unreadable database path: '+str(p)
        if p.is_file():
            with p.open('rb') as f: f.read(1)
        inventory.append(dict(path=str(p),kind='directory' if p.is_dir() else 'file',bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns))
for image in (run/'provenance/container_paths.txt').read_text().splitlines():
    p=pathlib.Path(image)
    assert p.is_file() and os.access(image,os.R_OK) and p.stat().st_size>0, image
assert os.access(str(run/'work'),os.W_OK)
out=run/'provenance'/('database_preflight_'+os.environ.get('SLURM_JOB_ID','login')+'.json')
with out.open('x') as f: json.dump(inventory,f,indent=2)
print('Input checksum, existing database files, containers, and work permissions verified.')
