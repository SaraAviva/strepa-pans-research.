#!/usr/bin/env python3
"""Run isolated topology supplements on a SLURM compute node and account for every input."""
import argparse, csv, datetime, hashlib, json, os, pathlib, re, subprocess, tarfile, urllib.request

EXPECTED_SHA = '3285311fbbe7d7088ff63a24ea515e680cabd1e0d3519803d47b24deae04aedb'
TMBED_IMAGE = 'docker://interpro/tmbed@sha256:1031a299e650a8d32e3e2a3e4417d63d7e1ae5f752100688dde58283902f32b4'

def sha(path):
    h = hashlib.sha256()
    with pathlib.Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''): h.update(block)
    return h.hexdigest()

def fasta(path):
    records = {}; ident = None
    for line in pathlib.Path(path).read_text().splitlines():
        if line.startswith('>'):
            ident = line[1:].split()[0]
            assert ident not in records, 'Duplicate ID'
            records[ident] = ''
        elif line.strip():
            assert ident is not None
            records[ident] += line.strip()
    return records

def exclusive_json(path, value):
    with pathlib.Path(path).open('x') as f:
        json.dump(value, f, indent=2); f.write('\n')

def execute(argv, **kwargs):
    print('COMMAND '+json.dumps([str(x) for x in argv]), flush=True)
    return subprocess.run([str(x) for x in argv], check=True, **kwargs)

def verify_phobius(path, seqs):
    records = {}; ident = None
    for line in path.read_text().splitlines():
        if line.startswith('ID '):
            assert ident is None
            ident = line.split()[1]
            assert ident in seqs and ident not in records
            records[ident] = []
        elif line.startswith('FT '):
            assert ident is not None
            fields = line.split(None, 4)
            kind, start, end = fields[1], int(fields[2]), int(fields[3])
            assert 1 <= start <= end <= len(seqs[ident]), (ident, line)
            records[ident].append((kind, start, end, fields[4] if len(fields)>4 else ''))
        elif line.strip() == '//':
            assert ident is not None
            ident = None
        elif line.strip(): raise AssertionError('Unexpected Phobius output: '+line)
    assert ident is None and set(records) == set(seqs), 'Incomplete Phobius output'
    return records

def verify_tmbed(path, seqs):
    lines = [x.strip() for x in path.read_text().splitlines() if x.strip()]
    assert len(lines) == 3 * len(seqs), 'Incomplete TMbed output'
    records = {}; raw_labels = {}
    classes = {'H':'TM helix inside-to-outside','h':'TM helix outside-to-inside',
               'B':'TM beta strand inside-to-outside','b':'TM beta strand outside-to-inside',
               'S':'Signal peptide','i':'Non-TM inside','o':'Non-TM outside'}
    for i in range(0, len(lines), 3):
        assert lines[i].startswith('>')
        ident = lines[i][1:].split()[0]; sequence, labels = lines[i+1:i+3]
        assert ident in seqs and ident not in records
        assert sequence == seqs[ident], 'Sequence changed: '+ident
        assert len(labels) == len(sequence) and set(labels) <= set(classes), (ident, set(labels))
        regions = []; start = 0
        for end in range(1, len(labels)+1):
            if end == len(labels) or labels[end] != labels[start]:
                regions.append((labels[start], start+1, end, classes[labels[start]])); start = end
        records[ident] = regions; raw_labels[ident] = labels
    assert set(records) == set(seqs)
    return records, raw_labels

def run(run, component):
    assert os.environ.get('SLURM_JOB_ID'), 'Compute-node SLURM allocation required'
    src = run/'input/all_significant_proteins.fasta'
    assert sha(src) == EXPECTED_SHA
    seqs = fasta(src); assert len(seqs) == 66 and all(seqs.values())
    output = run/'results'/component; output.mkdir()
    steps = []; labels = {}
    if component == 'phobius':
        package = run/'runtime/phobius101_linux.tgz'
        url = 'https://software.sbc.su.se/files/phobius101_linux.tgz'
        if not package.exists():
            with urllib.request.urlopen(url, timeout=90) as response, package.open('xb') as f:
                f.write(response.read())
        destination = run/'runtime/phobius_package'
        if not destination.exists():
            destination.mkdir()
            with tarfile.open(str(package)) as archive:
                for member in archive.getmembers():
                    target = (destination/member.name).resolve()
                    assert str(target).startswith(str(destination.resolve())+'/')
                    assert not member.issym() and not member.islnk(), 'Unexpected archive link'
                archive.extractall(str(destination))
        scripts = list(destination.rglob('phobius.pl')); assert len(scripts) == 1
        script = scripts[0]
        # Same 64-bit executable selection used by the installed InterProScan module.
        original = script.read_bytes()
        updated = original.replace(b'/decodeanhmm";', b'/decodeanhmm.64bit";')
        script.write_bytes(updated)
        raw = output/'phobius.out'
        execute(['perl','-MFindBin','-MCwd','-MFile::Basename','-e','print "Perl dependencies verified\\n"'])
        with raw.open('x') as f: execute(['perl', script, src], stdout=f, cwd=str(output))
        records = verify_phobius(raw, seqs)
        steps = [{'source':url,'package_sha256':sha(package),'version':'1.01',
                  'download_terms':'Private study, education or nonprofit research; UCSC research project',
                  'wrapper_64bit_selection_changed':original != updated}]
    elif component == 'tmbed':
        image = run/'runtime/tmbed-1.0.2.sif'
        execute(['singularity','pull','--disable-cache',image,TMBED_IMAGE])
        listing = subprocess.check_output(['singularity','sif','list',str(image)], universal_newlines=True)
        parts = [line.split('|') for line in listing.splitlines() if 'FS (Squashfs/*System/amd64)' in line]
        assert len(parts) == 1, listing
        offset = int(parts[0][3].strip().split('-')[0])
        old = pathlib.Path(json.loads((run/'provenance/supplement_setup.json').read_text())['previous_verified_run'])
        extraction = old/'recovery_20261004T163308Z/extraction_tools'
        sandbox = run/'runtime/tmbed-sandbox'
        env = dict(os.environ, LD_LIBRARY_PATH=str(extraction/'lib'))
        execute([extraction/'unsquashfs','-no-progress','-no-xattrs','-processors','2',
                 '-data-queue','64','-frag-queue','64','-offset',str(offset),'-dest',sandbox,image], env=env)
        raw = output/'tmbed.pred'
        command = ['singularity','exec','--cleanenv','--no-home','-B','/mnt/beegfs','-B','/hb',
                   sandbox,'tmbed','predict','-f',src,'-m','/opt/tmbed/models',
                   '-t',os.environ['SLURM_CPUS_PER_TASK'],'-p',raw,'--out-format','4',
                   '--batch-size','1024','--no-use-gpu']
        execute(command, cwd=str(output))
        records, labels = verify_tmbed(raw, seqs)
        steps = [{'image':TMBED_IMAGE,'sif_sha256':sha(image),'version':'1.0.2',
                  'output_format':4,'note':'Direct complete residue labels retain inside/outside topology and avoid the InterProScan 6.0.1 terminal-coordinate parser defect.'}]
    else: raise ValueError(component)
    rows = []; region_rows = []
    for ident, seq in seqs.items():
        features = records[ident]
        if component == 'phobius':
            tm = [r for r in features if r[0]=='TRANSMEM']
            signal = [r for r in features if r[0]=='SIGNAL']
        else:
            tm = [r for r in features if r[0] in 'hHbB']
            signal = [r for r in features if r[0]=='S']
        rows.append({'protein_id':ident,'length':len(seq),'status':'verified_complete',
                     'predictor':component,'tm_segments':len(tm),'tm_coordinates':'; '.join('%d-%d'%(r[1],r[2]) for r in tm),
                     'signal_peptide':bool(signal),'signal_coordinates':'; '.join('%d-%d'%(r[1],r[2]) for r in signal),
                     'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),'residue_labels':labels.get(ident,''),
                     'interpretation':'Sequence prediction; not experimental localization. No TM/signal does not establish a cytoplasmic location.'})
        for kind,start,end,description in features:
            region_rows.append({'protein_id':ident,'predictor':component,'feature':kind,'start':start,'end':end,'description':description})
    for name, values in [('per_protein.tsv',rows),('regions.tsv',region_rows)]:
        with (output/name).open('x',newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(values[0]), delimiter='\t');writer.writeheader();writer.writerows(values)
    verified = {'verified_complete':True,'analysis':component,'input_records':66,'processed_records':len(records),
                'input_sha256':EXPECTED_SHA,'slurm_job_id':os.environ['SLURM_JOB_ID'],
                'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
                'proteins_with_tm':sum(r['tm_segments']>0 for r in rows),'proteins_with_signal_peptide':sum(r['signal_peptide'] for r in rows),
                'provenance':steps,'output_sha256':{p.name:sha(p) for p in output.iterdir() if p.is_file()}}
    exclusive_json(output/'verification.json',verified)
    print(json.dumps(verified,indent=2),flush=True)

if __name__ == '__main__':
    ap=argparse.ArgumentParser();ap.add_argument('run',type=pathlib.Path);ap.add_argument('component',choices=['phobius','tmbed']);args=ap.parse_args()
    run(args.run,args.component)
