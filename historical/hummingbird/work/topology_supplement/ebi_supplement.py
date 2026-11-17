#!/usr/bin/env python3
"""Submit once, inspect status, and retrieve the explicitly authorized EBI supplement."""
import datetime, hashlib, json, pathlib, re, sys, urllib.error, urllib.parse, urllib.request, xml.etree.ElementTree as ET

ROOT=pathlib.Path(__file__).resolve().parent
REMOTE=(ROOT.parent/'provenance/supplement_setup.json').is_file()
if REMOTE:ROOT=ROOT.parent
API='https://www.ebi.ac.uk/Tools/services/rest/iprscan6'
RECEIPT=ROOT/('provenance/ebi_receipt.json' if REMOTE else 'ebi_receipt.json')
FASTA=(ROOT/'input/all_significant_proteins.fasta') if REMOTE else ROOT.parents[1]/'outputs/Hummingbird_results_20261004/input/all_significant_proteins.fasta'
EXPECTED='3285311fbbe7d7088ff63a24ea515e680cabd1e0d3519803d47b24deae04aedb'

def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(value):
    temp=RECEIPT.with_suffix('.tmp');temp.write_text(json.dumps(value,indent=2));temp.replace(RECEIPT)
def get(path):
    with urllib.request.urlopen(API+'/'+path,timeout=90) as r:return r.read()

if sys.argv[1]=='submit':
    assert not RECEIPT.exists(), 'Submission already attempted; inspect receipt, never duplicate blindly.'
    seq=FASTA.read_bytes();assert hashlib.sha256(seq).hexdigest()==EXPECTED
    params=[('email','researcher@example.org'),('title','StrepA 66 supplemental SignalP-Prok and InterPro-N 2026-10-06'),
            ('stype','p'),('appl','SignalP-Prok'),('appl','InterPro-N'),('goterms','false'),('pathways','false'),('sequence',seq.decode())]
    receipt={'authorized_by_user':True,'requested_at_utc':now(),'endpoint':API+'/run','applications':['SignalP-Prok','InterPro-N'],
             'input_sha256':EXPECTED,'input_records':66,'email':'researcher@example.org','state':'submission_started_no_acknowledgement'}
    save(receipt)
    try:
        request=urllib.request.Request(API+'/run',data=urllib.parse.urlencode(params).encode(),headers={'Content-Type':'application/x-www-form-urlencoded'},method='POST')
        with urllib.request.urlopen(request,timeout=120) as response: job=response.read().decode().strip()
        assert re.fullmatch(r'iprscan6-[A-Za-z0-9-]+',job), 'Unexpected submission response: '+job
        receipt.update(job_id=job,state='submitted',acknowledged_at_utc=now());save(receipt)
    except urllib.error.HTTPError as e:
        receipt.update(state='http_submission_error',http_status=e.code,response=e.read().decode(errors='replace'));save(receipt);raise
    print(json.dumps(receipt,indent=2))
elif sys.argv[1] in ('status','retrieve'):
    receipt=json.loads(RECEIPT.read_text());job=receipt['job_id'];status=get('status/'+job).decode().strip()
    receipt.update(last_status=status,last_checked_utc=now());save(receipt)
    if sys.argv[1]=='retrieve' and status=='FINISHED':
        types=get('resulttypes/'+job); out=ROOT/('results/ebi' if REMOTE else 'ebi_results');out.mkdir(exist_ok=True)
        (out/'resulttypes.xml').write_bytes(types)
        identifiers=[x.text for x in ET.fromstring(types).findall('.//identifier')]
        receipt['result_types']=identifiers; files={}
        for typ in identifiers:
            # Retain structured results, hits, input and complete execution output.
            if typ in ['json','xml','tsv','out','log','error','sequence','sequenceInput','input','jsonl','gff'] or any(s in typ.lower() for s in ['json','xml','tsv']):
                data=get('result/'+job+'/'+typ);path=out/(typ+'.txt' if '.' not in typ else typ)
                path.write_bytes(data);files[path.name]={'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
        receipt['retrieved_files']=files;save(receipt)
    print(json.dumps(receipt,indent=2))
else:raise ValueError('Expected submit, status or retrieve')
