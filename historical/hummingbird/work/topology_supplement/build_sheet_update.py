import pathlib,csv,json,hashlib,collections,datetime

ROOT=pathlib.Path(__file__).resolve().parents[2]
OUT=ROOT/'work/topology_supplement/sheet';OUT.mkdir(exist_ok=True)
DATA=ROOT/'outputs/Hummingbird_topology_20261006'
STATE=json.loads((ROOT/'work/topology_supplement/state.json').read_text())
ORIGINAL=ROOT/'outputs/Hummingbird_results_20261004'
records=list(csv.DictReader((ORIGINAL/'provenance/input_records.tsv').open(),delimiter='\t'))
ph={r['protein_id']:r for r in csv.DictReader((DATA/'phobius/per_protein.tsv').open(),delimiter='\t')}
eb={r['protein_id']:r for r in csv.DictReader((DATA/'ebi/per_protein.tsv').open(),delimiter='\t')}
ev=json.loads((DATA/'ebi/verification.json').read_text());pv=json.loads((DATA/'phobius/verification.json').read_text())
assert len(records)==len(ph)==len(eb)==66 and ev['verified_complete'] and pv['verified_complete']
structured=json.loads((DATA/'ebi/json.txt').read_text())
versions={(m['source'],m['signature']['signatureLibraryRelease']['library'],m['signature']['accession']):m['signature']['signatureLibraryRelease']['version'] for r in structured['results'] for m in r['matches']}
sequence={x['id']:r['sequence'] for r in structured['results'] for x in r['xref']}
tm_path=DATA/'tmbed/per_protein.tsv'
tm={r['protein_id']:r for r in csv.DictReader(tm_path.open(),delimiter='\t')} if tm_path.exists() else {}
if tm:assert len(tm)==66 and json.loads((DATA/'tmbed/verification.json').read_text())['verified_complete']
headers={};data={};titles={}
headers['Topology & signals']=['Protein ID','Length (aa)','Phobius TM count','Phobius TM coordinates','Phobius signal peptide','Phobius signal coordinates','SignalP-Prok type','SignalP coordinates','SignalP probability','Phobius vs SignalP signal prediction','TMbed status','TMbed TM coordinates','TMbed signal coordinates','InterPro-N match count','Verified analyses','Your follow-up notes']
data['Topology & signals']=[]
for record in records:
 ident=record['original_id'];p=ph[ident];e=eb[ident];assert p['sequence_sha256']==e['sequence_sha256']==record['sequence_sha256']
 ps=p['signal_peptide']=='True';es=int(e['signalp_hits'])>0
 comparison='Both predict a signal peptide' if ps and es else 'Neither predicts a signal peptide' if not ps and not es else 'REVIEW: Phobius only' if ps else 'REVIEW: SignalP-Prok only'
 t=tm.get(ident)
 data['Topology & signals'].append([ident,int(p['length']),int(p['tm_segments']),p['tm_coordinates'] or 'None predicted','Predicted' if ps else 'None predicted',p['signal_coordinates'] or 'None predicted',e['signalp_types'] or 'None predicted',e['signalp_regions'] or 'None predicted',e['signalp_scores'] or '—',comparison,'Verified complete' if t else 'Running — job 796204',t['tm_coordinates'] or 'None predicted' if t else 'Pending',t['signal_coordinates'] or 'None predicted' if t else 'Pending',int(e['interpro_n_matches']),'Phobius; SignalP-Prok; InterPro-N'+('; TMbed' if t else ''),''])
titles['Topology & signals']=['Membrane topology and signal peptides — all 66 proteins',
 'Phobius, SignalP-Prok and InterPro-N: 66/66 verified. TMbed '+('verified.' if tm else 'is running; pending cells are not negative results.'),
 'Phobius: 15 TM / 10 signal peptides. TMbed: 15 TM / 11 signal peptides. SignalP-Prok: 9 signal peptides.',
 'Predictions are not experimental localization. An outside segment or signal peptide alone does not prove surface exposure.']
headers['Supplement raw results']=['Protein ID','Sequence MD5','Length (aa)','Library','Signature / feature','Description (reported)','Start (aa)','End (aa)','Score (source scale)','Prediction source','Date (source)','InterPro accession','InterPro description','GO annotations','Pathways','Predictor / library version','Run / release','Source / interpretation']
raw=[]
for r in csv.DictReader((DATA/'phobius/regions.tsv').open(),delimiter='\t'):
 ident=r['protein_id'];s=sequence[ident]
 raw.append([ident,hashlib.md5(s.encode()).hexdigest().upper(),len(s),'Phobius',r['feature'],r['description'] or r['feature'],int(r['start']),int(r['end']),'-','Phobius','06-10-2026','-','-','-','-','1.01','Hummingbird supplement 2026-10-06','Direct Phobius output. Non-TM domain labels alone do not establish cellular localization.'])
for r in csv.reader((DATA/'ebi/tsv.txt').open(),delimiter='\t'):
 assert len(r)==15
 line=list(r)
 for c in (2,6,7):line[c]=int(line[c])
 raw.append(line+[versions[(r[9],r[3],r[4])],'EMBL-EBI InterProScan 6.0.2.2 / InterPro 110.0','InterPro-N prediction; distinct from a conventional signature search.' if r[9]=='InterPro-N' else 'SignalP-Prok fast model; predicted signal-peptide region.'])
if tm:
 for r in csv.DictReader((DATA/'tmbed/regions.tsv').open(),delimiter='\t'):
  s=sequence[r['protein_id']]
  raw.append([r['protein_id'],hashlib.md5(s.encode()).hexdigest().upper(),len(s),'TMbed',r['feature'],r['description'],int(r['start']),int(r['end']),'-','TMbed','06-10-2026','-','-','-','-','1.0.2','Hummingbird supplement','Direct residue labels. Inside/outside predictions do not establish experimental localization.'])
data['Supplement raw results']=raw
titles['Supplement raw results']=['Supplemental raw predictions — coordinates and provenance',
 str(len(raw))+' region rows: 200 direct Phobius regions + 337 EMBL-EBI TSV rows'+(' + 347 TMbed regions.' if tm else '. TMbed is pending.'),
 'All 15 EMBL-EBI TSV fields are preserved. Direct Phobius records include domain labels for proteins without TM/signal predictions.',
 'InterPro-N is explicitly labeled; probabilities are not E-values. Version 110.0 predictions are kept separate from the original 109.0 scan.']
headers['Supplement run details']=['Field','Verified value / current status']
data['Supplement run details']=[
 ['Supplement status','Phobius, SignalP-Prok and InterPro-N verified for all 66 proteins. '+('TMbed verified complete.' if tm else 'TMbed is running; this supplement is not yet complete.')],
 ['Phobius','COMPLETE — 66/66; version 1.01; SLURM 796170 COMPLETED 0:0.'],
 ['SignalP-Prok','COMPLETE — 66/66; SignalP 6.0i, fast, prokaryotic model; EMBL-EBI.'],
 ['InterPro-N','COMPLETE — 66/66; 302 signature matches across 52 proteins; EMBL-EBI. Neural predictions are explicitly labeled.'],
 ['TMbed','COMPLETE — 66/66; version 1.0.2; analysis 796204 and independent verifier 796208 both COMPLETED 0:0.' if tm else 'RUNNING — SLURM 796204; 4 CPUs, 32 GB, maximum 1 hour including image preparation. No results claimed yet.'],
 ['DeepTMHMM','NOT RUN — user confirmed that no separate licensed standalone package is available. Not present in the earlier 19-protein export either. Phobius and TMbed provide topology predictions, but are not identical methods.'],
 ['SignalP-Euk','NOT RUN — eukaryotic model is not appropriate for this bacterial dataset; SignalP-Prok was run.'],
 ['Original 18-analysis run','Preserved, verified, unchanged. First five tabs document the original InterProScan 6.0.1 / InterPro 109.0 results.'],
 ['Scientific scope','Topology and signal peptides are predictions, not proof of secretion, antibody accessibility, virulence, or disease relevance.'],
 ['Phobius vs SignalP agreement','8 proteins positive in both methods; Phobius-only: group_541, group_5591; SignalP-only: group_3684. Review these discrepancies.'],
 ['Input records',66],['Processed by Phobius',66],['Processed by SignalP-Prok',66],['Processed by InterPro-N',66],
 ['Phobius transmembrane proteins',15],['Phobius signal-peptide proteins',10],['SignalP-Prok signal-peptide proteins',9],
 ['InterPro-N annotated proteins',52],['InterPro-N signature matches',302],['EMBL-EBI TSV region rows',337],['Direct Phobius region rows',200],
 ['Supplement directory',STATE['run_dir']],['Input SHA-256',pv['input_sha256']],
 ['Input preservation','All 66 original IDs and exact sequences retained. Original input, earlier completed results and earlier attempts were not overwritten.'],
 ['EMBL-EBI job ID',ev['job_id']],['EMBL-EBI verification','FINISHED plus successful workflow log, 66 JSON/XML records with exact input sequences, matching MD5/SHA-256, valid coordinates and agreeing XML/JSON.'],
 ['EMBL-EBI version','InterProScan 6.0.2.2; InterPro 110.0; Nextflow 26.04.0. Hosted service version; shared Hummingbird installation was not upgraded.'],
 ['EBI result files',STATE['run_dir']+'/results/ebi/'],['Phobius result files',STATE['run_dir']+'/results/phobius/'],['TMbed result files',STATE['run_dir']+'/results/tmbed/'],
 ['Phobius successful job','796170'],['Preserved Phobius setup failures','796162: non-UTF-8 comment; 796165: minimal system Perl lacked FindBin. Corrected with byte-preserving wrapper handling and installed perl/5.40.0.'],
 ['Preserved TMbed setup failure','796163: cluster mksquashfs utility missing. Retry 796204 builds a directory container directly from the same pinned image, with a compute-node pilot before the full run.'],
 ['TMbed analysis job',STATE['tmbed_job_id']],
 ['Independent verification job',STATE.get('verification_job_id','Not submitted')],
 ['Verification requirements','All 66 input IDs and sequences, raw prediction coordinates, JSON/XML/TSV agreement, file checksums, and successful SLURM accounting.'],
 ['TMbed image digest','interpro/tmbed@sha256:1031a299e650a8d32e3e2a3e4417d63d7e1ae5f752100688dde58283902f32b4'],
 ['TMbed reporting choice','Direct TMbed 1.0.2 format 4 retains complete per-residue inside/outside states; avoids the installed InterProScan 6.0.1 terminal-coordinate parser defect. No upgrade to the shared pipeline.'],
 ['Original 19-protein reference','Saved export: 62 Phobius rows, 16 SignalP rows, and InterPro-N-derived member signatures. This confirmed the coverage gap.'],
 ['Analysis status reference','https://interproscan6.readthedocs.io/stable/analyses/'],['Licensed applications','https://interproscan6.readthedocs.io/stable/licensed_applications/'],
 ['Phobius method / software','https://software.sbc.su.se/phobius.html'],['TMbed method / output','https://github.com/BernhoferM/TMbed'],
 ['InterPro-N hosted availability','https://www.ebi.ac.uk/about/news/updates-from-data-resources/interpro-110-0-and-interproscan-6-a-re-engineered-pipeline-and-ai-ready-access-to-the-protein-universe/'],
 ['User authorization','2026-10-06: user approved submitting the same 66 IDs/sequences to EMBL-EBI for SignalP-Prok and InterPro-N; saving results to Hummingbird and the Google Sheet.'],
 ['TMbed transmembrane proteins',sum(int(r['tm_segments'])>0 for r in tm.values()) if tm else 'Pending'],
 ['TMbed signal-peptide proteins',sum(r['signal_peptide']=='True' for r in tm.values()) if tm else 'Pending'],
 ['Final server verification',STATE['run_dir']+'/summary/run_summary.json'],
 ['Combined per-protein TSV',STATE['run_dir']+'/summary/per_protein_summary.tsv'],
 ['Combined per-protein CSV',STATE['run_dir']+'/summary/per_protein_summary.csv'],
 ['Final transfer verification','22 result files matched server SHA-256 hashes; all 66 TMbed sequences and per-residue labels independently reconciled locally.'],
 ['Results verified at UTC',datetime.datetime.now(datetime.timezone.utc).isoformat()]
]
titles['Supplement run details']=['Supplement coverage and provenance — October 6, 2026',
 'Four supplemental analyses verified for all 66 proteins: Phobius, SignalP-Prok, InterPro-N and TMbed.' if tm else 'Three analyses verified; TMbed pending.',
 'The earlier Hummingbird results remain preserved. This tab separates completed, pending, unavailable and inapplicable analyses.',
 'Full raw outputs and verification files are stored in the separate supplement directory below.']
ids={'Topology & signals':660106,'Supplement raw results':660107,'Supplement run details':660108}
def cell(v):
 if v=='' or v is None:return {}
 return {'userEnteredValue':{('numberValue' if isinstance(v,(int,float)) else 'stringValue'):v}}
def grid(sid,r0,r1,c0,c1):return {'sheetId':sid,'startRowIndex':r0,'endRowIndex':r1,'startColumnIndex':c0,'endColumnIndex':c1}
plan={};structures=[];formats=[]
for name,sid in ids.items():
 n=len(headers[name]);rows=[['']*n]+[[t]+['']*(n-1) for t in titles[name]]+[['']*n]+[headers[name]]+data[name]
 assert all(len(r)==n for r in rows) and all(len(str(v))<50000 for r in rows for v in r)
 slug=name.lower().replace(' & ','_').replace(' ','_');(OUT/(slug+'_expected.json')).write_text(json.dumps(rows))
 structures.append({'addSheet':{'properties':{'sheetId':sid,'title':name,'gridProperties':{'rowCount':len(rows),'columnCount':n,'frozenRowCount':7,'frozenColumnCount':1,'hideGridlines':True}}}})
 files=[]
 for start in range(0,len(rows),100):
  filename=slug+'_'+str(start)+'.json';files.append(filename)
  (OUT/filename).write_text(json.dumps({'updateCells':{'range':grid(sid,start,min(start+100,len(rows)),0,n),'rows':[{'values':[cell(v) for v in r]} for r in rows[start:start+100]],'fields':'userEnteredValue'}}))
 formats.extend([
  {'repeatCell':{'range':grid(sid,0,len(rows),0,n),'cell':{'userEnteredFormat':{'verticalAlignment':'TOP','wrapStrategy':'WRAP','textFormat':{'fontFamily':'Arial','fontSize':11,'foregroundColor':{'red':.12,'green':.16,'blue':.22}}}},'fields':'userEnteredFormat'}},
  {'repeatCell':{'range':grid(sid,6,7,0,n),'cell':{'userEnteredFormat':{'backgroundColor':{'red':.15,'green':.23,'blue':.34},'textFormat':{'foregroundColor':{'red':1,'green':1,'blue':1},'bold':True},'verticalAlignment':'MIDDLE'}},'fields':'userEnteredFormat.backgroundColor,userEnteredFormat.textFormat,userEnteredFormat.verticalAlignment'}},
  {'repeatCell':{'range':grid(sid,1,5,0,n),'cell':{'userEnteredFormat':{'wrapStrategy':'OVERFLOW_CELL','verticalAlignment':'MIDDLE'}},'fields':'userEnteredFormat.wrapStrategy,userEnteredFormat.verticalAlignment'}},
  {'repeatCell':{'range':grid(sid,1,2,0,n),'cell':{'userEnteredFormat':{'textFormat':{'fontFamily':'Arial','fontSize':15,'bold':True}}},'fields':'userEnteredFormat.textFormat'}},
  {'addConditionalFormatRule':{'index':0,'rule':{'ranges':[grid(sid,7,len(rows),0,n)],'booleanRule':{'condition':{'type':'CUSTOM_FORMULA','values':[{'userEnteredValue':'=MOD(ROW(),2)=0'}]},'format':{'backgroundColor':{'red':.95,'green':.96,'blue':.98}}}}}},
  {'setBasicFilter':{'filter':{'range':grid(sid,6,len(rows),0,n)}}}
 ])
 widths=([150,95,125,220,150,185,185,175,160,270,180,190,190,170,250,320] if sid==660106 else [150,350,100,150,220,450,95,95,180,160,140,170,350,260,180,180,300,430] if sid==660107 else [300,900])
 for i,w in enumerate(widths):formats.append({'updateDimensionProperties':{'range':{'sheetId':sid,'dimension':'COLUMNS','startIndex':i,'endIndex':i+1},'properties':{'pixelSize':w},'fields':'pixelSize'}})
 for a,b,h in [(0,1,12),(1,2,34),(2,5,26),(5,6,24),(6,7,62)]:formats.append({'updateDimensionProperties':{'range':{'sheetId':sid,'dimension':'ROWS','startIndex':a,'endIndex':b},'properties':{'pixelSize':h},'fields':'pixelSize'}})
 formats.append({'autoResizeDimensions':{'dimensions':{'sheetId':sid,'dimension':'ROWS','startIndex':7,'endIndex':len(rows)}}})
 plan[name]={'sheet_id':sid,'rows':len(rows),'cols':n,'files':files,'expected':slug+'_expected.json'}
formats.append({'addConditionalFormatRule':{'index':0,'rule':{'ranges':[grid(660106,7,73,9,10)],'booleanRule':{'condition':{'type':'TEXT_CONTAINS','values':[{'userEnteredValue':'REVIEW:'}]},'format':{'backgroundColor':{'red':1,'green':.91,'blue':.69},'textFormat':{'bold':True}}}}}})
formats.append({'addConditionalFormatRule':{'index':0,'rule':{'ranges':[grid(660106,7,73,10,13)],'booleanRule':{'condition':{'type':'CUSTOM_FORMULA','values':[{'userEnteredValue':'=$K8<>"Verified complete"'}]},'format':{'backgroundColor':{'red':1,'green':.96,'blue':.82}}}}}})
(OUT/'plan.json').write_text(json.dumps(plan,indent=2));(OUT/'structure.json').write_text(json.dumps(structures));(OUT/'format.json').write_text(json.dumps(formats));print(json.dumps(plan,indent=2))
