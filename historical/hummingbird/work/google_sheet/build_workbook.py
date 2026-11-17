import csv,collections,datetime,hashlib,json,pathlib,re
from urllib.parse import urlparse,unquote,quote

root=pathlib.Path(__file__).resolve().parents[2]
out=root/'work/google_sheet'
run=root/'outputs/Hummingbird_results_20261004'
summary=json.loads((run/'summary/run_summary.json').read_text())
audit=json.loads((run/'provenance/input_audit.json').read_text())
db=json.loads((run/'provenance/database_versions.json').read_text())
proteins=list(csv.DictReader((run/'summary/per_protein_summary.csv').open()))
records=list(csv.DictReader((run/'provenance/input_records.tsv').open(),delimiter='\t'))
raw=list(csv.reader((run/'results/all_significant_proteins.tsv').open(),delimiter='\t'))
structured=json.loads((run/'results/all_significant_proteins.json').read_text())
byid={x['id']:r for r in structured['results'] for x in r['xref']}
assert len(proteins)==len(records)==len(byid)==66 and summary['verified_complete']
assert hashlib.sha256((run/'input/all_significant_proteins.fasta').read_bytes()).hexdigest()==audit['sha256']
ids={'All 66 proteins':1409174111,'Annotation review':1538520193,'Source results':2085426790,'Run details':660104,'Pathway mappings':660105}
pathway_strings=list(dict.fromkeys(r[14] for r in raw if r[14]!='-'))
pathway_ids={s:f'PW{i:03d}' for i,s in enumerate(pathway_strings,1)}
tableids={'All 66 proteins':'852420801','Annotation review':'969953018','Source results':'1990368'}
navy={'red':38/255,'green':60/255,'blue':88/255}
ink={'red':31/255,'green':41/255,'blue':55/255}
muted={'red':88/255,'green':101/255,'blue':117/255}
white={'red':1,'green':1,'blue':1}
light={'red':243/255,'green':246/255,'blue':250/255}
def cell(v):
    if isinstance(v,dict):return v
    if v is None or v=='':return {}
    key='boolValue' if isinstance(v,bool) else 'numberValue' if isinstance(v,(int,float)) else 'stringValue'
    return {'userEnteredValue':{key:v}}
def formula(s):return {'userEnteredValue':{'formulaValue':s}}
def richlinks(pairs):
    pairs=list(dict.fromkeys(pairs))
    if not pairs:return cell('No entry URL recorded')
    text='';runs=[]
    for label,url in pairs:
        if text:text+='\n'
        start=len(text.encode('utf-16-le'))//2
        runs.append({'startIndex':start,'format':{'link':{'uri':url},'foregroundColor':{'red':0.07,'green':0.33,'blue':0.8},'underline':True}})
        text+=label
        runs.append({'startIndex':len(text.encode('utf-16-le'))//2,'format':{}})
    runs=runs[:-1]
    return {'userEnteredValue':{'stringValue':text},'textFormatRuns':runs}
def plainlinks(text):
    pairs=[]
    for u in filter(None,text.split('; ')):
        path=urlparse(u).path
        label=unquote(path.rstrip('/').split('/')[-1])
        pairs.append((label,u))
    return richlinks(pairs)
def grid(sid,r0,r1,c0,c1):return dict(sheetId=sid,startRowIndex=r0,endRowIndex=r1,startColumnIndex=c0,endColumnIndex=c1)
def update(sid,r0,values):
    width=len(values[0]);assert all(len(row)==width for row in values)
    return {'updateCells':{'range':grid(sid,r0,r0+len(values),0,width),'rows':[{'values':[cell(v) for v in row]} for row in values],'fields':'userEnteredValue,textFormatRuns,note'}}
def repeat(r,fmt,fields='userEnteredFormat'):return {'repeatCell':{'range':r,'cell':{'userEnteredFormat':fmt},'fields':fields}}
def dim(sid,axis,start,end,pixels):return {'updateDimensionProperties':{'range':{'sheetId':sid,'dimension':axis,'startIndex':start,'endIndex':end},'properties':{'pixelSize':pixels},'fields':'pixelSize'}}
headers={}
rows={}
top={}
headers['All 66 proteins']=['Protein ID','Length (aa)','Processing status','Evidence category','Signature matches','Raw hit rows','Member databases','Member signatures / descriptions','Integrated InterPro entries','GO terms','Pathway sets (see mappings)','Entry links','Family/domain signature','AntiFam match','Sequence MD5','Sequence SHA-256','Original record #','Verification notes','Protein sequence','Your follow-up notes']
cat={'family_domain_or_functional_signature':'Family/domain/functional signature','disorder_only':'Disorder only','disorder_or_coiled_coil_only':'Disorder/coiled-coil only','no_matches':'No matches'}
raw_byid=collections.defaultdict(list)
for r in raw:raw_byid[r[0]].append(r)
rows['All 66 proteins']=[]
for i,(p,record) in enumerate(zip(proteins,records),8):
    ident=p['original_id'];assert ident==record['original_id']
    seq=byid[ident]['sequence']
    assert hashlib.sha256(seq.encode()).hexdigest()==record['sequence_sha256']
    pathways=sorted({pathway_ids[r[14]] for r in raw_byid[ident] if r[14]!='-'})
    status='Completed with annotations' if p['processing_status']=='completed_with_annotations' else 'Completed — no hits'
    notes='Exact sequence and ID verified in JSON/XML; all 51 workflow tasks completed.'
    if not byid[ident]['matches']:notes+=' Explicit empty match list and successful intermediate task evidence confirm no hits.'
    row=[ident,int(p['length']),status,cat[p['evidence_category']],int(p['match_count']),formula(f'=COUNTIF(\'Source results\'!$A$8:$A$377,A{i})'),p['member_databases'].replace('; ','\n'),p['member_signatures_descriptions'].replace('; ','\n'),p['integrated_interpro_entries'].replace('; ','\n'),p['go_terms'].replace('; ','\n'),'\n'.join(pathways),plainlinks(p['identifier_links']) if p['identifier_links'] else '',p['has_functional_or_domain_signature']=='True',p['has_antifam_match']=='True',p['sequence_md5'],record['sequence_sha256'],int(p['record_number']),notes,seq,'']
    rows['All 66 proteins'].append(row)
top['All 66 proteins']=[
 'All 66 significant proteins — verified InterProScan results',
 '66 processed • 51 with matches • 15 completed with no hits • 0 failed/unverified. Original FASTA order retained.',
 '43 proteins have integrated InterPro entries; 30 have GO terms. 44 have family/domain evidence; 7 have only disorder/coiled-coil predictions.',
 'Run completed 2026-10-04 • InterProScan 6.0.1 / InterPro 109.0. Scroll right for all annotations, sequence checks, and follow-up notes.']

# Group exactly the current raw library/accession/description tuples, retaining
# every source row in Source results and retaining all proteins in the first tab.
groups=collections.defaultdict(list)
for r in raw:groups[tuple(r[3:6])].append(r)
order={p['original_id']:i for i,p in enumerate(proteins)}
endpoints={'Pfam':'pfam','PANTHER':'panther','CDD':'cdd','CATH-Gene3D':'cathgene3d','CATH-FunFam':'cathfunfam','HAMAP':'hamap','NCBIFAM':'ncbifam','PIRSF':'pirsf','PRINTS':'prints','PROSITE patterns':'prosite','PROSITE profiles':'prosite','SFLD':'sfld','SMART':'smart','SUPERFAMILY':'ssf'}
headers['Annotation review']=['InterProScan hit','Library','Signature accession','Summary of what this hit means','Sources (database entries)','Protein IDs with this annotation','Your PANS/PANDAS notes','Description / evidence availability']
rows['Annotation review']=[]
for (lib,acc,desc),hits in sorted(groups.items(),key=lambda x:(x[0][1]=='mobidb-lite',x[0][2].lower(),x[0][0],x[0][1])):
    iprs=sorted({(r[11],r[12]) for r in hits if r[11]!='-'})
    proteinids=sorted({r[0] for r in hits},key=order.get)
    if lib=='MobiDB-lite':
        meaning='Predicted intrinsic disorder. This is a structural feature and does not assign a specific biochemical function.'
    elif lib=='COILS':
        meaning='Predicted coiled-coil region. This is a structural feature and does not assign a specific biochemical function.'
    elif re.search(r'unknown function|\bDUF\d+|uncharacteri[sz]ed',desc,re.I):
        meaning='Match to a family/domain labelled as unknown or uncharacterized in this run. A specific biochemical role is not established by this annotation.'
    elif desc in ('-',acc,''):
        meaning='A signature match was reported, but this raw export supplies no descriptive function. Use the linked entry, if available, for follow-up.'
    else:
        meaning='Sequence evidence for “'+desc+'”.'
        if iprs:meaning+=' Integrated InterPro: '+'; '.join(name for _,name in iprs)+'.'
        meaning+=' This annotation is not experimental confirmation of function.'
    pairs=[]
    if lib in endpoints:pairs.append((lib+': '+acc,'https://www.ebi.ac.uk/interpro/entry/'+endpoints[lib]+'/'+quote(acc,safe='')+'/'))
    if lib=='MobiDB-lite':pairs.append(('MobiDB-lite documentation','https://github.com/BioComputingUP/MobiDB-lite'))
    pairs.extend(('InterPro: '+acc,'https://www.ebi.ac.uk/interpro/entry/InterPro/'+acc+'/') for acc,_ in iprs)
    availability=('Description reported in current raw TSV' if desc not in ('-',acc,'') else 'Identifier only; no descriptive function in raw TSV')
    availability+='\n'+( '; '.join(acc for acc,_ in iprs) if iprs else 'No integrated InterPro entry reported')
    availability+=f'\n{len(proteinids)} protein(s); {len(hits)} raw location row(s)'
    rows['Annotation review'].append([desc if desc not in ('-','') else acc,lib,acc,meaning,richlinks(pairs),', '.join(proteinids),'',availability])
assert len(rows['Annotation review'])==228
top['Annotation review']=[
 'InterProScan annotation review — all significant proteins',
 '228 distinct library/accession/description combinations across 51 proteins with matches. All 66 proteins remain in the first tab.',
 'Summaries use the current run’s reported annotations. Source links are database entries; no new paper review was performed.',
 'Disorder/coiled-coil features are distinguished from family/domain evidence. Your PANS/PANDAS notes are left blank for interpretation.']

headers['Source results']=['Protein ID','Sequence MD5','Length (aa)','Library','Signature accession','Description (TSV)','Start (aa)','End (aa)','Score (source text)','Prediction source (TSV column 10)','Run date (DD-MM-YYYY)','InterPro accession','InterPro description','GO annotations','Pathway set (see mappings)']
rows['Source results']=[]
for r in raw:
    vals=list(r)
    if vals[14]!='-':vals[14]=pathway_ids[vals[14]]
    for i in [2,6,7]:vals[i]=int(vals[i])
    rows['Source results'].append(vals)
top['Source results']=[
 'Complete raw InterProScan results — 370 location rows',
 'Source: results/all_significant_proteins.tsv • all 370 rows retained in source order; all 15 fields accounted for.',
 'First 14 fields retain source content. Oversized pathway fields are losslessly stored by set/position in Pathway mappings; a dash remains a dash.',
 'This TSV contains hits only. The first tab confirms all 15 no-hit proteins using JSON/XML plus successful workflow evidence.']

meta=summary['software']
remote=meta['run_dir']
details=[
 ['Outcome','VERIFIED COMPLETE'],['Run completed','2026-10-04'],['Input records',66],['Processed records',66],['Failed or unverified records',0],
 ['Proteins with any match',51],['Completed with no hits',15],['With integrated InterPro entries',43],['With GO terms',30],
 ['Family/domain/functional evidence',44],['Disorder only',6],['Disorder/coiled-coil only',1],['Raw hit rows',370],['Distinct annotation groups',228],
 ['Analysis job','792566 — COMPLETED, exit 0:0, elapsed 00:13:30'],['Final verification job','792576 — COMPLETED, exit 0:0'],['Workflow tasks','51 completed; 0 unsuccessful'],
 ['Input FASTA',meta['source']],['Run directory',remote],['Input SHA-256',audit['sha256']],['Input preservation',audit['input_modifications']],['Unique protein IDs',audit['unique_ids']],
 ['Sequence lengths',f"{audit['minimum_length']}–{audit['maximum_length']} aa; {audit['total_residues']} total residues"],['Input validation','No empty records, duplicate IDs, duplicate sequences, gaps, or invalid records. Eight X residues were accepted by the installed version.'],
 ['InterProScan version',meta['software_version']],['InterPro release',meta['interpro_version']],['Nextflow',meta['nextflow_version']],['SingularityCE',meta['singularity_version']],['Pipeline commit',meta['pipeline_git_commit']],['Database location',meta['database_root']],
 ['Execution model',meta['execution_model']],['Allocation','8 CPUs, 32 GB, 2-hour bound; local task cap 6 CPUs / 24 GB / queue 6'],['Enabled analyses',', '.join(meta['enabled_analyses'])],
 ['Skipped analyses','\n'.join(k+': '+v for k,v in meta['skipped_analyses'].items())],
 ['Member database versions (enabled)','\n'.join(k+': '+db.get('MobiDB Lite' if k=='MobiDB-lite' else k,'see manifest') for k in meta['enabled_analyses'])],
 ['Verification method','All 66 IDs and exact sequences reconciled with structured JSON/XML; SHA-256 and MD5 checked; raw coordinates checked; all required analysis endpoints and 51 successful tasks verified; SLURM analysis and final verification COMPLETED 0:0.'],
 ['No-hit verification','15 explicit empty matches arrays, plus successful intermediate empty-map evidence. One CDD empty map was retained in COMBINE_MATCHES rather than REPORT_NO_MATCHES.'],
 ['Interpretation',summary['interpretation']],['Unresolved run issues','None'],['Annotation review scope','Grouped current-run annotations with database-entry links and cautious descriptions. New literature searches and disease relevance assessments were not performed.'],
 ['Original failed jobs','785450 failed during container startup; 785451 canceled after dependency failure. Preserved.'],
 ['Container recovery','Existing SIFs extracted into run-local sandboxes; all 12 containers passed compute-node tests; original software/database versions retained.'],
 ['Verifier recovery','792567 failed because of MD5 case comparison, success-log selection, and no-hit intermediate handling. Checks corrected and regression-tested; 792576 completed.'],
 ['Analysis command',(run/'provenance/retry_20261004_analysis_submission_command.txt').read_text().strip()],
 ['Nextflow command',(run/'provenance/nextflow_command.txt').read_text().strip()],
 ['Final verification command',(run/'provenance/final_postprocess_submission_command.txt').read_text().strip()],
 ['Status command','bash '+remote+'/scripts/status_final.sh'],
 ['Per-protein CSV',remote+'/summary/per_protein_summary.csv'],['Per-protein TSV',remote+'/summary/per_protein_summary.tsv'],
 ['Raw TSV',remote+'/results/all_significant_proteins.tsv'],['Raw JSON (full structured output)',remote+'/results/all_significant_proteins.json'],['Raw XML (full structured output)',remote+'/results/all_significant_proteins.xml'],
 ['Run summary JSON',remote+'/summary/run_summary.json'],['Run summary text',remote+'/summary/RUN_SUMMARY.txt'],['Full downloadable archive on Hummingbird',remote+'/final_results_20261004.tar.gz'],
 ['Archive SHA-256',json.loads((run/'LOCAL_VERIFICATION.json').read_text())['archive_sha256']],
 ['Transfer verification','107 files verified locally against server SHA-256 manifest; all 66 summary rows retained in original input order.'],
 ['Input audit',remote+'/provenance/input_audit.json'],['Database manifest',remote+'/provenance/database_versions.json'],['Scripts / logs',remote+'/scripts/ ; '+remote+'/logs/'],
 ['Style reference',richlinks([('Evan 19 Hypotheticals — original worksheet','https://docs.google.com/spreadsheets/d/SOURCE_SPREADSHEET_ID/edit')])],
 ['Data lineage','Copied the reference workbook’s native layout; replaced its old 19-protein scan content with the verified 66-protein run. Original reference unchanged.'],
 ['Source export status','Source results preserves the first 14 TSV fields and source row order. Column 15 references 49 pathway sets, whose 12,001 ordered tokens are fully retained in Pathway mappings. Join tokens with | in position order to reconstruct each original pathway field. Full nested JSON/XML, logs and provenance remain at the paths above.'],
 ['Pathway interpretation','Pathway mappings are reported database cross-references, including cross-species identifiers; they do not establish pathway activity in the scanned StrepA proteins.']
]
headers['Run details']=['Field','Verified value / provenance']
rows['Run details']=details
top['Run details']=['Run details and provenance — verified 2026-10-04','66/66 original proteins accounted for • analysis 792566 and final verification 792576 completed successfully.','This workbook presents the completed run; the original input, raw files, logs, scripts, and failed attempts remain saved.','Full nested JSON/XML and run archive paths are listed below.']
headers['Pathway mappings']=['Pathway set ID','Position in source field','Pathway mapping (source text)']
rows['Pathway mappings']=[[pathway_ids[s],i,t] for s in pathway_strings for i,t in enumerate(s.split('|'),1)]
assert len(rows['Pathway mappings'])==12001
top['Pathway mappings']=['Complete pathway fields — lossless normalization','49 distinct source fields • 12,001 ordered mappings. Source results column O points to the applicable set.','To reconstruct a source field, select its set, sort by position, and join mapping values with |. Nothing was dropped.','These database cross-references include cross-species IDs and are not evidence of pathway activity in the scanned proteins.']

structure=[]
for sid,n in [(1538520193,3),(1409174111,2)]:
    for index in reversed(range(n)):structure.append({'deleteConditionalFormatRule':{'sheetId':sid,'index':index}})
structure += [
 {'updateSheetProperties':{'properties':{'sheetId':1409174111,'title':'All 66 proteins','index':0},'fields':'title,index'}},
 {'updateSheetProperties':{'properties':{'sheetId':1538520193,'index':1},'fields':'index'}},
 {'addSheet':{'properties':{'sheetId':660104,'title':'Run details','index':3,'gridProperties':{'rowCount':1000,'columnCount':2,'frozenRowCount':7,'frozenColumnCount':1,'hideGridlines':True}}}},
 {'addSheet':{'properties':{'sheetId':660105,'title':'Pathway mappings','index':4,'gridProperties':{'rowCount':12008,'columnCount':3,'frozenRowCount':7,'frozenColumnCount':1,'hideGridlines':True}}}}
]
for name,sid in ids.items():
    width=len(headers[name]);end=7+len(rows[name])
    if name in tableids:
        structure.append({'updateTable':{'table':{'tableId':tableids[name],'range':grid(sid,6,end,0,width),'columnProperties':[{'columnIndex':i,'columnName':h} for i,h in enumerate(headers[name])]},'fields':'range,columnProperties'}})
    # Clear only copied row-specific values, hyperlinks and notes; native styles
    # and table identities remain available for the current result.
    if name not in ['Run details','Pathway mappings']:
        structure.append({'repeatCell':{'range':grid(sid,0,end,0,width),'cell':{},'fields':'userEnteredValue,textFormatRuns,note,dataValidation,chipRuns'}})
    if name in ['Annotation review','Source results']:
        structure.append({'copyPaste':{'source':grid(sid,7,8,0,width),'destination':grid(sid,7,end,0,width),'pasteType':'PASTE_FORMAT'}})
    # Trim only unused copied grid outside the exact result rectangle.
    if end<1000:structure.append({'deleteDimension':{'range':{'sheetId':sid,'dimension':'ROWS','startIndex':end,'endIndex':1000}}})
    if width<26 and name not in ['Run details','Pathway mappings']:
        structure.append({'deleteDimension':{'range':{'sheetId':sid,'dimension':'COLUMNS','startIndex':width,'endIndex':26}}})

(out/'01_structure.json').write_text(json.dumps(structure))
for name in ids:
    sid=ids[name];width=len(headers[name])
    title_rows=[['']*width for _ in range(7)]
    for i,t in enumerate(top[name],1):title_rows[i][0]=t
    title_rows[6]=headers[name]
    dataset=title_rows+rows[name]
    slug=name.lower().replace(' ','_')
    chunk=35 if name=='All 66 proteins' else 1000 if name=='Pathway mappings' else 80
    files=[]
    for start in range(0,len(dataset),chunk):
        path=out/(slug+f'_{start:03d}.json')
        path.write_text(json.dumps({'sheetId':sid,'startRow':start,'values':dataset[start:start+chunk]},separators=(',',':')))
        files.append(path.name)
    (out/(slug+'_expected.json')).write_text(json.dumps(dataset))
    rows[name]=dict(data_rows=len(rows[name]),total_rows=len(dataset),cols=width,files=files)

formats=[]
widths={'All 66 proteins':[165,100,200,260,125,115,175,540,450,165,300,245,150,130,275,470,115,460,600,340],
        'Annotation review':[290,155,185,490,310,410,360,390],
        'Source results':[160,280,100,150,175,380,100,100,190,220,170,165,380,290,330],
        'Run details':[300,1000],'Pathway mappings':[210,230,830]}
for name,sid in ids.items():
    spec=rows[name];end=spec['total_rows'];width=spec['cols']
    # Matching reference type, header palette, wrapping and frozen header.
    formats.append(repeat(grid(sid,7,end,0,width),{'verticalAlignment':'TOP','wrapStrategy':'WRAP','textFormat':{'fontFamily':'Arial','fontSize':11,'foregroundColor':ink},'numberFormat':{'type':'TEXT'},'backgroundColor':white}))
    formats.append(repeat(grid(sid,6,7,0,width),{'backgroundColor':navy,'verticalAlignment':'MIDDLE','wrapStrategy':'WRAP','textFormat':{'fontFamily':'Arial','fontSize':11,'bold':True,'foregroundColor':white}}))
    formats.append(repeat(grid(sid,1,2,0,width),{'textFormat':{'fontFamily':'Arial','fontSize':15,'bold':True,'foregroundColor':navy},'verticalAlignment':'MIDDLE','wrapStrategy':'WRAP'}))
    formats.append(repeat(grid(sid,2,5,0,width),{'textFormat':{'fontFamily':'Arial','fontSize':11,'foregroundColor':muted},'verticalAlignment':'MIDDLE','wrapStrategy':'WRAP'}))
    for rr in range(1,5):
        # Header text spans the first visible columns; data tables are unmerged.
        formats.append({'mergeCells':{'range':grid(sid,rr,rr+1,0, min(width, 7 if name=='All 66 proteins' else 4 if name!='Run details' else 2)),'mergeType':'MERGE_ALL'}})
    formats.append({'addConditionalFormatRule':{'index':0,'rule':{'ranges':[grid(sid,7,end,0,width)],'booleanRule':{'condition':{'type':'CUSTOM_FORMULA','values':[{'userEnteredValue':'=ISEVEN(ROW())'}]},'format':{'backgroundColor':light}}}}})
    formats += [dim(sid,'ROWS',0,1,12),dim(sid,'ROWS',1,2,34),dim(sid,'ROWS',2,5,38),dim(sid,'ROWS',5,6,12),dim(sid,'ROWS',6,7,52)]
    if name in ['All 66 proteins','Run details','Pathway mappings']:
        for i,pix in enumerate(widths[name]):formats.append(dim(sid,'COLUMNS',i,i+1,pix))
    if name=='All 66 proteins':
        for i in [1,4,5,16]:formats.append(repeat(grid(sid,7,end,i,i+1),{'numberFormat':{'type':'NUMBER','pattern':'0'}},'userEnteredFormat.numberFormat'))
        for status,color in [('Completed — no hits',{'red':1,'green':0.95,'blue':0.82}),('Completed with annotations',{'red':0.87,'green':0.95,'blue':0.89})]:
            formats.append({'addConditionalFormatRule':{'index':0,'rule':{'ranges':[grid(sid,7,end,2,3)],'booleanRule':{'condition':{'type':'TEXT_EQ','values':[{'userEnteredValue':status}]},'format':{'backgroundColor':color}}}}})
    if name=='Source results':
        for i in [2,6,7]:formats.append(repeat(grid(sid,7,end,i,i+1),{'numberFormat':{'type':'NUMBER','pattern':'0'}},'userEnteredFormat.numberFormat'))
    if name=='Pathway mappings':
        formats.append(dim(sid,'ROWS',7,end,26))
        formats.append(repeat(grid(sid,7,end,1,2),{'numberFormat':{'type':'NUMBER','pattern':'0'}},'userEnteredFormat.numberFormat'))
        formats.append({'setBasicFilter':{'filter':{'range':grid(sid,6,end,0,width)}}})
    else:formats.append({'autoResizeDimensions':{'dimensions':{'sheetId':sid,'dimension':'ROWS','startIndex':7,'endIndex':end}}})
(out/'03_format.json').write_text(json.dumps(formats))
plan={'spreadsheet_id':'OUTPUT_SPREADSHEET_ID','spreadsheet_url':'https://docs.google.com/spreadsheets/d/OUTPUT_SPREADSHEET_ID/edit?usp=drivesdk','source_id':'SOURCE_SPREADSHEET_ID','sheets':{name:{'sheet_id':ids[name],**spec} for name,spec in rows.items()}}
(out/'plan.json').write_text(json.dumps(plan,indent=2)+'\n')
for p in out.glob('*.json'):
    if p.name.endswith('_expected.json') or p.name=='plan.json':continue
    parsed=json.loads(p.read_text())
    if isinstance(parsed,dict):continue
    for req in parsed:
        assert len(req)==1
        if 'updateCells' in req:
            u=req['updateCells'];r=u['range']
            assert len(u['rows'])==r['endRowIndex']-r['startRowIndex']
            assert all(len(x['values'])==r['endColumnIndex']-r['startColumnIndex'] for x in u['rows'])
print(json.dumps(plan,indent=2))
