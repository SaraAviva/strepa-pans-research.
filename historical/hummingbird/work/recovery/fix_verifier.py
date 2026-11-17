import ast
import pathlib

base=pathlib.Path(__file__).resolve().parent
source=(base/'retry_scripts/postprocess_retry_20261004.py').read_text()
source=source.replace("run = pathlib.Path(__file__).resolve().parent.parent", '''def normalized_md5(value):
    # Hexadecimal digest spelling is case-insensitive; sequence SHA-256 stays exact.
    return value.lower() if isinstance(value,str) and re.fullmatch(r'[0-9A-Fa-f]{32}',value) else None

run = pathlib.Path(__file__).resolve().parent.parent''')
source=source.replace("if 'InterProScan completed successfully' not in logtext or 'all processes finished' not in logtext:\n    issues.append('Nextflow success and all-processes-finished log evidence missing')", '''analysis_job=completion.get('slurm_job_id','')
stdoutpath=run/'logs'/('analysis-'+str(analysis_job)+'.out')
stdouttext=stdoutpath.read_text(errors='replace') if stdoutpath.exists() else ''
if 'InterProScan completed successfully' not in stdouttext:
    issues.append('InterProScan success message missing from analysis stdout')
if 'all processes finished' not in logtext:
    issues.append('Nextflow all-processes-finished log evidence missing')''')
source=source.replace("md5!=original['sequence_md5']", "normalized_md5(md5)!=normalized_md5(original['sequence_md5'])")
source=source.replace("nohit_md5.update(read_json(p))", "nohit_md5.update(normalized_md5(value) for value in read_json(p))")
source=source.replace("cols[1]!=r['sequence_md5']", "normalized_md5(cols[1])!=normalized_md5(r['sequence_md5'])")
source=source.replace("result.get('md5')!=original['sequence_md5']", "normalized_md5(result.get('md5'))!=normalized_md5(original['sequence_md5'])")
source=source.replace("original['sequence_md5'] not in nohit_md5", "normalized_md5(original['sequence_md5']) not in nohit_md5")
source=source.replace("p=run/('provenance/retry_20261004_'+role+'_job_id.txt')", "p=run/('provenance/final_postprocess_job_id.txt' if role=='postprocess' else 'provenance/retry_20261004_analysis_job_id.txt')")
source=source.replace('provenance/retry_20261004_postprocess_submission_command.txt','provenance/final_postprocess_submission_command.txt')
source=source.replace("nohit_md5=set()\nfor p in (run/'work').glob('*/*/no_matches.json'):\n    nohit_md5.update(normalized_md5(value) for value in read_json(p))", '''# Use only explicit empty maps emitted by successful tasks in this run's trace.
# REPORT_NO_MATCHES omits a sequence if a parser emitted its key with an empty
# map (observed for CDD/group_2090). COMBINE_MATCHES retains that empty map.
nohit_md5=set()
nohit_evidence={}
for task in trace:
    if task.get('status') not in ('COMPLETED','CACHED') or task.get('exit') not in ('0','-',None,''):
        continue
    name=task.get('name','')
    pattern = ('no_matches.json' if ':REPORT_NO_MATCHES ' in name else
               'combined/*.json' if ':COMBINE_MATCHES ' in name else None)
    work_hash=task.get('hash','')
    if not pattern or not re.fullmatch(r'[0-9a-f]{2}/[0-9a-f]+',work_hash):
        continue
    for directory in (run/'work').glob(work_hash+'*'):
        for path in directory.glob(pattern):
            data=read_json(path)
            if not isinstance(data,dict):
                issues.append('Malformed no-match evidence: '+str(path)); continue
            for value,matches in data.items():
                digest=normalized_md5(value)
                if digest and matches == {}:
                    nohit_md5.add(digest)
                    nohit_evidence.setdefault(digest,[]).append(dict(
                        task=name, path=str(path.relative_to(run)), sha256=sha(path)))''')
source=source.replace("if result is not None and not isinstance(matches,list):", "if result is not None and ('matches' not in result or not isinstance(matches,list)):")
source=source.replace("if main_rows and (main_rows[0][1]!='COMPLETED' or main_rows[0][2]!='0:0'):", "if proc.returncode != 0 or len(main_rows) != 1 or main_rows[0][1]!='COMPLETED' or main_rows[0][2]!='0:0':")
source=source.replace("verified_at_utc=datetime.now(timezone.utc).isoformat(),issues=issues", "verified_at_utc=datetime.now(timezone.utc).isoformat(),no_match_evidence=nohit_evidence,issues=issues")
ast.parse(source)
(base/'postprocess_final_20261004.py').write_text(source)
print('Corrected digest comparisons and success-log source; sequence and task checks retained.')
