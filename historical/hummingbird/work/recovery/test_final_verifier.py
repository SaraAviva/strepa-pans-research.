import csv,json,pathlib,shutil,subprocess,sys

base=pathlib.Path(__file__).resolve().parent
source=(base.parent/'test_postprocess.py').read_text()
source=source.replace("script=here/'postprocess.py'", "script=here/'postprocess_final_20261004.py'")
source=source.replace("base=here/'verification-tests'", "base=here/'final-verification-tests'")
source=source.replace("['complete','missing_structured_record','failed_task','missing_nohit_evidence']", "['complete','empty_combined_map','missing_structured_record','failed_task','missing_nohit_evidence','wrong_digest','changed_sequence','missing_success_stdout','missing_matches_array','failed_combine_task']")
source=source.replace('provenance/run_metadata.json','provenance/retry_20261004_run_metadata.json')
source=source.replace('dict(exit_code=0)',"dict(exit_code=0,slurm_job_id='fixture')")
source=source.replace("(run/'logs/nextflow.log').write_text('InterProScan completed successfully\\nall processes finished\\n')", "(run/'logs/nextflow.log').write_text('all processes finished\\n')\n    (run/'logs/analysis-fixture.out').write_text('InterProScan completed successfully\\n' if case!='missing_success_stdout' else '')")
source=source.replace("fieldnames=['name','status','exit']", "fieldnames=['name','status','exit','hash']")
source=source.replace("    outputs=[]", "        w.writerow(dict(name='INTERPROSCAN:SCAN_SEQUENCES:REPORT_NO_MATCHES (1)',status='COMPLETED',exit='-',hash='aa/bb'))\n        w.writerow(dict(name='INTERPROSCAN:COMBINE:COMBINE_MATCHES (1)',status='FAILED' if case=='failed_combine_task' else 'COMPLETED',exit='1' if case=='failed_combine_task' else '-',hash='cc/dd'))\n    outputs=[]")
source=source.replace("    doc={'interproscan-version'", "    for result in outputs: result['md5']=result['md5'].upper()\n    if case=='wrong_digest': outputs[1]['md5']='0'*32\n    if case=='changed_sequence': outputs[1]['sequence']='MKVGA'\n    if case=='missing_matches_array': del outputs[1]['matches']\n    doc={'interproscan-version'")
source=source.replace("'md5':original[idx]['sequence_md5']", "'md5':original[idx]['sequence_md5'].upper()")
source=source.replace("['protein_hit',original[0]['sequence_md5']", "['protein_hit',original[0]['sequence_md5'].upper()")
source=source.replace("if case!='missing_nohit_evidence':", "if case not in ('missing_nohit_evidence','empty_combined_map','failed_combine_task'):")
source=source.replace("    proc=subprocess.run", "    if case in ('empty_combined_map','failed_combine_task'):\n        (run/'work/cc/dd/combined').mkdir(parents=True)\n        (run/'work/cc/dd/combined/1.json').write_text(json.dumps({original[1]['sequence_md5'].upper():{}}))\n    proc=subprocess.run")
source=source.replace("if case=='complete':", "if case in ('complete','empty_combined_map'):")
source=source.replace("here/'postprocessor_test_results.json'", "here/'final_verifier_test_results.json'")
exec(compile(source,str(base/'fixture_tests.py'),'exec'),dict(__file__=str(base/'fixture_tests.py')))

# Recheck the actual 66-record result with the same script. SLURM status is
# checked independently on Hummingbird; local fixtures do not call sacct.
run=base/'real_result_verification_final'
shutil.copytree(base/'real_result_verification_test',run,ignore=shutil.ignore_patterns('summary'))
(run/'summary').mkdir()
shutil.copy2(base/'postprocess_final_20261004.py',run/'scripts/postprocess.py')
combined=run/'work/95/75ef07b58b61ce737c5772581a2021/combined'
combined.mkdir(parents=True)
shutil.copy2(base/'combined_matches.json',combined/'1.json')
proc=subprocess.run([sys.executable,str(run/'scripts/postprocess.py')],capture_output=True,text=True)
assert proc.returncode==0,proc.stdout+proc.stderr
summary=json.loads((run/'summary/run_summary.json').read_text())
assert summary['verified_complete'] and summary['processed_records']==66
assert summary['completed_with_annotations']==51 and summary['completed_no_hits']==15
assert summary['workflow_tasks']==51 and summary['failed_workflow_tasks']==0
assert len(summary['no_match_evidence'])==15
print(json.dumps({'real_result_check':'passed','records':66,'with_matches':51,'no_hits':15,'successful_tasks':51},indent=2))
