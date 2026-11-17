#!/usr/bin/env python3
import datetime,json,pathlib,re,shlex,subprocess
run=pathlib.Path(__file__).resolve().parent.parent
def submit(role,args):
 command=['sbatch','--parsable']+args+[str(run/'scripts'/(role+'.sbatch'))]
 marker=run/'provenance'/(role+'_submission_command.txt')
 with marker.open('x') as f:
  f.write(shlex.join(command)+'\n'); f.flush()
 result=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
 with (run/'provenance'/(role+'_submission_receipt.json')).open('x') as f:
  json.dump(dict(command=command,returncode=result.returncode,stdout=result.stdout,stderr=result.stderr,time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f,indent=2)
 if result.returncode!=0:
  raise RuntimeError(role+' submission failed: '+result.stderr)
 jobid=result.stdout.strip().split(';')[0]
 if not re.fullmatch(r'[0-9]+',jobid):
  raise RuntimeError('Unclear submission acknowledgement; inspect squeue/sacct before retrying: '+result.stdout)
 with (run/'provenance'/(role+'_job_id.txt')).open('x') as f:
  f.write(jobid+'\n'); f.flush()
 print(role+' job ID: '+jobid,flush=True)
 return jobid
analysis=submit('analysis',[])
submit('postprocess',['--dependency=afterok:'+analysis,'--kill-on-invalid-dep=yes'])
