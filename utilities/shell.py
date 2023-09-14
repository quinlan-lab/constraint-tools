import subprocess

def shell(cmd):
  completed_process = subprocess.run(
    cmd,
    shell=True,
    executable='/usr/bin/bash',  # default shell is /bin/sh, but we need bash for <()
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
  ) 
  if completed_process.returncode == 0: # success
    print(completed_process.stderr.decode("utf-8").strip())
    return completed_process.stdout.decode("utf-8").strip()
  else: # failure
    print(completed_process.stdout.decode("utf-8").strip())
    raise Exception(completed_process.stderr.decode("utf-8").strip())
  
