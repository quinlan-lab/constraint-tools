import os
import subprocess
import multiprocessing

def get_hostname_process_cpu(): 
  hostname = os.uname()[1]
  pid = os.getpid()
  cpu = subprocess.run(
    ["ps", "-o", "psr=", "-p", f"{pid}"], 
    capture_output=True
  ).stdout.decode("utf-8").strip()
  return {
    'hostname': hostname,
    'process': pid, 
    'cpu': f'{cpu}/{multiprocessing.cpu_count()}'
  }

