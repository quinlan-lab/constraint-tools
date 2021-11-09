import yaml

with open('conda-environment.yml', 'r') as f:
  try:
    conda_environment = yaml.safe_load(f)
    print(conda_environment['name'])
  except yaml.YAMLError as yaml_error:
    print(yaml_error)
