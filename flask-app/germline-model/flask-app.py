import argparse 
import color_traceback
from flask import Flask, request, make_response
from flask_cors import CORS # required for development of vue app
import pyranges as pr

from compute_expected_observed_counts import (
  compute_expected_observed_counts,
  fetch_distribution_N
)
from colorize import print_string_as_error, print_string_as_info, print_json
from read_model import read_model
from pack_unpack import unpack 

def parse_arguments(): 
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('--model', type=str, help='')
  parser.add_argument('--port', type=int, help='')
  parser.add_argument('--region', type=str, help='')
  parser.add_argument('--window-stride', type=int, help='', dest='window_stride')
  return parser.parse_args()

app = Flask(__name__, static_folder='static', static_url_path="/static") # WSGI app

args = parse_arguments()
model = read_model(args.model)

# https://biocore-ntnu.github.io/pyranges/loadingcreating-pyranges.html
neutral_regions = pr.read_bed(model['neutralRegions'])

# https://developer.mozilla.org/en-US/docs/Web/HTTP/CORS
# this line of code can be removed once the vue.js app is served from the same port as the flask app: 
CORS(app) # required for development of vue app

# Notes on serving SPAs:
# https://flask.palletsprojects.com/en/2.0.x/patterns/singlepageapplications/
# https://github.com/petermchale/igv-grid/blob/98687612e25560caaffd3802905bd2bdcaefbb9d/launch-web-app.sh#L15

@app.route('/')
def serve_index_static_file():
  # https://flask.palletsprojects.com/en/2.0.x/api/#flask.Flask.send_static_file
  return app.send_static_file('index.html')

@app.route('/<path:path>')
def serve_other_static_file(path):
  # https://flask.palletsprojects.com/en/2.0.x/api/#flask.Flask.send_static_file
  return app.send_static_file(path)

def bad_request(): 
  response = make_response({ 
    'message': 'Please format the request payload in json format'
  })
  response.status_code = 400 # https://developer.mozilla.org/en-US/docs/Web/HTTP/Status/400
  return response

def internal_server_error(): 
  import traceback
  for line in traceback.format_exc().splitlines(): 
    print_string_as_error(line)
  response = make_response({ 
    'exceptionTraceback': traceback.format_exc().splitlines()
  })
  response.status_code = 500 # https://developer.mozilla.org/en-US/docs/Web/HTTP/Status/500
  return response

@app.route('/api/distribution-n', methods=['POST'])
def serve_api_distribution_N():
  if not request.is_json: 
    return bad_request()
  try: 
    # returning a dictionary makes flask respond with: 
    # "Content-Type: application/json" 
    # https://developer.mozilla.org/en-US/docs/Web/HTTP/Headers/Content-Type
    return fetch_distribution_N(window=request.json, model=model)
  except Exception: 
    return internal_server_error() 

@app.route('/api/expected-observed-counts', methods=['POST'])
def serve_api_expected_observed_counts():
  if not request.is_json: 
    return bad_request()
  try: 
    # returning a dictionary makes flask respond with: 
    # "Content-Type: application/json" 
    # https://developer.mozilla.org/en-US/docs/Web/HTTP/Headers/Content-Type
    return compute_expected_observed_counts(
      request.json['region'],
      model,
      int(request.json['windowStride']),
      number_chromosomes_min=model['numberChromosomesMin'],
      log = False
    )
  except Exception: 
    return internal_server_error() 

@app.route('/api/initial-plot-parameters', methods=['GET'])
def serve_api_initial_plot_parameters():
  # returning a dictionary makes flask respond with: 
  # "Content-Type: application/json" 
  # https://developer.mozilla.org/en-US/docs/Web/HTTP/Headers/Content-Type
  return { 
    'region': args.region,
    'windowStride': args.window_stride,
  }

@app.route('/api/model-parameters', methods=['GET'])
def serve_api_model_parameters():
  # returning a dictionary makes flask respond with: 
  # "Content-Type: application/json" 
  # https://developer.mozilla.org/en-US/docs/Web/HTTP/Headers/Content-Type
  return {
    'genomeBuild': model['build'], 
    'kmerSize': model['kmerSize'],
    'numberChromosomesMin': model['numberChromosomesMin'],
    'windowSize': model['windowSize']
  }

def get_neutral_regions(region): 
  chromosome, start, end = unpack(region)
  # https://biocore-ntnu.github.io/pyranges/manipulating-the-data-in-pyranges.html
  return neutral_regions[chromosome, start:end].df.to_dict(orient='records')

@app.route('/api/neutral-regions', methods=['POST'])
def serve_api_neutral_regions():
  if not request.is_json: 
    return bad_request()
  try: 
    # returning a dictionary makes flask respond with: 
    # "Content-Type: application/json" 
    # https://developer.mozilla.org/en-US/docs/Web/HTTP/Headers/Content-Type
    return { 
      'neutralRegions': get_neutral_regions(request.json['region'])
    }
  except Exception: 
    return internal_server_error() 

def print_app_info(): 
  # print_string_as_info(
  #   f'if index.html is located in the directory:\n'
  #   f'{app.static_folder}\n' 
  #   f'then it is available at the url:\n'
  #   f'{app.static_url_path}/index.html'
  # )
  print_string_as_info(
    f'Visit http://localhost:{args.port} in your web browser to view the web app.\n\n'
    f'If the software is run on a remote machine,\n'
    f'then to view the web-app on a local machine,\n'
    f'you may need to do the following on your local machine\n'
    f'(change username and host first):\n\n'
    f'ssh -N -L localhost:{args.port}:localhost:{args.port} username@host\n\n'
  )

def run_app():
  # pass `debug=True` to enable the debugger and reloader:
  app.run(debug=True, port=args.port) # WSGI server

if __name__ == "__main__":
  print_app_info() 
  run_app()



