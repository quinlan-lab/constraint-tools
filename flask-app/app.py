import argparse 
import color_traceback
from flask import Flask, request, make_response
from flask_cors import CORS # required for development of vue app

from compute_mutation_counts import compute_mutation_counts
from colorize import print_string_as_info, print_json

def parse_arguments(): 
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('--model', type=str, help='')
  parser.add_argument('--port', type=int, help='')
  return parser.parse_args()

app = Flask(__name__, static_folder='static', static_url_path="/static") # WSGI app

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

@app.route('/api', methods=['POST'])
def serve_api():
  if not request.is_json: 
    response = make_response({ 
      'message': 'Please format the request payload in json format'
    })
    response.status_code = 400 # https://developer.mozilla.org/en-US/docs/Web/HTTP/Status/400
    return response
  # returning a dictionary makes flask respond with: 
  # "Content-Type: application/json" 
  # https://developer.mozilla.org/en-US/docs/Web/HTTP/Headers/Content-Type
  return compute_mutation_counts(
    request.json['region'],
    parse_arguments().model,
    int(request.json['windowSize']),
    int(request.json['windowStride'])
  )

def print_app_info(): 
  print_string_as_info(
    f'if index.html is located in the directory:\n'
    f'{app.static_folder}\n' 
    f'then it is available at the url:\n'
    f'{app.static_url_path}/index.html'
  )

def run_app():
  # `debug=True` is passed to enable the debugger and reloader:
  app.run(debug=True, port=parse_arguments().port) # WSGI server

if __name__ == "__main__":
  print_app_info() 
  run_app()



