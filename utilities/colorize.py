import json

def print_unbuffered(string): 
  print(string, flush=True)

def print_json(data):
  formatted_json = json.dumps(data, indent=2)
  # https://stackoverflow.com/a/32166163/6674256
  from pygments import highlight, lexers, formatters
  # from pygments.styles import get_all_styles
  # print(list(get_all_styles()))
  colorful_json = highlight(formatted_json, lexers.JsonLexer(), formatters.Terminal256Formatter(style='monokai'))
  print_unbuffered(colorful_json)

from colorama import init as colorama_init
from colorama import Fore, Style
import sys

colorama_init(strip=False)

def print_string_as_error(text):
  print(Fore.RED + text + Style.RESET_ALL, file=sys.stderr)

def print_string_as_info(text):
  print(Fore.CYAN + text + Style.RESET_ALL, file=sys.stderr)

def print_string_as_info_dim(text):
  print(Fore.LIGHTBLACK_EX + text + Style.RESET_ALL, file=sys.stderr)
