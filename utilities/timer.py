import functools 
import time 
from colorize import print_string_as_info

# https://realpython.com/primer-on-python-decorators/
def timer(func):
  @functools.wraps(func)
  def wrapper(*args, **kwargs):
    start_time = time.time()
    result = func(*args, **kwargs)
    end_time = time.time()
    print_string_as_info(f'timer: {func.__name__} => {end_time-start_time:.2f} s')
    return result
  return wrapper

