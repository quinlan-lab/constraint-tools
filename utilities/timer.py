import functools 
import time 
from colorize import print_string_as_info, print_unbuffered

# https://realpython.com/primer-on-python-decorators/
def timer(func):
  @functools.wraps(func)
  def wrapper(*args, **kwargs):
    start_time = time.time()
    result = func(*args, **kwargs)
    end_time = time.time()
    print_string_as_info(f'timer: {func.__name__} => {end_time-start_time:.2f} s')
    print_unbuffered('')
    print_unbuffered('*************************************************')
    print_unbuffered('')
    return result
  return wrapper

