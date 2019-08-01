import gzip
import os
import requests
import sys
import urllib.request, urllib.parse, urllib.error

def log(message):
  print(message)
  sys.stdout.flush()

def download(url, path):
  if not os.path.exists(path):
    if requests.head(url).status_code != 200: return False
    urllib.request.urlretrieve(url, path)
  return True

def read(path):
  with open(path) as f:
    return f.read()

def read_gzip(path):
  with gzip.open(path) as f:
    return f.read()
