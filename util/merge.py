#! /usr/bin/env python

import gzip
import json
import sys

ret = {
  "samples": []
}

def opn(fn):
  if fn.endswith('.gz'):
    return gzip.open(fn)
  return open(fn)

for file_name in sys.argv[1:]:
  with opn(file_name) as f:
    file_content = json.load(f)
  ret["samples"].extend(file_content["samples"])

print(json.dumps(ret))
