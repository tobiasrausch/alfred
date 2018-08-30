#! /usr/bin/env python

import gzip
from itertools import chain
import json
import sys

ret = {
  "samples": []
}

def opn(fn):
  if fn.endswith('.gz'):
    return gzip.open(fn)
  return open(fn)

def unique(iterable):
  ret = []
  for value in iterable:
    if value not in ret:
      ret.append(value)
  return ret

for file_name in sys.argv[1:]:
  with opn(file_name) as f:
    file_content = json.load(f)
  ret["samples"].extend(file_content["samples"])

summary_columns = unique(
  chain(
    *[s["summary"]["data"]["columns"] for s in ret["samples"]]
  )
)

for sample in ret["samples"]:
  old_columns = sample["summary"]["data"]["columns"]
  old_rows = [
    dict(zip(old_columns, row)) for row in sample["summary"]["data"]["rows"]
  ]
  sample["summary"]["data"]["columns"] = summary_columns
  sample["summary"]["data"]["rows"] = []
  for row in old_rows:
    new_row = []
    for column in summary_columns:
      new_row.append(row.get(column))
    sample["summary"]["data"]["rows"].append(new_row)
  
print(json.dumps(ret))
