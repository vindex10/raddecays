#!/usr/bin/env python
from collections import OrderedDict
import json
import sys

out = OrderedDict()
for fname in sys.argv[1:]:
    with open(fname, "r") as f:
        out.update(\
                    json.load(f, object_pairs_hook=OrderedDict)\
                  )
print(json.dumps(out,indent=4))
