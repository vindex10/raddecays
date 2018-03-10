#!/usr/bin/env python
from collections import OrderedDict
import json
import sys

out = OrderedDict()
for fname in sys.argv[1:]:
    out.update(\
                json.load(open(fname, "r"), object_pairs_hook=OrderedDict)
            )
print(json.dumps(out,indent=4))
