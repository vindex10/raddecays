#!/usr/bin/env python

import json
import numpy as np
import os
from sys import argv
from collections import OrderedDict

try:
    path = argv[1]
except IndexError:
    print("setting dir to .")
    path = "."

pdata = json.load(open(os.path.join(path, "eigen_config"), "r"), object_pairs_hook=OrderedDict)

with open(os.path.join(path, "exclude"), "r") as f:
    exclude = f.read().split("\n")

for pname in pdata.keys():
    if pname in exclude:
        fin = os.path.join(path, "data", pname)

        qnums = np.array([pdata[pname]["eq"]["xJ"]\
                         ,pdata[pname]["eq"]["xL"]\
                         ,pdata[pname]["eq"]["xS"]], dtype=np.int)
        qnums = (qnums - 1)//2
        elevel = pname.split("_")[-1][:-1]
        fout = os.path.join(path, "data", "Rwf_{}_{}_{}({}).dat"\
                                              .format(*qnums, elevel))

        try:
            os.rename(fin, fout)
        except FileNotFoundError:
            print("no such file {}. skipping".format(fin))
