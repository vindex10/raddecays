import json
import pandas as pd
from collections import OrderedDict
import os

# particles = pd.read_csv(os.path.join(os.path.dirname(os.path.realpath(__file__)), "../../data/charmonia.csv"), index_col=0)
data = json.load(open("b-scr.cfg", "r"), object_pairs_hook=OrderedDict)

for name in data.keys():
    # data[name]["eq"]["env"]["alphaS"] = 0.5461
    # data[name]["eq"]["env"]["b"] = 0.1425
    # data[name]["eq"]["env"]["mC"] = 1.4830
    # data[name]["eq"]["env"]["sigma"] = 1.1384
    # data[name]["eq"]["env"]["rC"] = 0.202/0.19732697
    # data[name]["eq"]["E"] = particles["LP"].loc[name]/1000 - 2*data[name]["eq"]["env"]["mC"] 
    # data[name]["estep"] = 0.0001
    # data[name]["damping"] = 10
    data[name]["shrink"] = 0.4
    # data[name]["intstep"] = 0.001
    # data[name]["ewindow"] = 0.1
    # data[name]["cutscales"] = [37, 40, 43]
    # del data[name]["eq"]["env"]["mu"]
    # del data[name]["intabsTol"]
    # del data[name]["intrelTol"]
json.dump(data, open("b-scr.cfg", "w"), indent=4)
