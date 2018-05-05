import json
import pandas as pd
from collections import OrderedDict

particles = pd.read_csv(open("../../data/charmonia.csv", "r"), index_col=0)
data = OrderedDict()
for p in particles.index.values:
    data.update({p: {
                "limit": 1E-4,
                "rMax": 30
            }})
json.dump(data, open("c-scr.cfg", "w"), indent=4)
