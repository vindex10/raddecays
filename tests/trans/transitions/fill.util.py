import json
import pandas as pd
from collections import OrderedDict

data = pd.read_csv("../../data/c-lin.deng2017.csv")

out = list()
for i,row in data.iterrows():
    out.append({"instate": row["instate"], "outstate": row["outstate"]})

full = {"system": "c-scr", "alphaEM": 0.003244, "spec": out}

json.dump(full, open("c-lin.csv", "w"), indent=4)
