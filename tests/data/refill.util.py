import json
import pandas as pd
from collections import OrderedDict

df = pd.read_csv("c-lin.deng2017.csv")
df.drop(columns=["Gtot", "Gm1"], inplace=True)
df.columns = ["instate", "outstate", "frac_width"]

for i,row in df.iterrows():
    widths = input("{} -> {}: ".format(row["instate"], row["outstate"]))
    df.loc[i, "frac_width"] = widths

print(df)
df.to_csv("charm-trans-PDG.csv", index=False)

print("Done!")
