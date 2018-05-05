import json
import pandas as pd
from collections import OrderedDict

df = pd.read_csv("c.deng2017.csv")
data = list()

while True:
    instate = input("\ninstate: ")
    if instate == "":
        break

    while True:
        outstate = input("{} -> outstate:".format(instate))
        if outstate == "":
            break

        widths = input("widths(comma-sep):")
        widths = widths.split(",")
        data.append({"instate": instate, "outstate": outstate, "Gm1": widths[0], "Ge1": widths[1], "Gtot": widths[2]})

if len(data) > 0:
    df = df.append(data, ignore_index=True)
print(df)
df.to_csv("c.deng2017.new.csv", index=False)

print("Done!")
