import pandas as pd

data = pd.read_csv("c-lin.deng2017.csv", header=0, index_col=[0,1])/10**6

with open("c-lin.deng2017.csv", "w") as f:
    data.to_csv("c-lin.deng2017.csv")
