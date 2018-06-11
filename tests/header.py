from IPython import get_ipython
ipython = get_ipython()

ipython.magic("matplotlib inline")

import pandas as pd
import numpy as np
import scipy as sp
from scipy import integrate, special, interpolate
import json
import os
from collections import OrderedDict
from cycler import cycler
from IPython.display import display
import matplotlib.pyplot as plt
from functools import cmp_to_key
import re

SMALL_SIZE = 14
MEDIUM_SIZE = 18
BIGGER_SIZE = 20

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE, figsize=(12, 8))  # fontsize of the figure title

pd.set_option('precision', 10)
pd.options.display.float_format = "{:.3g}".format

def code2name(code):
    parts = code.split("_")
    
    vocab = {"yps": r"\Upsilon"
            ,"eta": r"\eta"
            ,"chi": r"\chi"
            ,"psi": r"\psi"}
    try:
        parts[0] = vocab[parts[0]]
    except KeyError:
        pass

    try:
        return "{0}_{{{1}}}({2})".format(*parts)
    except IndexError:
        return "{0}({1})".format(*parts)
    except:
        return code

def name2code(name):
    try:
        parts = [m \
                 for m in re.match("(\\\\\w+)(?:_\{(\w+)\})?\((\w+)\)", name).groups()\
                 if m is not None]
    except AttributeError:
        return name

    vocab = {r"\Upsilon":   "yps"
            ,r"\eta":       "eta"
            ,r"\chi":       "chi"
            ,r"\psi":       "psi"}
    try:
        parts[0] = vocab[parts[0]]
    except KeyError:
        return name

    return "_".join(parts)

def stateHash(prefix, cfgcode, pname):
    # hash [m, xL, xS, xJ]
    eigen_config = dict()
    eigen_config.update(json.load(open("../quarkEigen/output/{}.{}/config".format(prefix, cfgcode), "r")))
    eigen_config = eigen_config[pname]["eq"]
    return [int(pname.split("_")[-1][:-1])
           ,eigen_config["xL"]
           ,eigen_config["xS"]
           ,eigen_config["xJ"]]

def cmpStatesByCode(prefix, cfgcode, p1, p2):
    res = 0
    for v1, v2 in zip(stateHash(prefix, cfgcode, p1), stateHash(prefix, cfgcode, p2)):
        if abs(v2-v1) < 1E-5:
            continue
        if v1 < v2:
            res = -1
            break
        elif v1 > v2:
            res = 1
            break
    return res

def cmpMultyStatesByCode(prefix, cfgcode, p1, p2):
    cmp = 0
    try:
        i = 0
        while cmp == 0:
            cmp = cmpStatesByCode(prefix, cfgcode, p1[i], p2[i])
            i += 1
    except IndexError:
        pass
    return cmp

def cmpStatesByName(prefix, cfgcode, p1, p2):
    return cmpStatesByCode(prefix, cfgcode, name2code(p1), name2code(p2))

def cmpMultyStatesByName(prefix, cfgcode, p1, p2):
    return cmpMultyStatesByCode(prefix, cfgcode\
                               ,[name2code(p) for p in p1]
                               ,[name2code(p) for p in p2])

def dfsort(df, func):
    keys = sorted(df.index.values, key=cmp_to_key(func))
    return df.loc[keys]
