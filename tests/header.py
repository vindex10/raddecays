from IPython import get_ipython
from IPython.terminal.pt_inputhooks import UnknownBackend
ipython = get_ipython()

try:
    ipython.magic("matplotlib inline")
except UnknownBackend:
    ipython.magic("matplotlib gtk3")

import pandas as pd
import numpy as np
import scipy as sp
from scipy import integrate, special, interpolate
from sympy.physics.quantum import cg
import json
import os
from collections import OrderedDict
from cycler import cycler
from IPython.display import display
import matplotlib.pyplot as plt
from functools import cmp_to_key
import re
from pylatex import Tabular, LongTable, Tabularx, MultiColumn, MultiRow, NoEscape

SMALL_SIZE = 14
MEDIUM_SIZE = 18
BIGGER_SIZE = 20

FIG_HSIZE = 8
FIG_VSIZE = 6

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE, figsize=(8, 6), dpi=100)  # fontsize of the figure title

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
                 for m in re.match("(\\\\?[a-zA-Z]+)(?:_\{?(\w+)\}?)?\((\w+)\)", name).groups()\
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
        pass # leave as is. example: h_{b}(3P)

    return "_".join(parts)

def stateHash(prefix, cfgname, pname):
    # hash [m, xL, xS, xJ]
    eigen_config = dict()
    eigen_config.update(json.load(open("../quarkEigen/output/{}.{}/config".format(prefix, cfgname), "r")))
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

def cfgnameSplit(cfgname):
    cfgtag = cfgname.split(".")
    cfgcode = "" if len(cfgtag) == 1 else ".".join(cfgtag[:-1])
    cfgtag = cfgtag[-1]
    return cfgcode, cfgtag

def readSpec(fname):
    res = {}
    with open(fname) as f:
        in_columns = False
        in_data = False
        for line in f:
            line = line.strip()
            if not line.startswith("#"):
                break
            if line == "" and in_data == True:
                in_data = False
                break

            nospace = line.replace(" ", "")
            if nospace == "#Columns:":
                in_columns = True
                continue

            if in_columns and nospace.startswith("#*"):
                in_data = True
                parts = line.split("::")
                colname = parts[0].replace(" ", "")[2:]
                tags = parts[1].strip().split(" ")
                descr = parts[-1].strip()
                res.update({colname:
                               {  "tags": tags
                                , "descr": descr
                               }
                           })
    return res

def dimTrans(indim, outdim):
    dimre = re.compile("^\[?([A-Z])")
    codes = {
          "": 1
        , "K": 10**3
        , "M": 10**6
        , "G": 10**9
        , "T": 10**12
        , "P": 10**15
    }

    try:
        incode = dimre.match(indim).group(1)
    except AttributeError:
        incode = ""

    try:
        outcode = dimre.match(outdim).group(1)
    except AttributeError:
        outcode = ""

    return codes[incode]/codes[outcode]
