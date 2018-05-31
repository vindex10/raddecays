from IPython import get_ipython
ipython = get_ipython()

ipython.magic("matplotlib inline")

import pandas as pd
import numpy as np
import scipy as sp
import json
import os
from collections import OrderedDict
from cycler import cycler
import matplotlib.pyplot as plt

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
