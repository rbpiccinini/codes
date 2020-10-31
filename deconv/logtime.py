# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 22:10:16 2020

@author: rbpic
"""

import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt

# eventos
te = np.array([2., 4., 5., 6., 7., 8., 10., 40.])
logt = np.log10(te)
dlogt = np.diff(logt).min()

t=[]
for i in range(1, len(logt)):
    n = int((logt[i] - logt[i-1])/dlogt*120)
    t.append([np.logspace(logt[i-1], logt[i], n)])

t = np.ravel(np.concatenate(t, axis=1))
t = np.unique(t)

with open('wells.inc','w') as fid:
    for tt in t:
        fid.write('TIME\t {:5.4f}\n'.format(tt))