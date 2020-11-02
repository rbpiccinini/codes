#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 00:11:05 2020

@author: piccinini
"""

import numpy as np

ne = 10
nt = 100
nx = 100
nw = 5
nv = 5

c = np.ones([ne, nw, nw])
d = np.ones([nt, ne])
q = np.ones([nw, nt, nx])

r = np.einsum('wtx,xe,ewv -> tv', q, d, c)