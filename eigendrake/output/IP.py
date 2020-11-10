#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 21:52:39 2020

@author: piccinini
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

x = np.loadtxt('eigvalues_heterog.txt')

eig = x[1:,0]
psis = x[1:,1:].T

r = (psis/eig).T 

IP = 1e5/np.cumsum(r, axis=0)

# plt.semilogy(IP,'o', mfc='None', ms=3)
plt.plot(np.mean(IP, axis=1),'ko', mfc='None', ms=2)
plt.grid(1)