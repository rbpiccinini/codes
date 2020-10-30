# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 08:11:03 2020

@author: bfv8
"""


import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.signal as sg
from deconv_multiwell_v2 import *


def which_type(well):
    if well == 'P1':
        return 'producer'
    elif well == 'I1':
        return 'injector'
    elif well == 'OBS':
        return 'observer'
    else:
        return 'none'

# Lê dados
# df = pd.read_excel('multiwell.xlsx', sheet_name='Data')
# df.to_csv('Ex1.zip', compression='zip', encoding='utf-8', index=False)
df = pd.read_csv('Ex1.zip', compression='zip')

idx = (df['t']>1.5) & (df['t']<=25.)
df = df[idx]
df['t'] = df['t'] - df['t'].min()

#resample df
# wells = df['well'].drop_duplicates().tolist()
# dfs = []
# for well in wells:
#     idx = df['well'] == well 
#     dfs.append((df.loc[idx]).iloc[::20 ,:])
# df = pd.concat(dfs)

# initial pressure [kgf/cm2]
p0 = 421.839630126953

# Pore volume
ct = 79.08e-6 # cm2/kgf
Vpct = 274047.*ct
psi0 = np.sqrt(1.0/Vpct)

# Cria poços
eig = np.array([0.0000, 0.1414, 0.3391, 1.1721, 3.1149, 9.6688, 25.4431,
                76.1631, 190.4340, 312.3359, 400.])
psi = np.zeros(len(eig))
psi[0] = psi0




wells = []
for well in df['well'].drop_duplicates():
    # filter for well
    idx = df['well'] == well
    # identify well type
    well_type = which_type(well)
    # create well history
    hist = classHist(t=df.loc[idx,'t'].to_numpy(),
                     q=df.loc[idx,'q'].to_numpy(),
                     p=df.loc[idx,'p'].to_numpy())
    # append well to list of wells
    wells.append(classWell(name=well,
                           well_type=well_type,
                           hist=hist,
                           psi=psi,
                           eig=eig,
                           p0=p0))




# create time 
t = df['t'].drop_duplicates().sort_values().to_numpy()

# convolve P1
p1 = wells[2]
i1 = wells[0]
obs1 = wells[1]


# Create deconvolution class
ex1 = classConv(wells, t=t, ne=len(eig), p0=p0)
D = ex1.D(eig, ex1.t)
Q = ex1.Q()

psis = np.ravel([well.psi for well in wells])
C = ex1.C(psis)
p = ex1.convolve(ex1.t, eig, psis)

# Read x0 for deconvolution
eig, psis = ex1.x2conv(np.loadtxt('r_opt7.txt'))

# Set bounds for eingevalues
lb = [0.01]+[0.1]*(ex1.ne-1)+[-10]*(ex1.ne-1)*ex1.nw
ub = [1.0]+[1000.]*(ex1.ne-1)+[10]*(ex1.ne-1)*ex1.nw
sqbounds = [lb,ub]
print(np.array(sqbounds).shape)

# Run deconvolution
# r = ex1.deconvolve(eig, psis, bounds=sqbounds)
# eig, psis = ex1.x2conv(r.x)
# np.savetxt('r.txt', r.x)

# Update wells
ex1.set_wells(eig, psis)


#plotar
# plt.plot(ex1.t, ex1.p, 'ro', ms=3, mfc='None') 
# plt.plot(ex1.t, ex1.convolve(ex1.t, eig, psis), 'kx', mfc='None', ms=2)

fig, axes = plt.subplots(3,1,figsize=(6.3, 6.3))
colors = ['b', 'k', 'g']
for well, color in zip(wells, colors):
    axes[0].fill_between(well.hist.t, well.hist.q, alpha=0.5, step='pre', color=color, label=well.name)
    axes[1].plot(well.hist.t, well.hist.p, '--', ms=2, mfc='None', mec=color)
    axes[1].plot(well.hist.t, well.pconv, '-', ms=2, color=color, label=well.name)
    axes[2].loglog(well.hist.t, well.bourdet(), 'o', ms=2, mfc='None', label=well.name, color=color)
    print('PI of '+well.name+' = '+str(well.PI([40.])))
    
    axes[0].set_ylabel('Flow rate [m3/day]')
    axes[1].set_ylabel('BHP [kgf/cm2]')
    axes[2].set_ylabel('Burdet derivative [kgf/cm2]')

for ax in axes:
    ax.set_xlabel('Time [days]')
    ax.legend(loc='lower right')
    
fig.tight_layout()
plt.savefig('multiwell_imex.png', dpi=600)
plt.savefig('multiwell_imex.svg')