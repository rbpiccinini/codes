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
from deconv_multiwell import *


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
# df = pd.read_excel('multiwell.xlsx')
# df.to_csv('Ex1.zip', compression='zip', encoding='utf-8', index=False)
df = pd.read_csv('Ex1.zip', compression='zip')

idx = (df['t']>1.5) & (df['t']<=15.)
df = df[idx]
df['t'] = df['t'] - df['t'].min()

#resample df
# wells = df['well'].drop_duplicates().tolist()
# dfs = []
# for well in wells:
#     idx = df['well'] == well 
#     dfs.append((df.loc[idx]).iloc[::20 ,:])
# df = pd.concat(dfs)


# Cria poços
eig = np.array([0.0000, 0.1414, 0.3391, 1.1721, 3.1149, 9.6688, 25.4431,
                76.1631, 190.4340, 312.3359, 932.4149])

eig = np.linspace(0., 15, 20)**2
psi = 0.001*np.ones(len(eig))

# initial pressure [kgf/cm2]
p0 = 421.839630126953

# Pore volume
ct = 79.08e-6 # cm2/kgf
Vpct = 274047.*ct
psi0 = np.sqrt(1.0/Vpct)

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


# deconvolution
# from deconv_multiwell import *
spe10 = classDeconv(wells)
x0 = np.loadtxt('x0.txt')
spe10.x0_to_wells(x0)

# # # solve deconvolution
# r = spe10.deconv()
# np.savetxt('r.txt', r.x)
# spe10.x0_to_wells(r.x)


# plots convolucao
fig, axes = plt.subplots(3,1,figsize=(6.3, 6.3))
for well, ax in zip(wells, axes):
    ax.plot(well.hist.t, p0-well.convolve(wells), 'ko', ms=3, mfc='None', label='deconv')
    ax.plot(well.hist.t, well.hist.p, 'rx', ms=2, mfc='None', label='hist')

    # ax.plot(well.hist.t, p0-well.convolve(), 'k--', ms=1, mfc='None')
    ax.set_xlabel('Time [d]')
    ax.set_ylabel('Pressure [kgf/cm²]')
    ax.set_title(well.name)
    ax.legend(loc='lower right')

    twin = ax.twinx()
    twin.fill_between(well.hist.t, well.hist.q, alpha=0.35, step='pre')
    # twin.set_xlim([0., 500.])
    twin.set_ylim([-300., 300.])

    twin.set_ylabel('Flow rate [m³/days]')
    twin.grid(0)

fig.tight_layout()

# Bourdet
fig, ax = plt.subplots(1,1,figsize=(6.3, 6.3))
for well in wells:
    ax.loglog(well.hist.t, well.bourdet(), 'o', ms=2, mfc='None', label=well.name)
ax.set_xlabel('Time [d]')
ax.set_ylabel('Bourdet derivative [kgf/cm²]')
ax.legend(loc='lower right')
ax.set_xlim([1e-3, 1e2])
ax.set_aspect('equal', 'datalim')

fig.tight_layout()