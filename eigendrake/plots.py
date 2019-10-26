import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

x = pd.read_csv('solution.dat')
x.set_index('t (h)', inplace=True)

q = -x['q (m3/d)'].min()

x['dp'] = np.gradient(x['pwf (bar)'], x.index.values)

wf = x.loc[0:48].copy()
ws = x.loc[48:].copy()

wf.reset_index(inplace=True)
ws.reset_index(inplace=True)

ws['t (h)'] = ws['t (h)'] - 48.

pwf = wf['pwf (bar)'].values[-1]

# estática
plt.figure()
plt.loglog(ws['t (h)'], (ws['pwf (bar)']-pwf)/q, 'ob')
plt.loglog(ws['t (h)'], ws['t (h)']*ws['dp']/q, 'or')
plt.grid(1)

# fluxo 
plt.figure()
plt.loglog(wf['t (h)'], (400.-wf['pwf (bar)'])/q, 'ob')
plt.loglog(wf['t (h)'], -wf['t (h)']*wf['dp']/q, 'or')
plt.grid(1)

plt.show()
