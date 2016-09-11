#!/usr/bin/env python

from pylab import *


file='Sevap_centerline.csv'
xn=loadtxt(file,delimiter=',',skiprows=1)[:,-3]
un=loadtxt(file,delimiter=',',skiprows=1)[:,11]

fig_width = 7.  # Get this from LaTeX using \showthe\columnwidth
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
rcParams.update({'figure.figsize': fig_size})


D=0.0098
fig = figure()
# plot data
ax = fig.add_subplot(1,1,1)
subplots_adjust( left=0.15, bottom=0.15 )
#line1,=ax.plot(ls.yc,ls.uc,'--b',linewidth=2)
ax.plot(xn/D,un,'ko-',markerfacecolor='None',markeredgecolor='k',markeredgewidth=1.5,label='numerical')
ylim([0,22])
xlim([0,45])

legend().draw_frame(False)

xlabel('x/D')
ylabel('Spray Evaporation Rate - Sm (kg/m3/s)')
savefig('centerline_Sevap.png')

show()
