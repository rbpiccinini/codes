#!/usr/bin/env python

from pylab import *

def adjust_spines(ax,spines):
    for loc, spine in ax.spines.iteritems():
        if loc in spines:
            spine.set_position(('outward',10)) # outward by 10 points
            spine.set_smart_bounds(True)
        else:
            spine.set_color('none') # don't draw spine

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])


fig_width = 7.  # Get this from LaTeX using \showthe\columnwidth
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
params = {'backend': 'ps',
          'axes.labelsize': 14,
          'text.fontsize': 14,
          'legend.fontsize': 13,
          'xtick.labelsize': 18,
          'ytick.labelsize': 18,
          'text.usetex': False,
          'font.family': 'monospace',
          'figure.figsize': fig_size}
rcParams.update(params)

D=0.0098
yc=loadtxt('cfd_k.csv',skiprows=1,delimiter=',')[:,0]/D
kc=loadtxt('cfd_k.csv',skiprows=1,delimiter=',')[:,1]

ye=loadtxt('up.csv',skiprows=1,delimiter=',')[:,0]
uue=loadtxt('up.csv',skiprows=1,delimiter=',')[:,1]
vve=loadtxt('vp.csv',skiprows=1,delimiter=',')[:,1]

ke=0.5*(uue**2+2.*vve**2)

fig = figure()
# plot data
ax = fig.add_subplot(1,1,1)
subplots_adjust( left=0.2, bottom=0.2 )
#line1,=ax.plot(ls.yc,ls.uc,'--b',linewidth=2)
Um=24.
line2,=ax.plot(yc,kc,'k-o')
line3,=ax.plot(ye,ke,'s-r')

legend(('numerical','experimental')).draw_frame(False)


# adjust the spines
adjust_spines(ax,['left','bottom'])

# disable clipping of data points by axes range
for artist in (line2,line3):
    artist.set_clip_on(True)

# adjust spine to be within ticks
ax.spines['bottom'].set_bounds(0,1)
ax.set_xlim([0,1])
ax.set_xticks(linspace(0,1,5))

xlabel('y/D')
ylabel('Turbulent KE (m2/s2)')
title('Axial Position - x/D = 0.5')

savefig('k_bc.png')
show()

