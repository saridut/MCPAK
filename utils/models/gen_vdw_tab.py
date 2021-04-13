#!/usr/bin/env python

'''
Generates tabulated vdw potential.
'''

import math
import numpy as np
import matplotlib as mpl

#######################################################################
###FIGURE
fig_prop = {'figsize'        : (4, 3),
            'facecolor'      : '0.75',
            'edgecolor'      : 'white', 
            'subplot.left'   : 0.125,
            'subplot.right'  : 0.9, 
            'subplot.bottom' : 0.1,  
            'subplot.top'    : 0.9, 
            'autolayout'     : True
            }

###AXES
axes_prop = {'axisbelow'     : True, #'line',
             'unicode_minus' : False,
             'facecolor'     : 'white',
             'edgecolor'     : 'black',
             'linewidth'     : 0.75,
             'grid'          : False
             }

###LINES
lines_prop = {'linewidth'       : 0.8,
             'linestyle'       : '-',
             'color'           : 'black',
             'marker'          : None,
             'markeredgewidth' : 0.5,
             'markersize'      : 3,
             'dash_joinstyle'  : 'miter',
             'dash_capstyle'   : 'butt',
             'solid_joinstyle' : 'miter',
             'solid_capstyle'  : 'projecting',
             'antialiased'     : True
             }

###PATCH
patch_prop = {'linewidth'    : 1.0,
              'facecolor'    : 'blue',
              'edgecolor'    : 'black',
              'antialiased'  : True
              }

###TICK
xtick_prop = {'major.size' : 8, 
              'minor.size' : 4,
              'major.width': 0.4,
              'minor.width' : 0.4,
              'labelsize'  : 'small',
              'direction': 'in'
              }


ytick_prop = {'major.size' : 8, 
              'minor.size' : 4,
              'major.width': 0.4,
              'minor.width' : 0.4,
              'labelsize'  : 'small',
              'direction': 'in'
              }
###GRIDS
grid_prop = {'color'     : 'b0b0b0',
             'linestyle' : '-',
             'linewidth' : 0.4, # in points
             'alpha'     : 0.6  # transparency, between 0.0 and 1.0
             }

###LEGEND
leg_prop = {'fancybox'     : False,
            'numpoints'    : 1,    
            'fontsize'     : 'x-small',
            'borderpad'    : 0.5, 
            'markerscale'  : 1.0, 
            'labelspacing' : 0.5,
            'handlelength' : 2., 
            'handleheight' : 0.7,
            'handletextpad': 0.8,
            'borderaxespad': 0.5,
            'columnspacing': 2.,
            'frameon'      : False
            }

###FONT
#font_prop = {'family' : 'serif',
#              'serif' : 'Times',
#              'weight': 'medium',
#              'size'  : 10
#              }
              
###TEXT & LATEX
text_prop = {'usetex': True,
            'latex.preamble': r'\usepackage{amssymb}'+
                r'\usepackage{amsmath}'+
                r'\usepackage{sansmathfonts}'+
                r'\usepackage[T1]{fontenc}'
                 }

###PS & EPS BACKEND
ps_prop = {'useafm': True}

mpl.rc('figure', **fig_prop  )
mpl.rc('axes',   **axes_prop )
mpl.rc('lines',  **lines_prop)
mpl.rc('patch',  **patch_prop)
mpl.rc('xtick',  **xtick_prop)
mpl.rc('ytick',  **ytick_prop)
mpl.rc('grid',  **grid_prop)
mpl.rc('legend', **leg_prop  )
#mpl.rc('font',   **font_prop )
mpl.rc('text',   **text_prop )
mpl.rc('ps',     **ps_prop   )

#######################################################################
mpl.use('PDF')
import matplotlib.pyplot as plt

def vdw_gauss(r, A, B, rcut):
    '''
    Cut and shifted Gaussian potential.

    '''
    pot_rcut = A*math.exp(-B*rcut**2)

    pot = 0.0
    if (r < rcut):
        pot = A*math.exp(-B*r**2) - pot_rcut

    return pot


def vdw_mie(r, eps, sigma, n, m, rcut):
    '''
    Cut and shifted n-m Lennard Jones potential.

    '''
    assert n > m
    pot_pf = eps*(n/(n-m))*(n/m)**(m/(n-m))
    sir = sigma/rcut
    sirn = sir**int(n); sirm = sir**int(m)
    pot_rcut = pot_pf*(sirn - sirm)

    pot = 0.0
    if (r < rcut):
        sir = sigma/r
        sirn = sir**int(n); sirm = sir**int(m)
        pot = pot_pf*(sirn - sirm) - pot_rcut
    return pot


#LJ potential
#eps = 1.0; sigma = 2.0; n = 12; m = 6
#rcut = sigma*(n/m)**(1.0/(n-m))
#print(rcut)
#rmin = 0.8*sigma
#rmax = rcut
#
#rvals = np.logspace(math.log10(rmin), math.log10(rmax), 100)
#q = len(rvals)
#Uvals = np.zeros((q,))
#
#for i in range(q):
#    r = rvals[i]
#    Uvals[i] = vdw_mie(r, eps, sigma, n, m, rcut)

#Gaussian potential
A = 2.0; B = 2.0; rcut = 2.245
rmin = 0.0
rmax = rcut

rvals = np.linspace(rmin, rmax, 100)
q = len(rvals)
Uvals = np.zeros((q,))

for i in range(q):
    r = rvals[i]
    Uvals[i] = vdw_gauss(r, A, B, rcut)

#Plotting
figh, axh = plt.subplots(nrows=1, ncols=1)
axh.plot(rvals, Uvals, ls='-', marker='None', label=r'$U$')
#axh.axvline(x=rcut)
axh.grid(True)
axh.tick_params(which='both', top=True, right=True, color='0.5')
axh.set_yscale('linear')
axh.set_xscale('linear')
legend = axh.legend(loc='best', ncol=1)
#axh.set_xlim(1, 1.2)
axh.set_ylim(-0.5, 5.25)

axh.set_xlabel(r'$r$', labelpad=5.0)
axh.set_ylabel(r"$U, U'$", labelpad=5.0)
plt.savefig('fig-pot.pdf', transparent=True)
plt.close()

#Writing out the potential
#with open('vdw_tab_lj.txt', 'w') as fh:
#    fh.write('#Tabulated LJ potential\n')
#    fh.write('\n')
#    fh.write('%d\n'%q)
#    for i in range(q):
#        fh.write('%f  %f\n'%(rvals[i], Uvals[i]))

with open('vdw_tab_gauss.txt', 'w') as fh:
    fh.write('#Tabulated Gaussian potential\n')
    fh.write('\n')
    fh.write('%d\n'%q)
    for i in range(q):
        fh.write('%f  %f\n'%(rvals[i], Uvals[i]))
