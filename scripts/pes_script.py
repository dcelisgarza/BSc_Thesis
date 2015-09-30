"""
PES Plotting Script.
Author: Daniel Celis Garza
"""

import numpy as np
import matplotlib
matplotlib.use('Agg') # Lets us save in .svg format.
import matplotlib.pyplot as plt
from matplotlib import ticker

plt.close('all') # Close all figures when replotting.
# Plotting parameters for use at scale = 0.5 with the thesis template.
params = {'backend': 'ps',
          'axes.labelsize': 24,
          'axes.titlesize': 24,
          'text.fontsize': 24,
          'legend.fontsize': 24,
          'xtick.labelsize': 24,
          'ytick.labelsize': 24,
          'text.usetex': True, # LaTeX Masterrace Huuuh haaaah!
          'font.family': 'serif',
          'mathtext.fontset': 'custom'}
plt.rcParams.update(params)
# LaTeX font activated by " r'' ". Accepts all environments.

dpt = '../data/pes/'
ipt = '../images/'
ext = 'svg'

# Plot diabatic PES.
def plotd(data,legend,b,xlim,name,ext):
    x, y0, y1, y2 = data[:,0], data[:,1], data[:,2], data[:,4]
    # Creating figure.
    fig, ax = plt.subplots(figsize=(12,7.42)) # Golden ratio ;)
    # Plotting PES.
    ax.plot(x, y0,label=legend[0], ls = '-', lw = 2.5) # Continuous line
    ax.plot(x, y1,label=legend[1], ls = '--', dashes=[8, 4, 2, 4, 2, 4], lw = 2.5) # Custom dashes.
    ax.plot(x, y2,label=legend[2], ls = '--', lw = 2.5) # Dashed line.
    # Adding grid.
    ax.grid(b=True, which='major', color='0.5',linestyle=':') # Dotted line.
    # Axis labels.
    ax.set_xlabel(r'Position ($\mathbf{R}$)')
    ax.set_ylabel(r'Energy (A.U.)')
    # Scientific notation on y-axis for the win.
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True) 
    formatter.set_powerlimits((-1,1)) 
    ax.yaxis.set_major_formatter(formatter)
    # Minor and major tick locators.
    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=b[0])) # Major x-tick every 1.
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(n=4)) # 4 Minor x-tick per major tick.
    ax.yaxis.set_major_locator(ticker.MultipleLocator(base=b[1])) # Major y-tick every 0.002.
    # x-range
    ax.set_xlim(xlim)
    # Legend options (smart legend positioning, deciding whether to remove frame or not).
    ax.legend(loc=0,frameon=False,fancybox=True,framealpha=1)
    # Adjusting the plot so the bottom legend fits into the figure.
    plt.subplots_adjust(bottom=0.12)
    # Saving into a .svg so it can be converted into anything else with Inkscape.
    plt.savefig(name+'.'+ext, format=ext, dpi=1200)

# Plot adiabatic PES.
def plota(data,legend,b,xlim,name,ext):
    x, y0, y1 = data[:,0], data[:,5], data[:,6]
    # Creating figure.
    fig, ax = plt.subplots(figsize=(12,7.42)) # Golden ratio ;)
    # Plotting PES.
    ax.plot(x, y0,label=legend[0], ls = '-', c = 'b', lw = 2.5) # Continuous line
    # Custom dashes.
    ax.plot(x, y1,label=legend[1], ls = '--', c = 'r', dashes=[8, 4, 2, 4, 2, 4], lw = 2.5)
    # Adding grid.
    ax.grid(b=True, which='major', color='0.5',linestyle=':') # Dotted line.
    # Axis labels.
    ax.set_xlabel(r'Position ($\mathbf{R}$)')
    ax.set_ylabel(r'Energy (A.U.)')
    # Scientific notation on y-axis for the win.
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True) 
    formatter.set_powerlimits((-1,1)) 
    ax.yaxis.set_major_formatter(formatter)
    # Minor and major tick locators.
    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=b[0])) # Major x-tick every 1.
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(n=4)) # 4 Minor x-tick per major tick.
    ax.yaxis.set_major_locator(ticker.MultipleLocator(base=b[1])) # Major y-tick every 0.002.
    # x-range
    ax.set_xlim(xlim)
    # Legend options (smart legend positioning, deciding whether to remove frame or not).
    ax.legend(loc=0,frameon=False,fancybox=True,framealpha=1)
    # Adjusting the plot so the bottom legend fits into the figure.
    plt.subplots_adjust(bottom=0.12)
    # Saving into a .svg so it can be converted into anything else with Inkscape.
    plt.savefig(name+'.'+ext, format=ext, dpi=1200)

# Plot E2 - E1
def plotad(data,legend,b,xlim,name,ext):
    x, y0, y1 = data[:,0], data[:,5], data[:,6]
    # Creating figure.
    fig, ax = plt.subplots(figsize=(12,7.42)) # Golden ratio ;)
    # Plotting PES.
    ax.plot(x, y1-y0,label=legend[0], ls = '-', c = 'b', lw = 2.5) # Continuous line
    # Adding grid.
    ax.grid(b=True, which='major', color='0.5',linestyle=':') # Dotted line.
    # Axis labels.
    ax.set_xlabel(r'Position ($\mathbf{R}$)')
    ax.set_ylabel(r'$\Delta$ Energy (A.U.)')
    # Scientific notation on y-axis for the win.
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True) 
    formatter.set_powerlimits((-1,1)) 
    ax.yaxis.set_major_formatter(formatter)
    # Minor and major tick locators.
    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=b[0])) # Major x-tick every 1.
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(n=4)) # 4 Minor x-tick per major tick.
    ax.yaxis.set_major_locator(ticker.MultipleLocator(base=b[1])) # Major y-tick every 0.002.
    # x-range
    ax.set_xlim(xlim)
    # Legend options (smart legend positioning, deciding whether to remove frame or not).
    ax.legend(loc=0,frameon=False,fancybox=True,framealpha=1)
    # Adjusting the plot so the bottom legend fits into the figure.
    plt.subplots_adjust(bottom=0.12)
    # Saving into a .svg so it can be converted into anything else with Inkscape.
    plt.savefig(name+'.'+ext, format=ext, dpi=1200)

# Single Avoided Crossing Diabatic
data   = np.loadtxt(dpt+'sac.dat')
legend = [r'$ H_{11} $',r'$ H_{12} = H_{21} $',r'$ H_{22} $']
b      = [1,0.002]
xlim   = [-4,4]
name   = ipt+'scpes'
plotd(data,legend,b,xlim,name,ext)

# Single Avoided Crossing Adiabatic
legend = [r'$ E_{1} $',r'$ E_{2}']
b      = [1,0.002]
xlim   = [-4,4]
name   = ipt+'ascpes'
plota(data,legend,b,xlim,name,ext)

# Single Avoided Crossing E2 - E1
legend = [r'$ E_{2} - E_{1} $']
b      = [1,0.002]
xlim   = [-4,4]
name   = ipt+'del_ascpes'
plotad(data,legend,b,xlim,name,ext)

# Single Avoided Crossing Diabatic Derivatives
data = np.loadtxt(dpt+'dsac.dat')
legend = [r'$ \frac{\partial H_{11}}{\partial R} $',
          r'$ \frac{\partial H_{12}}{\partial R} = \frac{\partial H_{21}}{\partial R} $',
          r'$ \frac{\partial H_{22}}{\partial R} $']
b      = [1,0.005]
xlim   = [-4,4]
name   = ipt+'dscpes'
plotd(data,legend,b,xlim,name,ext)

# ----------------------------------------- #

# Double Avoided Crossing Diabatic
data = np.loadtxt(dpt+'dac.dat')
legend = [r'$ H_{11} $',r'$ H_{12} = H_{21} $',r'$ H_{22} $']
b      = [2,0.01]
xlim   = [-8,8]
name   = ipt+'dcpes'
plotd(data,legend,b,xlim,name,ext)

# Double Avoided Crossing Adiabatic
legend = [r'$ E_{1} $',r'$ E_{2}']
b      = [2,0.01]
xlim   = [-8,8]
name   = ipt+'adcpes'
plota(data,legend,b,xlim,name,ext)

# Double Avoided Crossing E2 - E1
legend = [r'$ E_{2} - E_{1} $']
b      = [1,0.005]
xlim   = [-8,8]
name   = ipt+'del_adcpes'
plotad(data,legend,b,xlim,name,ext)

# Double Avoided Crossing Diabatic Derivatives
data = np.loadtxt(dpt+'ddac.dat')
legend = [r'$ \frac{\partial H_{11}}{\partial R} $',
          r'$ \frac{\partial H_{12}}{\partial R} = \frac{\partial H_{21}}{\partial R} $',
          r'$ \frac{\partial H_{22}}{\partial R} $']
b      = [2,0.01]
xlim   = [-8,8]
name   = ipt+'ddcpes'
plotd(data,legend,b,xlim,name,ext)

# ----------------------------------------- #

# Extended Coupling Diabatic
data = np.loadtxt(dpt+'ec.dat')
legend = [r'$ H_{11} \times 50 $',r'$ H_{12} = H_{21} $',r'$ H_{22} \times 50 $']
b      = [2,0.025]
xlim   = [-10,10]
name   = ipt+'ecpes'
plotd(data,legend,b,xlim,name,ext)

# Extended Coupling Adiabatic
legend = [r'$ E_{1} $',r'$ E_{2}']
b      = [2,0.05]
xlim   = [-10,10]
name   = ipt+'aecpes'
plota(data,legend,b,xlim,name,ext)

# Double Avoided Crossing E2 - E1
legend = [r'$ E_{2} - E_{1} $']
b      = [2,0.05]
xlim   = [-10,10]
name   = ipt+'del_aecpes'
plotad(data,legend,b,xlim,name,ext)

# Extended Coupling Diabatic Derivatives
data = np.loadtxt(dpt+'dec.dat')
legend = [r'$ \frac{\partial H_{11}}{\partial R} = \frac{\partial H_{22}}{\partial R} $',
          r'$ \frac{\partial H_{12}}{\partial R} = \frac{\partial H_{21}}{\partial R} $',
          r'$ \frac{\partial H_{22}}{\partial R} = \frac{\partial H_{11}}{\partial R} $']
b      = [2,0.01]
xlim   = [-10,10]
name   = ipt+'decpes'
plotd(data,legend,b,xlim,name,ext)

# Extended Coupling Adiabatic Derivatives
legend = [r'$ \frac{\partial E_{1}}{\partial R} $',r'$ \frac{\partial E_{2}}{\partial R} $']
b      = [2,0.02]
xlim   = [-10,10]
name   = ipt+'daecpes'
plota(data,legend,b,xlim,name,ext)