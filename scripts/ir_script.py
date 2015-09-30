"""
IR Plotting Script.
Author: Daniel Celis Garza
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import ticker

plt.close('all')
params = {'backend': 'ps',
          'axes.labelsize': 24,
          'axes.titlesize': 24,
          'text.fontsize': 24,
          'legend.fontsize': 24,
          'xtick.labelsize': 24,
          'ytick.labelsize': 24,
          'text.usetex': True,
          'font.family': 'serif',
          'mathtext.fontset': 'custom'}
plt.rcParams.update(params)

dpt = '../data/no3/spectra/'
ipt = '../images/'

data1 = np.loadtxt(dpt+'no3_f_orca.txt')
x1, y1 = data1[:,0], data1[:,1]/np.max(data1[:,1]) # Normalising absorbance.
fig, ax = plt.subplots(figsize=(12,7.42))
ax.plot(x1, y1, lw=2, ls='-',c='b',label=None)
ax.set_xlabel(r'Wavenumber (cm$^{-1}$)')
ax.set_ylabel(r'Relative Absorbance')
ax.legend(loc=0)
formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True) 
formatter.set_powerlimits((-1,1)) 
ax.yaxis.set_major_formatter(formatter)
ax.legend(loc=0,frameon=False,fancybox=True,framealpha=1)
ax.xaxis.set_major_locator(ticker.MultipleLocator(base=200))
ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(n=5))
ax.set_ylim([0,1])
ax.set_xlim([np.min(data1[:,0]),np.max(data1[:,0])]) # x-range adjusted to its min and max values.
ax.invert_yaxis()
ax.invert_xaxis()
plt.subplots_adjust(bottom=0.12,left=0.07,right=0.97) # Adjusting so everything fits nicely.
plt.savefig(ipt+'irtno3.eps', format='eps', dpi=1200)

data1 = np.loadtxt(dpt+'form_f_exp.txt')
# Averaging absorbance data and normalising the average.
x1, y1 = data1[:,0], 1-np.mean(data1[:,1:],axis=1)/np.max(np.mean(data1[:,1:],axis=1))
data2 = np.loadtxt(dpt+'form_f_teo.txt')
x2, y2 = data2[:,0], data2[:,1]/np.max(data2[:,1])
fig, ax = plt.subplots(figsize=(12,7.42))
ax.plot(x1, y1, lw=2, ls = '-', c='r',label=r'Experimental')
ax.plot(x2, y2, lw=2, ls='--',c='b',label=r'Simulated')
ax.set_xlabel(r'Wavenumber (cm$^{-1}$)')
ax.set_ylabel(r'Relative Absorbance')
formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True) 
formatter.set_powerlimits((-1,1)) 
ax.yaxis.set_major_formatter(formatter)
ax.legend(loc=0,frameon=False,fancybox=True,framealpha=1)
ax.xaxis.set_major_locator(ticker.MultipleLocator(base=300))
ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(n=3))
ax.set_ylim([0,1])
ax.set_xlim([np.min(data1[:,0]),3900])
ax.invert_yaxis()
ax.invert_xaxis()
plt.subplots_adjust(bottom=0.12,left=0.07,right=0.97)
plt.savefig(ipt+'irform.eps', format='eps', dpi=1200)