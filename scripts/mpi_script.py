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

dpt = '../mpi_tests/'
ipt = '../images/'

data = np.loadtxt(dpt+'mpi_dat.dat')
x, y0, y1, y2, y3 = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4]
fig, ax = plt.subplots(figsize=(9.19,14.87))
ax.plot(x, y3,label=r'500 MC reps', ls = '--', c = 'b', lw = 2.5, dashes=[8, 4, 2, 4], 
        marker='o', markersize=8)
ax.plot(x, y2,label=r'100 MC reps', ls = '--', c = 'g', lw = 2.5, marker='s', markersize=8)
ax.plot(x, y1,label=r'50 MC reps', ls = '--', c = 'r', lw = 2.5, dashes=[8, 4, 2, 4, 2, 4], 
        marker='D', markersize=8)
ax.plot(x, y0,label=r'10 MC reps', ls = '-', c = 'k', lw = 2.5, marker='^', markersize=8)
ax.grid(b=True, which='major', color='0.5',linestyle=':')
ax.set_xlabel(r'Number of Cores')
ax.set_ylabel(r'Computational Time (s)')
formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True) 
formatter.set_powerlimits((-1,1)) 
ax.yaxis.set_major_formatter(formatter)
ax.legend(loc=0)
ax.set_yscale('log')
ax.set_yticks([10, 15, 20, 50, 60, 70, 80, 90, 110, 130,150,170,195,520,600,680,770,870,1000])
ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.xaxis.set_major_locator(ticker.MultipleLocator(base=1))
ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(n=4))
ax.set_xlim([1,8])
plt.legend(bbox_to_anchor=(0.0, 0.32),loc=2,borderaxespad=0.,frameon=False,fancybox=True,framealpha=1)
plt.subplots_adjust(left=0.15)
plt.savefig(ipt+'mpi.eps', format='eps', dpi=1200)