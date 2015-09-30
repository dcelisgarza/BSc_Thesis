"""
Probabilities Plotting Script.
Author: Daniel Celis Garza
"""

import numpy as np
import matplotlib
matplotlib.use('Agg') # Use .svg format.
import matplotlib.pyplot as plt
from matplotlib import ticker

plt.close('all')
# Plotting parameters
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

scp = '../data/single_avoided_crossing/'
dcp = '../data/double_avoided_crossing/'
ecp = '../data/extended_coupling/'
ipt = '../images/'
ext = 'eps'

def plot(data,xlabel,ylabel,lpos,b,xlim,ylim,name,ext):
    x, y0, y1, y2, y3 = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4]
    fig, ax = plt.subplots(figsize=(12,7.42))
    ax.plot(x, y0,label=r'$ R_{1\leftarrow 1} $', ls = '-', c = 'b', lw = 2.5, 
            marker = 'o', markersize = 10)
    ax.plot(x, y1,label=r'$ R_{2\leftarrow 1} $', ls = '-.', c = 'g', lw = 2.5, 
            marker = 's', markersize = 10)
    ax.plot(x, y2,label=r'$ T_{1\leftarrow 1} $', ls = '--', c = 'r', lw = 2.5, 
            marker = 'd', markersize = 10)
    ax.plot(x, y3,label=r'$ T_{2\leftarrow 1} $', ls = '--', c = 'k', lw = 2.5, 
            dashes=[8, 4, 2, 4, 2, 4], marker = '*', markersize = 10)
    ax.grid(b=True, which='major', color='0.5',linestyle=':')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True) 
    formatter.set_powerlimits((-1,1)) 
    ax.yaxis.set_major_formatter(formatter)
    plt.legend(bbox_to_anchor=lpos,
           bbox_transform=plt.gcf().transFigure, frameon=False)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=b))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(n=5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(base=0.1))
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    plt.subplots_adjust(bottom=0.12)
    plt.savefig(name+'.'+ext, format=ext, dpi=1200)

# Generating the file with the average value of the two parallel runs for the Single Avoided Crossing
data1 = np.loadtxt(scp+'sc_prob1_par.dat')
data2 = np.loadtxt(scp+'sc_prob2_par.dat')

data = np.empty([30,5])

data[:,0]=np.mean([data1[:,0],data2[:,0]],axis=0)
data[:,1]=np.mean([data1[:,1],data2[:,1]],axis=0)
data[:,2]=np.mean([data1[:,2],data2[:,2]],axis=0)
data[:,3]=np.mean([data1[:,3],data2[:,3]],axis=0)
data[:,4]=np.mean([data1[:,4],data2[:,4]],axis=0)

text_file = open(scp+'sc_prob_parallel.dat', 'w')
text_file.write('# 1/5ip, two 7500 MC step runs')

for i in range(30):
    text_file.write(str(data[i,0])+' '+str(data[i,1])+' '+str(data[i,2])+' '+
    str(data[i,3])+' '+str(data[i,4])+'\n')
text_file.close()

xlabel = r'Initial Nuclear Momentum ($P_{i}$)'
ylabel = r'Probability'
lpos = (0.52, 0.7)
b = 5
xlim = [0,30]
ylim = [0,1]
name = ipt+'sc_prob_parallel'
plot(data,xlabel,ylabel,lpos,b,xlim,ylim,name,ext)

# Generating the file with the average value of the two parallel runs for the Double Avoided Crossing
data1 = np.loadtxt(dcp+'dc_prob1_par.dat')
data2 = np.loadtxt(dcp+'dc_prob2_par.dat')

data = np.empty([26,5])

data[:,0]=np.mean([data1[:,0],data2[:,0]],axis=0)
data[:,1]=np.mean([data1[:,1],data2[:,1]],axis=0)
data[:,2]=np.mean([data1[:,2],data2[:,2]],axis=0)
data[:,3]=np.mean([data1[:,3],data2[:,3]],axis=0)
data[:,4]=np.mean([data1[:,4],data2[:,4]],axis=0)

text_file = open(dcp+'dc_prob_parallel.dat', 'w')
text_file.write('# 1/5ip, two 7500 MC step runs')

for i in range(26):
    text_file.write(str(data[i,0])+' '+str(data[i,1])+' '+str(data[i,2])+' '+
    str(data[i,3])+' '+str(data[i,4])+'\n')
text_file.close()

xlabel = r'$\ln{(E_{i})}$'
ylabel = r'Probability'
lpos = (0.4, 0.7)
b = 1
xlim = [-4.2,1]
ylim = [0,1]
name = ipt+'dc_prob_parallel'
plot(data,xlabel,ylabel,lpos,b,xlim,ylim,name,ext)

# Generating the file with the average value of the two parallel runs for the Extended Coupling
data1 = np.loadtxt(ecp+'ec_prob1_par.dat')
data2 = np.loadtxt(ecp+'ec_prob2_par.dat')
data3 = np.loadtxt(ecp+'ec_prob3_par.dat')
data4 = np.loadtxt(ecp+'ec_prob4_par.dat')

data = np.empty([29,5])

data[:,0]=np.mean([data1[:,0],data2[:,0],data3[:,0],data4[:,0]],axis=0)
data[:,1]=np.mean([data1[:,1],data2[:,1],data3[:,1],data4[:,1]],axis=0)
data[:,2]=np.mean([data1[:,2],data2[:,2],data3[:,2],data4[:,2]],axis=0)
data[:,3]=np.mean([data1[:,3],data2[:,3],data3[:,3],data4[:,3]],axis=0)
data[:,4]=np.mean([data1[:,4],data2[:,4],data3[:,4],data4[:,4]],axis=0)

text_file = open(ecp+'ec_prob_parallel.dat', 'w')
text_file.write('# 1/5ip, four 3750 MC step runs')

for i in range(26):
    text_file.write(str(data[i,0])+' '+str(data[i,1])+' '+str(data[i,2])+' '+
    str(data[i,3])+' '+str(data[i,4])+'\n')
text_file.close()

xlabel = r'Initial Nuclear Momentum ($P_{i}$)'
ylabel = r'Probability'
lpos = (0.35, 0.6)
b = 2
xlim = [2,30]
ylim = [0,1]
name = ipt+'ec_prob_parallel'
plot(data,xlabel,ylabel,lpos,b,xlim,ylim,name,ext)