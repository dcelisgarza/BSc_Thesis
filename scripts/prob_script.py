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
sbp = '../data/spin_boson/'
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

def plotsb(data1,data2,b,xlim,ylim,name,ext):
    x0, y0, x1, y1 = data1[:,0], data1[:,1], data2[:,0], data2[:,1]
    fig, ax = plt.subplots(figsize=(12,7.42))
    ax.plot(x0, y0,label=r'$ \epsilon = 0 $', ls = '-', c = 'b', lw = 2.5, 
            marker = 'o', markersize = 10)
    ax.plot(x1, y1,label=r'$ \epsilon = 1 $', ls = '--', c = 'r', lw = 2.5, 
            marker = '*', markersize = 10)
    ax.grid(b=True, which='major', color='0.5',linestyle=':')
    ax.set_xlabel(r'Time (A.U.)')
    ax.set_ylabel(r'$ D(t) $')
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True) 
    formatter.set_powerlimits((-1,1)) 
    ax.yaxis.set_major_formatter(formatter)
    plt.legend(loc=0,#bbox_to_anchor=lpos
           bbox_transform=plt.gcf().transFigure, frameon=False)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=b[0]))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(n=5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(base=b[1]))
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    plt.subplots_adjust(bottom=0.12)
    plt.savefig(name+'.'+ext, format=ext, dpi=1200)

def plott(data,b,xlim,ylim,name,ext):
    t, r, p, n1, n2 = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4]
    fig, ax = plt.subplots(figsize=(12,7.42))
    ax.plot(t, r,label=r'$ R $', ls = '-.', c = 'b', lw = 2.5)
    ax.plot(t, p,label=r'$ P $', ls = '--', dashes=[8, 4, 2, 4, 2, 4], c = 'g', lw = 2.5)
    ax.plot(t, n1,label=r'$ n_{1} + \gamma $', ls = '--', c = 'r', lw = 2.5)
    ax.plot(t, n2,label=r'$ n_{2} + \gamma$', ls = '-', c = 'k', lw = 2.5)
    ax.grid(b=True, which='major', color='0.5',linestyle=':')
    ax.set_xlabel(r'Time (A.U.)')
    ax.set_ylabel(r'Value (A.U.)')
#    formatter = ticker.ScalarFormatter(useMathText=True)
#    formatter.set_scientific(True) 
#    formatter.set_powerlimits((-1,1))
#    ax.yaxis.set_major_formatter(formatter)
    plt.legend(loc=0,#bbox_to_anchor=lpos
           bbox_transform=plt.gcf().transFigure, frameon=False)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=b[0]))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(n=5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(base=b[1]))
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
#    plt.subplots_adjust(bottom=0.12)
    plt.savefig(name+'.'+ext, format=ext, dpi=1200)

def plott2(data,b,lpos,xlim,ylim,name,ext):
    t, r, p, n1, n2 = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4]
    fig, ax = plt.subplots(figsize=(12,7.42))
    ax.plot(t, r,label=r'$ R $', ls = '-.', c = 'b', lw = 2.5)
    ax.plot(t, p,label=r'$ P $', ls = '--', dashes=[8, 4, 2, 4, 2, 4], c = 'g', lw = 2.5)
    ax.plot(t, n1,label=r'$ n_{1} + \gamma$', ls = '--', c = 'r', lw = 2.5)
    ax.plot(t, n2,label=r'$ n_{2} + \gamma$', ls = '-', c = 'k', lw = 2.5)
    ax.grid(b=True, which='major', color='0.5',linestyle=':')
    ax.set_xlabel(r'Time (A.U.)')
    ax.set_ylabel(r'Value (A.U.)')
#    formatter = ticker.ScalarFormatter(useMathText=True)
#    formatter.set_scientific(True) 
#    formatter.set_powerlimits((-1,1))
#    ax.yaxis.set_major_formatter(formatter)
    plt.legend(bbox_to_anchor=lpos,
           bbox_transform=plt.gcf().transFigure, frameon=False)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=b[0]))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(n=5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(base=b[1]))
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
#    plt.subplots_adjust(bottom=0.12)
    plt.savefig(name+'.'+ext, format=ext, dpi=1200)

def plote(data,b,xlim,ylim,name,ext):
    t, n1, n2 = data[:,0], data[:,3], data[:,4]
    fig, ax = plt.subplots(figsize=(12,7.42))
    ax.plot(t, n1,label=r'$ n_{1} + \gamma$', ls = '--', c = 'r', lw = 2.5)
    ax.plot(t, n2,label=r'$ n_{2} + \gamma$', ls = '-', c = 'k', lw = 2.5)
    ax.grid(b=True, which='major', color='0.5',linestyle=':')
    ax.set_xlabel(r'Time (A.U.)')
    ax.set_ylabel(r'Value (A.U.)')
#    formatter = ticker.ScalarFormatter(useMathText=True)
#    formatter.set_scientific(True) 
#    formatter.set_powerlimits((-1,1))
#    ax.yaxis.set_major_formatter(formatter)
    plt.legend(loc=0,#bbox_to_anchor=lpos
           bbox_transform=plt.gcf().transFigure, frameon=False)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=b[0]))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(n=5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(base=b[1]))
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
#    plt.subplots_adjust(bottom=0.12)
    plt.savefig(name+'.'+ext, format=ext, dpi=1200)

# Trajectories
data = np.loadtxt(scp+'sc_traj_R11.dat')
name = ipt+'sc_traj_r11'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [min(data[:,1]),max(data[:,3])+0.1]
b = [1000,0.5]
plott(data,b,xlim,ylim,name,ext)
name = ipt+'sc_traj_r11_e'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [0.3,1.7]
b = [1000,0.1]
plote(data,b,xlim,ylim,name,ext)

data = np.loadtxt(scp+'sc_traj_T12.dat')
name = ipt+'sc_traj_t12'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [min(data[:,1]),max(data[:,2])+0.1]
b = [500,1]
plott(data,b,xlim,ylim,name,ext)
name = ipt+'sc_traj_t12_e'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [0,1.15]
b = [500,0.1]
plote(data,b,xlim,ylim,name,ext)

data = np.loadtxt(scp+'sc_traj_T11.dat')
name = ipt+'sc_traj_t11'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [min(data[:,1]),max(data[:,2])+1]
b = [100,5]
plott(data,b,xlim,ylim,name,ext)
name = ipt+'sc_traj_t11_e'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [0.4,1.3]
b = [100,0.1]
plote(data,b,xlim,ylim,name,ext)

data = np.loadtxt(scp+'sc_traj_R22.dat')
name = ipt+'sc_traj_r22'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [-5,5.5]
b = [5000,1]
plott(data,b,xlim,ylim,name,ext)
name = ipt+'sc_traj_r22_e'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [0,2]
b = [5000,0.2]
plote(data,b,xlim,ylim,name,ext)

data = np.loadtxt(scp+'sc_traj_T21.dat')
name = ipt+'sc_traj_t21'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [min(data[:,1]),max(data[:,2])+0.1]
b = [200,1]
plott(data,b,xlim,ylim,name,ext)
name = ipt+'sc_traj_t21_e'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [0,1.7]
b = [200,0.2]
plote(data,b,xlim,ylim,name,ext)

data = np.loadtxt(scp+'sc_traj_T22.dat')
name = ipt+'sc_traj_t22'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [min(data[:,1]),max(data[:,2])+1]
b = [100,3]
plott(data,b,xlim,ylim,name,ext)
name = ipt+'sc_traj_t22_e'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [0.55,1.25]
b = [100,0.2]
plote(data,b,xlim,ylim,name,ext)

data = np.loadtxt(dcp+'dc_traj_T11.dat')
name = ipt+'dc_traj_t11'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [min(data[:,1]),20]
b = [200,2]
lpos = (0.3,0.67)
plott2(data,b,lpos,xlim,ylim,name,ext)
name = ipt+'dc_traj_t11_e'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [0,1.8]
b = [200,0.2]
plote(data,b,xlim,ylim,name,ext)

data = np.loadtxt(dcp+'dc_traj_T12.dat')
name = ipt+'dc_traj_t12'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [-8,28]
b = [200,3]
lpos = (0.3,0.67)
plott2(data,b,lpos,xlim,ylim,name,ext)
name = ipt+'dc_traj_t12_e'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [0,2.2]
b = [200,0.2]
plote(data,b,xlim,ylim,name,ext)

data = np.loadtxt(dcp+'dc_traj_T21.dat')
name = ipt+'dc_traj_t21'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [min(data[:,1]),20]
b = [200,2]
lpos = (0.3,0.7)
plott2(data,b,lpos,xlim,ylim,name,ext)
name = ipt+'dc_traj_t21_e'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [0,2.2]
b = [200,0.2]
plote(data,b,xlim,ylim,name,ext)

data = np.loadtxt(dcp+'dc_traj_T22.dat')
name = ipt+'dc_traj_t22'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [-8,27]
b = [200,3]
lpos = (0.3,0.67)
plott2(data,b,lpos,xlim,ylim,name,ext)
name = ipt+'dc_traj_t22_e'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [0,1.6]
b = [200,0.2]
plote(data,b,xlim,ylim,name,ext)

data = np.loadtxt(dcp+'dc_traj_low_momentum.dat')
name = ipt+'dc_traj_low_momentum'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [min(data[:,1]),22]
b = [5000,2]
plott(data,b,xlim,ylim,name,ext)
name = ipt+'dc_traj_low_momentum_e'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [0,1.6]
b = [5000,0.2]
plote(data,b,xlim,ylim,name,ext)

data = np.loadtxt(dcp+'dc_traj_low_momentum_2.dat')
name = ipt+'dc_traj_low_momentum_2'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [-8,2]
b = [2500,1]
plott(data,b,xlim,ylim,name,ext)
name = ipt+'dc_traj_low_momentum_e_2'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [0.1,1.65]
b = [2500,0.2]
plote(data,b,xlim,ylim,name,ext)

data = np.loadtxt(ecp+'ec_traj_R11.dat')
name = ipt+'ec_traj_r11'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [-16,16]
b = [250,3]
plott(data,b,xlim,ylim,name,ext)
name = ipt+'ec_traj_r11_e'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [0,1.6]
b = [250,0.2]
plote(data,b,xlim,ylim,name,ext)

data = np.loadtxt(ecp+'ec_traj_R12.dat')
name = ipt+'ec_traj_r12'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [-13.5,13.5]
b = [500,3]
plott(data,b,xlim,ylim,name,ext)
name = ipt+'ec_traj_r12_e'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [0,1.8]
b = [500,0.2]
plote(data,b,xlim,ylim,name,ext)

data = np.loadtxt(ecp+'ec_traj_R21.dat')
name = ipt+'ec_traj_r21'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [-16,16]
b = [250,3]
plott(data,b,xlim,ylim,name,ext)
name = ipt+'ec_traj_r21_e'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [0,1.7]
b = [250,0.2]
plote(data,b,xlim,ylim,name,ext)

data = np.loadtxt(ecp+'ec_traj_R22.dat')
name = ipt+'ec_traj_r22'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [-13.5,13.5]
b = [250,3]
plott(data,b,xlim,ylim,name,ext)
name = ipt+'ec_traj_r22_e'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [0,1.75]
b = [250,0.2]
plote(data,b,xlim,ylim,name,ext)

data = np.loadtxt(ecp+'ec_traj_T11.dat')
name = ipt+'ec_traj_t11'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [min(data[:,1]),27]
b = [250,3]
plott(data,b,xlim,ylim,name,ext)
name = ipt+'ec_traj_t11_e'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [0,1.7]
b = [250,0.2]
plote(data,b,xlim,ylim,name,ext)

data = np.loadtxt(ecp+'ec_traj_T12.dat')
name = ipt+'ec_traj_t12'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [-9,max(data[:,2])+1]
b = [250,3]
plott(data,b,xlim,ylim,name,ext)
name = ipt+'ec_traj_t12_e'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [0,1.9]
b = [250,0.2]
plote(data,b,xlim,ylim,name,ext)
#
data = np.loadtxt(ecp+'ec_traj_T21.dat')
name = ipt+'ec_traj_t21'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [min(data[:,1]),31]
b = [250,3]
plott(data,b,xlim,ylim,name,ext)
name = ipt+'ec_traj_t21_e'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [0,2]
b = [250,0.2]
plote(data,b,xlim,ylim,name,ext)

data = np.loadtxt(ecp+'ec_traj_T22.dat')
name = ipt+'ec_traj_t22'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [-9,max(data[:,2])+1]
b = [250,3]
plott(data,b,xlim,ylim,name,ext)
name = ipt+'ec_traj_t22_e'
xlim = [min(data[:,0]),max(data[:,0])]
ylim = [0,1.6]
b = [250,0.2]
plote(data,b,xlim,ylim,name,ext)

# Single Avoided Crossing
data = np.loadtxt(scp+'sc_prob_1o5ip_15000.dat')
xlabel = r'Initial Nuclear Momentum ($P_{i}$)'
ylabel = r'Probability'
lpos = (0.52, 0.7)
b = 5
xlim = [0,30]
ylim = [0,1]
name = ipt+'sc_prob_1o5ip'
plot(data,xlabel,ylabel,lpos,b,xlim,ylim,name,ext)

data = np.loadtxt(scp+'sc_prob_1o5ip_15000.dat')
data1 = np.loadtxt(scp+'sc_prob_1oip_15000.dat')

x, y0, y1, y2, y3 = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4]
x1, y01, y11, y21, y31 = data1[:,0], data1[:,1], data1[:,2], data1[:,3], data1[:,4]
fig, ax = plt.subplots(figsize=(12,7.42))
ax.plot(x, y0,label=r'$ R_{1\leftarrow 1} $', ls = '-', lw = 2.5, marker = 'o', markersize = 10)
ax.plot(x1, y01,label=r'$ R_{1\leftarrow 1} $', ls = '-', lw = 2.5, marker = 'o', markersize = 10)
ax.plot(x, y1,label=r'$ R_{2\leftarrow 1} $', ls = '-.', lw = 2.5, marker = 's', markersize = 10)
ax.plot(x1, y11,label=r'$ R_{2\leftarrow 1} $', ls = '-.', lw = 2.5, marker = 's', markersize = 10)
ax.plot(x, y2,label=r'$ T_{1\leftarrow 1} $', ls = '--', lw = 2.5, marker = 'd', markersize = 10)
ax.plot(x1, y21,label=r'$ T_{1\leftarrow 1} $', ls = '--', lw = 2.5, marker = 'd', markersize = 10)
ax.plot(x, y3,label=r'$ T_{2\leftarrow 1} $', ls = '--', dashes=[8, 4, 2, 4, 2, 4], lw = 2.5, 
        marker = '*', markersize = 10)
ax.plot(x1, y31,label=r'$ T_{2\leftarrow 1} $', ls = '--', dashes=[8, 4, 2, 4, 2, 4], lw = 2.5, 
        marker = '*', markersize = 10)
ax.grid(b=True, which='major', color='0.5',linestyle=':')
ax.set_xlabel(r'Initial Nuclear Momentum ($P_{i}$)')
ax.set_ylabel(r'Probability')
formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True) 
formatter.set_powerlimits((-1,1)) 
ax.yaxis.set_major_formatter(formatter)
plt.legend(bbox_to_anchor=(0.49, 0.8),
           bbox_transform=plt.gcf().transFigure, frameon=False)
ax.xaxis.set_major_locator(ticker.MultipleLocator(base=5))
ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(n=5))
ax.yaxis.set_major_locator(ticker.MultipleLocator(base=0.1))
ax.set_xlim([0,30])
ax.set_ylim([0,1])
plt.subplots_adjust(bottom=0.12)
plt.savefig(ipt+'sc_prob_1oip_vs_1o5ip.eps', format='eps', dpi=1200)

# Double Avoided Crossing
data = np.loadtxt(dcp+'dc_prob_1o5ip.dat')
xlabel = r'$\ln{(E_{i})}$'
ylabel = r'Probability'
lpos = (0.4, 0.7)
b = 1
xlim = [-4.2,1]
ylim = [0,1]
name = ipt+'dc_prob_1o5ip'
plot(data,xlabel,ylabel,lpos,b,xlim,ylim,name,ext)

# Extended Coupling
data = np.loadtxt(ecp+'ec_prob_1o5ip.dat')
xlabel = r'Initial Nuclear Momentum ($P_{i}$)'
ylabel = r'Probability'
lpos = (0.35, 0.6)
b = 2
xlim = [2,30]
ylim = [0,1]
name = ipt+'ec_prob_1o5ip'
plot(data,xlabel,ylabel,lpos,b,xlim,ylim,name,ext)

# Spin-Boson
data1 = np.loadtxt(sbp+'sb_prob_exp1A1.dat')
data2 = np.loadtxt(sbp+'sb_prob_exp1B1.dat')
b = [5,0.2]
xlim = [0,40]
ylim = [-1,1]
name = ipt+'spin_boson_e11'
plotsb(data1,data2,b,xlim,ylim,name,ext)

data1 = np.loadtxt(sbp+'sb_prob_exp1A2.dat')
data2 = np.loadtxt(sbp+'sb_prob_exp1B2.dat')
b = [5,0.2]
xlim = [0,40]
ylim = [-1,1]
name = ipt+'spin_boson_e12'
plotsb(data1,data2,b,xlim,ylim,name,ext)

data1 = np.loadtxt(sbp+'sb_prob_exp2A1.dat')
data2 = np.loadtxt(sbp+'sb_prob_exp2B1.dat')
b = [5,0.2]
xlim = [0,50]
ylim = [-1,1]
name = ipt+'spin_boson_e21'
plotsb(data1,data2,b,xlim,ylim,name,ext)

data1 = np.loadtxt(sbp+'sb_prob_exp2A2.dat')
data2 = np.loadtxt(sbp+'sb_prob_exp2B2.dat')
b = [5,0.2]
xlim = [0,50]
ylim = [-1,1]
name = ipt+'spin_boson_e22'
plotsb(data1,data2,b,xlim,ylim,name,ext)