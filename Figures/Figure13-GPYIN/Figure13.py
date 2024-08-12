# -*- coding: utf-8 -*-

from brian2 import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#Make sure you are working in the good directory
VPY_1=loadtxt('../../Data/PY_v_disconnected_gPYIN_02.txt')
VPY_2=loadtxt('../../Data/PY_v_disconnected_gPYIN_07.txt')
runtime=10
N=100
#see file "propagation_speeds.ods
gPYPY=[0.12,0.13,0.14,0.15,0.16,0.17,0.18]
speed_gPYPY=[48.255,153.439,236.683,296.382,394.173,445.183,503.808]
std_gPYPY=[39.044,111.567,125.536,154.122,135.885,135.915,183.156]

gPYIN=[0.02,0.03,0.04,0.05,0.06,0.07,0.08]
speed_gPYIN=[423.773,420.824,346.384,296.382,266.485,256.199,152.877]
std_gPYIN=[237.496,159.579,110.732,154.122,152.759,125.132,78.666]

#Top pannel
fig = plt.figure(figsize=(15, 12))
gs = gridspec.GridSpec(2, 2, width_ratios=[20, 0.5], height_ratios=[1, 1], wspace=0.05, hspace=0.3)
ax1 = plt.subplot(gs[0, 0])
im1 = ax1.imshow(VPY_1, aspect='auto', cmap='YlGnBu', vmax=-60, vmin=-75, 
                 extent=[0, runtime, N-1, 0], interpolation='bicubic')
ax1.set_title('Conductance PY-IN = 0.02$\mu$S', size=30, loc='left')
ax1.set_ylabel('Neuron index', size=30, labelpad=30)
ax1.set_xlim(0.1, 5.1)
ax1.yaxis.set_major_locator(MultipleLocator(base=25))
ax1.tick_params(axis='both', which='major', labelsize=25, width=2)
ax2 = plt.subplot(gs[1, 0])
im2 = ax2.imshow(VPY_2, aspect='auto', cmap='YlGnBu', vmax=-60, vmin=-75, 
                 extent=[0, runtime, N-1, 0], interpolation='bicubic')
ax2.set_title('Conductance PY-IN = 0.07$\mu$S', size=30, loc='left')
ax2.set_xlabel('Time (s)', size=30, labelpad=10)
ax2.set_ylabel('Neuron index', size=30, labelpad=30)
ax2.set_xlim(0.1, 5.1)
ax2.yaxis.set_major_locator(MultipleLocator(base=25))
ax2.tick_params(axis='both', which='major', labelsize=25, width=2)
cbar_ax = plt.subplot(gs[:, 1])
cbar = fig.colorbar(im1, cax=cbar_ax, orientation='vertical')
cbar.set_label('Membrane potential (mV)', size=30, labelpad=30)
cbar.ax.tick_params(labelsize=25, width=2)
plt.savefig('Figure13cmap.png', dpi=300,bbox_inches='tight')
plt.show()

#Bottom pannel
plt.figure(figsize=(15, 6))
ax1 = plt.subplot(121)
ax1.errorbar(gPYPY, speed_gPYPY, std_gPYPY, marker='^', color="black")
ax1.set_xlabel('Conductance PY-PY $(μS)$', size=30, labelpad=20)
ax1.set_ylabel('Velocity (cell/s)', size=30, labelpad=20)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_visible(True)
ax1.spines['left'].set_visible(True)
ax1.tick_params(axis='both', which='major', labelsize=30, width=2)
ax2 = plt.subplot(122)
ax2.errorbar(gPYIN, speed_gPYIN, std_gPYIN, marker='^', color="black")
ax2.set_xlabel('Conductance PY-IN $(μS)$', size=30, labelpad=20)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_visible(True)
ax2.spines['left'].set_visible(True)
ax2.tick_params(axis='both', which='major', labelsize=30, width=2)
plt.savefig('Figure13velocity.png', dpi=300, bbox_inches='tight')
plt.show()