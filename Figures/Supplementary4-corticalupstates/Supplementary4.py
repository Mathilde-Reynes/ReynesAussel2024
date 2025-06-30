from brian2 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator

VPY20=loadtxt(r'C:\Users\Mathilde\Desktop\Reynes_Aussel\Model_bazhenov\results_bazhenov_disconnected_cortex\results_bazhenov_disconnected_cortex\v_SOMA_20')
VPY60=loadtxt(r'C:\Users\Mathilde\Desktop\Reynes_Aussel\Model_bazhenov\results_bazhenov_disconnected_cortex\results_bazhenov_disconnected_cortex\v_SOMA_60')
VPY100=loadtxt(r'C:\Users\Mathilde\Desktop\Reynes_Aussel\Model_bazhenov\results_bazhenov_disconnected_cortex\results_bazhenov_disconnected_cortex\v_SOMA_100')
time=np.arange(0, 20000, 0.02)
time_s=time/1000

###
fig,ax = subplots(3,1, sharex = True,figsize=(15,15))
ax[0].plot(time_s, VPY20[:,10],color="tab:blue",linewidth=1.5)
ax[0].set_title('PY, N=20', size=35, loc='left')
ax[0].set_ylabel('mV',size=30, labelpad=25)
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(True)
ax[0].spines['left'].set_visible(True)
ax[0].yaxis.set_major_locator(MultipleLocator(base=25))
ax[0].set_xlim([7.5,15.5])
ax[0].tick_params(axis='both', which='major', labelsize=30, width=2)
ax[1].plot(time_s, VPY60[:,30],color="tab:blue",linewidth=1.5)
ax[1].set_title('PY, N=60', size=35, loc='left')
ax[1].set_ylabel('mV',size=30, labelpad=25)
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(True)
ax[1].spines['left'].set_visible(True)
ax[1].yaxis.set_major_locator(MultipleLocator(base=25))
ax[1].set_xlim([7.5,15.5])
ax[1].tick_params(axis='both', which='major', labelsize=30, width=2)
ax[2].plot(time_s, VPY100[:,50],color="tab:blue",linewidth=1.5)
ax[2].set_title('PY, N=100', size=35, loc='left')
ax[2].set_ylabel('mV',size=30, labelpad=25)
ax[2].set_xlabel('Time (s)',size=30, labelpad=25)
ax[2].spines['top'].set_visible(False)
ax[2].spines['right'].set_visible(False)
ax[2].spines['bottom'].set_visible(True)
ax[2].spines['left'].set_visible(True)
ax[2].yaxis.set_major_locator(MultipleLocator(base=25))
ax[2].set_xlim([7.5,15.5])
ax[2].tick_params(axis='both', which='major', labelsize=30, width=2)
plt.savefig('Supp_cell.png', dpi=300, bbox_inches='tight')


dt = 0.02 / 1000
time_steps = VPY20.shape[0]
time = np.linspace(0, time_steps * dt, time_steps)

fig = plt.figure(figsize=(15, 15))
gs = gridspec.GridSpec(3, 2, width_ratios=[20, 0.5], height_ratios=[1, 1, 1], hspace=0.2, wspace=0.05)

ax1 = plt.subplot(gs[0, 0])
im1 = ax1.imshow(VPY20.T, aspect='auto', cmap='Greys', vmin=-75, vmax=-60,
                 extent=[0, 20, VPY20.shape[1], 0], interpolation='bicubic')
ax1.set_title('PY, N=20', size=35, loc='left')
#ax1.set_ylabel('Neuron index', size=25, labelpad=30)
ax1.set_xlim(0, 20)
ax1.set_ylim(VPY20.shape[1], 0)
ax1.yaxis.set_major_locator(MultipleLocator(base=5))
ax1.set_xticklabels([])
ax1.tick_params(axis='y', labelsize=30)

cbar1_ax = plt.subplot(gs[0, 1])
cbar1 = fig.colorbar(im1, cax=cbar1_ax, orientation='vertical')
cbar1.locator = MultipleLocator(5)
cbar1.update_ticks()
cbar1.ax.tick_params(labelsize=30, width=2)

ax2 = plt.subplot(gs[1, 0])
im2 = ax2.imshow(VPY60.T, aspect='auto', cmap='Greys', vmin=-75, vmax=-60,
                 extent=[0, 20, VPY60.shape[1], 0], interpolation='bicubic')
ax2.set_title('PY, N=60', size=35, loc='left')
ax2.set_ylabel('Neuron index', size=35, labelpad=30)
ax2.set_xlim(0, 20)
ax2.set_ylim(VPY60.shape[1], 0)
ax2.yaxis.set_major_locator(MultipleLocator(base=20))
ax2.set_xticklabels([]) 
ax2.tick_params(axis='y', labelsize=30)

cbar2_ax = plt.subplot(gs[1, 1])
cbar2 = fig.colorbar(im2, cax=cbar2_ax, orientation='vertical')
cbar2.locator = MultipleLocator(5)
cbar2.update_ticks()
cbar2.set_label('Membrane potential (mV)', size=30, labelpad=30)
cbar2.ax.tick_params(labelsize=30, width=2)

ax3 = plt.subplot(gs[2, 0])
im3 = ax3.imshow(VPY100.T, aspect='auto', cmap='Greys', vmin=-75, vmax=-60,
                 extent=[0, 20, VPY100.shape[1], 0], interpolation='bicubic')
ax3.set_title('PY, N=100', size=35, loc='left')
ax3.set_xlabel('Time (s)', size=30, labelpad=10)
#ax3.set_ylabel('Neuron index', size=25, labelpad=30)
ax3.set_xlim(0, 20)
ax3.set_ylim(VPY100.shape[1], 0)
ax3.yaxis.set_major_locator(MultipleLocator(base=25))
ax3.xaxis.set_major_locator(MultipleLocator(base=5))
ax3.tick_params(axis='both', which='major', labelsize=30, width=2)

cbar3_ax = plt.subplot(gs[2, 1])
cbar3 = fig.colorbar(im3, cax=cbar3_ax, orientation='vertical')
cbar3.locator = MultipleLocator(5)
cbar3.update_ticks()
cbar3.ax.tick_params(labelsize=30, width=2)
plt.savefig('Supp_raster.png', dpi=300, bbox_inches='tight')
plt.show()
