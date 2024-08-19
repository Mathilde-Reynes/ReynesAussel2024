# -*- coding: utf-8 -*-

from brian2 import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

VPY_1 = np.loadtxt('../../Data/PY_v_c_disconnected.txt')
VPY_2 = np.loadtxt('../../Data/PY_v_c_disconnected_mini50.txt')
VPY_3 = np.loadtxt('../../Data/PY_v_c_disconnected_synapses50.txt')
time20=np.arange(0, 20000, 0.02)
time20_s=time20/1000

###Figure 11
fig,ax = subplots(3,1, sharex = True,figsize=(20,15))
ax[0].plot(time20_s, VPY_1,color="tab:blue",linewidth=1.5)
ax[0].set_title('PY, 100% synaptic conductance, 100% mini conductance', size=25, loc='left')
ax[0].set_ylabel('mV',size=20, labelpad=25)
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(True)
ax[0].spines['left'].set_visible(True)
ax[0].yaxis.set_major_locator(MultipleLocator(base=25))
ax[0].set_xlim([0,20])
ax[0].tick_params(axis='both', which='major', labelsize=20, width=2)
ax[1].plot(time20_s, VPY_2,color="tab:blue",linewidth=1.5)
ax[1].set_title('PY, 50% synaptic conductance, 100% mini conductance', size=25, loc='left')
ax[1].set_ylabel('mV',size=20, labelpad=25)
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(True)
ax[1].spines['left'].set_visible(True)
ax[1].yaxis.set_major_locator(MultipleLocator(base=25))
ax[1].set_xlim([0,20])
ax[1].tick_params(axis='both', which='major', labelsize=20, width=2)
ax[2].plot(time20_s, VPY_3,color="tab:blue",linewidth=1.5)
ax[2].set_title('PY, 100% synaptic conductance, 50% mini conductance', size=25, loc='left')
ax[2].set_ylabel('mV',size=20, labelpad=25)
ax[2].set_xlabel('Time (s)',size=20, labelpad=25)
ax[2].spines['top'].set_visible(False)
ax[2].spines['right'].set_visible(False)
ax[2].spines['bottom'].set_visible(True)
ax[2].spines['left'].set_visible(True)
ax[2].yaxis.set_major_locator(MultipleLocator(base=25))
ax[2].set_xlim([0,20])
ax[2].tick_params(axis='both', which='major', labelsize=20, width=2)
plt.savefig('Figure11.png', dpi=300, bbox_inches='tight')
