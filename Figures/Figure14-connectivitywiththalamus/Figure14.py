# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 17:38:07 2024

@author: Mathilde
"""

from brian2 import *
import numpy as np
import matplotlib.pyplot as plt
RPY_1 = np.loadtxt('../../Data/PY_raster.txt')
RPY_2 = np.loadtxt('../../Data/PY_raster_disconnected.txt')
time20=np.arange(0, 20000, 0.02)
time20_s=time20/1000
time30=np.arange(0, 30000, 0.02)
time30_s=time20/1000

###Figure 6
fig,ax = subplots(2,1, sharex = True,figsize=(19,15))
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(True)
ax[0].spines['left'].set_visible(True)
ax[0].plot(time30_s, RPY_1,'.',markersize=2,alpha=0.5,color="tab:blue")
ax[0].set_ylabel('Neuron index',fontsize=30)
ax[0].set_xlabel('Time (s)',size=30)
ax[0].tick_params(axis='both', which='major', labelsize=25, width=2)
ax[0].set_title('PY', size=30, loc='left')
ax[1].plot(time20_s, RPY_2, '.',markersize=2,alpha=0.5,color="tab:green")
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(True)
ax[1].spines['left'].set_visible(True)
ax[1].set_ylabel('Neuron index',fontsize=30)
ax[1].set_xlabel('Time (s)',size=30)
ax[1].tick_params(axis='both', which='major', labelsize=25, width=2)
ax[1].set_title('IN', size=30, loc='left')
ax[1].set_xlim([0,10])
fig.suptitle('Cortical cells raster plot',fontsize=30)
fig.tight_layout()
#plt.savefig('Figure6Disconnected', dpi=300)