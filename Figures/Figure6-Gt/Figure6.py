# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 11:53:05 2024

@author: Mathilde
"""

from brian2 import *
import numpy as np
import matplotlib.pyplot as plt

###Figure 6
fig,ax = subplots(2,1, sharex = True,figsize=(19,15))
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(True)
ax[0].spines['left'].set_visible(True)
ax[0].plot(R2_PYs.t/second, R2_PYs.i,'.',markersize=2,alpha=0.5,color="tab:blue")
ax[0].set_ylabel('Neuron index',fontsize=30)
ax[0].set_xlabel('Time (s)',size=30)
ax[0].tick_params(axis='both', which='major', labelsize=25, width=2)
ax[0].set_title('PY', size=30, loc='left')
ax[1].plot(R4_INs.t/second, R4_INs.i, '.',markersize=2,alpha=0.5,color="tab:green")
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