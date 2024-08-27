# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 11:53:05 2024

@author: Mathilde
"""

from brian2 import *
import numpy as np
import matplotlib.pyplot as plt
RPY_1 = np.loadtxt('../../Data/PY_raster.txt')
RPY_2 = np.loadtxt('../../Data/PY_raster_disconnected.txt')
timeRPY_1=np.loadtxt('../../Data/PY_raster_time.txt')
timeRPY_2=np.loadtxt('../../Data/PY_raster_time_disconnected.txt')


fig1, ax1 = plt.subplots(figsize=(19, 9))

ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_visible(True)
ax1.spines['left'].set_visible(True)
ax1.plot(timeRPY_1, RPY_1, '.', markersize=2, alpha=0.5, color="tab:blue")
ax1.set_ylabel('Neuron index', fontsize=30)
ax1.set_xlabel('Time (s)', size=30)
ax1.set_xlim([0,10])
ax1.tick_params(axis='both', which='major', labelsize=25, width=2)
ax1.set_title('Raster plot, PY', size=30, loc='left')
fig1.tight_layout()
#fig1.savefig('Figure13A.png', dpi=300) 
#
fig2, ax2 = plt.subplots(figsize=(19, 9))
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_visible(True)
ax2.spines['left'].set_visible(True)
ax2.plot(timeRPY_2, RPY_2, '.', markersize=2, alpha=0.5, color="tab:blue")
ax2.set_ylabel('Neuron index', fontsize=30)
ax2.set_xlabel('Time (s)', size=30)
ax2.tick_params(axis='both', which='major', labelsize=25, width=2)
ax2.set_title('Raster plot, PY', size=30, loc='left')
ax2.set_xlim([0, 10])
fig2.tight_layout()
#fig2.savefig('Figure13B.png', dpi=300)

plt.show()