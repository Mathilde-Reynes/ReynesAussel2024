# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 14:46:53 2024

@author: Mathilde
"""

from brian2 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

VPY=loadtxt('../../Data/PY_v.txt')
VTC=loadtxt('../../Data/TC_v.txt')
runtime=10
N=100

###Figure 5
#Top panel
fig = plt.figure(figsize=(20, 12))
gs = gridspec.GridSpec(2, 2, width_ratios=[20, 0.5], height_ratios=[2, 1], wspace=0.05, hspace=0.3)
ax1 = plt.subplot(gs[0, 0])
im1 = ax1.imshow(VPY, aspect='auto', cmap='YlGnBu', vmax=-60, vmin=-75, 
                 extent=[0, runtime, N-1, 0], interpolation='bicubic')
ax1.set_title('PY', size=30, loc='left')
ax1.set_ylabel('Neuron index', size=30, labelpad=30)
ax1.set_xlim(12, 30)
ax1.yaxis.set_major_locator(MultipleLocator(base=25))
ax1.tick_params(axis='both', which='major', labelsize=25, width=2)
ax2 = plt.subplot(gs[1, 0])
im2 = ax2.imshow(VTC, aspect='auto', cmap='YlGnBu', vmax=-60, vmin=-75, 
                 extent=[0, runtime, 50, 0], interpolation='bicubic')
ax2.set_title('TC', size=30, loc='left')
ax2.set_xlabel('Time (s)', size=30, labelpad=10)
ax2.set_ylabel('Neuron index', size=30, labelpad=30)
ax2.set_xlim(12, 30)
ax2.yaxis.set_major_locator(MultipleLocator(base=25))
ax2.tick_params(axis='both', which='major', labelsize=25, width=2)
cbar_ax = plt.subplot(gs[:, 1])
cbar = fig.colorbar(im1, cax=cbar_ax, orientation='vertical')
cbar.set_label('Membrane potential (mV)', size=30, labelpad=30)
cbar.ax.tick_params(labelsize=25, width=2)
#plt.savefig('Figure5top.png', dpi=300)
plt.show()
#Bottom panel
fig = plt.figure(figsize=(20, 12))
gs = gridspec.GridSpec(2, 2, width_ratios=[20, 0.5], height_ratios=[2, 1], wspace=0.05, hspace=0.3)
ax1 = plt.subplot(gs[0, 0])
im1 = ax1.imshow(VPY, aspect='auto', cmap='YlGnBu', vmax=-60, vmin=-75, 
                 extent=[0, runtime, N-1, 0], interpolation='bicubic')
ax1.set_title('PY', size=30, loc='left')
ax1.set_ylabel('Neuron index', size=30, labelpad=30)
ax1.set_xlim(15.2,18.7)
ax1.yaxis.set_major_locator(MultipleLocator(base=25))
ax1.tick_params(axis='both', which='major', labelsize=25, width=2)
ax2 = plt.subplot(gs[1, 0])
im2 = ax2.imshow(VTC, aspect='auto', cmap='YlGnBu', vmax=-60, vmin=-75, 
                 extent=[0, runtime, 50, 0], interpolation='bicubic')
ax2.set_title('TC', size=30, loc='left')
ax2.set_xlabel('Time (s)',size=30, labelpad=10)
ax2.set_ylabel('Neuron index', size=30, labelpad=30)
ax2.set_xlim(15.2,18.7)
ax2.yaxis.set_major_locator(MultipleLocator(base=25))
ax2.tick_params(axis='both', which='major', labelsize=25, width=2)
cbar_ax = plt.subplot(gs[:, 1])
cbar = fig.colorbar(im1, cax=cbar_ax, orientation='vertical')
cbar.set_label('Membrane potential (mV)', size=30, labelpad=30)
cbar.ax.tick_params(labelsize=25, width=2)
#plt.savefig('Figure5bottom.png', dpi=300)
plt.show()

