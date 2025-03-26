# -*- coding: utf-8 -*-
"""
Created on Mon Mar 24 18:00:09 2025

@author: Mathilde
"""
from brian2 import *
import numpy as np
import matplotlib.pyplot as plt
from brian2.units.constants import *
import matplotlib.gridspec as gridspec
from scipy import signal

v_PY_bazhenov=genfromtxt(r'C:\Users\Mathilde\Desktop\Reynes_Aussel\Data_bazhenov\v_SOMA')
v_TC_bazhenov=genfromtxt(r'C:\Users\Mathilde\Desktop\Reynes_Aussel\Data_bazhenov\v_TC')

# Number of time steps (1500000) and time step in seconds (0.02 ms = 0.00002 s)
dt = 0.02 / 1000
time_steps = v_PY_bazhenov.shape[0]
time = np.linspace(0, time_steps * dt, time_steps)  # Time in seconds (0 to 30 seconds)


fig = plt.figure(figsize=(12, 4))
gs = gridspec.GridSpec(1, 2, width_ratios=[20, 0.5], height_ratios=[1], wspace=0.05)

ax1 = plt.subplot(gs[0])
im1 = ax1.imshow(v_PY_bazhenov.T, aspect='auto', cmap='YlGnBu', vmin=-75, vmax=-60, 
                 extent=[0, 25, 100, 0], interpolation='bicubic')

ax1.set_title('PY Neurons', size=30, loc='left', pad=30)
ax1.set_xlabel('Time (s)', size=25, labelpad=10)
ax1.set_ylabel('Neuron index', size=25, labelpad=30)
ax1.set_xlim(0, 25)  
ax1.yaxis.set_major_locator(MultipleLocator(base=25)) 
ax1.tick_params(axis='both', which='major', labelsize=25, width=2)

cbar_ax = plt.subplot(gs[1])
cbar = fig.colorbar(im1, cax=cbar_ax, orientation='vertical')
cbar.set_label('Membrane potential (mV)', size=25, labelpad=30)
cbar.ax.tick_params(labelsize=25, width=2)
plt.savefig(r"C:\Users\Mathilde\Desktop\Reynes_Aussel\Figures\Figure6-cmapfullmodel\raster_bazhenov.png", dpi=300, bbox_inches='tight')
plt.show()
