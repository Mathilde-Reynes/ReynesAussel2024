# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 11:56:29 2024

@author: Mathilde
"""

from brian2 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

VPY=loadtxt('../../Data/PY_v_c.txt')
VIN=loadtxt('../../Data/IN_v_c.txt')
VTC=loadtxt('../../Data/TC_v_c.txt')
VRE=loadtxt('../../Data/RE_v_c.txt')
time=np.arange(0, 30000, 0.02)
time_s=time/1000

###Figure 7
#top
fig,ax = subplots(4,1, sharex = True,figsize=(19,18))
ax[0].plot(time_s, VPY,color="tab:blue",linewidth=1.5)
ax[0].set_title('PY', size=30, loc='left')
ax[0].set_ylabel('mV',size=30, labelpad=30)
ax[0].set_xlim([15,30])
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(True)
ax[0].spines['left'].set_visible(True)
ax[0].yaxis.set_major_locator(MultipleLocator(base=25))
ax[0].tick_params(axis='both', which='major', labelsize=25, width=2)
#
ax[1].plot(time_s, VIN,color="tab:green",linewidth=1.5)
ax[1].set_title('IN', size=30, loc='left')
ax[1].set_ylabel('mV',size=30, labelpad=30)
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(True)
ax[1].spines['left'].set_visible(True)
ax[1].yaxis.set_major_locator(MultipleLocator(base=25))
ax[1].tick_params(axis='both', which='major', labelsize=25, width=2)
#
ax[2].plot(time_s, VRE,color="tab:red",linewidth=1.5)
ax[2].set_title('RE', size=30, loc='left')
ax[2].set_ylabel('mV',size=30, labelpad=30)
ax[2].spines['top'].set_visible(False)
ax[2].spines['right'].set_visible(False)
ax[2].spines['bottom'].set_visible(True)
ax[2].spines['left'].set_visible(True)
ax[2].yaxis.set_major_locator(MultipleLocator(base=25))
ax[2].tick_params(axis='both', which='major', labelsize=25, width=2)
#
ax[3].plot(time_s, VTC,color="tab:orange",linewidth=1.5)
ax[3].set_title('TC', size=30, loc='left')
ax[3].set_ylabel('mV', size=30, labelpad=30)
ax[3].set_xlabel('second', size=30)
ax[3].tick_params(axis='both',labelbottom=True)
ax[3].spines['top'].set_visible(False)
ax[3].spines['right'].set_visible(False)
ax[3].spines['bottom'].set_visible(True)
ax[3].spines['left'].set_visible(True)
ax[3].yaxis.set_major_locator(MultipleLocator(base=25))
ax[3].tick_params(axis='both', which='major', labelsize=25, width=2)
#plt.savefig('Figure7top.png', dpi=300)

#bottom
fig,ax = subplots(4,1, sharex = True,figsize=(19,18))
#
ax[0].plot(time, VPY,color="tab:blue",linewidth=1.5)
ax[0].axhline(y=-68, color='gray', linestyle='--', linewidth=1.5)
ax[0].text(6500, -66, '-68 mV', color='gray', fontsize=30,verticalalignment='bottom', horizontalalignment='right')
ax[0].set_title('PY', size=30, loc='left')
ax[0].set_ylabel('mV',size=30, labelpad=30)
ax[0].set_xlim([6100,8600])
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(True)
ax[0].spines['left'].set_visible(True)
ax[0].yaxis.set_major_locator(MultipleLocator(base=25))
ax[0].tick_params(axis='both', which='major', labelsize=25, width=2)
#
ax[1].plot(time, VIN,color="tab:green",linewidth=1.5)
ax[1].text(6500, -64, '-66 mV', color='gray', fontsize=30,verticalalignment='bottom', horizontalalignment='right')
ax[1].axhline(y=-66, color='gray', linestyle='--', linewidth=1.5)
ax[1].set_title('IN', size=30, loc='left')
ax[1].set_ylabel('mV',size=30, labelpad=30)
ax[1].set_xlim([6100,8600])
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(True)
ax[1].spines['left'].set_visible(True)
ax[1].yaxis.set_major_locator(MultipleLocator(base=25))
ax[1].tick_params(axis='both', which='major', labelsize=25, width=2)
#
ax[2].plot(time, VRE,color="tab:red",linewidth=1.5)
ax[2].text(6500, -71, '-73 mV', color='gray', fontsize=30,verticalalignment='bottom', horizontalalignment='right')
ax[2].axhline(y=-73, color='gray', linestyle='--', linewidth=1.5)
ax[2].set_title('RE', size=30, loc='left')
ax[2].set_ylabel('mV',size=30, labelpad=30)
ax[2].set_xlim([6100,8600])
ax[2].spines['top'].set_visible(False)
ax[2].spines['right'].set_visible(False)
ax[2].spines['bottom'].set_visible(True)
ax[2].spines['left'].set_visible(True)
ax[2].yaxis.set_major_locator(MultipleLocator(base=25))
ax[2].tick_params(axis='both', which='major', labelsize=25, width=2)
#
ax[3].plot(time, VTC,color="tab:orange",linewidth=1.5)
ax[3].text(6500, -70, '-72 mV', color='gray', fontsize=30,verticalalignment='bottom', horizontalalignment='right')
ax[3].axhline(y=-72, color='gray', linestyle='--', linewidth=1.5)
ax[3].set_title('TC', size=30, loc='left')
ax[3].set_ylabel('mV',size=30, labelpad=30)
ax[3].set_xlabel('ms', size=30)
ax[3].set_xlim([6100,8600])
ax[3].spines['top'].set_visible(False)
ax[3].spines['right'].set_visible(False)
ax[3].spines['bottom'].set_visible(True)
ax[3].spines['left'].set_visible(True)
ax[3].yaxis.set_major_locator(MultipleLocator(base=25))
ax[3].tick_params(axis='both', which='major', labelsize=25, width=2)
#plt.savefig('Figure7bottom.png', dpi=300)