# -*- coding: utf-8 -*-

from brian2 import *
import matplotlib.pyplot as plt

#Make sure your working directory is properly set to access the .txt
VPY_1 = np.loadtxt('../../Data/PY_v.txt')
VPY_2 = np.loadtxt('../../Data/PY_v_Miniexp.txt')
time_10 = np.loadtxt('../../Data/time_10.txt')
time_array = np.arange(0, 30000, 0.02)

###Figure 12
fig,ax = subplots(2,1, sharex = True,figsize=(10,6))
ax[0].plot(time_array/1000, VPY_1, color="tab:blue", linewidth=1.5)
ax[0].set_title('PY, logarithmic based mini rate', size=25, loc='left')
ax[0].set_ylabel('mV', size=20, labelpad=25)
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(True)
ax[0].spines['left'].set_visible(True)
ax[0].yaxis.set_major_locator(MultipleLocator(base=25))
ax[0].set_xlim([0, 8])
ax[0].tick_params(axis='both', which='major', labelsize=20, width=2)
ax[1].plot(time_10/1000, VPY_2, color="tab:blue", linewidth=1.5)
ax[1].set_title('PY, exponential based mini rate', size=25, loc='left')
ax[1].set_ylabel('mV', size=20, labelpad=25)
ax[1].set_xlabel('Time (s)', size=20, labelpad=15)
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(True)
ax[1].spines['left'].set_visible(True)
ax[1].yaxis.set_major_locator(MultipleLocator(base=25))
ax[1].set_xlim([0, 8])
ax[1].tick_params(axis='both', which='major', labelsize=20, width=2)
plt.savefig('Figure12.jpg', dpi=300, bbox_inches='tight')

