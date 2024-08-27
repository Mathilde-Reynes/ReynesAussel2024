# -*- coding: utf-8 -*-

from brian2 import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

VPY_1 = np.loadtxt('../../Data/PY_v_c_nokl_strong.txt')
VRE_1 = np.loadtxt('../../Data/RE_v_c_nokl_strong.txt')
VTC_1 = np.loadtxt('../../Data/TC_v_c_nokl_strong.txt')

VPY_2 = np.loadtxt('../../Data/PY_v_c_nokl_weak.txt')
VRE_2 = np.loadtxt('../../Data/RE_v_c_nokl_weak.txt')
VTC_2 = np.loadtxt('../../Data/TC_v_c_nokl_weak.txt')

VPY_3 = np.loadtxt('../../Data/PY_v_c_nokl_weakweak.txt')
VRE_3 = np.loadtxt('../../Data/RE_v_c_nokl_weakweak.txt')
VTC_3 = np.loadtxt('../../Data/TC_v_c_nokl_weakweak.txt')

time=np.arange(0, 10000, 0.02)
time_s=time/1000
time_20=np.arange(0, 20000, 0.02)
time20_s=time_20/1000

#For firing frequency measures see the file FiringFrequency.ods
gPYPY=[0.08,0.10,0.12,0.14,0.15,0.16]
freqPY=[14.488,18.9555,22.2055,27.294,29.961,32.2]

###Figure 15 top
fig,ax = subplots(3,1, sharex = True,figsize=(12,15))
ax[0].plot(time20_s, VPY_1,color="tab:blue",linewidth=1.5)
ax[0].set_title('PY, strong PY-PY', size=35, loc='left')
ax[0].set_ylabel('mV',size=30, labelpad=25)
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(True)
ax[0].spines['left'].set_visible(True)
ax[0].yaxis.set_major_locator(MultipleLocator(base=25))
ax[0].tick_params(axis='both', which='major', labelsize=30, width=2)
ax[1].plot(time20_s, VRE_1,color="tab:red",linewidth=1.5)
ax[1].set_title('RE, strong PY-PY', size=35, loc='left')
ax[1].set_ylabel('mV',size=30, labelpad=25)
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(True)
ax[1].spines['left'].set_visible(True)
ax[1].yaxis.set_major_locator(MultipleLocator(base=25))
ax[1].tick_params(axis='both', which='major', labelsize=30, width=2)
ax[2].plot(time20_s, VTC_1,color="tab:orange",linewidth=1.5)
ax[2].set_title('TC, strong PY-PY', size=35, loc='left')
ax[2].set_ylabel('mV',size=30, labelpad=25)
ax[2].set_xlabel('Time (s)',size=30, labelpad=25)
ax[2].spines['top'].set_visible(False)
ax[2].spines['right'].set_visible(False)
ax[2].spines['bottom'].set_visible(True)
ax[2].spines['left'].set_visible(True)
ax[2].yaxis.set_major_locator(MultipleLocator(base=25))
ax[2].set_xlim([15,18])
ax[2].tick_params(axis='both', which='major', labelsize=30, width=2)
#plt.savefig('Figure15strongPYPY.png', dpi=300, bbox_inches='tight')
#
fig,ax = subplots(3,1, sharex = True,figsize=(12,15))
ax[0].plot(time20_s, VPY_2,color="tab:blue",linewidth=1.5)
ax[0].set_title('PY, weak PY-PY', size=35, loc='left')
ax[0].set_ylabel('mV',size=30, labelpad=25)
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(True)
ax[0].spines['left'].set_visible(True)
ax[0].yaxis.set_major_locator(MultipleLocator(base=25))
ax[0].tick_params(axis='both', which='major', labelsize=30, width=2)
ax[1].plot(time20_s, VRE_2,color="tab:red",linewidth=1.5)
ax[1].set_title('RE, weak PY-PY', size=35, loc='left')
ax[1].set_ylabel('mV',size=30, labelpad=25)
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(True)
ax[1].spines['left'].set_visible(True)
ax[1].yaxis.set_major_locator(MultipleLocator(base=25))
ax[1].tick_params(axis='both', which='major', labelsize=30, width=2)
ax[2].plot(time20_s, VTC_2,color="tab:orange",linewidth=1.5)
ax[2].set_title('TC, weak PY-PY', size=35, loc='left')
ax[2].set_ylabel('mV',size=30, labelpad=25)
ax[2].set_xlabel('Time (s)',size=30, labelpad=25)
ax[2].spines['top'].set_visible(False)
ax[2].spines['right'].set_visible(False)
ax[2].spines['bottom'].set_visible(True)
ax[2].spines['left'].set_visible(True)
ax[2].yaxis.set_major_locator(MultipleLocator(base=25))
ax[2].set_xlim([15,18])
ax[2].tick_params(axis='both', which='major', labelsize=30, width=2)
#plt.savefig('Figure15weakPYPY.png', dpi=300, bbox_inches='tight')
#
fig,ax = subplots(3,1, sharex = True,figsize=(12,15))
ax[0].plot(time20_s, VPY_3,color="tab:blue",linewidth=1.5)
ax[0].set_title('PY, weak PY-PY & RE-TC-RE', size=35, loc='left')
ax[0].set_ylabel('mV',size=30, labelpad=25)
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(True)
ax[0].spines['left'].set_visible(True)
ax[0].yaxis.set_major_locator(MultipleLocator(base=25))
ax[0].tick_params(axis='both', which='major', labelsize=30, width=2)
ax[1].plot(time20_s, VRE_3,color="tab:red",linewidth=1.5)
ax[1].set_title('RE, weak PY-PY & RE-TC-RE', size=35, loc='left')
ax[1].set_ylabel('mV',size=30, labelpad=25)
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(True)
ax[1].spines['left'].set_visible(True)
ax[1].yaxis.set_major_locator(MultipleLocator(base=25))
ax[1].tick_params(axis='both', which='major', labelsize=30, width=2)
ax[2].plot(time20_s, VTC_3,color="tab:orange",linewidth=1.5)
ax[2].set_title('TC, weak PY-PY & RE-TC-RE', size=35, loc='left')
ax[2].set_ylabel('mV',size=30, labelpad=25)
ax[2].set_xlabel('Time (s)',size=30, labelpad=25)
ax[2].spines['top'].set_visible(False)
ax[2].spines['right'].set_visible(False)
ax[2].spines['bottom'].set_visible(True)
ax[2].spines['left'].set_visible(True)
ax[2].yaxis.set_major_locator(MultipleLocator(base=25))
ax[2].set_xlim([15,18])
ax[2].tick_params(axis='both', which='major', labelsize=30, width=2)
#plt.savefig('Figure15weakPYPYRETC.png', dpi=300, bbox_inches='tight')

###Figure 15 bottom
plt.figure(figsize=(7, 6))
plt.plot(gPYPY, freqPY, marker='^', color="black")
plt.xlabel('Conductance PY-PY $(μS)$', size=30, labelpad=20)
plt.ylabel('Frequency (Hz)', size=30, labelpad=20)
ax = plt.gca()
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(True)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(True)
plt.tick_params(axis='both', which='major', labelsize=30, width=2)
#plt.savefig('Figure15bottom.png', dpi=300, bbox_inches='tight')
plt.show()