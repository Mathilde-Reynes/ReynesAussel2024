# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 12:30:05 2024

@author: Mathilde
"""

###Figure S1
fig,ax = subplots(4,1, sharex = True,figsize=(19,24))
ax[0].plot(V2_PYs.t[::50]/second, V2_PYs.v[N//2][::50]/mV,color="tab:blue",linewidth=1.5)
ax[0].set_title('PY', size=30, loc='left')
ax[0].set_ylabel('mV',size=30, labelpad=30)
#ax[0].set_xlim([15000,30000])
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(True)
ax[0].spines['left'].set_visible(True)
ax[0].yaxis.set_major_locator(MultipleLocator(base=25))
ax[0].tick_params(axis='both', which='major', labelsize=25, width=2)
#
ax[1].plot(V4_INs.t[::50]/second, V4_INs.v[int(N_IN//2)][::50]/mV,color="tab:green",linewidth=1.5)
ax[1].set_title('IN', size=30, loc='left')
ax[1].set_ylabel('mV',size=30, labelpad=30)
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(True)
ax[1].spines['left'].set_visible(True)
ax[1].yaxis.set_major_locator(MultipleLocator(base=25))
ax[1].tick_params(axis='both', which='major', labelsize=25, width=2)
#
ax[2].plot(V1_RE_TCo.t[::50]/second, V1_RE_TCo.v[20][::50]/mV,color="tab:red",linewidth=1.5)
ax[2].set_title('RE', size=30, loc='left')
ax[2].set_ylabel('mV',size=30, labelpad=30)
ax[2].spines['top'].set_visible(False)
ax[2].spines['right'].set_visible(False)
ax[2].spines['bottom'].set_visible(True)
ax[2].spines['left'].set_visible(True)
ax[2].yaxis.set_major_locator(MultipleLocator(base=25))
ax[2].tick_params(axis='both', which='major', labelsize=25, width=2)
#
ax[3].plot(V2_TC_TCo.t[::50]/second, V2_TC_TCo.v[20][::50]/mV,color="tab:orange",linewidth=1.5)
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

###Figure S2
fig,ax = subplots(4,1, sharex = True,figsize=(19,24))
ax[0].plot(V2_PYs.t/second, V2_PYs.v[N//2]/mV,color="tab:blue",linewidth=1.5)
ax[0].set_title('PY axosomatic compartment', size=30, loc='left')
ax[0].set_ylabel('mV',size=30, labelpad=30)
#ax[0].set_xlim([15000,30000])
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(True)
ax[0].spines['left'].set_visible(True)
ax[0].set_ylim(-80, 50) 
ax[0].yaxis.set_major_locator(MultipleLocator(base=25))
ax[0].tick_params(axis='both', which='major', labelsize=25, width=2)
#
ax[1].plot(V1_PYd.t/second, V1_PYd.v[N//2]/mV,color="tab:blue",linewidth=1.5)
ax[1].set_title('PY dendritic compartment', size=30, loc='left')
ax[1].set_ylabel('mV',size=30, labelpad=30)
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(True)
ax[1].spines['left'].set_visible(True)
ax[1].set_ylim(-80, 50) 
ax[1].yaxis.set_major_locator(MultipleLocator(base=25))
ax[1].tick_params(axis='both', which='major', labelsize=25, width=2)
#
ax[2].plot(V4_INs.t/second, V4_INs.v[int(N_IN//2)]/mV,color="tab:green",linewidth=1.5)
ax[2].set_title('IN axosomatic compartment', size=30, loc='left')
ax[2].set_ylabel('mV',size=30, labelpad=30)
ax[2].spines['top'].set_visible(False)
ax[2].spines['right'].set_visible(False)
ax[2].spines['bottom'].set_visible(True)
ax[2].spines['left'].set_visible(True)
ax[2].set_ylim(-100, 50) 
ax[2].yaxis.set_major_locator(MultipleLocator(base=25))
ax[2].tick_params(axis='both', which='major', labelsize=25, width=2)
#
ax[3].plot(V3_INd.t/second, V3_INd.v[int(N_IN//2)]/mV,color="tab:green",linewidth=1.5)
ax[3].set_title('IN dendritic compartment', size=30, loc='left')
ax[3].set_ylabel('mV', size=30, labelpad=30)
ax[3].set_xlabel('second', size=30)
ax[3].tick_params(axis='both',labelbottom=True)
ax[3].spines['top'].set_visible(False)
ax[3].spines['right'].set_visible(False)
ax[3].spines['bottom'].set_visible(True)
ax[3].spines['left'].set_visible(True)
ax[3].set_ylim(-100, 50) 
ax[3].yaxis.set_major_locator(MultipleLocator(base=25))
ax[3].tick_params(axis='both', which='major', labelsize=25, width=2)
fig.tight_layout()