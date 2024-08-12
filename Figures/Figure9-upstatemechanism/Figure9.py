# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 12:18:37 2024

@author: Mathilde
"""

###Figure 9
n_start = 75
fig,ax = subplots(7,1, sharex = True,figsize=(19,34))
ax[0].plot(V2_PYs.t/ms, V2_PYs.v[n_start-4]/mV,linewidth=0.4, color='tab:blue', alpha=0.6)
ax[0].plot(V2_PYs.t/ms, V2_PYs.v[n_start-3]/mV,linewidth=0.4, color='tab:blue', alpha=0.6)
ax[0].plot(V2_PYs.t/ms, V2_PYs.v[n_start-2]/mV,linewidth=0.4, color='tab:blue', alpha=0.6)
ax[0].plot(V2_PYs.t/ms, V2_PYs.v[n_start-1]/mV,linewidth=0.4, color='tab:blue', alpha=0.6)
ax[0].plot(V2_PYs.t/ms, V2_PYs.v[n_start]/mV,linewidth=2, color='tab:blue')
ax[0].plot(V2_PYs.t/ms, V2_PYs.v[n_start+1]/mV,linewidth=0.4, color='tab:blue', alpha=0.6)
ax[0].plot(V2_PYs.t/ms, V2_PYs.v[n_start+2]/mV,linewidth=0.4, color='tab:blue', alpha=0.6)
ax[0].plot(V2_PYs.t/ms, V2_PYs.v[n_start+3]/mV,linewidth=0.4, color='tab:blue', alpha=0.6)
ax[0].plot(V2_PYs.t/ms, V2_PYs.v[n_start+4]/mV,linewidth=0.4, color='tab:blue', alpha=0.6)
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['bottom'].set_visible(True)
ax[0].spines['left'].set_visible(True)
ax[0].set_xlim([3250,3700])
ax[0].yaxis.set_major_locator(MultipleLocator(base=25))
ax[0].tick_params(axis='both', which='major', labelsize=25, width=2)
ax[0].set_title('Pyramidal cells', size=30, loc='left')
ax[0].set_ylabel(r'$\mathrm{mV}$',size=30,  labelpad=37)
#
ax[1].plot(I1_PYd.t/ms, I1_PYd.IEPSPs_PY_PY[n_start]/(0.001*amp*meter**-2),linewidth=1.5, color="black", label='EPSPs from PY')
ax[1].plot(I1_PYd.t/ms, I1_PYd.IEPSPs_IN_PY[n_start]/(0.001*amp*meter**-2),linewidth=1.5, color="gray", label='IPSPs from IN')
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['bottom'].set_visible(True)
ax[1].spines['left'].set_visible(True)
ax[1].tick_params(axis='both', which='major', labelsize=25, width=2)
ax[1].set_title('miniPSPs', size=30, loc='left')
ax[1].set_ylabel(r'$1 \, \mathrm{mA} \cdot \mathrm{m}^{-2}$',size=30,  labelpad=25)
ax[1].legend(fontsize=25, loc='lower right')
#
ax[2].plot(I1_PYd.t/ms, I1_PYd.I_nap[n_start]/(0.001*amp*meter**-2),color="black",linewidth=1.5, label='10*I_Na(p)')
ax[2].plot(I1_PYd.t/ms, I1_PYd.I_na[n_start]/(0.01*amp*meter**-2),color="gray",linewidth=1.5, label='I_Na')
ax[2].spines['top'].set_visible(False)
ax[2].spines['right'].set_visible(False)
ax[2].spines['bottom'].set_visible(True)
ax[2].spines['left'].set_visible(True)
ax[2].tick_params(axis='both', which='major', labelsize=25, width=2)
ax[2].set_title('Sodium currents', size=30, loc='left')
ax[2].set_ylabel(r'$10 \, \mathrm{mA} \cdot \mathrm{m}^{-2}$',size=30,  labelpad=30)
ax[2].legend(fontsize=25, loc='lower right')
#
ax[3].plot(I1_PYd.t/ms, I1_PYd.IsynAMPA_PY_PY[n_start]/(0.01*amp*meter**-2),color="black",linewidth=1.5, label='I_AMPAs')
ax[3].plot(I1_PYd.t/ms, I1_PYd.IsynNMDA_PY_PY[n_start]/(0.001*amp*meter**-2),color="gray",linewidth=1.5, label='10*I_NMDAs')
ax[3].spines['top'].set_visible(False)
ax[3].spines['right'].set_visible(False)
ax[3].spines['bottom'].set_visible(True)
ax[3].spines['left'].set_visible(True)
ax[3].tick_params(axis='both', which='major', labelsize=25, width=2)
ax[3].set_title('Synaptic currents PY-PY', size=30, loc='left')
ax[3].set_ylabel(r'$10 \, \mathrm{mA} \cdot \mathrm{m}^{-2}$',size=30,  labelpad=44)
ax[3].legend(fontsize=25, loc='lower right')
#
ax[4].plot(I1_PYd.t/ms, I1_PYd.I_kca[n_start]/(0.001*amp*meter**-2),color="black",linewidth=1.5)
ax[4].spines['top'].set_visible(False)
ax[4].spines['right'].set_visible(False)
ax[4].spines['bottom'].set_visible(True)
ax[4].spines['left'].set_visible(True)
ax[4].tick_params(axis='both', which='major', labelsize=25, width=2)
ax[4].set_title('Calcium-dependent potassium current', size=30, loc='left')
ax[4].set_ylabel(r'$1 \, \mathrm{mA} \cdot \mathrm{m}^{-2}$',size=30,  labelpad=42)
#
ax[5].plot(I1_PYd.t/ms, I1_PYd.IsynGABAA_IN_PY[n_start]/(0.01*amp*meter**-2),color="black",linewidth=1.5)
ax[5].spines['top'].set_visible(False)
ax[5].spines['right'].set_visible(False)
ax[5].spines['bottom'].set_visible(True)
ax[5].spines['left'].set_visible(True)
ax[5].tick_params(axis='both', which='major', labelsize=25, width=2)
ax[5].set_title('Synaptic GABAA current IN-PY', size=30, loc='left')
ax[5].set_ylabel(r'$10 \, \mathrm{mA} \cdot \mathrm{m}^{-2}$',size=30,  labelpad=65)
#
ax[6].plot(S1.t/ms, S1.D[(S_AMPA_PY_PY.j[:] == n_start).nonzero()[0][0]],color="black",linewidth=1.5) #see the first synapse that target the neuron of interest
ax[6].plot(S1.t/ms, S1.D[(S_AMPA_PY_PY.j[:] == n_start).nonzero()[0][1]],color="black",linewidth=0.4, alpha=0.6)
ax[6].plot(S1.t/ms, S1.D[(S_AMPA_PY_PY.j[:] == n_start).nonzero()[0][2]],color="black",linewidth=0.4, alpha=0.6)
ax[6].plot(S1.t/ms, S1.D[(S_AMPA_PY_PY.j[:] == n_start).nonzero()[0][3]],color="black",linewidth=0.4, alpha=0.6)
ax[6].spines['top'].set_visible(False)
ax[6].spines['right'].set_visible(False)
ax[6].spines['bottom'].set_visible(True)
ax[6].spines['left'].set_visible(True)
ax[6].tick_params(axis='both', which='major', labelsize=25, width=2)
ax[6].set_title('AMPA synaptic depression from input synapses', size=30, loc='left')
ax[6].set_xlabel('ms',size=30,  labelpad=30)
fig.tight_layout()
#plt.savefig('Figure5.jpg', dpi=300)