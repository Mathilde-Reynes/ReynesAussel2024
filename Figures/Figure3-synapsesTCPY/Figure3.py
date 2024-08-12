#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from brian2 import *
from Synapses_propre import *
from Layer_Thalamus import *
from Soma_eqs import *


close('all')    
start_scope()
prefs.codegen.target = "numpy"
defaultclock.dt = 0.02*ms

simulation_time=1000*ms

v_PY_soma_bazhenov=genfromtxt('../../Data/v_SOMA')
v_PY_dend_bazhenov=genfromtxt('../../Data/time_cx_dend')
array_v_PY_soma_bazhenov=TimedArray(v_PY_soma_bazhenov[:,50],0.02*ms)
array_v_PY_dend_bazhenov=TimedArray(v_PY_dend_bazhenov[:,50],0.02*ms)

v_TC_bazhenov=genfromtxt('../../Data/v_TC')
array_v_TC_bazhenov=TimedArray(v_TC_bazhenov[:,25],0.02*ms)


time_bazhenov=arange(0,int(simulation_time/ms),0.02)

a_CX_TC_bazhenov=genfromtxt('../../Data/a_cx_tc')[:,10] #synapse from neuron 50 to neuron 25 should be k=10 since radius PY_TC = 5 and N_PY = 2*N_TC
a_TC_CX_bazhenov=genfromtxt('../../Data/a_tc_cx')[:,5] #synapse from neuron 25 to neuron 50 should be k=5 since radius TC_PY = 10 and N_TC = 0.5*N_PY

eq_PY_soma='''
    v = array_v_PY_soma_bazhenov(t)*mV  : volt  (constant over dt)
    t_last_spike_PY : second
    
    t_last_spike_IN : second
    
    mean_rate_PY = (log((t-t_last_spike_PY + 50*ms)/(50*ms))/400*kHz)*int((t-t_last_spike_PY)>70*ms) : Hz
    
    mean_rate_IN = (log((t-t_last_spike_PY + 50*ms)/(50*ms))/400*kHz)*int((t-t_last_spike_PY)>70*ms) : Hz
    
    mean_rate_GABAA = (log((t-t_last_spike_IN + 50*ms)/(50*ms))/400*kHz)*int((t-t_last_spike_IN)>70*ms) : Hz

'''

eq_PY_dend='''
    v = array_v_PY_dend_bazhenov(t)*mV  : volt  (constant over dt)
    IsynAMPA_TC_PY : amp * meter**-2
    IsynGABAA_IN_PY: amp * meter**-2
    IEPSPs_IN_PY : amp * meter**-2
'''


eq_TC='''
    v = array_v_TC_bazhenov(t)*mV  : volt  (constant over dt)
    IsynAMPA_PY_TC :  amp * meter**-2
    IsynGABAA_RE_TC :  amp * meter**-2
    IsynGABAB_RE_TC :  amp * meter**-2
    
'''


g_syn_ampa_tcpy = 0.0001*msiemens #0.0001
g_syn_ampa_tcin = 0.0001*msiemens #0.0001
g_syn_ampa_pytc = 0.000025*msiemens #0.000025
g_syn_ampa_pyre = 0.00005*msiemens #0.00005

g_syn_gabaa_retc = 0.0002*msiemens #0.0002*msiemens

g_syn_ampa_pypy = 0.00015*msiemens
g_syn_nmda_pypy = 0.00001*msiemens
g_syn_ampa_pyin = 0.00005*msiemens
g_syn_nmda_pyin = 0.000008*msiemens
g_syn_gabaa_inpy = 0.00005*msiemens

A_PY_PY = 0.00006*msiemens
A_PY_IN = (0.000025*msiemens)

alpha_gabaa = 10*ms**-1*mM**-1
beta_gabaa = 0.25*kHz

E_gabaa = -70*mV
E_gabaa_RE_TC = -83*mV

alpha_gabaa_thal = 10.5*ms**-1*mM**-1
beta_gabaa_thal = 0.166*kHz


s_Soma_PYIN = 10**-6*cm**2
s_Dend_PY = 165*s_Soma_PYIN
s_Dend_IN = 50*s_Soma_PYIN
s_TC = 2.9E-4*cm**2

Bazhenov_PY_soma = NeuronGroup(1,eq_PY_soma,method='rk4',threshold='v>0*mV',refractory=3*ms,events={'custom_poisson_PY':'rand()<mean_rate_PY*dt','custom_poisson_IN':'rand()<mean_rate_IN*dt'})
Bazhenov_PY_dendrite = NeuronGroup(1,eq_PY_dend,method='rk4',threshold='v>0*mV',refractory=3*ms)
Bazhenov_TC = NeuronGroup(1,eq_TC,method='rk4',threshold='v>0*mV',refractory=3*ms)

#Check AMPA "thal" synapses
S_AMPA_PYTC = syn_ampa_thal(Bazhenov_PY_soma,Bazhenov_TC,'IsynAMPA_PY_TC',s_TC,'',g_syn_ampa_pytc) 
S_AMPA_PYTC.t_last_spike = -100*ms
        
S_AMPA_TCPY = syn_ampa_thal(Bazhenov_TC,Bazhenov_PY_dendrite,'IsynAMPA_TC_PY',s_Dend_PY,'',g_syn_ampa_tcpy) 
S_AMPA_TCPY.t_last_spike = -1000*ms

V_PY_dend = StateMonitor(Bazhenov_PY_dendrite,('v','IsynAMPA_TC_PY','IsynGABAA_IN_PY','IEPSPs_IN_PY'),record=True)
V_TC = StateMonitor(Bazhenov_TC,('v','IsynAMPA_PY_TC','IsynGABAA_RE_TC','IsynGABAB_RE_TC'),record=True)

run(simulation_time, report='text',report_period=60*second)

figure()
ax1=subplot(311)
plot(time_bazhenov,v_PY_soma_bazhenov[:int((simulation_time/0.02)/ms),50])
title('v_soma')
subplot(312, sharex=ax1)
plot(time_bazhenov,v_TC_bazhenov[:int((simulation_time/0.02)/ms),25])
title('v_TC')
subplot(313, sharex=ax1)
plot(time_bazhenov,a_CX_TC_bazhenov[:int((simulation_time/0.02)/ms)], label='Bazhenov 2002')
plot(V_TC.t/ms,V_TC.IsynAMPA_PY_TC[0]/(0.01*amp*meter**-2), label='Replicated model')
legend()
title('AMPA current')
xlabel('Time (ms)')
suptitle('AMPA CX-TC')

a, b = 'Present model', 'Original model'
alphaa, linea, lineb = 0.60, 1.5, 1.3
fig, ax = plt.subplots(3, 1, sharex=True, figsize=(15, 15))

def configure_axis(axis, title, ylabel, y_major_locator_base, beg, end, show_legend=False):
    axis.set_title(title, size=35, loc='left')
    axis.set_ylabel(ylabel, size=30, labelpad=25)
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)
    axis.spines['bottom'].set_visible(True)
    axis.spines['left'].set_visible(True)
    axis.yaxis.set_major_locator(MultipleLocator(base=y_major_locator_base))
    axis.set_xlim([beg, end])
    axis.tick_params(axis='both', which='major', labelsize=30, width=2)
    if show_legend:
        axis.legend(loc="upper right", fontsize=23)

ax[0].plot(time_bazhenov,v_PY_soma_bazhenov[:int((simulation_time/0.02)/ms),50], color="black", linewidth=lineb)
configure_axis(ax[0], 'Membrane potential used for computations', 'mV', 25, 0, 1000)
#
ax[1].plot(time_bazhenov,v_TC_bazhenov[:int((simulation_time/0.02)/ms),25], color="black", linewidth=lineb)
configure_axis(ax[1], 'Fast sodium current, $I_{Na}$', 'mA/cm²', 10, 0, 1000, show_legend=True)
#
ax[2].plot(V_TC.t/ms,V_TC.IsynAMPA_PY_TC[0]/(0.01*amp*meter**-2), label=a, color="#2A52BE", linewidth=linea)
ax[2].plot(time_bazhenov,a_CX_TC_bazhenov[:int((simulation_time/0.02)/ms)], label=b, color="#4B9CD3", linewidth=lineb, alpha=alphaa)
configure_axis(ax[2], 'Fast potassium current, $I_{K}$', 'mA/cm²', 5, 0, 1000, show_legend=True)
ax[2].set_xlabel('Time (ms)', size=30, labelpad=25)
#
#plt.savefig('Figure3A.png', dpi=300, bbox_inches='tight')
plt.show()
