#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from brian2 import *
import sys
import os

# Add the directory containing the model .py
module_path = os.path.abspath('../../Model/')
if module_path not in sys.path:
    sys.path.append(module_path)
from Synapses import *

close('all')    
start_scope()
prefs.codegen.target = "numpy"
defaultclock.dt = 0.02*ms
simulation_time=10000*ms

v_PY_soma_bazhenov=genfromtxt('../../Data_bazhenov/v_SOMA')
array_v_PY_soma_bazhenov=TimedArray(v_PY_soma_bazhenov[:,50],0.02*ms)
v_TC_bazhenov=genfromtxt('../../Data_bazhenov/v_TC')
array_v_TC_bazhenov=TimedArray(v_TC_bazhenov[:,25],0.02*ms)
time_bazhenov=arange(0,int(simulation_time/ms),0.02)
time_bazhenov_s=time_bazhenov/1000
a_CX_TC_bazhenov=genfromtxt('../../Data_bazhenov/a_cx_tc')[:,10] #synapse from neuron 50 to neuron 25 should be k=10 since radius PY_TC = 5 and N_PY = 2*N_TC

eq_PY_soma='''
    v = array_v_PY_soma_bazhenov(t)*mV  : volt  (constant over dt)
    t_last_spike_PY : second
    t_last_spike_IN : second
    mean_rate_PY = (log((t-t_last_spike_PY + 50*ms)/(50*ms))/400*kHz)*int((t-t_last_spike_PY)>70*ms) : Hz
    mean_rate_IN = (log((t-t_last_spike_PY + 50*ms)/(50*ms))/400*kHz)*int((t-t_last_spike_PY)>70*ms) : Hz
    mean_rate_GABAA = (log((t-t_last_spike_IN + 50*ms)/(50*ms))/400*kHz)*int((t-t_last_spike_IN)>70*ms) : Hz
'''

eq_TC='''
    v = array_v_TC_bazhenov(t)*mV  : volt  (constant over dt)
    IsynAMPA_PY_TC :  amp * meter**-2
    IsynGABAA_RE_TC :  amp * meter**-2
    IsynGABAB_RE_TC :  amp * meter**-2
    
'''

g_syn_ampa_pytc = 0.000025*msiemens

A_PY_PY = 0.00006*msiemens
A_PY_IN = (0.000025*msiemens)

s_Soma_PYIN = 10**-6*cm**2
s_Dend_PY = 165*s_Soma_PYIN
s_TC = 2.9E-4*cm**2

Bazhenov_PY_soma = NeuronGroup(1,eq_PY_soma,method='rk4',threshold='v>40*mV',refractory=3*ms,events={'custom_poisson_PY':'rand()<mean_rate_PY*dt','custom_poisson_IN':'rand()<mean_rate_IN*dt'})
Bazhenov_TC = NeuronGroup(1,eq_TC,method='rk4',threshold='v>20*mV',refractory=3*ms)

#Check AMPA "thal" synapses
S_AMPA_PYTC = syn_ampa_thal(Bazhenov_PY_soma,Bazhenov_TC,'IsynAMPA_PY_TC',s_TC,'',g_syn_ampa_pytc) 
S_AMPA_PYTC.t_last_spike = -100*ms

V_TC = StateMonitor(Bazhenov_TC,('v','IsynAMPA_PY_TC'),record=True)

run(simulation_time, report='text',report_period=60*second)

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

ax[0].plot(time_bazhenov_s,v_PY_soma_bazhenov[:int((simulation_time/0.02)/ms),50], color="black", linewidth=lineb)
configure_axis(ax[0], 'Membrane potential of axosomatic compartment', 'mV', 25, 5, 6)
#
ax[1].plot(time_bazhenov_s,v_TC_bazhenov[:int((simulation_time/0.02)/ms),25], color="black", linewidth=lineb)
configure_axis(ax[1], 'Membrane potential of thalamic relay cell', 'mV', 25, 5, 6, show_legend=True)
#
ax[2].plot(V_TC.t/second,V_TC.IsynAMPA_PY_TC[0]/(20*0.01*amp*meter**-2), label=a, color="#2A52BE", linewidth=linea) #Normalize by 20 (which is the number of incoming synapses in the model of Bazhenov but is not taking into account here as brian consider N_income = 1 in that test scenario)
ax[2].plot(time_bazhenov_s,a_CX_TC_bazhenov[:int((simulation_time/0.02)/ms)], label=b, color="#4B9CD3", linewidth=lineb, alpha=alphaa)
configure_axis(ax[2], 'AMPA mediated synapse from PY to TC', 'mA/cm²', 0.02, 5, 6, show_legend=True)
ax[2].set_xlabel('Time (s)', size=30, labelpad=25)
#
plt.savefig('Figure3A.png', dpi=300, bbox_inches='tight')
plt.show()


fig, ax = plt.subplots(3, 1, sharex=True, figsize=(15, 15))
ax[0].plot(time_bazhenov_s,v_PY_soma_bazhenov[:int((simulation_time/0.02)/ms),50], color="black", linewidth=lineb)
configure_axis(ax[0], 'Membrane potential of axosomatic compartment', 'mV', 25, 5.230, 5.236)
#
ax[1].plot(time_bazhenov_s,v_TC_bazhenov[:int((simulation_time/0.02)/ms),25], color="black", linewidth=lineb)
configure_axis(ax[1], 'Membrane potential of thalamic relay cell','mV', 25, 5.230, 5.236, show_legend=True)
#
ax[2].plot(V_TC.t/second,V_TC.IsynAMPA_PY_TC[0]/(20*0.01*amp*meter**-2), label=a, color="#2A52BE", linewidth=linea)
ax[2].plot(time_bazhenov_s,a_CX_TC_bazhenov[:int((simulation_time/0.02)/ms)], label=b, color="#4B9CD3", linewidth=lineb, alpha=alphaa)
configure_axis(ax[2], 'AMPA mediated synapse from PY to TC', 'mA/cm²', 0.02, 5.230, 5.236, show_legend=True)
ax[2].set_xlabel('Time (s)', size=30, labelpad=25)
#
plt.savefig('Figure3B.png', dpi=300, bbox_inches='tight')
plt.show()