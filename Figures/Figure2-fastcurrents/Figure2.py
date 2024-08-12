# -*- coding: utf-8 -*-

from brian2 import *
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

simulation_time=6000*ms

v_bazhenov=genfromtxt('../../Data/v_SOMA_IN')
Ik_bazhenov=genfromtxt('../../Data/I_k_IN')[:,15]
Ina_bazhenov=genfromtxt('../../Data/I_na_IN_soma')[:,15]
hna_bazhenov=genfromtxt('../../Data/h_na_IN')[:,15]
mna_bazhenov=genfromtxt('../../Data/m_na_IN')[:,15]
v_Ena=genfromtxt('../../Data/v_m_E_Na_in')[:,15]
time_bazhenov=arange(0,int(simulation_time/ms),0.02)
array_v_bazhenov=TimedArray(v_bazhenov[:,15],0.02*ms)

###Parameters

#Area
s_soma = 10**-6*cm**2

#Time constants
tau_m_nap = 0.1991*ms

#Ion-specific reversal potential
E_na = 50*mV
E_k = -90*mV

#Other constants 
Qt = pow(2.3,((36-23)/10))
Tad = pow(2.3,((36-23)/10))
Phi_m = pow(2.3,((36-23)/10))


###Equations

Soma_eqs = '''

    v = array_v_bazhenov(t)*mV  : volt  (constant over dt)
    
    truev = (vgap + (R * s_soma) * g2) / (1 + (R * s_soma) * g1) : volt    
                            
    g1 = g_na *(m_na ** 3) * h_na * Phi_m + g_k * n_k * Tad + g_nap *  m_nap : siemens * meter**-2
    g2 = g_na *(m_na ** 3) * h_na * Phi_m * E_na + g_k * n_k * Tad * E_k + g_nap * m_nap * E_na + 6.74172*uA*cm**-2 + Iext : amp*meter**-2
   
    I_na = Phi_m * g_na * (m_na ** 3) * h_na * (v - E_na) : amp * meter**-2
        dm_na/dt = -(m_na - m_nainf) / tau_m_na : 1
        dh_na/dt = -(h_na - h_nainf) / tau_h_na : 1
        alpham_na = (0.182/ms * (v + 25*mV)/mV / (1 - exp(-(v + 25*mV)/9/mV))) * int(abs((v-10*mV)/(-35*mV)) > 1e-6) + (0.182/ms * 9) * int(abs((v-10*mV)/(-35*mV)) < 1e-6)  : Hz
        betam_na = (-0.124/ms * (v + 25*mV)/mV / (1 - exp((v + 25*mV)/9/mV))) * int(abs((-v+10*mV)/(35*mV)) > 1e-6) + (0.124/ms * 9) * int(abs((-v+10*mV)/(35*mV)) < 1e-6)  : Hz
        alphah_na = (0.024/ms * (v + 40*mV)/mV / (1 - exp(-(v + 40*mV)/5/mV))) * int(abs((v-10*mV)/(-50*mV)) > 1e-6) + (0.024/ms * 5) * int(abs((v-10*mV)/(-50*mV)) < 1e-6)  : Hz
        betah_na = (-0.0091/ms * (v + 65*mV)/mV / (1 - exp((v + 65*mV)/5/mV))) * int(abs((-v+10*mV)/(75*mV)) > 1e-6) + (0.0091/ms * 5) * int(abs((-v+10*mV)/(75*mV)) < 1e-6)  : Hz
        h_nainf = 1 / (1 + exp((v + 55*mV)/6.2/mV)) : 1
        tau_h_na = (1 / (alphah_na + betah_na)) / Qt : second   
        m_nainf = alpham_na / (alpham_na + betam_na) : 1
        tau_m_na = (1 / (alpham_na + betam_na)) / Qt : second
        
    I_nap = g_nap * m_nap * (v - E_na) : amp * meter**-2 
        dm_nap/dt = -(m_nap -  m_napinf)/tau_m_nap : 1
        m_napinf = 0.02 / (1 + exp(-(v + 42*mV)/5/mV)) : 1 
    
    I_k = Tad * g_k * n_k * (v - E_k) : amp * meter**-2
        dn_k/dt = -(n_k - n_kinf) / tau_n_k : 1
        alphan_k = (0.02/mV) * (v - 25*mV) / (1 - exp(-(v - 25*mV)/(9*mV))) : 1
        betan_k = (-0.002/mV) * (v - 25*mV) / (1 - exp((v - 25*mV)/(9*mV))) : 1
        n_kinf = alphan_k / (alphan_k + betan_k) : 1
        tau_n_k = (1*msecond / (alphan_k  + betan_k)) / Qt : second
        
    Iext : amp * meter**-2
        
    vgap : volt
     
    g_na : siemens * meter**-2
    g_k : siemens * meter**-2
    g_nap : siemens * meter**-2

    t_last_spike_PY : second
    
    t_last_spike_IN : second
    
    mean_rate_PY = (log((t-t_last_spike_PY + 50*ms)/(50*ms))/400*kHz)*int((t-t_last_spike_PY)>70*ms) : Hz
    
    mean_rate_IN = (log((t-t_last_spike_PY + 50*ms)/(50*ms))/400*kHz)*int((t-t_last_spike_PY)>70*ms) : Hz
    
    mean_rate_GABAA = (log((t-t_last_spike_IN + 50*ms)/(50*ms))/400*kHz)*int((t-t_last_spike_IN)>70*ms) : Hz
     

    '''
    
# Simulation parameters
close('all')    
start_scope()
prefs.codegen.target = "numpy"
defaultclock.dt = 0.02*ms

# Initialization
Bazhenov_IN_soma = NeuronGroup(1,Soma_eqs,method='rk4',threshold='v>0*mV',refractory=3*ms)
Bazhenov_IN_soma.h_na = 0.95
Bazhenov_IN_soma.m_na = 0.05
Bazhenov_IN_soma.m_nap = 0.00
Bazhenov_IN_soma.n_k = 0.05
Bazhenov_IN_soma.g_na = 3000*msiemens*cm**-2
Bazhenov_IN_soma.g_nap = 15*msiemens*cm**-2
Bazhenov_IN_soma.g_k = 200*msiemens*cm**-2

# Monitoring
V1 = StateMonitor(Bazhenov_IN_soma,('v','I_na','I_k','m_na','h_na'),record=True)

# Run simulation
run(simulation_time, report='text',report_period=30*second)


# Plotting
# Figure 2A

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

ax[0].plot(time_bazhenov, v_bazhenov[:int((simulation_time / 0.02) / ms), 15], color="black", linewidth=lineb)
configure_axis(ax[0], 'Membrane potential used for computations', 'mV', 25, 0, 1000)
#
ax[1].plot(V1.t / ms, V1.I_na[0] / (10 * amp * meter**-2), label=a, color="#2A52BE", linewidth=linea)
ax[1].plot(time_bazhenov, Ina_bazhenov[:int((simulation_time / 0.02) / ms)] / 1000, label=b, color="#4B9CD3", alpha=alphaa, linewidth=lineb)
configure_axis(ax[1], 'Fast sodium current, $I_{Na}$', 'mA/cm²', 10, 0, 1000, show_legend=True)
#
ax[2].plot(V1.t / ms, V1.I_k[0] / (10 * amp * meter**-2), label=a, color="#2A52BE", linewidth=linea)
ax[2].plot(time_bazhenov, Ik_bazhenov[:int((simulation_time / 0.02) / ms)] / 1000, label=b, color="#4B9CD3", linewidth=lineb, alpha=alphaa)
configure_axis(ax[2], 'Fast potassium current, $I_{K}$', 'mA/cm²', 5, 0, 1000, show_legend=True)
ax[2].set_xlabel('Time (ms)', size=30, labelpad=25)
#
plt.savefig('Figure2A.png', dpi=300, bbox_inches='tight')
plt.show()

# Figure 2B
v_Ena_us = V1.v[0] - E_na
fig, ax = plt.subplots(5, 1, sharex=True, figsize=(15, 25))
ax[0].plot(time_bazhenov, v_bazhenov[:int((simulation_time / 0.02) / ms), 15], color="black", linewidth=lineb)
configure_axis(ax[0], 'Membrane potential used for computations', 'mV', 25, 414.25, 415.5)
#
ax[1].plot(V1.t / ms, v_Ena_us/mV, label=a, color="#2A52BE", linewidth=linea)
ax[1].plot(time_bazhenov, v_Ena[:int((simulation_time / 0.02) / ms)], label=b, color="#4B9CD3", alpha=alphaa, linewidth=lineb)
configure_axis(ax[1], 'Fast sodium electrochemical driving force (v-E)', 'mV', 25, 414.25, 415.5, show_legend=True)
#
ax[2].plot(V1.t / ms, V1.I_na[0] / (10 * amp * meter**-2), label=a, color="#2A52BE", linewidth=linea)
ax[2].plot(time_bazhenov, Ina_bazhenov[:int((simulation_time / 0.02) / ms)] / 1000, label=b, color="#4B9CD3", alpha=alphaa, linewidth=lineb)
configure_axis(ax[2], 'Fast sodium current, $I_{Na}$', 'mA/cm²', 10, 414.25, 415.5, show_legend=True)
#
ax[3].plot(V1.t / ms, V1.m_na[0], label=a, color="#2A52BE", linewidth=linea)
ax[3].plot(time_bazhenov, mna_bazhenov[:int((simulation_time / 0.02) / ms)], label=b, color="#4B9CD3", linewidth=lineb, alpha=alphaa)
configure_axis(ax[3], 'Activation gating variable, $m$', '', 0.5, 414.25, 415.5, show_legend=True)
#
ax[4].plot(V1.t / ms, V1.h_na[0], label=a, color="#2A52BE", linewidth=linea)
ax[4].plot(time_bazhenov, hna_bazhenov[:int((simulation_time / 0.02) / ms)], label=b, color="#4B9CD3", linewidth=lineb, alpha=alphaa)
configure_axis(ax[4], 'Inactivation gating variable, $h$', '', 0.5, 414.25, 415.5, show_legend=True)
ax[4].set_xlabel('Time (ms)', size=30, labelpad=25)
#
#plt.savefig('Figure2B.png', dpi=300, bbox_inches='tight')
plt.show()