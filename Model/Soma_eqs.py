#!/usr/bin/env python
# coding: utf-8

from brian2 import *

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

#Somatic and dendritic compartments
R = 10 * Mohm
g_ds = 1 / (R * s_soma)

###Equations

Soma_eqs = '''

    v = -68*mV*init_timestep + truev*(1-init_timestep) : volt
    init_timestep = init_timedarray(t) : 1
    
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
    
    #exp and log used for figure 12
    
    #log
    # mean_rate_PY = (log((t-t_last_spike_PY + 50*ms)/(50*ms))/400*kHz)*int((t-t_last_spike_PY)>70*ms) : Hz
    # mean_rate_IN = (log((t-t_last_spike_PY + 50*ms)/(50*ms))/400*kHz)*int((t-t_last_spike_PY)>70*ms) : Hz
    # mean_rate_GABAA = (log((t-t_last_spike_IN + 50*ms)/(50*ms))/400*kHz)*int((t-t_last_spike_IN)>70*ms) : Hz
    
    #exp
    # mean_rate_PY = ((2/(1 +exp(-(t-t_last_spike_PY)/(400*ms))) - 1)/100/ms)*int((t-t_last_spike_PY)>70*ms) : Hz
    # mean_rate_IN = ((2/(1 +exp(-(t-t_last_spike_PY)/(400*ms))) - 1)/100/ms)*int((t-t_last_spike_PY)>70*ms) : Hz
    # mean_rate_GABAA = ((2/(1 +exp(-(t-t_last_spike_IN)/(400*ms))) - 1)/100/ms)*int((t-t_last_spike_PY)>70*ms) : Hz
