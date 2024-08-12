#!/usr/bin/env python
# coding: utf-8

#Thalamus.ipy provides a function to create the thalamus subpart of the model.

from TC_eqs import *
from RE_eqs import *
from Synapses_propre import *
from brian2 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


def create_thalamic_subparts(N):
    
    ###Preferences
    prefs.codegen.target = 'cython'
    defaultclock.dt = 0.02*ms
    
    ###Parameters
    
    #Channel-specific conductances per unit of surface
    g_syn_gabaa_retc = 0.0002*msiemens #0.0002*msiemens 
    g_syn_gabaa_rere = 0.0002*msiemens
    g_syn_ampa_tcre = 0.0004*msiemens #0.0004*msiemens
    
    #Reverse potentials
    E_gabaa = -70*mV
    E_gabaa_RE_TC = -83*mV
    
    #Rate constants
    alpha_gabaa_thal = 10.5*ms**-1*mM**-1
    beta_gabaa_thal = 0.166*kHz
    
    #Number of neurons
    N_RE,N_TC = N,N
    
    #Areas of the different neurons
    s_RE = 1.43e-4*cm**2
    s_TC = 2.9E-4*cm**2
    
    #Radius of connection
    TC_RE = 5
    RE_TC = 5
    RE_RE = 5

    
    ###Instantiate neurons
    #Reticular neurons
    RE = NeuronGroup(N,RE_eqs,method='rk4',threshold='v>0*mV',refractory=3*ms,events={'custom_poisson_PY':'v>1000*mV','custom_poisson_IN':'v>1000*mV','custom_poisson_inter':'v>1000*mV'})
    RE.v = -61*mV
    RE.m_na = 0.01
    RE.h_na = 0.99
    RE.m_t = 0.01
    RE.h_t = 0.01
    RE.n_k = 0.01
    RE.CA_i_RE = 1e-4*uM
    RE.g_kl_RE = 0.005*msiemens*cm**-2 #0.005*msiemens*cm**-2
    
    #Thalamic relay cells
    TC = NeuronGroup(N,TC_eqs,method='rk4',threshold='v>20*mV',refractory=3*ms,events={'custom_poisson_PY':'v>1000*mV','custom_poisson_IN':'v>1000*mV','custom_poisson_inter':'v>1000*mV'})
    TC.v = -68*mV
    TC.m_na = 0.01
    TC.h_na = 0.99
    TC.m_t = 0.01
    TC.h_t = 0.01
    TC.n_k = 0.01
    TC.CA_i_TC = 1e-4*uM
    TC.P1 = 0.0
    TC.Op = 0.5
    TC.Op_L = 0.0
    TC.g_kl_TC = 0.03*msiemens*cm**-2 #0.03*msiemens*cm**-2
    
    #Introduce variability in parameters to ensure robustness
    r = [2.0 * rand() - 1.0, 2.0 * rand() - 1.0]  # Initialize random numbers
    
    for w in range(int(N_TC)):
        RA = 2.0 * rand() - 1.0
        while (r[0] * r[1] > 0) and (r[1] * RA > 0):
            RA = 2.0 * rand() - 1.0
        r[1] = r[0]
        r[0] = RA
        # Adjust the parameters for each TC cell
        TC.g_kl_TC[w] += RA * 0.001*msiemens*cm**-2
        
    for w in range(int(N_RE)):
        RA = 2.0 * rand() - 1.0
        while (r[0] * r[1] > 0) and (r[1] * RA > 0):
            RA = 2.0 * rand() - 1.0
        r[1] = r[0]
        r[0] = RA
        # Adjust the parameters for each RE cell
        RE.g_kl_RE[w] += RA * 0.001*msiemens*cm**-2
    
    ###Instantiate synapses
    S_GABAA_RE_TC = syn_gabaa_thal(RE,TC,'IsynGABAA_RE_TC',s_TC,'abs(i-j)<='+str(RE_TC)+'',g_syn_gabaa_retc,E_gabaa_RE_TC,alpha_gabaa_thal,beta_gabaa_thal)
    S_GABAA_RE_TC.t_last_spike = -1000*ms
    S_GABAB_RE_TC = syn_gabab(RE,TC,'IsynGABAB_RE_TC',s_TC,'abs(i-j)<='+str(RE_TC)+'') 
    S_GABAB_RE_TC.t_last_spike = -1000*ms
    S_GABAA_RE_RE = syn_gabaa_thal(RE,RE,'IsynGABAA_RE_RE',s_RE,'abs(i-j)<='+str(RE_RE)+' and i!=j',g_syn_gabaa_rere,E_gabaa,alpha_gabaa_thal,beta_gabaa_thal) 
    S_GABAA_RE_RE.t_last_spike = -1000*ms
    S_AMPA_TC_RE = syn_ampa_thal(TC,RE,'IsynAMPA_TC_RE',s_RE,'abs(i-j)<='+str(TC_RE)+'',g_syn_ampa_tcre) 
    S_AMPA_TC_RE.t_last_spike = -1000*ms
    
    ###Define monitors
    #General neuron monitoring
    V1=StateMonitor(RE,('v'),record=True)
    V2=StateMonitor(TC,('v'),record=True)
    #Spike monitoring
    R1=SpikeMonitor(RE,record=True)
    R2=SpikeMonitor(TC,record=True)
    #Synaptic currents monitoring
    I1=StateMonitor(RE,('IsynAMPA_PY_RE'),record=True)
    I2=StateMonitor(TC,('IsynGABAB_RE_TC','IsynGABAA_RE_TC','I_h','I_t','ratio','E_ca_Thal'),record=True)
    #Synapses monitoring
    #S1=StateMonitor(S_GABAA_RE_TC,('W','g_syn'),record=True)
    #S2=StateMonitor(S_GABAB_RE_TC,('R_gab','G'),record=True)
    #S3=StateMonitor(S_GABAA_RE_RE,('W'),record=True)
    #S4=StateMonitor(S_AMPA_TC_RE,('W'),record=True)
    #P0=PopulationRateMonitor(TC)
    
    all_neurons=RE,TC
    all_synapses=S_GABAA_RE_TC,S_GABAB_RE_TC,S_GABAA_RE_RE,S_AMPA_TC_RE
    all_monitors=V1,V2,R1,R2,I1,I2
    
    return all_neurons,all_synapses,all_monitors  
