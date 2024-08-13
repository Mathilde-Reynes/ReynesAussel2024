#!/usr/bin/env python
# coding: utf-8

#Cortical_layer.ipy provides a function to create the cortical subpart of the model.

from Soma_eqs import *
from Dendritic_eqs import *
from Synapses import *
from brian2 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

###Defining a cortical layer: instantiating neurons, synapses and Poisson inputs.
#The function create_cortical_layer permits the generation of a single cortical layer of the prefrontal cortex. A cortical layer is composed of pyramidal (PY) and interneurons/inhibitory cells (IN).
#N is considered to be the number of PY cells needed. The number of IN is logically deduced from it (N/4).

def create_cortical_layer(N):
    
    ###PREFERENCES
    prefs.codegen.target = 'cython'
    defaultclock.dt = 0.02*ms
    
    ###PARAMETERS
    #Channel-specific conductances per unit of surface
    g_syn_ampa_pypy = 0.00015*msiemens #0.00015*msiemens
    g_syn_nmda_pypy = 0.00001*msiemens
    g_syn_ampa_pyin = 0.00005*msiemens #0.00005
    g_syn_nmda_pyin = 0.000008*msiemens
    g_syn_gabaa_inpy = 0.00005*msiemens
    
    #For comparison with Figure4
    # # #Channel-specific conductances per unit of surface
    # g_syn_ampa_pypy = 0.00015*msiemens*0.6
    # g_syn_nmda_pypy = 0.00001*msiemens*0.6
    # g_syn_ampa_pyin = 0.00005*msiemens*0.6
    # g_syn_nmda_pyin = 0.000008*msiemens*0.6
    # g_syn_gabaa_inpy = 0.00005*msiemens*0.6
    
    # #For comparison with Figure5
    # g_syn_ampa_pypy = 0.00015*msiemens
    # g_syn_nmda_pypy = 0.00001*msiemens
    # g_syn_ampa_pyin = 0.00005*msiemens
    # g_syn_nmda_pyin = 0.000008*msiemens
    # g_syn_gabaa_inpy = 0.00005*msiemens   
    
    
    #Rate constants
    alpha_gabaa = 10*ms**-1*mM**-1
    beta_gabaa = 0.25*kHz
    #Number of neurons
    N_PY,N_IN = N,N/4
    #Fictional & arbitrary position of neurons for proper synapses.connect() condition
    neuron_spacing = 1*um
    #Amplitudes for minis
    A_PY_PY = 0.00006*msiemens
    A_PY_IN = 0.000025*msiemens
    
    # #For reproduction of Figure4
    # A_PY_PY = 0.00006*msiemens*1.5
    # A_PY_IN = (0.000025*msiemens)*1.5
    # A_PY_IN = (0.000025*msiemens)*0.1
    
    #Areas of the different neurons
    s_Soma_PYIN = 10**-6*cm**2
    s_Dend_PY = 165*s_Soma_PYIN
    s_Dend_IN = 50*s_Soma_PYIN
    
    #Radii
    PY_IN = 1
    IN_PY = 5
    PY_PY = 5
    
    print('gPYPY/sPYdend')
    print(g_syn_ampa_pypy/s_Dend_PY)
    print('gPYPY_mini/sPYdend')
    print(A_PY_PY/s_Dend_PY)
    print('gPYIN/sINdend')
    print(g_syn_ampa_pyin/s_Dend_IN)
    print('gPYIN_mini/sINdend')
    print(A_PY_IN/s_Dend_IN)
    
    ###Instantiate neurons
    
    #Pyramidal dendrite
    PY_dendrite = NeuronGroup(N_PY,Dendritic_eqs,method='rk4',threshold='v>-20*mV',refractory=3*ms)
    PY_dendrite.v = -68*mvolt
    PY_dendrite.m_na = 0.05
    PY_dendrite.h_na = 0.95
    PY_dendrite.m_nap = 0.00
    PY_dendrite.m_km = 0.05
    PY_dendrite.m_hva = 0.05
    PY_dendrite.h_hva = 0.6
    PY_dendrite.CA_i = 0.0001*mM
    PY_dendrite.g_na = 0.8*msiemens*cm**-2 
    PY_dendrite.g_nap = 3.5*msiemens*cm**-2 
    PY_dendrite.g_km = 0.01*msiemens*cm**-2 
    PY_dendrite.g_kl = 0.0025*msiemens*cm**-2 #0.0025*msiemens*cm**-2
    PY_dendrite.E_l = -68*mV
    PY_dendrite.rho = 165
    
    #Pyramidal axosomatic
    PY_soma = NeuronGroup(N_PY,Soma_eqs,method='rk4',threshold='v>40*mV',refractory=3*ms,events={'custom_poisson_PY':'rand()<mean_rate_PY*dt','custom_poisson_IN':'rand()<mean_rate_IN*dt'})
    PY_soma.h_na = 0.95
    PY_soma.m_na = 0.05
    PY_soma.m_nap = 0.00
    PY_soma.n_k = 0.05
    PY_soma.g_na = 3000*msiemens*cm**-2 #*1.2
    PY_soma.g_nap = 15*msiemens*cm**-2
    PY_soma.g_k = 200*msiemens*cm**-2
    # PY_soma.Iext = 1*mamp*cm**-2
    
    #Interneurons dendrite
    IN_dendrite = NeuronGroup(N_IN,Dendritic_eqs,method='rk4',threshold='v>-25*mV',refractory=3*ms)
    IN_dendrite.v = -68*mvolt
    IN_dendrite.m_na = 0.05
    IN_dendrite.h_na = 0.95
    IN_dendrite.m_nap = 0.00
    IN_dendrite.m_km = 0.05
    IN_dendrite.m_hva = 0.05
    IN_dendrite.h_hva = 0.6
    IN_dendrite.CA_i = 0.0001*mM
    IN_dendrite.g_na = 0.8*msiemens*cm**-2 
    IN_dendrite.g_nap = 0.0*msiemens*cm**-2
    IN_dendrite.g_km = 0.01*msiemens*cm**-2
    IN_dendrite.g_kl = 0.00*msiemens*cm**-2
    IN_dendrite.E_l = -70*mV
    IN_dendrite.rho = 50 
    
    #Interneurons axosomatic
    IN_soma = NeuronGroup(N_IN,Soma_eqs,method='rk4',threshold='v>20*mV',refractory=3*ms,events={'custom_poisson_gabaa':'rand()<mean_rate_GABAA*dt'})
    IN_soma.h_na = 0.95
    IN_soma.m_na = 0.05
    IN_soma.m_nap = 0.00
    IN_soma.n_k = 0.05
    IN_soma.g_na = 2500*msiemens*cm**-2
    IN_soma.g_nap = 0.00*msiemens*cm**-2
    IN_soma.g_k = 200*msiemens*cm**-2
    
    #Introduce variability in parameters to ensure robustness
    r = [2.0 * rand() - 1.0, 2.0 * rand() - 1.0] 
    
    for w in range(int(N_IN)):
        RA = 2.0 * rand() - 1.0
        while (r[0] * r[1] > 0) and (r[1] * RA > 0):
            RA = 2.0 * rand() - 1.0
        r[1] = r[0]
        r[0] = RA
        # Adjust the parameters for each IN cell
        IN_dendrite.E_l[w] += RA * 0.5*mV
        
    for w in range(int(N_IN)):
        RA = 2.0 * rand() - 1.0
        while (r[0] * r[1] > 0) and (r[1] * RA > 0):
            RA = 2.0 * rand() - 1.0
        r[1] = r[0]
        r[0] = RA
        # Adjust the parameters for each IN cell
        IN_soma.g_na[w] += RA * 500*msiemens*cm**-2
        
    for w in range(int(N_IN)):
        RA = 2.0 * rand() - 1.0
        while (r[0] * r[1] > 0) and (r[1] * RA > 0):
            RA = 2.0 * rand() - 1.0
        r[1] = r[0]
        r[0] = RA
        # Adjust the parameters for each IN cell
        IN_dendrite.g_na[w] += RA * 0.5*msiemens*cm**-2
    
    for w in range(int(N_IN)):
        RA = 2.0 * rand() - 1.0
        while (r[0] * r[1] > 0) and (r[1] * RA > 0):
            RA = 2.0 * rand() - 1.0
        r[1] = r[0]
        r[0] = RA
        # Adjust the parameters for each IN cell
        IN_soma.g_k[w] += RA * 50*msiemens*cm**-2
    
    
    ###Define junctions between soma and dendritic compartments for both PY and IN
    eq_gap_dendrite ='''
        Igap_post = g * (v_post - v_pre) : amp * meter**-2 (summed)
        g : siemens * meter**-2
    '''
    eq_gap_soma ='''
        vgap_post = v_pre : volt (summed)
    '''
    
    #From soma to dendrites: pyramidal cells
    gapPY_som_den = Synapses(PY_soma,PY_dendrite,model=eq_gap_dendrite)
    gapPY_som_den.connect(j = 'i')
    gapPY_som_den.g = 1 / (R * s_Dend_PY)
    #From soma to dendrites: interneurons
    gapIN_som_den = Synapses(IN_soma,IN_dendrite,model=eq_gap_dendrite)
    gapIN_som_den.connect(j = 'i')
    gapIN_som_den.g = 1 / (R * s_Dend_IN)
    
    #From dendrites to soma: pyramidal cells
    gapPY_den_som = Synapses(PY_dendrite,PY_soma,model=eq_gap_soma)
    gapPY_den_som.connect(j = 'i')
    #From dendrites to soma: interneurons
    gapIN_den_som = Synapses(IN_dendrite,IN_soma,model=eq_gap_soma)
    gapIN_den_som.connect(j = 'i')
    
    ###Instantiate synapses
    #AMPA
    S_AMPA_PY_PY = syn_ampa(PY_soma,PY_dendrite,'IsynAMPA_PY_PY',s_Dend_PY,'abs(i-j)<='+str(PY_PY)+' and i!=j',g_syn_ampa_pypy,A_PY_PY,'IEPSPs_PY_PY',0)
    S_AMPA_PY_PY.t_last_spike = -1000*ms
    S_AMPA_PY_PY.t_last_spike_Poisson_PY = -100*ms
    S_AMPA_PY_PY.t_last_spike_Poisson_IN = -100*ms
    S_AMPA_PY_PY.D = 1
    S_AMPA_PY_PY.A_mEPSP = A_PY_PY
    
    S_AMPA_PY_IN = syn_ampa(PY_soma,IN_dendrite,'IsynAMPA_PY_IN',s_Dend_IN,'abs(floor(i*'+str(N_IN)+'/'+str(N_PY)+') -j)<='+str(PY_IN)+'',g_syn_ampa_pyin,A_PY_IN,'IEPSPs_PY_IN',1)
    S_AMPA_PY_IN.t_last_spike = -1000*ms
    S_AMPA_PY_IN.t_last_spike_Poisson_PY = -100*ms
    S_AMPA_PY_IN.t_last_spike_Poisson_IN = -100*ms
    S_AMPA_PY_IN.D = 1
    S_AMPA_PY_IN.A_mEPSP = A_PY_IN
    
    #NMDA
    S_NMDA_PY_PY = syn_nmda(PY_soma,PY_dendrite,'IsynNMDA_PY_PY',s_Dend_PY,'abs(i-j)<='+str(PY_PY)+' and i!=j',g_syn_nmda_pypy) 
    S_NMDA_PY_PY.t_last_spike = -1000*ms
    S_NMDA_PY_PY.D = 1
    S_NMDA_PY_IN = syn_nmda(PY_soma,IN_dendrite,'IsynNMDA_PY_IN',s_Dend_IN,'abs(floor(i*'+str(N_IN)+'/'+str(N_PY)+') -j)<='+str(PY_IN)+'',g_syn_nmda_pyin) 
    S_NMDA_PY_IN.t_last_spike = -1000*ms
    S_NMDA_PY_IN.D = 1
    
    #GABAA
    S_GABAA_IN_PY = syn_gabaa(IN_soma,PY_dendrite,'IsynGABAA_IN_PY',s_Dend_PY,'abs(floor(i*'+str(N_PY)+'/'+str(N_IN)+') -j)<='+str(IN_PY)+'',g_syn_gabaa_inpy,alpha_gabaa,beta_gabaa)
    S_GABAA_IN_PY.t_last_spike = -1000*ms
    S_GABAA_IN_PY.t_last_spike_Poisson_fromIN = -100*ms
    S_GABAA_IN_PY.D = 1
    S_GABAA_IN_PY.D_Poisson = 1
    
    ###Define monitors
    
    #General neuron monitoring 
    V1=StateMonitor(PY_dendrite,('v'),record=True)
    # V2=StateMonitor(PY_soma,('v','I_na'),record=True)
    V2=StateMonitor(PY_soma,('v'),record=True)
    V3=StateMonitor(IN_dendrite,('v'),record=True)
    # V4=StateMonitor(IN_soma,('v','I_na'),record=True) 
    V4=StateMonitor(IN_soma,('v'),record=True) 

    #Spike monitoring
    #R1=SpikeMonitor(PY_dendrite,record=True)
    R2=SpikeMonitor(PY_soma,record=True)
    #R3=SpikeMonitor(IN_dendrite,record=True)
    R4=SpikeMonitor(IN_soma,record=True)
    
    #Synaptic currents monitoring
    I1=StateMonitor(PY_dendrite,('I_na','I_nap','I_kca','IEPSPs_PY_PY','IEPSPs_IN_PY','IsynAMPA_PY_PY','IsynNMDA_PY_PY','IsynGABAA_IN_PY'),record=True)
    # I2=StateMonitor(IN_dendrite,('IEPSPs_PY_IN','IsynAMPA_PY_IN','IsynNMDA_PY_IN'),record=True)
    I2=StateMonitor(IN_dendrite,('I_na','I_nap','I_kca'),record=False)
    
    #Synapses monitoring
    S1=StateMonitor(S_AMPA_PY_PY,('D'),record=True)
    S2=StateMonitor(S_AMPA_PY_IN,('D','t_last_spike_Poisson_IN','W'),record=False)
    #S3=StateMonitor(S_NMDA_PY_PY,('W'),record=True)
    #S4=StateMonitor(S_NMDA_PY_IN,('W'),record=True)
    #S5=StateMonitor(S_GABAA_IN_PY,('W'),record=True)
    
    #miniEPSPs monitoring
    M0=EventMonitor(PY_soma,'custom_poisson_PY',record=False)
    M1=EventMonitor(PY_soma,'custom_poisson_IN',record=False)
    
    all_neurons=PY_dendrite,PY_soma,IN_dendrite,IN_soma
    all_gap_junctions=gapPY_som_den,gapPY_den_som,gapIN_som_den,gapIN_den_som
    all_synapses=S_AMPA_PY_PY,S_AMPA_PY_IN,S_NMDA_PY_PY,S_NMDA_PY_IN,S_GABAA_IN_PY
    all_monitors=V1,V2,V3,V4,R2,R4,I1,I2,S1,S2,M0,M1
    
    return all_neurons,all_synapses,all_gap_junctions,all_monitors  