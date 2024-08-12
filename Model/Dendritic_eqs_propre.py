#!/usr/bin/env python
# coding: utf-8


from brian2 import *
import numpy as np

###Parameters

#Area
s_soma = 10**-6*cm**2

#Membrane capacitance per unit of surface
Cm = 0.75*ufarad/cm**2 

#Ion-specific conductances per unit of surface
g_hva = 0.01*msiemens*cm**-2
g_kca = 0.3*msiemens*cm**-2 
g_l = 0.033*msiemens*cm**-2

#Ion-specific reversal potential
E_na = 50*mV 
E_k = -90*mV
E_ca = 140*mV
E_kl = -95*mV
E_kca = -90*mV

#Calcium constants
CAinf = 2.4E-4*mM

#Time constants
tauCA = 165*ms
tau_m_nap = 0.1991*ms

#Other constants 
Qt = pow(2.3,((36-23)/10))
Tad = pow(2.3,((36-23)/10))
Phi_m = pow(2.3,((36-23)/10))


###Equations

Dendritic_eqs = '''

    dv/dt = (- I_kl - I_na - I_nap - I_km - I_kca - I_hva - I_l - Isyn - Igap - IEPSPs + Iext) * (1/Cm)  : volt 
   
    I_kl = g_kl * (v - E_kl) : amp * meter**-2
        
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

    I_km = Tad * g_km * m_km * (v - E_k) : amp * meter**-2
        dm_km/dt = -(m_km - m_kminf) / tau_m_km : 1
        alpham_km = 0.001/ms * (v + 30*mV)/mV / (1 - exp(-(v + 30*mV)/9/mV)) : Hz
        betam_km = -0.001/ms * (v + 30*mV)/mV / (1 - exp((v + 30*mV)/9/mV)) : Hz
        m_kminf = alpham_km / (alpham_km + betam_km) : 1
        tau_m_km = (1 / (alpham_km + betam_km)) / Qt : second
        
    I_hva = Phi_m * g_hva * (m_hva ** 2) * h_hva * (v - E_ca) : amp * meter **-2
        dm_hva/dt = -(m_hva - m_hvainf) / tau_m_hva : 1
        dh_hva/dt = -(h_hva - h_hvainf) / tau_h_hva : 1
        alpham_hva = 0.055*(mV**-1)*3.8*mV/exprel((-v-27*mV)/(3.8*mV))/ms : Hz
        betam_hva = 0.94 * exp((-v - 75*mV)/17/mV)/ms : Hz
        alphah_hva = 0.000457 * exp((-13*mV - v)/50/mV)/ms : Hz
        betah_hva = 0.0065 / (1 + exp(-(v + 15*mV)/28/mV))/ms : Hz
        m_hvainf = alpham_hva / (alpham_hva + betam_hva) : 1
        tau_m_hva = (1 / (alpham_hva + betam_hva)) / Qt : second 
        h_hvainf = alphah_hva / (alphah_hva + betah_hva) : 1
        tau_h_hva = (1 / (alphah_hva + betah_hva)) / Qt : second
    drive = -A * I_hva : katal * meter**-3
    A = (5.1819E-5*mM*cm**2)/(ms*uA) : meter**-1 * second**-1 * amp**-1 * mol
    dCA_i/dt = (drive + (CAinf - CA_i) / tauCA) * int(drive > 0*katal*meter**-3) + (0*katal*meter**-3 + (CAinf - CA_i) / tauCA) * int(drive <= 0*katal*meter**-3) : mM 
        
    I_kca = Tad * g_kca * m_kca * (v - E_kca) : amp * meter**-2
        dm_kca/dt = -(m_kca - m_kcainf) / tau_m_kca : 1
        alpham_kca = (0.01/ms) * (CA_i/mM) :  Hz
        betam_kca = (0.02/ms) :  Hz
        m_kcainf = alpham_kca / (alpham_kca + betam_kca) : 1
        tau_m_kca = (1 / (alpham_kca + betam_kca)) / Qt: second 
        
    I_l = g_l * (v - E_l) : amp * meter**-2
    
    Isyn = (IsynAMPA_PY_PY + IsynNMDA_PY_PY + IsynAMPA_PY_IN + IsynNMDA_PY_IN + IsynGABAA_IN_PY + IsynAMPA_TC_PY + IsynAMPA_TC_IN): amp * meter**-2
        IsynAMPA_PY_PY : amp * meter**-2
        IsynNMDA_PY_PY :  amp * meter**-2
        IsynAMPA_PY_IN :  amp * meter**-2
        IsynNMDA_PY_IN : amp * meter**-2
        IsynGABAA_IN_PY :  amp * meter**-2 
        IsynAMPA_TC_PY : amp * meter**-2
        IsynAMPA_TC_IN : amp * meter**-2
        
    IEPSPs = IEPSPs_intracortical + IEPSPs_thalcort : amp * meter**-2
        IEPSPs_intracortical = IEPSPs_PY_PY + IEPSPs_PY_IN + IEPSPs_IN_PY : amp * meter**-2
            IEPSPs_PY_PY : amp * meter**-2
            IEPSPs_PY_IN : amp * meter**-2
            IEPSPs_IN_PY : amp * meter**-2
        IEPSPs_thalcort = IEPSPs_TC_PY + IEPSPs_TC_IN : amp * meter**-2
            IEPSPs_TC_PY : amp * meter**-2
            IEPSPs_TC_IN : amp * meter**-2

        
    Igap : amp * meter**-2 
    
    Iext : amp * meter**-2
    
    rho : 1
    
    g_na : siemens * meter**-2
    
    g_nap : siemens * meter**-2
    
    g_km : siemens * meter**-2
    
    g_kl : siemens * meter**-2
    
    E_l : volt 

    '''

