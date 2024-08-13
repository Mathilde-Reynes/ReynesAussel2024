#!/usr/bin/env python
# coding: utf-8

#TC_eqs.ipy provides the equations to create the reticular cells of the thalamus.

from brian2 import *
import numpy as np

###Parameters

#Membrane capacitance per unit of surface
Cm_RE = 1*ufarad/cm**2

#Conductances
g_na_RE = 100*msiemens*cm**-2
g_t_RE = 2.3*msiemens*cm**-2
g_l_RE = 0.05*msiemens*cm**-2
g_k_RE = 10*msiemens*cm**-2

#Reversal potentials
E_kl = -95*mV
E_l_TC = -70*mV
E_l_RE = -77*mV
E_na = 50*mV
E_k_RE = -95*mV
E_ca0 = 1000*8.31441*(273.15 + 36)/(2*96489)*mV #13.31*mV

#Calcium parameters
tau_CA_RE = 5*ms
A_RE = 5.1819E-5*(mM*cm**2)/(ms*uA)
CA_inf = 2.4E-4*mM 
CA_0 = 2*mM

#Temperature-dependent variables 
T = 36
Qm_RE = 5**((T-24)/10)
Qh_RE = 3**((T-24)/10)
Q = 2.3
Qhyp = pow(3,((T-36)/10))

#Shift in voltage
Vtr_RE = -50*mV
VtrK_RE = -50*mV


###Equations

RE_eqs = '''

    dv/dt = (- I_kl - I_na - I_k - I_t - I_l - Isyn_RE + Iext) * (1/Cm_RE) : volt
    v2 = v - Vtr_RE : volt
    v2K = v - VtrK_RE : volt
    
    I_kl = g_kl_RE * (v - E_kl) : amp * meter**-2
    
    I_l = g_l_RE * (v - E_l_RE) : amp * meter**-2

    I_na = g_na_RE * (m_na ** 3) * h_na * (v - E_na) : amp * meter**-2
        dm_na/dt = Qhyp*(alpham_na*(1-m_na)-betam_na*m_na) : 1
        dh_na/dt = Qhyp*(alphah_na*(1-h_na)-betah_na*h_na) : 1
        
        alpham_na = 0.32/ms * (13*mV - v2)/mV / (exp((13*mV - v2)/4/mV) - 1) : Hz
        betam_na = 0.28/ms * (v2 - 40*mV)/mV / (exp((v2 - 40*mV)/5/mV) - 1) : Hz
        
        alphah_na = 0.128 * exp((17*mV - v2)/18/mV)/ms : Hz
        betah_na = 4/(exp((40*mV - v2)/5/mV) + 1)/ms  : Hz 

    I_k = g_k_RE * (n_k ** 4) * (v - E_k_RE) : amp * meter**-2 
        dn_k/dt = Qhyp*(alphan_k*(1-n_k)-betan_k*n_k) : 1
        
        alphan_k = 0.032/ms * (15*mV - v2K)/mV / (exp((15*mV - v2K)/5/mV) - 1) : Hz
        betan_k = 0.5/ms * exp((10*mV - v2K)/40/mV) : Hz

    I_t = g_t_RE * (m_t ** 2) * h_t * (v - E_ca_Thal) : amp * meter**-2
        dm_t/dt = -(m_t - m_tinf) / tau_m_t : 1
        dh_t/dt = -(h_t - h_tinf) / tau_h_t : 1
        
        tau_m_t = (3*ms + 1*ms/(exp((v + 27*mV)/10/mV) + exp(-(v + 102*mV)/15/mV))) / Qm_RE : second
        m_tinf = 1 / (1 + exp(-(v + 52*mV)/7.4/mV)) : 1
        
        tau_h_t = (85*ms + 1*ms/(exp((v + 48*mV)/4/mV) + exp(-(v + 407*mV)/50/mV))) / Qh_RE : second
        h_tinf = 1 / (1 + exp((v + 80*mV)/5/mV)) : 1
        
        drive = -A_RE * I_t : katal * meter**-3
        dCA_i_RE/dt = (drive + (CA_inf - CA_i_RE)/tau_CA_RE) * int(drive > 0*katal*meter**-3) + (0*katal*meter**-3 + (CA_inf - CA_i_RE)/tau_CA_RE) * int(drive < 0*katal*meter**-3) : mM 
        ratio = CA_0/CA_i_RE : 1
        E_ca_Thal = E_ca0 * log(ratio) : volt
        
    Isyn_RE = IsynGABAA_RE_RE + IsynAMPA_TC_RE + IsynAMPA_PY_RE : amp * meter**-2 
        IsynGABAA_RE_RE : amp * meter**-2
        IsynAMPA_TC_RE : amp * meter**-2
        IsynAMPA_PY_RE :  amp * meter**-2 
    
    Iext : amp * meter**-2
        
    g_kl_RE : siemens * meter**-2
        
    '''