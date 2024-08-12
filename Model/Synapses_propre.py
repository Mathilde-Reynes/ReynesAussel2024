#!/usr/bin/env python
# coding: utf-8

from brian2 import *

###GABAa with short-term depression

def syn_gabaa(source,target,syntype,surfacetarget,connection_pattern,g_syn,alpha,beta):
    eq_syn_GABAA = '''_post = (D * g_syn * (1/surfacetarget) * (W/mM) * (v_post - E_syn)) / N_incoming : amp * meter**-2 (summed)
        dW/dt = alpha * (1*mM - W) * T - beta * W : mM (clock-driven)
        T = A_syn * int((t_last_spike + tmax - t) > 0*ms)* int((t - t_last_spike) > 0*ms) : mM
        D : 1
        D_Poisson : 1
        g_syn : siemens
        E_syn = -70*mV : volt
        t_last_spike : second
        alpha : second**-1*mM**-1
        beta : Hz
        q1 : second
        U = 0.07 : 1
        tau = 700*ms : second
        tmax = 0.3*ms : second
        A_syn = 0.5*mM : mM
        surfacetarget : meter**2
        IEPSPs_IN_PY_post = ((D * A_mEPSP * (W_Poisson/mM) * (v_post - E_syn) * (1/surfacetarget)))/ N_incoming : amp * meter**-2 (summed)
        dW_Poisson/dt = alpha * (1*mM - W_Poisson) * T_Poisson - beta * W_Poisson : mM (clock-driven)
        T_Poisson = (A_syn * int((t_last_spike_Poisson_fromIN + tmax - t) > 0*ms)* int((t - t_last_spike_Poisson_fromIN) > 0*ms)) : mM
        t_last_spike_Poisson_fromIN : second
        A_mEPSP = g_syn/10 : siemens
    '''
    pre_code_gabaa = '''q1 = ((t - t_last_spike) - tmax) 
        D = 1 - (1 - D*(1 - U)) * exp(-q1/tau)*int((-q1/tau)>-10)*int((-q1/tau)<10)
        D_Poisson = 1 - (1 - D_Poisson*(1 - 0)) * exp(-q1/tau)*int((-q1/tau)>-10)*int((-q1/tau)<10)
        t_last_spike = t
        t_last_spike_IN_pre = t
    '''
    poisson_gabaa = '''
        q1 = ((t - t_last_spike) - tmax) 
        t_last_spike_Poisson_fromIN = t
        D = 1 - (1 - D*(1 - U)) * exp(-q1/tau)*int((-q1/tau)>-10)*int((-q1/tau)<10)
    '''
    
    S=Synapses(source,target,model=syntype+eq_syn_GABAA,method='rk4',on_pre={'pre': pre_code_gabaa,'mini_gabaa': poisson_gabaa}, on_event={'mini_gabaa': 'custom_poisson_gabaa'})
    if connection_pattern=='':
        S.connect()
    else :
        S.connect(condition=connection_pattern, skip_if_invalid=True)
    S.surfacetarget = surfacetarget
    S.g_syn = g_syn
    S.alpha = alpha
    S.beta = beta
    # print("GABAA synapses from "+str(source)+" to "+str(target))
    # print(S.N_outgoing_pre)
    # print(S.N_incoming_post)
    return S

###GABAa without short-term depression (thalamic synapses)

def syn_gabaa_thal(source,target,syntype,surfacetarget,connection_pattern,g_syn,E_syn,alpha,beta):
    eq_syn_GABAA = '''_post = (g_syn * (1/surfacetarget) * (W/mM) * (v_post - E_syn)) / N_incoming : amp * meter**-2 (summed)
        dW/dt = alpha * (1*mM - W) * T - beta * W : mM (clock-driven)
        T = A_syn * int((t_last_spike + tmax - t) > 0*ms)* int((t - t_last_spike) > 0*ms) : mM
        g_syn : siemens
        E_syn : volt
        t_last_spike : second
        alpha : second**-1*mM**-1
        beta : Hz
        U = 0.07 : 1
        tmax = 0.3*ms : second
        A_syn = 0.5*mM : mM
        surfacetarget : meter**2
    '''
    pre_code_gabaa = '''
        t_last_spike = t
    '''
    S=Synapses(source,target,model=syntype+eq_syn_GABAA,method='rk4',on_pre=pre_code_gabaa)
    if connection_pattern=='':
        S.connect()
    else :
        S.connect(condition=connection_pattern, skip_if_invalid=True)
    S.surfacetarget = surfacetarget
    S.g_syn = g_syn
    S.E_syn = E_syn
    S.alpha = alpha
    S.beta = beta
    # print("GABAA thal synapses from "+str(source)+" to "+str(target))
    # print(S.N_outgoing_pre)
    # print(S.N_incoming_post)
    return S


###AMPA with short term depression, Poisson inputs and STDP rule

def syn_ampa(source,target,syntype,surfacetarget,connection_pattern,g_syn,A_mEPSP,poissontype,connection_pattern_mini):
    eq_syn_AMPA = '''_post = (D * g_syn * (1/surfacetarget) * (W/mM) * (v_post - E_syn)) / N_incoming : amp * meter**-2 (summed)
        dW/dt = alpha * (1*mM - W) * T - beta * W : mM (clock-driven)
        T = A_syn * int((t_last_spike + tmax - t) > 0*ms)* int((t - t_last_spike) > 0*ms) : mM
        D : 1
        q1 : second
        g_syn : siemens
        E_syn = 0*mV : volt
        t_last_spike : second
        alpha = 0.94*ms**-1*mM**-1 : second**-1*mM**-1
        beta = 0.18*kHz : Hz
        U = 0.073 : 1
        tau = 700*ms : second
        tmax = 0.3*ms : second
        A_syn = 0.5*mM : mM
        connection_pattern_mini : 1
        surfacetarget : meter**2
    '''
    eq_syn_AMPA2 = '''_post = (D * A_mEPSP * (1/surfacetarget) * (W_Poisson/mM) * (v_post - E_syn)) / N_incoming : amp * meter**-2 (summed)
        dW_Poisson/dt = alpha * (1*mM - W_Poisson) * T_Poisson - beta * W_Poisson : mM (clock-driven)
        T_Poisson = (A_syn * int((t_last_spike_Poisson_PY + tmax - t) > 0*ms)* int((t - t_last_spike_Poisson_PY) > 0*ms))*int(connection_pattern_mini == 0)  + (A_syn * int((t_last_spike_Poisson_IN + tmax - t) > 0*ms)* int((t - t_last_spike_Poisson_IN) > 0*ms))*int(connection_pattern_mini == 1) : mM
        t_last_spike_Poisson_PY : second
        t_last_spike_Poisson_IN : second
        A_mEPSP : siemens
    '''
    pre_code_ampa = '''q1 = ((t - t_last_spike) - tmax) 
        D = 1 - (1 - D*(1 - U)) * exp(-q1/tau)*int((-q1/tau)>-10)*int((-q1/tau)<10)
        t_last_spike = t
        t_last_spike_PY_pre = t
    '''
    poisson_PY = '''
        q1 = ((t - t_last_spike) - tmax) 
        t_last_spike_Poisson_PY = t
        D = (1 - (1 - D*(1 - 0)) * exp(-q1/tau)*int((-q1/tau)>-10)*int((-q1/tau)<10))*int(connection_pattern_mini == 0) + D*int(connection_pattern_mini == 1)
        W0 = W
    '''
    poisson_IN = '''
        q1 = ((t - t_last_spike) - tmax) 
        t_last_spike_Poisson_IN = t
        D = (1 - (1 - D*(1 - 0)) * exp(-q1/tau)*int((-q1/tau)>-10)*int((-q1/tau)<10))*int(connection_pattern_mini == 1) + D*int(connection_pattern_mini == 0)
        W0 = W
    '''
    S=Synapses(source,target,model=(syntype+eq_syn_AMPA)+(poissontype+eq_syn_AMPA2),method='rk4',on_pre={'pre': pre_code_ampa,'mini_PY': poisson_PY,'mini_IN': poisson_IN},on_event={'mini_PY':'custom_poisson_PY','mini_IN':'custom_poisson_IN'})
    if connection_pattern=='':
        S.connect()
    else :
        S.connect(condition=connection_pattern, skip_if_invalid=True)
    S.surfacetarget = surfacetarget
    S.g_syn = g_syn
    S.A_mEPSP = A_mEPSP
    S.connection_pattern_mini = connection_pattern_mini
    # print("AMPA synapses from "+str(source)+" to "+str(target))
    # print(S.N_outgoing_pre)
    # print(S.N_incoming_post)
    return S

###AMPA without short term depression, Poisson inputs for thalamic synapses

def syn_ampa_thal(source,target,syntype,surfacetarget,connection_pattern,g_syn):
    eq_syn_AMPA = '''_post = (g_syn * (1/surfacetarget) * (W/mM) * (v_post - E_syn)) / N_incoming : amp * meter**-2 (summed)
        dW/dt = alpha * (1*mM - W) * T - beta * W : mM (clock-driven)
        T = A_syn * int((t_last_spike + tmax - t) > 0*ms)* int((t - t_last_spike) > 0*ms) : mM
        g_syn : siemens
        E_syn = 0*mV : volt
        t_last_spike : second
        alpha = 0.94*ms**-1*mM**-1 : second**-1*mM**-1
        beta = 0.18*kHz : Hz
        U = 0.073 : 1
        tau = 700*ms : second
        tmax = 0.3*ms : second
        A_syn = 0.5*mM : mM
        surfacetarget : meter**2
    '''
    pre_code_ampa = '''
        t_last_spike = t
    '''

    S=Synapses(source,target,model=syntype+eq_syn_AMPA,method='rk4',on_pre=pre_code_ampa)
    if connection_pattern=='':
        S.connect()
    else :
        S.connect(condition=connection_pattern, skip_if_invalid=True)
    S.surfacetarget = surfacetarget
    S.g_syn = g_syn
    # print("AMPA thal synapses from "+str(source)+" to "+str(target))
    # print(S.N_outgoing_pre)
    # print(S.N_incoming_post)
    return S


###NMDA with short-term depression and v_post relationship

def syn_nmda(source,target,syntype,surfacetarget,connection_pattern,g_syn):
    eq_syn_NMDA = '''_post = (g_syn * (1/surfacetarget) * (W/mM) * f * (v_post - E_syn)) / N_incoming : amp * meter**-2 (summed)
        dW/dt = alpha * (1*mM - W) * T - beta * W : mM (clock-driven)
        T = A_syn * int((t_last_spike + tmax - t) > 0*ms)* int((t - t_last_spike) > 0*ms) : mM
        f = 1/(1 + exp(-(v_post - v_th)/(sigma))) : 1
        g_syn : siemens
        E_syn = 0*mV : volt
        t_last_spike : second
        alpha = 1*ms**-1*mM**-1 : second**-1*mM**-1
        beta = 0.0067*kHz : Hz
        A_syn = 0.5*mM : mM
        v_th = -25*mV : volt
        sigma = 12.5*mV : volt
        surfacetarget : meter**2
        D : 1
        tau = 700*ms : second
        tmax = 0.3*ms : second
        U = 0 : 1
    '''
    pre_code_nmda = '''q1 = ((t - t_last_spike) - tmax) 
        D = 1 - (1 - D*(1 - U)) * exp(-q1/tau)*int((-q1/tau)>-10)*int((-q1/tau)<10)
        t_last_spike = t
    '''
    S=Synapses(source,target,model=syntype+eq_syn_NMDA,method='rk4',on_pre=pre_code_nmda)
    if connection_pattern=='':
        S.connect()
    else :
        S.connect(condition=connection_pattern, skip_if_invalid=True)
    S.surfacetarget = surfacetarget
    S.g_syn = g_syn
    # print("NMDA synapses from "+str(source)+" to "+str(target))
    # print(S.N_outgoing_pre)
    # print(S.N_incoming_post)
    return S


###GABAb with fraction of activated receptors and G-proteins concentration
    
def syn_gabab(source,target,syntype,surfacetarget,connection_pattern):

    eq_syn_GABAB = '''_post = (g_gabab * (1/surfacetarget) * (G**4 / (G**4 + K)) * (v - E_k_gab)) / N_incoming : amp * meter**-2 (summed)
        dR_gab/dt = K1 * (1*uM - R_gab)*T - K2*R_gab : mM (clock-driven)
        dG/dt = K3*R_gab - K4*G : mM (clock-driven)
        T = A_syn * int((t_last_spike + tmax - t) >= 0*ms)* int((t - t_last_spike) >= 0*ms) : mM
        g_gabab = 0.00004*msiemens : siemens
        E_k_gab = -95*mV : volt
        K = 100*uM**4 : mM**4
        K1 = 0.5*ms**-1*mM**-1 : second**-1*mM**-1
        K2 = 0.0012*ms**-1 : Hz
        K3 = 0.1*ms**-1 : Hz
        K4 = 0.034*ms**-1 : Hz
        A_syn = 0.5*mM : mM
        t_last_spike : second
        tmax = 0.3*ms : second
        surfacetarget : meter**2
    '''
    
    pre_code_gabab = '''
        t_last_spike = t
    '''
    S=Synapses(source,target,model=syntype+eq_syn_GABAB,method='rk4',on_pre=pre_code_gabab)
    if connection_pattern=='':
        S.connect()
    else :
        S.connect(condition=connection_pattern, skip_if_invalid=True)
    S.surfacetarget = surfacetarget
    print("GABAB synapses from "+str(source)+" to "+str(target))
    print(S.N_incoming_post)
    return S

