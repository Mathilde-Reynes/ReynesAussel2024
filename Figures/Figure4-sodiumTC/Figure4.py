# -*- coding: utf-8 -*-

from brian2 import *
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

simulation_time=6000*ms

v_bazhenov=genfromtxt('v_TC')
Ik_bazhenov=genfromtxt('I_k_TC')[:,15]
Ina_bazhenov=genfromtxt('I_na_TC')[:,15]
It_bazhenov=genfromtxt('I_t_TC')[:,15]
Ih_bazhenov=genfromtxt('I_h_TC')[:,15]
Il_bazhenov=genfromtxt('I_l_TC')[:,15]
Ikl_bazhenov=genfromtxt('I_kl_TC')[:,15]
hna_bazhenov=genfromtxt('h_na_TC')[:,15]
mna_bazhenov=genfromtxt('m_na_TC')[:,15]
mnainf_bazhenov=genfromtxt('m_inf_na_tc')[:,15]
hnainf_bazhenov=genfromtxt('h_inf_na_tc')[:,15]
taumna_bazhenov=genfromtxt('tau_m_na_tc')[:,15]
tauhna_bazhenov=genfromtxt('tau_h_na_tc')[:,15]
alphamna_bazhenov=genfromtxt('alpha_m_TC')[:,15]
alphahna_bazhenov=genfromtxt('alpha_h_TC')[:,15]
betamna_bazhenov=genfromtxt('beta_m_TC')[:,15]
betahna_bazhenov=genfromtxt('beta_h_TC')[:,15]

time_bazhenov=arange(0,int(simulation_time/ms),0.02)
array_v_bazhenov=TimedArray(v_bazhenov[:,15],0.02*ms)

#Membrane capacitance per 2unit of surface
Cm_TC = 1*ufarad/cm**2

#Conductances
g_na_TC = 90*msiemens*cm**-2 
g_t_TC = 2.2*msiemens*cm**-2 
g_l_TC = 0.01*msiemens*cm**-2
g_k_TC = 12*msiemens*cm**-2
g_h = 0.017*msiemens*cm**-2 
g_a = 0*msiemens*cm**-2 

#Reversal potentials
E_kl = -95*mV
E_l_TC = -70*mV
E_h = -40*mV
E_na = 50*mV
E_k = -95*mV
E_ca0 = 1000*8.31441*(273.15 + 36)/(2*96489)*mV #13.31*mV approx

#Calcium parameters
tau_CA_TC = 5*ms
A_TC = 5.1819E-5*(mM*cm**2)/(ms*uA)
CA_inf = 2.4E-4*mM 
CA_0 = 2*mM #unit was found in Vijayan and Kopell 10.1073/pnas.1215385109

#Temperature-dependent variables 
T = 36
Qm_TC = 3.55**((T-24)/10)
Qh_TC = 3**((T-24)/10)
Qhyp = pow(3,((T-36)/10))
Q = 2.3
Tad = pow(3,((T-23.5)/10))

#Rates for open and close channels dynamics
k = 2
k1 = 7.9012E7*(mM**-4)*ms**-1
k2 = 0.0004*ms**-1 
k3 = 0.1*ms**-1
k4 = 0.001*ms**-1

#ModelDB parameters
cac=0.0015*mM
pc=0.007

#Parameters to match modelDB
Vtr_TC= -40*mV
VtrK_TC= -25*mV
tau_m_nap = 0.1991*ms


TC_eqs_modifs = '''
    
    v = array_v_bazhenov(t)*mV  : volt  (constant over dt)
    v2 = v - Vtr_TC : volt
    v2K = v - VtrK_TC : volt
    
    I_kl = g_kl_TC * (v - E_kl) : amp * meter**-2
    
    I_l = g_l_TC * (v - E_l_TC) : amp * meter**-2

    I_na = g_na_TC * (m_na ** 3) * h_na * (v - E_na) : amp * meter**-2 
        dm_na/dt = (-(m_na - mna_inf) / tau_m_na) : 1
        dh_na/dt = (-(h_na - hna_inf) / tau_h_na) : 1
        
        alpham_na = 0.32*(mV**-1)*4*mV/exprel((13*mV-v2)/(4*mV))/ms : Hz
        betam_na = 0.28*(mV**-1)*5*mV/exprel((v2 - 40*mV)/(5*mV))/ms : Hz
        tau_m_na = (1/(alpham_na + betam_na) / Qhyp): second
        mna_inf = alpham_na/(alpham_na + betam_na) : 1
        
        alphah_na = 0.128 * exp((17*mV - v2)/18/mV)/ms : Hz
        betah_na = 4/(exp((40*mV - v2)/5/mV) + 1)/ms  : Hz     
        tau_h_na =  1/(alphah_na + betah_na) / Qhyp : second
        hna_inf = alphah_na/(alphah_na + betah_na) : 1
        
    I_k = g_k_TC * (n_k ** 4) * (v - E_k) : amp * meter**-2
        dn_k/dt = -(n_k - nk_inf) / tau_n_k : 1
        
        alphan_k = 0.032*(mV**-1)*5*mV/exprel((15*mV - v2K)/(5*mV))/ms : Hz
        betan_k = 0.5/ms * exp((10*mV - v2K)/40/mV) : Hz
        tau_n_k =  (1/(alphan_k + betan_k) / Qhyp) : second
        nk_inf = alphan_k/(alphan_k + betan_k) : 1
        
    I_t = g_t_TC * (m_t ** 2) * h_t * (v - E_ca_Thal) : amp * meter**-2
        dm_t/dt = -(m_t - m_tinf) / tau_m_t : 1
        dh_t/dt = -(h_t - h_tinf) / tau_h_t : 1
        
        tau_m_t = ((1 / (exp(-(v + 131.6*mV)/16.7/mV) + exp((v + 16.8*mV)/18.2/mV)) + 0.612)/Qm_TC)*ms : second
        m_tinf = 1 / (1 + exp(-(v + 59*mV)/6.2/mV)) : 1
        
        tau_h_t = ((30.8 + (211.4 + exp((v + 115.2*mV)/5/mV)) / (1 + exp((v + 86*mV)/3.2/mV)))/Qh_TC)*ms : second
        h_tinf = 1 / (1 + exp((v + 83*mV)/4/mV)) : 1
                      
        drive = -A_TC * I_t/2 : katal * meter**-3
        dCA_i_TC/dt = (drive + (CA_inf - CA_i_TC)/tau_CA_TC) * int(drive > 0*katal*meter**-3) + (0*katal*meter**-3 + (CA_inf - CA_i_TC)/tau_CA_TC) * int(drive < 0*katal*meter**-3) : mM 
        ratio = CA_0/CA_i_TC : 1
        E_ca_Thal = E_ca0 * log(ratio) : volt
   
   I_h = g_h * (Op + k*Op_L) * (v - E_h) : amp * meter**-2
        dOp/dt = alpha_Op * (1 - Op - Op_L) - beta_Op*Op : 1
        dP1/dt = k1ca*(1 - P1) - k2*P1 : 1
        dOp_L/dt = k3p*Op*int((Op_L+Op)<1) - k4*Op_L : 1
       
        tau_Op = (20*ms + 1000*ms / (exp((v + 71.5*mV)/(14.2*mV)) + exp(-(v + 89*mV)/(11.6*mV))))/Qhyp : second
        m_Opinf = 1 / (1 + exp((v + 75*mV)/5.5/mV)) : 1
        
        alpha_Op = m_Opinf / tau_Op : Hz
        beta_Op = (1 - m_Opinf) / tau_Op : Hz
        
        k1ca = k2 * (CA_i_TC/cac)**4 : Hz
        k3p = k4 * (P1/pc) : Hz 
        
    I_a = g_a * (m_a **4) * h_a * (v - E_k) : amp * meter**-2
        dm_a/dt = -(1/tau_m_a)*(m_a - m_ainf) : 1
        dh_a/dt = -(1/tau_h_a)*(h_a - h_ainf) : 1
        
        tau_m_a = (1*ms / (exp((v + 35.82*mV)/19.69/mV) + exp(-(v + 79.69*mV)/12.7/mV)) + 0.37*ms) / Tad : second
        m_ainf = 1 / (1 + exp(-(v + 60*mV)/8.5/mV)) : 1
        
        tau_h_a = (19*ms/Tad) * int(v >= -63*mV) + (1*ms / ((exp((v + 46.05*mV)/5/mV) + exp(-(v + 238.4*mV)/37.45/mV))) / Tad) * int(v < -63*mV) : second
        h_ainf = 1 / (1 + exp((v + 78*mV)/6/mV)) : 1

    Isyn_TC = IsynGABAA_RE_TC + IsynGABAB_RE_TC + IsynAMPA_PY_TC + IsynAMPA_stim_TC : amp * meter**-2
        IsynGABAA_RE_TC : amp * meter**-2
        IsynGABAB_RE_TC : amp * meter**-2
        IsynAMPA_PY_TC :  amp * meter**-2
        IsynAMPA_stim_TC :  amp * meter**-2
    
    Iext : amp * meter**-2
    
    Idummy = IEPSPs_L_VI_TCTCo + IEPSPs_L_V_TCTCo : amp * meter**-2
        IEPSPs_L_VI_TCTCo : amp * meter**-2
        IEPSPs_L_V_TCTCo : amp * meter**-2
    
    g_kl_TC : siemens * meter**-2
    
    '''

    
### SIMULATION and PLOT
close('all')    
start_scope()
prefs.codegen.target = "numpy"
defaultclock.dt = 0.02*ms

# Initialization
Bazhenov_TC = NeuronGroup(1,TC_eqs_modifs,method='rk4',threshold='v>20*mV',refractory=3*ms,events={'custom_poisson_PY':'v>1000*mV','custom_poisson_IN':'v>1000*mV','custom_poisson_inter':'v>1000*mV'})
Bazhenov_TC.v = -68*mV
Bazhenov_TC.m_na = 0.01
Bazhenov_TC.h_na = 0.99
Bazhenov_TC.m_t = 0.01
Bazhenov_TC.h_t = 0.01
Bazhenov_TC.n_k = 0.01
Bazhenov_TC.CA_i_TC = 1e-4*uM
Bazhenov_TC.P1 = 0.0
Bazhenov_TC.Op = 0.5
Bazhenov_TC.Op_L = 0.0
Bazhenov_TC.g_kl_TC = 0.03*msiemens*cm**-2 #0.03*msiemens*cm**-2

# Monitoring
V1 = StateMonitor(Bazhenov_TC,('v','I_na','I_k','I_t','I_kl','I_l','I_a','I_h','m_na','h_na','alpham_na','alphah_na','betam_na','betah_na','tau_m_na','mna_inf','tau_h_na','hna_inf'),record=True)

# Simulation
run(simulation_time, report='text',report_period=180*second)

# Plotting
# Figure 3A

a, b = 'Present model', 'Original model'
alphaa, linea, lineb = 0.60, 1.5, 1.3

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

fig, ax = plt.subplots(7, 1, sharex=True, figsize=(15, 30))
ax[0].plot(time_bazhenov, v_bazhenov[:int((simulation_time / 0.02) / ms), 15], color="black", linewidth=lineb)
configure_axis(ax[0], 'Membrane potential used for computations', 'mV', 25, 0, 3000)
#
ax[1].plot(V1.t / ms, V1.I_na[0] / (10 * amp * meter**-2), label=a, color="#2A52BE", linewidth=linea)
ax[1].plot(time_bazhenov, Ina_bazhenov[:int((simulation_time / 0.02) / ms)] / 1000, label=b, color="#4B9CD3", alpha=alphaa, linewidth=lineb)
configure_axis(ax[1], 'Fast sodium current, $I_{Na}$', 'mA/cm²', 0.25, 0, 3000, show_legend=True)
#
ax[2].plot(V1.t / ms, V1.I_k[0] / (10 * amp * meter**-2), label=a, color="#2A52BE", linewidth=linea)
ax[2].plot(time_bazhenov, Ik_bazhenov[:int((simulation_time / 0.02) / ms)] / 1000, label=b, color="#4B9CD3", linewidth=lineb, alpha=alphaa)
configure_axis(ax[2], 'Fast potassium current, $I_{K}$', 'mA/cm²', 0.10, 0, 3000, show_legend=True)
#
ax[3].plot(V1.t / ms, V1.I_t[0] / (10 * amp * meter**-2), label=a, color="#2A52BE", linewidth=linea)
ax[3].plot(time_bazhenov, It_bazhenov[:int((simulation_time / 0.02) / ms)] / 1000, label=b, color="#4B9CD3", linewidth=lineb, alpha=alphaa)
configure_axis(ax[3], 'Low threshold calcium current, $I_{t}$', 'mA/cm²', 0.005, 0, 3000, show_legend=True)
#
ax[4].plot(V1.t / ms, V1.I_h[0] / (10 * amp * meter**-2), label=a, color="#2A52BE", linewidth=linea)
ax[4].plot(time_bazhenov, Ih_bazhenov[:int((simulation_time / 0.02) / ms)] / 1000, label=b, color="#4B9CD3", linewidth=lineb, alpha=alphaa)
configure_axis(ax[4], 'Hyperpolarization activated cation current, $I_{h}$', 'mA/cm²', 0.0005, 0, 3000, show_legend=True)
#
ax[5].plot(V1.t / ms, V1.I_kl[0] / (10 * amp * meter**-2), label=a, color="#2A52BE", linewidth=linea)
ax[5].plot(time_bazhenov, -Ikl_bazhenov[:int((simulation_time / 0.02) / ms)] / 1000, label=b, color="#4B9CD3", linewidth=lineb, alpha=alphaa)
configure_axis(ax[5], 'Potassium leak current, $I_{Kl}$', 'mA/cm²', 0.002, 0, 3000, show_legend=True)
#
ax[6].plot(V1.t / ms, V1.I_l[0] / (10 * amp * meter**-2), label=a, color="#2A52BE", linewidth=linea)
ax[6].plot(time_bazhenov, -Il_bazhenov[:int((simulation_time / 0.02) / ms)] / 1000, label=b, color="#4B9CD3", linewidth=lineb, alpha=alphaa)
configure_axis(ax[6], 'Leak current, $I_{L}$', 'mA/cm²', 0.0005, 0, 3000, show_legend=True)
ax[6].set_xlabel('Time (ms)', size=30, labelpad=25)
#
plt.savefig('Figure4A.png', dpi=300, bbox_inches='tight')
plt.show()

# Figure 3B
fig, ax = plt.subplots(4, 1, sharex=True, figsize=(15, 25))
ax[0].plot(time_bazhenov, v_bazhenov[:int((simulation_time / 0.02) / ms), 15], color="black", linewidth=lineb)
configure_axis(ax[0], 'Membrane potential used for computations', 'mV', 25, 0, 500)
#
ax[1].plot(V1.t / ms, V1.I_na[0] / (10 * amp * meter**-2), label=a, color="#2A52BE", linewidth=linea)
ax[1].plot(time_bazhenov, Ina_bazhenov[:int((simulation_time / 0.02) / ms)] / 1000, label=b, color="#4B9CD3", alpha=alphaa, linewidth=lineb)
configure_axis(ax[1], 'Fast sodium current, $I_{Na}$', 'nA/nm²', 0.25, 1312, 1314, show_legend=True)
#
ax[2].plot(V1.t / ms, V1.m_na[0], label=a, color="#2A52BE", linewidth=linea)
ax[2].plot(time_bazhenov, mna_bazhenov[:int((simulation_time / 0.02) / ms)], label=b, color="#4B9CD3", linewidth=lineb, alpha=alphaa)
configure_axis(ax[2], 'Activation gating variable, $m$', '', 0.5, 1312, 1314, show_legend=True)
#
ax[3].plot(V1.t / ms, V1.h_na[0], label=a, color="#2A52BE", linewidth=linea)
ax[3].plot(time_bazhenov, hna_bazhenov[:int((simulation_time / 0.02) / ms)], label=b, color="#4B9CD3", linewidth=lineb, alpha=alphaa)
configure_axis(ax[3], 'Inactivation gating variable, $h$', '', 0.5, 1312, 1314, show_legend=True)
ax[3].set_xlabel('Time (ms)', size=30, labelpad=25)
ax[3].xaxis.set_major_locator(MultipleLocator(0.5))
#
plt.savefig('Figure4B.png', dpi=300, bbox_inches='tight')
plt.show()

# Figure 3C
fig, ax = plt.subplots(4, 1, sharex=True, figsize=(15, 25))
#
ax[0].plot(V1.t / ms, (V1.alpham_na[0])/1000, label=a, color="#2A52BE", linewidth=linea)
ax[0].plot(time_bazhenov, alphamna_bazhenov[:int((simulation_time / 0.02) / ms)], label=b, color="#4B9CD3", linewidth=lineb, alpha=alphaa)
configure_axis(ax[0], r'$\alpha_m$', 'kHz', 5, 1312, 1314, show_legend=True)
#
ax[1].plot(V1.t / ms, (V1.betam_na[0])/1000, label=a, color="#2A52BE", linewidth=linea)
ax[1].plot(time_bazhenov, betamna_bazhenov[:int((simulation_time / 0.02) / ms)], label=b, color="#4B9CD3", linewidth=lineb, alpha=alphaa)
configure_axis(ax[1], r'$\beta_m$', 'kHz', 5, 1312, 1314, show_legend=True)
#
ax[2].plot(V1.t / ms, (V1.tau_m_na[0])/0.001, label=a, color="#2A52BE", linewidth=linea)
ax[2].plot(time_bazhenov, taumna_bazhenov[:int((simulation_time / 0.02) / ms)], label=b, color="#4B9CD3", linewidth=lineb, alpha=alphaa)
configure_axis(ax[2], r'$\tau_m$', 'ms', 0.02, 1312, 1314, show_legend=True)
#
ax[3].plot(V1.t / ms, V1.mna_inf[0], label=a, color="#2A52BE", linewidth=linea)
ax[3].plot(time_bazhenov, mnainf_bazhenov[:int((simulation_time / 0.02) / ms)], label=b, color="#4B9CD3", linewidth=lineb, alpha=alphaa)
configure_axis(ax[3], r'$m \infty$', '', 0.5, 1312, 1314, show_legend=True)
ax[3].set_xlabel('Time (ms)', size=30, labelpad=25)
ax[3].xaxis.set_major_locator(MultipleLocator(0.5))
#
plt.savefig('Figure4C.png', dpi=300, bbox_inches='tight')
plt.show()

# Figure 3D
fig, ax = plt.subplots(4, 1, sharex=True, figsize=(15, 25))
#
ax[0].plot(V1.t / ms, (V1.alphah_na[0])/1000, label=a, color="#2A52BE", linewidth=linea)
ax[0].plot(time_bazhenov, alphahna_bazhenov[:int((simulation_time / 0.02) / ms)], label=b, color="#4B9CD3", linewidth=lineb, alpha=alphaa)
configure_axis(ax[0], r'$\alpha_h$', 'kHz', 1, 1312, 1314, show_legend=True)
#
ax[1].plot(V1.t / ms, (V1.betah_na[0])/1000, label=a, color="#2A52BE", linewidth=linea)
ax[1].plot(time_bazhenov, betahna_bazhenov[:int((simulation_time / 0.02) / ms)], label=b, color="#4B9CD3", linewidth=lineb, alpha=alphaa)
configure_axis(ax[1], r'$\beta_h$', 'kHz', 1, 1312, 1314, show_legend=True)
#
ax[2].plot(V1.t / ms, (V1.tau_h_na[0])/0.001, label=a, color="#2A52BE", linewidth=linea)
ax[2].plot(time_bazhenov, tauhna_bazhenov[:int((simulation_time / 0.02) / ms)], label=b, color="#4B9CD3", linewidth=lineb, alpha=alphaa)
configure_axis(ax[2], r'$\tau_h$', 'ms', 2, 1312, 1314, show_legend=True)
#
ax[3].plot(V1.t / ms, V1.hna_inf[0], label=a, color="#2A52BE", linewidth=linea)
ax[3].plot(time_bazhenov, hnainf_bazhenov[:int((simulation_time / 0.02) / ms)], label=b, color="#4B9CD3", linewidth=lineb, alpha=alphaa)
configure_axis(ax[3], r'$h \infty$', '', 0.5, 1312, 1314, show_legend=True)
ax[3].set_xlabel('Time (ms)', size=30, labelpad=25)
ax[3].xaxis.set_major_locator(MultipleLocator(0.5))
#
plt.savefig('Figure4D.png', dpi=300, bbox_inches='tight')
plt.show()