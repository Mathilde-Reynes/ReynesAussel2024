#!/usr/bin/env python
# coding: utf-8

from brian2 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from brian2.units.constants import *
import matplotlib.gridspec as gridspec
seed(4068)
from Soma_eqs import *
from Dendritic_eqs import *
from TC_eqs import *
from RE_eqs import *
from Synapses import *
from Cortical_layer import *
from Thalamus import *

close('all')
#clear_cache('cython')

###Parameters

#Number of pyramidal neurons
N = 100
N_PY = N
N_TC = N_PY/2
N_RE = N_PY/2
N_IN = N_PY/4 
#Conductances
g_syn_ampa_tcpy = 0.0001*msiemens
g_syn_ampa_tcin = 0.0001*msiemens
g_syn_ampa_pytc = 0.000025*msiemens
g_syn_ampa_pyre = 0.00005*msiemens
#Areas of the different neurons
s_Soma_PYIN = 10**-6*cm**2
s_Dend_PY = 165*s_Soma_PYIN
s_Dend_IN = 50*s_Soma_PYIN
s_TC = 2.9E-4*cm**2
s_RE = 1.43e-4*cm**2
#Radius of connection
TC_PY = 10 
TC_IN = 2 
PY_RE = 5 
PY_TC = 5 

###Creation of the substructures

net=Network(collect())
#Thalamus ("T")
all_neurons_T, all_synapses_T, all_monitors_T = create_thalamic_subparts(N/2)
RE,TC = all_neurons_T  
V1_RE,V2_TC,R1_RE,R2_TC,I1_RE,I2_TC = all_monitors_T
net.add(all_neurons_T)
net.add(all_synapses_T)
net.add(all_monitors_T) 
#Layer Cortex
all_neurons, all_synapses, all_gap_junctions, all_monitors = create_cortical_layer(N)
PY_dendrite, PY_soma, IN_dendrite, IN_soma = all_neurons
V1_PYd,V2_PYs,V3_INd,V4_INs,R2_PYs,R4_INs,I1_PYd,I2_INd,S1,S2,M0,M1 = all_monitors
S_AMPA_PY_PY,S_AMPA_PY_IN,S_NMDA_PY_PY,S_NMDA_PY_IN,S_GABAA_IN_PY = all_synapses
net.add(all_neurons)
net.add(all_synapses)
net.add(all_gap_junctions)
net.add(all_monitors)
        
###Creation of the synapses
        
#Cortico-thalamic synapses
S_AMPA_PY_TC = syn_ampa_thal(PY_soma,TC,'IsynAMPA_PY_TC',s_TC,'abs(floor(i*'+str(N_TC)+'/'+str(N_PY)+') -j)<='+str(PY_TC)+'',g_syn_ampa_pytc) 
S_AMPA_PY_TC.t_last_spike = -100*ms
net.add(S_AMPA_PY_TC)
S_AMPA_PY_RE = syn_ampa_thal(PY_soma,RE,'IsynAMPA_PY_RE',s_RE,'abs(floor(i*'+str(N_RE)+'/'+str(N_PY)+') -j)<='+str(PY_RE)+'',g_syn_ampa_pyre) 
S_AMPA_PY_RE.t_last_spike = -100*ms
net.add(S_AMPA_PY_RE) 
#Thalamo-cortical synapses
S_AMPA_TC_PY = syn_ampa_thal(TC,PY_dendrite,'IsynAMPA_TC_PY',s_Dend_PY,'abs(floor(i*'+str(N_PY)+'/'+str(N_TC)+') -j)<='+str(TC_PY)+'',g_syn_ampa_tcpy) 
S_AMPA_TC_PY.t_last_spike = -1000*ms
net.add(S_AMPA_TC_PY)
S_AMPA_TC_IN = syn_ampa_thal(TC,IN_dendrite,'IsynAMPA_TC_IN',s_Dend_IN,'abs(floor(i*'+str(N_IN)+'/'+str(N_TC)+') -j)<='+str(TC_IN)+'',g_syn_ampa_tcin) 
S_AMPA_TC_IN.t_last_spike = -1000*ms
net.add(S_AMPA_TC_IN)

###Change parameters as needed

#Figure 6
# TC.g_t_TC = 2*msiemens*cm**-2 

#Figure 10 (N=20,60,100 or minis=120%,140%,170%)
# N = 20
# g_syn_ampa_tcpy = 0*msiemens 
# g_syn_ampa_tcin = 0*msiemens 
# g_syn_ampa_pytc = 0*msiemens 
# g_syn_ampa_pyre = 0*msiemens
#to implement in Cortical_layer.py
# A_PY_PY = 0.00006*msiemens*1.4
# A_PY_IN = 0.000025*msiemens*1.4

#Figure 11 (minis=50% or cortical synapses 50%)
# g_syn_ampa_tcpy = 0*msiemens 
# g_syn_ampa_tcin = 0*msiemens 
# g_syn_ampa_pytc = 0*msiemens 
# g_syn_ampa_pyre = 0*msiemens
#to implement in Cortical_layer.py
# A_PY_PY = 0.00006*msiemens*0.5
# A_PY_IN = 0.000025*msiemens*0.5
# g_syn_ampa_pypy = 0.00015*msiemens*0.5
# g_syn_nmda_pypy = 0.00001*msiemens*0.5
# g_syn_ampa_pyin = 0.00005*msiemens*0.5
# g_syn_nmda_pyin = 0.000008*msiemens*0.5
# g_syn_gabaa_inpy = 0.00005*msiemens*0.5

#Figure 12 (exp or log mini arrival)
# changes to implement in Soma_eqs.py

#Figure 13 (g PY-IN)
# g_syn_ampa_tcpy = 0*msiemens 
# g_syn_ampa_tcin = 0*msiemens 
# g_syn_ampa_pytc = 0*msiemens 
# g_syn_ampa_pyre = 0*msiemens
# g_syn_ampa_pyin = 0.00007*msiemens #or 0.00002*msiemens

#Figure 14 (disconnected versus connected)
# g_syn_ampa_tcpy = 0*msiemens 
# g_syn_ampa_tcin = 0*msiemens 
# g_syn_ampa_pytc = 0*msiemens 
# g_syn_ampa_pyre = 0*msiemens

#Figure 16
#Strong PYPY: gPYPY=0.15uS, gRETC=0.2, gTC-RE=0.4
# PY_dendrite.g_kl=0*msiemens*cm**-2
# TC.g_kl_TC = 0*msiemens*cm**-2
# RE.g_kl_RE = 0*msiemens*cm**-2
# syn_PYPY=all_synapses[0]
# syn_PYPY.g_syn=0.00015*msiemens
# syn_RETC=all_synapses_T[0]
# syn_RETC.g_syn=0.0002*msiemens
# syn_TCRE=all_synapses_T[-1]
# syn_TCRE.g_syn=0.0004*msiemens
#Weak PYPY : gPYPY=0.09uS, gRETC=0.2, gTC-RE=0.4
# PY_dendrite.g_kl=0*msiemens*cm**-2
# TC.g_kl_TC = 0*msiemens*cm**-2
# RE.g_kl_RE = 0*msiemens*cm**-2
# syn_PYPY=all_synapses[0]
# syn_PYPY.g_syn=0.00009*msiemens
# syn_RETC=all_synapses_T[0]
# syn_RETC.g_syn=0.0002*msiemens
# syn_TCRE=all_synapses_T[-1]
# syn_TCRE.g_syn=0.0004*msiemens
#Weak PYRETC : gPYPY=0.09uS, gRETC=0.1, gTC-RE=0.2
# PY_dendrite.g_kl=0*msiemens*cm**-2
# TC.g_kl_TC = 0*msiemens*cm**-2
# RE.g_kl_RE = 0*msiemens*cm**-2
# syn_PYPY=all_synapses[0]
# syn_PYPY.g_syn=0.00009*msiemens
# syn_RETC=all_synapses_T[0]
# syn_RETC.g_syn=0.0001*msiemens
# syn_TCRE=all_synapses_T[-1]
# syn_TCRE.g_syn=0.0002*msiemens

#Figure 17 (transition from sleep to (wakefulness)
#A: gPYPY=0.15uS, gRETC=0.2, gTC-RE=0.4, gKL=0.3 #v2 : gKL=0.25 in cortex
# PY_dendrite.g_kl=0.0025*msiemens*cm**-2 #0.003*msiemens*cm**-2
# TC.g_kl_TC = 0.003*msiemens*cm**-2
# syn_PYPY=all_synapses[0]
# syn_PYPY.g_syn=0.00015*msiemens
# syn_RETC=all_synapses_T[0]
# syn_RETC.g_syn=0.0002*msiemens
# syn_TCRE=all_synapses_T[-1]
# syn_TCRE.g_syn=0.0004*msiemens
#B:
# PY_dendrite.g_kl=0.0025*5/6*msiemens*cm**-2
# TC.g_kl_TC = 0.0025*msiemens*cm**-2
# syn_PYPY=all_synapses[0]
# syn_PYPY.g_syn=0.00013833*msiemens
# syn_RETC=all_synapses_T[0]
# syn_RETC.g_syn=0.00018333*msiemens
# syn_TCRE=all_synapses_T[-1]
# syn_TCRE.g_syn=0.0003666*msiemens
#C:
# PY_dendrite.g_kl=0.0025*4/6*msiemens*cm**-2
# TC.g_kl_TC = 0.002*msiemens*cm**-2
# syn_PYPY=all_synapses[0]
# syn_PYPY.g_syn=0.0001266*msiemens
# syn_RETC=all_synapses_T[0]
# syn_RETC.g_syn=0.0001666*msiemens
# syn_TCRE=all_synapses_T[-1]
# syn_TCRE.g_syn=0.000333*msiemens
#D:
# PY_dendrite.g_kl=0.0025*3/6*msiemens*cm**-2
# TC.g_kl_TC = 0.0015*msiemens*cm**-2
# syn_PYPY=all_synapses[0]
# syn_PYPY.g_syn=0.000115*msiemens
# syn_RETC=all_synapses_T[0]
# syn_RETC.g_syn=0.00015*msiemens
# syn_TCRE=all_synapses_T[-1]
# syn_TCRE.g_syn=0.0003*msiemens
#E:
# PY_dendrite.g_kl=0.0025*2/6*msiemens*cm**-2
# TC.g_kl_TC = 0.001*msiemens*cm**-2
# syn_PYPY=all_synapses[0]
# syn_PYPY.g_syn=0.0001033*msiemens
# syn_RETC=all_synapses_T[0]
# syn_RETC.g_syn=0.0001333*msiemens
# syn_TCRE=all_synapses_T[-1]
# syn_TCRE.g_syn=0.0002666*msiemens
#F:
# PY_dendrite.g_kl=0.00025*1/6*msiemens*cm**-2
# TC.g_kl_TC = 0.0005*msiemens*cm**-2
# syn_PYPY=all_synapses[0]
# syn_PYPY.g_syn=0.00009166*msiemens
# syn_RETC=all_synapses_T[0]
# syn_RETC.g_syn=0.0001166*msiemens
# syn_TCRE=all_synapses_T[-1]
# syn_TCRE.g_syn=0.0002333*msiemens
#G:"Activated state"
# PY_dendrite.g_kl=0*msiemens*cm**-2
# TC.g_kl_TC = 0*msiemens*cm**-2
# syn_PYPY=all_synapses[0]
# syn_PYPY.g_syn=0.00008*msiemens
# syn_RETC=all_synapses_T[0]
# syn_RETC.g_syn=0.0001*msiemens
# syn_TCRE=all_synapses_T[-1]
# syn_TCRE.g_syn=0.0002*msiemens

#Figures 18-20: with stimulation to 25% of TC cells, which is about 12 cells
# base_rate=25*Hz
# modulation=2.5*Hz #0.4 or 1 or 2.5 Hz
# g_syn_ampa_stim = 0.0004*msiemens #g_syn_ampa_tcre=0.0004*msiemens
# Poisson_stim=PoissonGroup(N_TC//4, rates='base_rate+0.9*base_rate*sin(2*pi*modulation*t)')
# S_AMPA_stim_TC = syn_ampa_thal(Poisson_stim,TC_TCo,'IsynAMPA_stim_TC',s_TC,'j>18 and j<31',g_syn_ampa_stim,10) 
# S_AMPA_stim_TC.t_last_spike = -1000*ms
# monitor_poisson=SpikeMonitor(Poisson_stim)
# net.add([Poisson_stim,S_AMPA_stim_TC,monitor_poisson])


###Simulation

#Define the parameters of the simulation
runtime=30*second
np.seterr(all='raise')
prefs.codegen.target = 'cython'
#
num_samples = int(runtime/defaultclock.dt)
init_arr = zeros(num_samples)
init_arr[0]=1
init_timedarray = TimedArray(init_arr, dt=defaultclock.dt)

#Run the simulation
net.run(runtime,report='text',report_period=120*second)

#Function to estimate the propagation speed of cortical up-states. 
def analyze_propagation_speed(raster_PY):
    #we need to : 1-detect the up states and 2-compute the propagation speed in each up state
    list_up_states_beginning=[]
    min_silence=30*msecond #up-states will be defined as groups of spikes more than min_silence apart
    
    list_up_states_speed=[] #we'll store the up-states propagation speed here
    
    for k in range(len(raster_PY.t[1:])):
        if (raster_PY.t[k]-raster_PY.t[k-1])>min_silence or list_up_states_beginning==[]:
            list_up_states_beginning.append(k) #we found an up state !
            try : #if the simulations ends in an up-state there can be bugs
                
                first_neuron_spiking=raster_PY.i[k] #the index of the first neuron spiking
                time_first_neuron_spiking=raster_PY.t[k] #the time at which it spikes
                
                #let's compute the time at which each neuron spikes for the first time in the up-state (we assume they all spike in the up state)
                #for each neuron N, the first time it spikes in the up-state is raster_PY.t[k+argwhere(raster_PY.i[k:]==N)[0]]
                last_neuron_spiking=argmax([raster_PY.t[k+argwhere(raster_PY.i[k:]==N)[0]] for N in range(100)]) #the array index of the first spike of the last neuron
                time_last_neuron_spiking=raster_PY.t[k+argwhere(raster_PY.i[k:]==last_neuron_spiking)[0]] #the time at which the last neuron spikes
                speed_first_to_last= abs(last_neuron_spiking-first_neuron_spiking)/((time_last_neuron_spiking-time_first_neuron_spiking)/second) #propagation speed between 1st spiking neuron and neuron 0
                
                list_up_states_speed.append(speed_first_to_last)
                
                #some prints to check if everything looks fine
                print("Up state detected beginning at time t="+str(time_first_neuron_spiking))
                print("The first neuron to spike is N="+str(first_neuron_spiking))
                print("The last neuron to spike in this up-state is neuron N="+str(last_neuron_spiking)+" at time t="+str(time_last_neuron_spiking))
                print("The propagation speed is: "+str(speed_first_to_last)+" cells/s")
                print(" ")
            except : 
                pass
            
    print(str(len(list_up_states_speed))+" upstates were found.")
    print("Mean propagation speed: "+str(mean(list_up_states_speed)))
    print("Std propagation speed: "+str(std(list_up_states_speed)))
    return

analyze_propagation_speed(R2_PYs)

#Mean firing rate of PY neurons
print('Mean firing rate of PY neurons:')
print(len(R2_PYs.t)/N_PY/runtime)


###Export raw data

#Uncomment to save the raw data in .txt as needed

# savetxt('PY_v.txt',V2_PYs.v/mV) #All PY soma membrane potentials
# savetxt('IN_v.txt',V4_INs.v/mV) #All IN soma membrane potentials
# savetxt('PYdend_v.txt',V1_PYd.v/mV) #All PY dendrite membrane potentials
# savetxt('INdend_v.txt',V3_INd.v/mV) #All IN dendrite membrane potentials
# savetxt('TC_v.txt',V2_TC.v/mV) #All TC membrane potentials
# savetxt('RE_v.txt',V1_RE.v/mV) #All RE membrane potentials
# savetxt('time30.txt',V2_PYs.t/ms) #Time array

# savetxt('PY_v_c.txt',V2_PYs.v[N_PY/2]/mV) #One PY soma membrane potential
# savetxt('IN_v_c.txt',V4_INs.v[N_IN/2]/mV) #One IN soma membrane potential
# savetxt('PYdend_v_c.txt',V1_PYd.v[N_PY/2]/mV) #One PY dendrite membrane potential
# savetxt('INdend_v_c.txt',V3_INd.v[N_IN/2]/mV) #One IN dendrite membrane potential
# savetxt('TC_v_c.txt',V2_TC.v[N_TC/2]/mV) #One TC membrane potential
# savetxt('RE_v_c.txt',V1_RE.v[N_RE/2]/mV) #One RE membrane potential
# savetxt('time_c.txt',V2_PYs.t[N_PY/2]/ms) #Time array

# savetxt('PY_raster.txt',R2_PYs.i) #PY soma spikes
# savetxt('IN_v.txt',R4_INs.i) #IN soma spikes
# savetxt('TC_v.txt',R2_TC.i) #TC spikes
# savetxt('RE_v.txt',R1_RE.i) #RE spikes

# ###Figure 6
# fig,ax = subplots(2,1, sharex = True,figsize=(19,15))
# ax[0].spines['top'].set_visible(False)
# ax[0].spines['right'].set_visible(False)
# ax[0].spines['bottom'].set_visible(True)
# ax[0].spines['left'].set_visible(True)
# ax[0].plot(R2_PYs.t/second, R2_PYs.i,'.',markersize=2,alpha=0.5,color="tab:blue")
# ax[0].set_ylabel('Neuron index',fontsize=30)
# ax[0].set_xlabel('Time (s)',size=30)
# ax[0].tick_params(axis='both', which='major', labelsize=25, width=2)
# ax[0].set_title('PY', size=30, loc='left')
# ax[1].plot(R4_INs.t/second, R4_INs.i, '.',markersize=2,alpha=0.5,color="tab:green")
# ax[1].spines['top'].set_visible(False)
# ax[1].spines['right'].set_visible(False)
# ax[1].spines['bottom'].set_visible(True)
# ax[1].spines['left'].set_visible(True)
# ax[1].set_ylabel('Neuron index',fontsize=30)
# ax[1].set_xlabel('Time (s)',size=30)
# ax[1].tick_params(axis='both', which='major', labelsize=25, width=2)
# ax[1].set_title('IN', size=30, loc='left')
# ax[1].set_xlim([0,10])
# fig.suptitle('Cortical cells raster plot',fontsize=30)
# fig.tight_layout()
# #plt.savefig('Figure6Disconnected', dpi=300)