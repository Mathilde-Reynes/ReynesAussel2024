#!/usr/bin/env python
# coding: utf-8

from brian2 import *
import numpy as np
import matplotlib.pyplot as plt
from brian2.units.constants import *
import matplotlib.gridspec as gridspec
import gc
from Soma_eqs import *
from Soma_eqs_exp import *
from Dendritic_eqs import *
from TC_eqs import *
from RE_eqs import *
from Synapses import *
from Cortical_layer import *
from Thalamus import *
from Figure_conditions import *
from Figure_rawdata import *
from Figure_plotting import *

### Analyze propagation speed
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


###Create the complete model and run it
def thalamocortical_network(seed_value, analyze_speed, fig_number, raw_data, plot_figure):
    # Set seed for reproducibility
    seed(seed_value)
    
    # Close all plots
    close('all')
    start_scope()
    
    ### Standard parameters
    N = 100
    N_PY = N
    N_TC = N_PY // 2
    N_RE = N_PY // 2
    N_IN = N_PY // 4 
    
    # Conductances
    g_syn_ampa_tcpy = 0.0001 * msiemens
    g_syn_ampa_tcin = 0.0001 * msiemens
    g_syn_ampa_pytc = 0.000025 * msiemens
    g_syn_ampa_pyre = 0.00005 * msiemens 
    
    # Areas of the different neurons
    s_Soma_PYIN = 10**-6 * cm**2
    s_Dend_PY = 165 * s_Soma_PYIN
    s_Dend_IN = 50 * s_Soma_PYIN
    s_TC = 2.9E-4 * cm**2
    s_RE = 1.43e-4 * cm**2
    
    # Radius of connection
    TC_PY = 10 
    TC_IN = 2 
    PY_RE = 5 
    PY_TC = 5
    
    #Amplitudes for minis
    A_PY_PY = 0.00006*msiemens
    A_PY_IN = 0.000025*msiemens
    
    ### Creation of the substructures
    net = Network()
        
    print("Updating N and minis amplitudes before instantiating Thalamus & Cortical later for Figure n°: "+str(fig_number))
    A_PY_PY, A_PY_IN, N = figure_conditions_pre(str(fig_number),A_PY_PY, A_PY_IN, N)
    
    # Thalamus ("T")
    all_neurons_T, all_synapses_T, all_monitors_T = create_thalamic_subparts(N // 2)
    RE, TC = all_neurons_T  
    V1_RE, V2_TC, R1_RE, R2_TC, I1_RE, I2_TC = all_monitors_T
    
    net.add(all_neurons_T)
    net.add(all_synapses_T)
    net.add(all_monitors_T)
    
    # Layer Cortex
    if fig_number == "11-1":
        all_neurons, all_synapses, all_gap_junctions, all_monitors = create_cortical_layer(N,log,A_PY_PY,A_PY_IN)
    else:
        all_neurons, all_synapses, all_gap_junctions, all_monitors = create_cortical_layer(N,exp,A_PY_PY,A_PY_IN)
    PY_dendrite, PY_soma, IN_dendrite, IN_soma = all_neurons
    V1_PYd, V2_PYs, V3_INd, V4_INs, R2_PYs, R4_INs, I1_PYd, I2_INd, S1, S2, M0, M1 = all_monitors
    S_AMPA_PY_PY, S_AMPA_PY_IN, S_NMDA_PY_PY, S_NMDA_PY_IN, S_GABAA_IN_PY = all_synapses
    
    net.add(all_neurons)
    net.add(all_synapses)
    net.add(all_gap_junctions)
    net.add(all_monitors)
    print("Thalamus & Cortical layer initialized")

    print("Updating other parameters to plot Figure n°: "+str(fig_number))
    runtime, g_syn_ampa_tcpy, g_syn_ampa_tcin, g_syn_ampa_pytc, g_syn_ampa_pyre, modulation, base_rate, g_syn_ampa_stim, monitor_poisson = figure_conditions(
    fig_number, all_synapses, all_synapses_T, all_neurons_T, all_neurons, g_syn_ampa_tcpy, g_syn_ampa_tcin, g_syn_ampa_pytc, g_syn_ampa_pyre
    )
    
    ### Creation of the synapses
    # Cortico-thalamic synapses
    S_AMPA_PY_TC = syn_ampa_thal(PY_soma, TC, 'IsynAMPA_PY_TC', s_TC,
                                 'abs(floor(i*'+str(N_TC)+'/'+str(N_PY)+') -j)<='+str(PY_TC)+'',
                                 g_syn_ampa_pytc) 
    #In Brian2, i are the presynaptic neuron indices (PY) and j are the postsynaptic indices (TC).
    #Here, we scale i by (N_TC / N_PY) using floor() to map PY neurons to TC neurons,
    #ensuring each TC neuron receives input from nearby PY neurons within a defined range (PY_TC).
    S_AMPA_PY_TC.t_last_spike = -100 * ms
    net.add(S_AMPA_PY_TC)
    
    S_AMPA_PY_RE = syn_ampa_thal(PY_soma, RE, 'IsynAMPA_PY_RE', s_RE,
                                 'abs(floor(i*'+str(N_RE)+'/'+str(N_PY)+') -j)<='+str(PY_RE)+'',
                                 g_syn_ampa_pyre) 
    S_AMPA_PY_RE.t_last_spike = -100 * ms
    net.add(S_AMPA_PY_RE) 
    
    # Thalamo-cortical synapses
    S_AMPA_TC_PY = syn_ampa_thal(TC, PY_dendrite, 'IsynAMPA_TC_PY', s_Dend_PY,
                                 'abs(floor(i*'+str(N_PY)+'/'+str(N_TC)+') -j)<='+str(TC_PY)+'',
                                 g_syn_ampa_tcpy) 
    S_AMPA_TC_PY.t_last_spike = -1000 * ms
    net.add(S_AMPA_TC_PY)
    
    S_AMPA_TC_IN = syn_ampa_thal(TC, IN_dendrite, 'IsynAMPA_TC_IN', s_Dend_IN,
                                 'abs(floor(i*'+str(N_IN)+'/'+str(N_TC)+') -j)<='+str(TC_IN)+'',
                                 g_syn_ampa_tcin) 
    S_AMPA_TC_IN.t_last_spike = -1000 * ms
    net.add(S_AMPA_TC_IN)
    
    if fig_number in ["17-A1", "18-A1", "17-A2", "18-A2", "17-A3", "18-A3", "17-B1", "18-B1", "17-B2", "18-B2", "17-B3", "18-B3", "19"]: 
        Poisson_stim = PoissonGroup(N_TC // 4, rates='base_rate + 0.9 * base_rate * sin(2 * pi * modulation * t)')
        S_AMPA_stim_TC = syn_ampa_thal(Poisson_stim, TC, 'IsynAMPA_stim_TC', s_TC, 'j > 18 and j < 31', g_syn_ampa_stim)
        S_AMPA_stim_TC.t_last_spike = -1000 * ms
        monitor_poisson = SpikeMonitor(Poisson_stim)
        net.add([Poisson_stim, S_AMPA_stim_TC, monitor_poisson])
    print("Synapses initialized")

    ### Simulation
    # Define the parameters of the simulation
    np.seterr(all='raise')
    prefs.codegen.target = 'cython'
    
    # Create a TimedArray for initialization
    num_samples = int(runtime / defaultclock.dt)
    init_arr = zeros(num_samples)
    init_arr[0] = 1
    init_timedarray = TimedArray(init_arr, dt=defaultclock.dt)
       
    # Run the simulation
    print("Simulation ready to run")
       
    net.run(runtime, report='text', report_period=120*second)
    
    # If analyze speed is needed
    if analyze_speed:
        analyze_propagation_speed(R2_PYs)
    # Mean firing rate of PY neurons
    print('Mean firing rate of PY neurons:')
    print(len(R2_PYs.t) / N_PY / runtime)
    
    ###plot figures here
    if plot_figure:
         figure_plotting(str(fig_number),all_monitors,all_monitors_T,all_synapses,monitor_poisson,runtime, N)
    
    ### Raw data here
    if raw_data:
        figure_rawdata(str(fig_number),all_monitors,all_monitors_T,monitor_poisson)

    return net

### Call the function to run the simulation
if __name__ == "__main__":
    gc.collect()
    print("Forced garbage collection to free up memory")
    print(f"Objects in memory: {len(gc.get_objects())}")
    print("Start")
    seed_value = 4168
    analyze_speed = False #True of False
    fig_number = "5"
    raw_data = False #True or False
    plot_figure = True #True or False
    thalamocortical_network(seed_value,analyze_speed,fig_number,raw_data,plot_figure)
