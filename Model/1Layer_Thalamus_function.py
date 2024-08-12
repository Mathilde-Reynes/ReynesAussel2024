#!/usr/bin/env python
# coding: utf-8

from brian2 import *
import numpy as np
import matplotlib.pyplot as plt
from brian2.units.constants import *
import matplotlib.gridspec as gridspec
from Soma_eqs import *
from Dendritic_eqs import *
from TC_eqs import *
from RE_eqs import *
from Synapses import *
from Cortical_layer import *
from Thalamus import *

### Analyze propagation speed
def analyze_propagation_speed(raster_PY):
    list_up_states_beginning = []
    min_silence = 30 * ms  # up-states will be defined as groups of spikes more than min_silence apart
        
    list_up_states_speed = []  # we'll store the up-states propagation speed here
        
    for k in range(1, len(raster_PY.t)):
        if (raster_PY.t[k] - raster_PY.t[k - 1]) > min_silence or not list_up_states_beginning:
            list_up_states_beginning.append(k)  # we found an up state !
            try:
                first_neuron_spiking = raster_PY.i[k]  # the index of the first neuron spiking
                time_first_neuron_spiking = raster_PY.t[k]  # the time at which it spikes
                    
                last_neuron_spiking = argmax(
                    [raster_PY.t[k + argwhere(raster_PY.i[k:] == N)[0]] for N in range(N_PY)]
                )  # the array index of the first spike of the last neuron
                time_last_neuron_spiking = raster_PY.t[
                    k + argwhere(raster_PY.i[k:] == last_neuron_spiking)[0]
                ]  # the time at which the last neuron spikes
                    
                speed_first_to_last = abs(last_neuron_spiking - first_neuron_spiking) / (
                    (time_last_neuron_spiking - time_first_neuron_spiking) / second
                )  # propagation speed between 1st spiking neuron and neuron 0
                    
                list_up_states_speed.append(speed_first_to_last)
                    
                # Print for debugging
                print("Up state detected beginning at time t=" + str(time_first_neuron_spiking))
                print("The first neuron to spike is N=" + str(first_neuron_spiking))
                print("The last neuron to spike in this up-state is neuron N=" + str(last_neuron_spiking) +
                         " at time t=" + str(time_last_neuron_spiking))
                print("The propagation speed is: " + str(speed_first_to_last) + " cells/s\n")
                    
            except:
                pass
                
    print(str(len(list_up_states_speed)) + " upstates were found.")
    print("Mean propagation speed: " + str(mean(list_up_states_speed)))
    print("Std propagation speed: " + str(std(list_up_states_speed)))
    return list_up_states_speed


###Create the complete model and run it
def thalamocortical_network(seed_value, N, runtime, analyze_speed):
    # Set seed for reproducibility
    seed(seed_value)
    
    # Close all plots
    close('all')

    ### Parameters
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

    ### Creation of the substructures
    net = Network(collect())

    # Thalamus ("T")
    all_neurons_T, all_synapses_T, all_monitors_T = create_thalamic_subparts(N // 2)
    RE, TC = all_neurons_T  
    V1_RE, V2_TC, R1_RE, R2_TC, I1_RE, I2_TC = all_monitors_T
    
    net.add(all_neurons_T)
    net.add(all_synapses_T)
    net.add(all_monitors_T)
    
    # Layer Cortex
    all_neurons, all_synapses, all_gap_junctions, all_monitors = create_cortical_layer(N, 5)
    PY_dendrite, PY_soma, IN_dendrite, IN_soma = all_neurons
    V1_PYd, V2_PYs, V3_INd, V4_INs, R2_PYs, R4_INs, I1_PYd, I2_INd, S1, S2, M0, M1 = all_monitors
    S_AMPA_PY_PY, S_AMPA_PY_IN, S_NMDA_PY_PY, S_NMDA_PY_IN, S_GABAA_IN_PY = all_synapses
    
    net.add(all_neurons)
    net.add(all_synapses)
    net.add(all_gap_junctions)
    net.add(all_monitors)
        
    ### Creation of the synapses
    # Cortico-thalamic synapses
    S_AMPA_PY_TC = syn_ampa_thal(PY_soma, TC, 'IsynAMPA_PY_TC', s_TC,
                                 'abs(floor(i*'+str(N_TC)+'/'+str(N_PY)+') -j)<='+str(PY_TC)+'',
                                 g_syn_ampa_pytc) 
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
    net.run(runtime, report='text', report_period=120*second)

    # Analyze propagation speed if needed
    if analyze_speed:
        analyze_propagation_speed(R2_PYs)

    # Mean firing rate of PY neurons
    print('Mean firing rate of PY neurons:')
    print(len(R2_PYs.t) / N_PY / runtime)

    return net

### Call the function to run the simulation
if __name__ == "__main__":
    thalamocortical_network(4168,100,2000*ms,True)
    
###Export raw data
#Uncomment to save the raw data in .txt as needed

# savetxt('PY_v.txt',V2_PYs.v/mV) #All PY soma membrane potentials
# savetxt('PY_v.txt',V4_INs.v/mV) #All IN soma membrane potentials
# savetxt('PY_v.txt',V1_PYd.v/mV) #All PY dendrite membrane potentials
# savetxt('PY_v.txt',V3_INd.v/mV) #All IN dendrite membrane potentials
# savetxt('PY_v.txt',V2_TC.v/mV) #All TC membrane potentials
# savetxt('PY_v.txt',V1_RE.v/mV) #All RE membrane potentials
# savetxt('time.txt',V2_PYs.t/ms) #Time array

# savetxt('PY_v.txt',V2_PYs.v[N_PY//2]/mV) #One PY soma membrane potential
# savetxt('PY_v.txt',V4_INs.v[N_IN//2]/mV) #One IN soma membrane potential
# savetxt('PY_v.txt',V1_PYd.v[N_PY//2]/mV) #One PY dendrite membrane potential
# savetxt('PY_v.txt',V3_INd.v[N_IN//2]/mV) #One IN dendrite membrane potential
# savetxt('PY_v.txt',V2_TC.v[N_TC//2]/mV) #One TC membrane potential
# savetxt('PY_v.txt',V1_RE.v[N_RE//2]/mV) #One RE membrane potential
# savetxt('time.txt',V2_PYs.t[N_PY//2]/ms) #Time array

# savetxt('PY_v.txt',R2_PYs.i) #PY soma spikes
# savetxt('PY_v.txt',R4_INs.i) #IN soma spikes
# savetxt('PY_v.txt',R2_TC.i) #TC spikes
# savetxt('PY_v.txt',R1_RE.i) #RE spikes
