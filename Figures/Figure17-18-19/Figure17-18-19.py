# -*- coding: utf-8 -*-


from brian2 import *
from scipy import signal
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

#Open files
#Input modulated at 0.4Hz, SWS
raster_PY_i_04Hz_sws=loadtxt('../../Data/stimfreq_04_sws/raster_PY_i.txt')
raster_PY_t_04Hz_sws=loadtxt('../../Data/stimfreq_04_sws/raster_PY_t.txt')
raster_TC_i_04Hz_sws=loadtxt('../../Data/stimfreq_04_sws/raster_TC_i.txt')
raster_TC_t_04Hz_sws=loadtxt('../../Data/stimfreq_04_sws/raster_TC_t.txt')
raster_input_i_04Hz_sws=loadtxt('../../Data/stimfreq_04_sws/input_i.txt')
raster_input_t_04Hz_sws=loadtxt('../../Data/stimfreq_04_sws/input_t.txt')

#Input modulated at 1Hz, SWS
raster_PY_i_1Hz_sws=loadtxt('../../Data/stimfreq_1_sws/raster_PY_i.txt')
raster_PY_t_1Hz_sws=loadtxt('../../Data/stimfreq_1_sws/raster_PY_t.txt')
raster_TC_i_1Hz_sws=loadtxt('../../Data/stimfreq_1_sws/raster_TC_i.txt')
raster_TC_t_1Hz_sws=loadtxt('../../Data/stimfreq_1_sws/raster_TC_t.txt')
raster_input_i_1Hz_sws=loadtxt('../../Data/stimfreq_1_sws/input_i.txt')
raster_input_t_1Hz_sws=loadtxt('../../Data/stimfreq_1_sws/input_t.txt')

#Input modulated at 2.5Hz, SWS
raster_PY_i_25Hz_sws=loadtxt('../../Data/stimfreq_25_sws/raster_PY_i.txt')
raster_PY_t_25Hz_sws=loadtxt('../../Data/stimfreq_25_sws/raster_PY_t.txt')
raster_TC_i_25Hz_sws=loadtxt('../../Data/stimfreq_25_sws/raster_TC_i.txt')
raster_TC_t_25Hz_sws=loadtxt('../../Data/stimfreq_25_sws/raster_TC_t.txt')
raster_input_i_25Hz_sws=loadtxt('../../Data/stimfreq_25_sws/input_i.txt')
raster_input_t_25Hz_sws=loadtxt('../../Data/stimfreq_25_sws/input_t.txt')

#Input modulated at 0.4Hz, activated state
raster_PY_i_04Hz_wake=loadtxt('../../Data/stimfreq_04_wake/raster_PY_i.txt')
raster_PY_t_04Hz_wake=loadtxt('../../Data/stimfreq_04_wake/raster_PY_t.txt')
raster_TC_i_04Hz_wake=loadtxt('../../Data/stimfreq_04_wake/raster_TC_i.txt')
raster_TC_t_04Hz_wake=loadtxt('../../Data/stimfreq_04_wake/raster_TC_t.txt')
raster_input_i_04Hz_wake=loadtxt('../../Data/stimfreq_04_wake/input_i.txt')
raster_input_t_04Hz_wake=loadtxt('../../Data/stimfreq_04_wake/input_t.txt')

#Input modulated at 1Hz, activated state
raster_PY_i_1Hz_wake=loadtxt('../../Data/stimfreq_1_wake/raster_PY_i.txt')
raster_PY_t_1Hz_wake=loadtxt('../../Data/stimfreq_1_wake/raster_PY_t.txt')
raster_TC_i_1Hz_wake=loadtxt('../../Data/stimfreq_1_wake/raster_TC_i.txt')
raster_TC_t_1Hz_wake=loadtxt('../../Data/stimfreq_1_wake/raster_TC_t.txt')
raster_input_i_1Hz_wake=loadtxt('../../Data/stimfreq_1_wake/input_i.txt')
raster_input_t_1Hz_wake=loadtxt('../../Data/stimfreq_1_wake/input_t.txt')

#Input modulated at 2.5Hz, activated state
raster_PY_i_25Hz_wake=loadtxt('../../Data/stimfreq_25_wake/raster_PY_i.txt')
raster_PY_t_25Hz_wake=loadtxt('../../Data/stimfreq_25_wake/raster_PY_t.txt')
raster_TC_i_25Hz_wake=loadtxt('../../Data/stimfreq_25_wake/raster_TC_i.txt')
raster_TC_t_25Hz_wake=loadtxt('../../Data/stimfreq_25_wake/raster_TC_t.txt')
raster_input_i_25Hz_wake=loadtxt('../../Data/stimfreq_25_wake/input_i.txt')
raster_input_t_25Hz_wake=loadtxt('../../Data/stimfreq_25_wake/input_t.txt')

#Compute "Running Spike Histograms" (RHS)
#The time window used to compute the RSH is not provided, nor the overlap between these windows
#to get something similar to Figure10 in the paper, let's use
#time_window=100ms
#overlap=50ms
OVERLAP=50

def RHS(raster_i,raster_t,neuron_range,time_window=100,overlap=OVERLAP):
    #first, remove all spikes that are not emitted by the neurons in neuron_range from raster
    range_raster=where(logical_and(raster_i>=neuron_range[0],raster_i<=neuron_range[-1]))[0]
    small_raster_i=raster_i[range_raster]
    small_raster_t=raster_t[range_raster]
    
    N_neurons=len(neuron_range)
    
    #then, compute histogram
    rhs=[]
    for t in arange(0,15000,overlap): #simulation time goes from 0ms to 15000ms
        rhs.append(len(where(logical_and(small_raster_t>=t, small_raster_t<t+time_window))[0])/N_neurons)
    return rhs

rhs_PY_04Hz_sws=RHS(raster_PY_i_04Hz_sws,raster_PY_t_04Hz_sws,list(range(38,63)))
rhs_TC_04Hz_sws=RHS(raster_TC_i_04Hz_sws,raster_TC_t_04Hz_sws,list(range(19,32)))
rhs_input_04Hz_sws=RHS(raster_input_i_04Hz_sws,raster_input_t_04Hz_sws,list(range(12)))
rhs_PY_1Hz_sws=RHS(raster_PY_i_1Hz_sws,raster_PY_t_1Hz_sws,list(range(38,63)))
rhs_TC_1Hz_sws=RHS(raster_TC_i_1Hz_sws,raster_TC_t_1Hz_sws,list(range(19,32)))
rhs_input_1Hz_sws=RHS(raster_input_i_1Hz_sws,raster_input_t_1Hz_sws,list(range(12)))
rhs_PY_25Hz_sws=RHS(raster_PY_i_25Hz_sws,raster_PY_t_25Hz_sws,list(range(38,63)))
rhs_TC_25Hz_sws=RHS(raster_TC_i_25Hz_sws,raster_TC_t_25Hz_sws,list(range(19,32)))
rhs_input_25Hz_sws=RHS(raster_input_i_25Hz_sws,raster_input_t_25Hz_sws,list(range(12)))

rhs_PY_04Hz_wake=RHS(raster_PY_i_04Hz_wake,raster_PY_t_04Hz_wake,list(range(38,63)))
rhs_TC_04Hz_wake=RHS(raster_TC_i_04Hz_wake,raster_TC_t_04Hz_wake,list(range(19,32)))
rhs_input_04Hz_wake=RHS(raster_input_i_04Hz_wake,raster_input_t_04Hz_wake,list(range(12)))
rhs_PY_1Hz_wake=RHS(raster_PY_i_1Hz_wake,raster_PY_t_1Hz_wake,list(range(38,63)))
rhs_TC_1Hz_wake=RHS(raster_TC_i_1Hz_wake,raster_TC_t_1Hz_wake,list(range(19,32)))
rhs_input_1Hz_wake=RHS(raster_input_i_1Hz_wake,raster_input_t_1Hz_wake,list(range(12)))
rhs_PY_25Hz_wake=RHS(raster_PY_i_25Hz_wake,raster_PY_t_25Hz_wake,list(range(38,63)))
rhs_TC_25Hz_wake=RHS(raster_TC_i_25Hz_wake,raster_TC_t_25Hz_wake,list(range(19,32)))
rhs_input_25Hz_wake=RHS(raster_input_i_25Hz_wake,raster_input_t_25Hz_wake,list(range(12)))

close('all')

#Figure17
time_rhs=arange(0,15000,OVERLAP)/1000
plt.figure(figsize=(20, 20))
plt.subplots_adjust(hspace=0.5)

# Create a helper function to handle plotting with fill_between
def plot_with_fill(ax, x, y, label, color='black', show_xlabel=False, show_xticks=True):
    ax.plot(x, y, color=color)
    ax.fill_between(x, y, color=color, alpha=0.8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='x', labelsize=25 if show_xticks else 0)  # x-axis tick labels
    ax.tick_params(axis='y', labelsize=25)  # y-axis tick labels
    if show_xlabel:
        ax.set_xlabel('Time (s)', fontsize=30)
    else:
        ax.set_xlabel('')  # Remove xlabel if not showing

ax1 = plt.subplot(631)
plot_with_fill(ax1, time_rhs, rhs_PY_04Hz_wake, 'PY', show_xlabel=False, show_xticks=False)
ax1.set_ylabel('PY', fontsize=30)
ax1.set_title('0.4Hz \n Activated', fontsize=35)

plt.subplot(634, sharex=ax1)
plot_with_fill(plt.gca(), time_rhs, rhs_TC_04Hz_wake, 'TC', show_xlabel=False, show_xticks=False)
plt.gca().set_ylabel('TC', fontsize=30)

plt.subplot(637, sharex=ax1)
plot_with_fill(plt.gca(), time_rhs, rhs_input_04Hz_wake, 'Input', show_xlabel=False, show_xticks=False)
plt.gca().set_ylabel('Input', fontsize=30)

plt.subplot(632, sharex=ax1)
plot_with_fill(plt.gca(), time_rhs, rhs_PY_1Hz_wake, 'PY', show_xlabel=False, show_xticks=False)
plt.gca().set_title('1Hz \n Activated', fontsize=35)

plt.subplot(635, sharex=ax1)
plot_with_fill(plt.gca(), time_rhs, rhs_TC_1Hz_wake, 'TC', show_xlabel=False, show_xticks=False)

plt.subplot(638, sharex=ax1)
plot_with_fill(plt.gca(), time_rhs, rhs_input_1Hz_wake, 'Input', show_xlabel=False, show_xticks=False)

plt.subplot(633, sharex=ax1)
plot_with_fill(plt.gca(), time_rhs, rhs_PY_25Hz_wake, 'PY', show_xlabel=False, show_xticks=False)
plt.gca().set_title('2.5Hz \n Activated', fontsize=35)

plt.subplot(636, sharex=ax1)
plot_with_fill(plt.gca(), time_rhs, rhs_TC_25Hz_wake, 'TC', show_xlabel=False, show_xticks=False)

plt.subplot(639, sharex=ax1)
plot_with_fill(plt.gca(), time_rhs, rhs_input_25Hz_wake, 'Input', show_xlabel=False, show_xticks=False)

plt.subplot(6, 3, 10, sharex=ax1)
plot_with_fill(plt.gca(), time_rhs, rhs_PY_04Hz_sws, 'PY', show_xlabel=False, show_xticks=False)
plt.gca().set_ylabel('PY', fontsize=30)
plt.gca().set_title('SWS', fontsize=35)

plt.subplot(6, 3, 13, sharex=ax1)
plot_with_fill(plt.gca(), time_rhs, rhs_TC_04Hz_sws, 'TC', show_xlabel=False, show_xticks=False)
plt.gca().set_ylabel('TC', fontsize=30)

plt.subplot(6, 3, 16, sharex=ax1)
plot_with_fill(plt.gca(), time_rhs, rhs_input_04Hz_sws, 'Input', show_xlabel=True, show_xticks=True)
plt.gca().set_ylabel('Input', fontsize=30)

plt.subplot(6, 3, 11, sharex=ax1)
plot_with_fill(plt.gca(), time_rhs, rhs_PY_1Hz_sws, 'PY', show_xlabel=False, show_xticks=False)
plt.gca().set_title('SWS', fontsize=35)

plt.subplot(6, 3, 14, sharex=ax1)
plot_with_fill(plt.gca(), time_rhs, rhs_TC_1Hz_sws, 'TC', show_xlabel=False, show_xticks=False)

plt.subplot(6, 3, 17, sharex=ax1)
plot_with_fill(plt.gca(), time_rhs, rhs_input_1Hz_sws, 'Input', show_xlabel=True, show_xticks=True)

plt.subplot(6, 3, 12, sharex=ax1)
plot_with_fill(plt.gca(), time_rhs, rhs_PY_25Hz_sws, 'PY', show_xlabel=False, show_xticks=False)
plt.gca().set_title('SWS', fontsize=35)

plt.subplot(6, 3, 15, sharex=ax1)
plot_with_fill(plt.gca(), time_rhs, rhs_TC_25Hz_sws, 'TC', show_xlabel=False, show_xticks=False)

plt.subplot(6, 3, 18, sharex=ax1)
plot_with_fill(plt.gca(), time_rhs, rhs_input_25Hz_sws, 'Input', show_xlabel=True, show_xticks=True)

plt.xlim(5, 15)

plt.savefig('Figure17.png', dpi=300, bbox_inches='tight')
plt.show()


#For figure 18, we're going to need the crosscorrelation of PY and input RSHs
#and the power spectra of all RSHs

crosscorr_04Hz_sws=signal.correlate(rhs_PY_04Hz_sws,rhs_input_04Hz_sws)
crosscorr_1Hz_sws=signal.correlate(rhs_PY_1Hz_sws,rhs_input_1Hz_sws)
crosscorr_25Hz_sws=signal.correlate(rhs_PY_25Hz_sws,rhs_input_25Hz_sws)
crosscorr_04Hz_wake=signal.correlate(rhs_PY_04Hz_wake,rhs_input_04Hz_wake)
crosscorr_1Hz_wake=signal.correlate(rhs_PY_1Hz_wake,rhs_input_1Hz_wake)
crosscorr_25Hz_wake=signal.correlate(rhs_PY_25Hz_wake,rhs_input_25Hz_wake)

lags = signal.correlation_lags(len(rhs_PY_04Hz_sws), len(rhs_input_04Hz_sws))
# print(lags[0])

rhs_sampling_freq=1/(time_rhs[1]) #sampling frequency of the running spike histogram in Hz
lags= lags/rhs_sampling_freq #lags in second
# print(lags[0])

freqs,Spectrum_PY_04_sws=signal.periodogram(rhs_PY_04Hz_sws, rhs_sampling_freq, scaling='spectrum')
_,Spectrum_TC_04_sws=signal.periodogram(rhs_TC_04Hz_sws, rhs_sampling_freq, scaling='spectrum')
_,Spectrum_input_04_sws=signal.periodogram(rhs_input_04Hz_sws, rhs_sampling_freq, scaling='spectrum')

_,Spectrum_PY_1_sws=signal.periodogram(rhs_PY_1Hz_sws, rhs_sampling_freq, scaling='spectrum')
_,Spectrum_TC_1_sws=signal.periodogram(rhs_TC_1Hz_sws, rhs_sampling_freq, scaling='spectrum')
_,Spectrum_input_1_sws=signal.periodogram(rhs_input_1Hz_sws, rhs_sampling_freq, scaling='spectrum')

_,Spectrum_PY_25_sws=signal.periodogram(rhs_PY_25Hz_sws, rhs_sampling_freq, scaling='spectrum')
_,Spectrum_TC_25_sws=signal.periodogram(rhs_TC_25Hz_sws, rhs_sampling_freq, scaling='spectrum')
_,Spectrum_input_25_sws=signal.periodogram(rhs_input_25Hz_sws, rhs_sampling_freq, scaling='spectrum')

_,Spectrum_PY_04_wake=signal.periodogram(rhs_PY_04Hz_wake, rhs_sampling_freq, scaling='spectrum')
_,Spectrum_TC_04_wake=signal.periodogram(rhs_TC_04Hz_wake, rhs_sampling_freq, scaling='spectrum')
_,Spectrum_input_04_wake=signal.periodogram(rhs_input_04Hz_wake, rhs_sampling_freq, scaling='spectrum')

_,Spectrum_PY_1_wake=signal.periodogram(rhs_PY_1Hz_wake, rhs_sampling_freq, scaling='spectrum')
_,Spectrum_TC_1_wake=signal.periodogram(rhs_TC_1Hz_wake, rhs_sampling_freq, scaling='spectrum')
_,Spectrum_input_1_wake=signal.periodogram(rhs_input_1Hz_wake, rhs_sampling_freq, scaling='spectrum')

_,Spectrum_PY_25_wake=signal.periodogram(rhs_PY_25Hz_wake, rhs_sampling_freq, scaling='spectrum')
_,Spectrum_TC_25_wake=signal.periodogram(rhs_TC_25Hz_wake, rhs_sampling_freq, scaling='spectrum')
_,Spectrum_input_25_wake=signal.periodogram(rhs_input_25Hz_wake, rhs_sampling_freq, scaling='spectrum')


#Figure18
# Define the figure and adjust spacing
plt.figure(figsize=(22, 25))
plt.subplots_adjust(hspace=1.4, wspace=0.3)

# Compute range values for cross-correlation data
rangey_04_wake = [
    min(crosscorr_04Hz_wake[int(len(crosscorr_04Hz_wake)//2 - 5*rhs_sampling_freq):
                            int(len(crosscorr_04Hz_wake)//2 + 5*rhs_sampling_freq)]),
    max(crosscorr_04Hz_wake[int(len(crosscorr_04Hz_wake)//2 - 5*rhs_sampling_freq):
                            int(len(crosscorr_04Hz_wake)//2 + 5*rhs_sampling_freq)])
]
rangey_1_wake = [
    min(crosscorr_1Hz_wake[int(len(crosscorr_1Hz_wake)//2 - 2*rhs_sampling_freq):
                            int(len(crosscorr_1Hz_wake)//2 + 2*rhs_sampling_freq)]),
    max(crosscorr_1Hz_wake[int(len(crosscorr_1Hz_wake)//2 - 2*rhs_sampling_freq):
                            int(len(crosscorr_1Hz_wake)//2 + 2*rhs_sampling_freq)])
]
rangey_25_wake = [
    min(crosscorr_25Hz_wake[int(len(crosscorr_25Hz_wake)//2 - 1*rhs_sampling_freq):
                             int(len(crosscorr_25Hz_wake)//2 + 1*rhs_sampling_freq)]),
    max(crosscorr_25Hz_wake[int(len(crosscorr_25Hz_wake)//2 - 1*rhs_sampling_freq):
                             int(len(crosscorr_25Hz_wake)//2 + 1*rhs_sampling_freq)])
]

rangey_04_sleep = [
    min(crosscorr_04Hz_sws[int(len(crosscorr_04Hz_sws)//2 - 5*rhs_sampling_freq):
                            int(len(crosscorr_04Hz_sws)//2 + 5*rhs_sampling_freq)]),
    max(crosscorr_04Hz_sws[int(len(crosscorr_04Hz_sws)//2 - 5*rhs_sampling_freq):
                            int(len(crosscorr_04Hz_sws)//2 + 5*rhs_sampling_freq)])
]
rangey_1_sleep = [
    min(crosscorr_1Hz_sws[int(len(crosscorr_1Hz_sws)//2 - 2*rhs_sampling_freq):
                           int(len(crosscorr_1Hz_sws)//2 + 2*rhs_sampling_freq)]),
    max(crosscorr_1Hz_sws[int(len(crosscorr_1Hz_sws)//2 - 2*rhs_sampling_freq):
                           int(len(crosscorr_1Hz_sws)//2 + 2*rhs_sampling_freq)])
]
rangey_25_sleep = [
    min(crosscorr_25Hz_sws[int(len(crosscorr_25Hz_sws)//2 - 1*rhs_sampling_freq):
                            int(len(crosscorr_25Hz_sws)//2 + 1*rhs_sampling_freq)]),
    max(crosscorr_25Hz_sws[int(len(crosscorr_25Hz_sws)//2 - 1*rhs_sampling_freq):
                            int(len(crosscorr_25Hz_sws)//2 + 1*rhs_sampling_freq)])
]

def plot_without_fill(ax, x, y, xlabel, ylabel, title, show_xlabel=False, x_major_locator_base=None):
    ax.plot(x, y, color='black')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Configure x-axis ticks if specified
    if x_major_locator_base is not None:
        ax.xaxis.set_major_locator(MultipleLocator(base=x_major_locator_base))
    
    # Always display x-axis ticks
    ax.tick_params(axis='x', width=2)  # Ensure x-axis ticks are visible
    
    # Configure x-axis labels and xlabel conditionally
    if show_xlabel:
        ax.tick_params(axis='x', labelsize=25)  # Display x-axis tick labels
        ax.set_xlabel(xlabel, fontsize=30)
    else:
        ax.tick_params(axis='x', labelsize=0)  # Hide x-axis tick labels
        ax.set_xlabel('')  # Ensure xlabel is empty
    
    # Always display y-axis labels
    ax.tick_params(axis='y', labelsize=25, width=2)
    ax.set_ylabel(ylabel, fontsize=30)
    ax.set_title(title, fontsize=35)

# Plot cross-correlation data
plt.subplot(8, 3, 1)
plot_without_fill(plt.gca(), lags, crosscorr_04Hz_wake, 'Lag (s)', 'Crosscorr', '0.4Hz \n Activated', show_xlabel=True)
plt.xlim(-5, 5)
plt.ylim(rangey_04_wake[0] * 0.9, rangey_04_wake[1] * 1.1)

plt.subplot(8, 3, 2)
plot_without_fill(plt.gca(), lags, crosscorr_1Hz_wake, 'Lag (s)', '', '1Hz \n Activated', show_xlabel=True)
plt.xlim(-2, 2)
plt.ylim(rangey_1_wake[0] * 0.9, rangey_1_wake[1] * 1.1)

plt.subplot(8, 3, 3)
plot_without_fill(plt.gca(), lags, crosscorr_25Hz_wake, 'Lag (s)', '', '2.5Hz \n Activated', show_xlabel=True)
plt.xlim(-1, 1)
plt.ylim(rangey_25_wake[0] * 0.9, rangey_25_wake[1] * 1.1)

plt.subplot(8, 3, 13)
plot_without_fill(plt.gca(), lags, crosscorr_04Hz_sws, 'Lag', 'Crosscorr', 'SWS', show_xlabel=True)
plt.xlim(-5, 5)
plt.ylim(rangey_04_sleep[0] * 0.9, rangey_04_sleep[1] * 1.1)

plt.subplot(8, 3, 14)
plot_without_fill(plt.gca(), lags, crosscorr_1Hz_sws, 'Lag', '', 'SWS', show_xlabel=True)
plt.xlim(-2, 2)
plt.ylim(rangey_1_sleep[0] * 0.9, rangey_1_sleep[1] * 1.1)

plt.subplot(8, 3, 15)
plot_without_fill(plt.gca(), lags, crosscorr_25Hz_sws, 'Lag', '', 'SWS', show_xlabel=True)
plt.xlim(-1, 1)
plt.ylim(rangey_25_sleep[0] * 0.9, rangey_25_sleep[1] * 1.1)

# Plot spectrum data
ax1 = plt.subplot(8, 3, 4)
plot_without_fill(ax1, freqs, Spectrum_PY_04_wake, 'Frequency (Hz)', 'PY', '', show_xlabel=False, x_major_locator_base=1)
ax1.xaxis.set_major_locator(MultipleLocator(0.5))
ax1.set_xticks(MultipleLocator(0.5).tick_values(freqs.min(), freqs.max()))
ax1.set_ylabel('PY', fontsize=30)

plt.subplot(8, 3, 7, sharex=ax1)
plot_without_fill(plt.gca(), freqs, Spectrum_TC_04_wake, 'Frequency (Hz)', 'TC', '', show_xlabel=False, x_major_locator_base=1)
plt.gca().set_xticks(ax1.get_xticks())

plt.subplot(8, 3, 10, sharex=ax1)
plot_without_fill(plt.gca(), freqs, Spectrum_input_04_wake, 'Frequency (Hz)', 'Input', '', show_xlabel=True, x_major_locator_base=1)

plt.subplot(8, 3, 5, sharex=ax1)
plot_without_fill(plt.gca(), freqs, Spectrum_PY_1_wake, 'Frequency (Hz)', '', '', show_xlabel=False, x_major_locator_base=1)
plt.gca().set_xticks(ax1.get_xticks())

plt.subplot(8, 3, 8, sharex=ax1)
plot_without_fill(plt.gca(), freqs, Spectrum_TC_1_wake, 'Frequency (Hz)', '', '', show_xlabel=False, x_major_locator_base=1)
plt.gca().set_xticks(ax1.get_xticks())

plt.subplot(8, 3, 11, sharex=ax1)
plot_without_fill(plt.gca(), freqs, Spectrum_input_1_wake, 'Frequency (Hz)', '', '', show_xlabel=True, x_major_locator_base=1)

plt.subplot(8, 3, 6, sharex=ax1)
plot_without_fill(plt.gca(), freqs, Spectrum_PY_25_wake, 'Frequency (Hz)', '', '', show_xlabel=False, x_major_locator_base=1)
plt.gca().set_xticks(ax1.get_xticks())

plt.subplot(8, 3, 9, sharex=ax1)
plot_without_fill(plt.gca(), freqs, Spectrum_TC_25_wake, 'Frequency (Hz)', '', '', show_xlabel=False, x_major_locator_base=1)
plt.gca().set_xticks(ax1.get_xticks())

plt.subplot(8, 3, 12, sharex=ax1)
plot_without_fill(plt.gca(), freqs, Spectrum_input_25_wake, 'Frequency (Hz)', '', '', show_xlabel=True, x_major_locator_base=1)

plt.subplot(8, 3, 16, sharex=ax1)
plot_without_fill(plt.gca(), freqs, Spectrum_PY_04_sws, 'Frequency (Hz)', 'PY', '', show_xlabel=False, x_major_locator_base=1)
plt.gca().set_ylabel('PY', fontsize=30)

plt.subplot(8, 3, 19, sharex=ax1)
plot_without_fill(plt.gca(), freqs, Spectrum_TC_04_sws, 'Frequency (Hz)', 'TC', '', show_xlabel=False, x_major_locator_base=1)

plt.subplot(8, 3, 22, sharex=ax1)
plot_without_fill(plt.gca(), freqs, Spectrum_input_04_sws, 'Frequency (Hz)', 'Input', '', show_xlabel=True, x_major_locator_base=1)

plt.subplot(8, 3, 17, sharex=ax1)
plot_without_fill(plt.gca(), freqs, Spectrum_PY_1_sws, 'Frequency (Hz)', '', '', show_xlabel=False, x_major_locator_base=1)

plt.subplot(8, 3, 20, sharex=ax1)
plot_without_fill(plt.gca(), freqs, Spectrum_TC_1_sws, 'Frequency (Hz)', '', '', show_xlabel=False, x_major_locator_base=1)

plt.subplot(8, 3, 23, sharex=ax1)
plot_without_fill(plt.gca(), freqs, Spectrum_input_1_sws, 'Frequency (Hz)', '', '', show_xlabel=True, x_major_locator_base=1)

plt.subplot(8, 3, 18, sharex=ax1)
plot_without_fill(plt.gca(), freqs, Spectrum_PY_25_sws, 'Frequency (Hz)', '', '', show_xlabel=False, x_major_locator_base=1)

plt.subplot(8, 3, 21, sharex=ax1)
plot_without_fill(plt.gca(), freqs, Spectrum_TC_25_sws, 'Frequency (Hz)', '', '', show_xlabel=False, x_major_locator_base=1)

plt.subplot(8, 3, 24, sharex=ax1)
plot_without_fill(plt.gca(), freqs, Spectrum_input_25_sws, 'Frequency (Hz)', '', '', show_xlabel=True, x_major_locator_base=1)

# Set x-limits for the last row
plt.xlim(0, 4)

# Save and show the figure
plt.savefig('Figure18.png', dpi=300, bbox_inches='tight')
plt.show()


#Figure 19
PY_v_close=loadtxt('../../Data/stimfreq_04_sws/PY_v_close.txt')
PY_v_away=loadtxt('../../Data/stimfreq_04_sws/PY_v_away.txt')
time=loadtxt('../../Data/stimfreq_04_sws/time.txt')
time_s=time/1000

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

plt.figure(figsize=(12,12))
#
ax1 = plt.subplot(211)
ax1.plot(time_s, PY_v_away, color='tab:blue')
configure_axis(ax1, 'Away from stimulus', 'mV', 50, 0, 5)
#
ax2 = plt.subplot(212)
ax2.plot(time_s, PY_v_close, color='tab:blue')
ax2.set_xlabel('Time (ms)', size=30, labelpad=25)
configure_axis(ax2, 'Close to stimulus', 'mV', 50, 0, 5)
plt.tight_layout()
plt.savefig('Figure19.png', dpi=300, bbox_inches='tight')
plt.show()

