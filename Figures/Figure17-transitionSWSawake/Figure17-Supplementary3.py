# -*- coding: utf-8 -*-

from brian2 import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def normalize(array, min_val, max_val):
    return (array - min_val) / (max_val - min_val)

#Make sure you are working in the proper depository
#Data v2
VPY_A1_v2=loadtxt('FigureA1v2/PY_v.txt')
VPY_A2_v2=loadtxt('FigureA2v2/PY_v.txt')
VPY_A3_v2=loadtxt('FigureA3v2/PY_v.txt')
VPY_A4_v2=loadtxt('FigureA4v2/PY_v.txt')
time=loadtxt('FigureA-1/time.txt')
time_s = time/1000
#See the file spectral_powers
freq1_v2=[2.89E-05,2.55E-05,1.84E-05,1.63E-05,1.34E-05,6.91E-06,3.44E-06]
freq2_v2=[5.79E-06,6.01E-06,9.06E-06,7.70E-06,1.06E-05,1.23E-05,7.34E-06]
freq3_v2=[8.13E-06,9.75E-06,1.07E-05,7.95E-06,1.11E-05,1.17E-05,6.86E-06]
freq4_v2=[9.27E-06,8.02E-06,1.34E-05,1.03E-05,1.09E-05,9.66E-06,5.21E-06]
freq1_v2_array = np.array(freq1_v2)
freq2_v2_array = np.array(freq2_v2)
freq3_v2_array = np.array(freq3_v2)
freq4_v2_array = np.array(freq4_v2)
combined_freqs_v2 = np.concatenate([freq1_v2_array, freq2_v2_array, freq3_v2_array, freq4_v2_array])
min_freq_v2 = np.min(combined_freqs_v2)
max_freq_v2 = np.max(combined_freqs_v2)
normalized_v2_freq1 = normalize(freq1_v2_array, min_freq_v2, max_freq_v2)
normalized_v2_freq2 = normalize(freq2_v2_array, min_freq_v2, max_freq_v2)
normalized_v2_freq3 = normalize(freq3_v2_array, min_freq_v2, max_freq_v2)
normalized_v2_freq4 = normalize(freq4_v2_array, min_freq_v2, max_freq_v2)

#Data v2
VPY_A1_v3=loadtxt('FigureA1v3/PY_v.txt')
VPY_A2_v3=loadtxt('FigureA2v3/PY_v.txt')
VPY_A3_v3=loadtxt('FigureA3v3/PY_v.txt')
VPY_A4_v3=loadtxt('FigureA4v3/PY_v.txt')
time=loadtxt('FigureA-1/time.txt')
time_s = time/1000
#See the file spectral_powers
freq1_v3=[3.20E-05,1.99E-05,1.26E-05,4.24E-06,4.27E-06,4.24E-06,2.30E-06]
freq2_v3=[7.44E-06,1.02E-05,1.29E-05,1.71E-05,1.22E-05,1.80E-05,1.65E-05]
freq3_v3=[8.38E-06,1.29E-05,1.41E-05,1.81E-05,1.66E-05,1.55E-05,1.44E-05]
freq4_v3=[1.10E-05,1.47E-05,1.62E-05,1.11E-05,8.69E-06,8.47E-06,9.40E-06]
freq1_v3_array = np.array(freq1_v3)
freq2_v3_array = np.array(freq2_v3)
freq3_v3_array = np.array(freq3_v3)
freq4_v3_array = np.array(freq4_v3)
combined_freqs_v3 = np.concatenate([freq1_v3_array, freq2_v3_array, freq3_v3_array, freq4_v3_array])
min_freq_v3 = np.min(combined_freqs_v3)
max_freq_v3 = np.max(combined_freqs_v3)
normalized_v3_freq1 = normalize(freq1_v3_array, min_freq_v3, max_freq_v3)
normalized_v3_freq2 = normalize(freq2_v3_array, min_freq_v3, max_freq_v3)
normalized_v3_freq3 = normalize(freq3_v3_array, min_freq_v3, max_freq_v3)
normalized_v3_freq4 = normalize(freq4_v3_array, min_freq_v3, max_freq_v3)



#Figure17Top
fig, ax = plt.subplots(4, 1, sharex=True, figsize=(15, 20))

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

ax[0].plot(time_s, VPY_A1_v2, color="tab:blue")
configure_axis(ax[0], 'PY, slow wave sleep', 'mV', 25, 0, 30)
#
ax[1].plot(time_s, VPY_A2_v2, color="tab:blue")
configure_axis(ax[1], 'PY', 'mV', 25, 0, 30, show_legend=True)
#
ax[2].plot(time_s, VPY_A3_v2, color="tab:blue")
configure_axis(ax[2], 'PY', 'mV', 25, 0, 30, show_legend=True)
#
ax[3].plot(time_s, VPY_A4_v2, color="tab:blue")
configure_axis(ax[3], 'PY, activated', 'mV', 25, 0, 30, show_legend=True)
ax[3].set_xlabel('Time (ms)', size=30, labelpad=25)
#
#plt.savefig('Figure17Top.png', dpi=300, bbox_inches='tight')
plt.show()

#Supplementary3Top
fig, ax = plt.subplots(4, 1, sharex=True, figsize=(15, 20))
#
ax[0].plot(time_s, VPY_A1_v3, color="tab:blue")
configure_axis(ax[0], 'PY, slow wave sleep', 'mV', 25, 0, 30)
#
ax[1].plot(time_s, VPY_A2_v3, color="tab:blue")
configure_axis(ax[1], 'PY', 'mV', 25, 0, 30, show_legend=True)
#
ax[2].plot(time_s, VPY_A3_v3, color="tab:blue")
configure_axis(ax[2], 'PY', 'mV', 25, 0, 30, show_legend=True)
#
ax[3].plot(time_s, VPY_A4_v3, color="tab:blue")
configure_axis(ax[3], 'PY, activated', 'mV', 25, 0, 30, show_legend=True)
ax[3].set_xlabel('Time (ms)', size=30, labelpad=25)
#
plt.savefig('Supplementary3Top.png', dpi=300, bbox_inches='tight')
plt.show()

#Figure17Bottom
x_labels = ['SWS', '', '', '', '', '', 'Activated']  # X-axis labels

def plot_spectral_data(ax, data, label, color):
    ax.plot(data, color=color, label=label)

plt.figure(figsize=(8, 8))
ax = plt.gca()
plot_spectral_data(ax, normalized_v2_freq1, '0-2Hz', '#2A52BE')
plot_spectral_data(ax, normalized_v2_freq2, '10-20Hz', '#4B9CD3')
plot_spectral_data(ax, normalized_v2_freq3, '20-30Hz', '#009E60')
plot_spectral_data(ax, normalized_v2_freq4, '30-40Hz', 'green')

# Customize x-ticks and labels
ax.set_xticks(range(len(x_labels)))
ax.set_xticklabels(x_labels)
ax.set_xlabel('Condition', size=25, labelpad=20)
ax.set_ylabel('Power', size=25, labelpad=20)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)
ax.tick_params(axis='both', which='major', labelsize=20, width=2)

# Adding legend
ax.legend(fontsize=20)

plt.tight_layout()
plt.savefig('Figure17Bottom.png', dpi=300, bbox_inches='tight')
plt.show()

#Supplementary3Bottom
plt.figure(figsize=(8, 8))
ax = plt.gca()
plot_spectral_data(ax, normalized_v3_freq1, '0-2Hz', '#2A52BE')
plot_spectral_data(ax, normalized_v3_freq2, '10-20Hz', '#4B9CD3')
plot_spectral_data(ax, normalized_v3_freq3, '20-30Hz', '#009E60')
plot_spectral_data(ax, normalized_v3_freq4, '30-40Hz', 'green')

# Customize x-ticks and labels
ax.set_xticks(range(len(x_labels)))
ax.set_xticklabels(x_labels)
ax.set_xlabel('Condition', size=25, labelpad=20)
ax.set_ylabel('Power', size=25, labelpad=20)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)
ax.tick_params(axis='both', which='major', labelsize=20, width=2)

# Adding legend
ax.legend(fontsize=20)

plt.tight_layout()
plt.savefig('Supplementary3Bottom.png', dpi=300, bbox_inches='tight')
plt.show()