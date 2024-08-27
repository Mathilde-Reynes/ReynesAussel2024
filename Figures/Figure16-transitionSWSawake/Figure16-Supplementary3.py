# -*- coding: utf-8 -*-

from brian2 import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def normalize(array, min_val, max_val):
    return (array - min_val) / (max_val - min_val)

#Make sure you are working in the proper depository
#Data v2
VPY_A1_v2=loadtxt('../../Data/PY_v_c_A.txt')
VPY_A2_v2=loadtxt('../../Data/PY_v_c_C.txt')
VPY_A3_v2=loadtxt('../../Data/PY_v_c_E.txt')
VPY_A4_v2=loadtxt('../../Data/PY_v_c_G.txt')
time=np.arange(0, 30000, 0.02)
time_s = time/1000
#See the file spectral_powers
freq1_v2=[4.32E-05,3.77E-05,3.18E-05,2.49E-05,1.82E-05,5.03E-06,4.02E-06]
freq2_v2=[7.77E-06,7.57E-06,5.98E-061,8.55E-06,1.06E-05,1.91E-05,1.24E-05]
freq3_v2=[9.25E-06,1.12E-05,9.28E-06,8.32E-06,1.19E-05,1.40E-05,8.25E-06]
freq4_v2=[1.26E-05,8.55E-06,8.19E-06,1.21E-05,1.04E-05,9.48E-06,5.93E-06]
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
VPY_A1_v3=loadtxt('../../Data/PY_v_c_A_supp3.txt')
VPY_A2_v3=loadtxt('../../Data/PY_v_c_C_supp3.txt')
VPY_A3_v3=loadtxt('../../Data/PY_v_c_G_supp3.txt')
VPY_A4_v3=loadtxt('../../Data/PY_v_c_E_supp3.txt')
#See the file spectral_powers
freq1_v3=[4.32E-05,3.47E-05,2.10E-05,5.18E-06,3.60E-06,1.68E-06,2.84E-06]
freq2_v3=[7.77E-06,7.65E-06,9.63E-06,1.46E-05,1.81E-05,1.86E-05,2.02E-05]
freq3_v3=[9.25E-06,1.09E-05,1.46E-05,2.09E-05,1.85E-05,1.94E-05,1.56E-05]
freq4_v3=[1.26E-05,1.21E-05,1.54E-05,1.44E-05,1.52E-05,1.41E-05,9.94E-06]
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



#Figure16Top
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
configure_axis(ax[1], 'PY', 'mV', 25, 0, 30)
#
ax[2].plot(time_s, VPY_A3_v2, color="tab:blue")
configure_axis(ax[2], 'PY', 'mV', 25, 0, 30)
#
ax[3].plot(time_s, VPY_A4_v2, color="tab:blue")
configure_axis(ax[3], 'PY, activated', 'mV', 25, 0, 30)
ax[3].set_xlabel('Time (ms)', size=30, labelpad=25)
#
#plt.savefig('Figure16Top.png', dpi=300, bbox_inches='tight')
plt.show()

#Supplementary3Top
fig, ax = plt.subplots(4, 1, sharex=True, figsize=(15, 20))
#
ax[0].plot(time_s, VPY_A1_v3, color="tab:blue")
configure_axis(ax[0], 'PY, slow wave sleep', 'mV', 25, 0, 30)
#
ax[1].plot(time_s, VPY_A2_v3, color="tab:blue")
configure_axis(ax[1], 'PY', 'mV', 25, 0, 30)
#
ax[2].plot(time_s, VPY_A3_v3, color="tab:blue")
configure_axis(ax[2], 'PY', 'mV', 25, 0, 30)
#
ax[3].plot(time_s, VPY_A4_v3, color="tab:blue")
configure_axis(ax[3], 'PY, activated', 'mV', 25, 0, 30)
ax[3].set_xlabel('Time (ms)', size=30, labelpad=25)
#
#plt.savefig('Supplementary3Top.png', dpi=300, bbox_inches='tight')
plt.show()

#Figure16Bottom
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
plt.savefig('Figure16Bottom.png', dpi=300, bbox_inches='tight')
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