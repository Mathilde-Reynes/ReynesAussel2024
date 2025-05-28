#!/usr/bin/env python
# coding: utf-8

from brian2 import *
import numpy as np
import matplotlib.pyplot as plt
from brian2.units.constants import *
import matplotlib.gridspec as gridspec
from scipy import signal


def figure_plotting(fig_number, all_monitors, all_monitors_T, all_synapses, monitor_poisson, runtime, N):
    V1_PYd, V2_PYs, V3_INd, V4_INs, R2_PYs, R4_INs, I1_PYd, I2_INd, S1, S2, M0, M1 = all_monitors
    V1_RE, V2_TC, R1_RE, R2_TC, I1_RE, I2_TC = all_monitors_T
    S_AMPA_PY_PY, S_AMPA_PY_IN, S_NMDA_PY_PY, S_NMDA_PY_IN, S_GABAA_IN_PY = all_synapses

    # Figure 5
    if fig_number == "5":
        fig = plt.figure(figsize=(20, 12))
        gs = gridspec.GridSpec(2, 2, width_ratios=[20, 0.5], height_ratios=[2, 1], wspace=0.05, hspace=0.3)
        ax1 = plt.subplot(gs[0, 0])
        im1 = ax1.imshow(
            V2_PYs.v / mV,
            aspect="auto",
            cmap="Greys",
            vmax=-60,
            vmin=-75,
            extent=[0, runtime, N - 1, 0],
            interpolation="bicubic",
        )
        ax1.set_title("PY", size=30, loc="left")
        ax1.set_ylabel("Neuron index", size=30, labelpad=30)
        ax1.set_xlim(0, 30)
        ax1.yaxis.set_major_locator(MultipleLocator(base=25))
        ax1.tick_params(axis="both", which="major", labelsize=25, width=2)
        ax2 = plt.subplot(gs[1, 0])
        im2 = ax2.imshow(
            V2_TC.v / mV,
            aspect="auto",
            cmap="Greys",
            vmax=-60,
            vmin=-75,
            extent=[0, runtime, 50, 0],
            interpolation="bicubic",
        )
        ax2.set_title("TC", size=30, loc="left")
        ax2.set_xlabel("Time (s)", size=30, labelpad=10)
        ax2.set_ylabel("Neuron index", size=30, labelpad=30)
        ax2.set_xlim(0, 30)
        ax2.yaxis.set_major_locator(MultipleLocator(base=25))
        ax2.tick_params(axis="both", which="major", labelsize=25, width=2)

        cbar_ax = plt.subplot(gs[:, 1])
        cbar = fig.colorbar(im1, cax=cbar_ax, orientation="vertical")
        cbar.set_label("Membrane potential (mV)", size=30, labelpad=30)
        cbar.ax.tick_params(labelsize=25, width=2)
        plt.show()

    # Figure 7
    elif fig_number == "7":
        fig, ax = subplots(4, 1, sharex=True, figsize=(19, 18))
        ax[0].plot(V2_PYs.t, V2_PYs.v[N // 2] / mV, color="tab:blue", linewidth=1.5)
        ax[0].set_title("PY", size=30, loc="left")
        ax[0].set_ylabel("mV", size=30, labelpad=30)
        ax[0].set_xlim([15, 30])
        ax[0].spines["top"].set_visible(False)
        ax[0].spines["right"].set_visible(False)
        ax[0].spines["bottom"].set_visible(True)
        ax[0].spines["left"].set_visible(True)
        ax[0].yaxis.set_major_locator(MultipleLocator(base=25))
        ax[0].tick_params(axis="both", which="major", labelsize=25, width=2)
        #
        ax[1].plot(V4_INs.t, V4_INs.v[N // 8] / mV, color="tab:green", linewidth=1.5)
        ax[1].set_title("IN", size=30, loc="left")
        ax[1].set_ylabel("mV", size=30, labelpad=30)
        ax[1].spines["top"].set_visible(False)
        ax[1].spines["right"].set_visible(False)
        ax[1].spines["bottom"].set_visible(True)
        ax[1].spines["left"].set_visible(True)
        ax[1].yaxis.set_major_locator(MultipleLocator(base=25))
        ax[1].tick_params(axis="both", which="major", labelsize=25, width=2)
        #
        ax[2].plot(V1_RE.t, V1_RE.v[N // 4] / mV, color="tab:red", linewidth=1.5)
        ax[2].set_title("RE", size=30, loc="left")
        ax[2].set_ylabel("mV", size=30, labelpad=30)
        ax[2].spines["top"].set_visible(False)
        ax[2].spines["right"].set_visible(False)
        ax[2].spines["bottom"].set_visible(True)
        ax[2].spines["left"].set_visible(True)
        ax[2].yaxis.set_major_locator(MultipleLocator(base=25))
        ax[2].tick_params(axis="both", which="major", labelsize=25, width=2)
        #
        ax[3].plot(V2_TC.t, V2_TC.v[N // 4] / mV, color="tab:orange", linewidth=1.5)
        ax[3].set_title("TC", size=30, loc="left")
        ax[3].set_ylabel("mV", size=30, labelpad=30)
        ax[3].set_xlabel("second", size=30)
        ax[3].tick_params(axis="both", labelbottom=True)
        ax[3].spines["top"].set_visible(False)
        ax[3].spines["right"].set_visible(False)
        ax[3].spines["bottom"].set_visible(True)
        ax[3].spines["left"].set_visible(True)
        ax[3].yaxis.set_major_locator(MultipleLocator(base=25))
        ax[3].tick_params(axis="both", which="major", labelsize=25, width=2)
        plt.show()

    # Figure 8
    elif fig_number == "8":
        n_start = 25
        fig, ax = subplots(7, 1, sharex=True, figsize=(19, 34))
        ax[0].plot(V2_PYs.t / ms, V2_PYs.v[n_start - 4] / mV, linewidth=0.4, color="tab:blue", alpha=0.6)
        ax[0].plot(V2_PYs.t / ms, V2_PYs.v[n_start - 3] / mV, linewidth=0.4, color="tab:blue", alpha=0.6)
        ax[0].plot(V2_PYs.t / ms, V2_PYs.v[n_start - 2] / mV, linewidth=0.4, color="tab:blue", alpha=0.6)
        ax[0].plot(V2_PYs.t / ms, V2_PYs.v[n_start - 1] / mV, linewidth=0.4, color="tab:blue", alpha=0.6)
        ax[0].plot(V2_PYs.t / ms, V2_PYs.v[n_start] / mV, linewidth=2, color="tab:blue")
        ax[0].plot(V2_PYs.t / ms, V2_PYs.v[n_start + 1] / mV, linewidth=0.4, color="tab:blue", alpha=0.6)
        ax[0].plot(V2_PYs.t / ms, V2_PYs.v[n_start + 2] / mV, linewidth=0.4, color="tab:blue", alpha=0.6)
        ax[0].plot(V2_PYs.t / ms, V2_PYs.v[n_start + 3] / mV, linewidth=0.4, color="tab:blue", alpha=0.6)
        ax[0].plot(V2_PYs.t / ms, V2_PYs.v[n_start + 4] / mV, linewidth=0.4, color="tab:blue", alpha=0.6)
        ax[0].spines["top"].set_visible(False)
        ax[0].spines["right"].set_visible(False)
        ax[0].spines["bottom"].set_visible(True)
        ax[0].spines["left"].set_visible(True)
        ax[0].set_xlim([3200, 3650])
        ax[0].yaxis.set_major_locator(MultipleLocator(base=25))
        ax[0].tick_params(axis="both", which="major", labelsize=25, width=2)
        ax[0].set_title("Pyramidal cells", size=30, loc="left")
        ax[0].set_ylabel(r"$\mathrm{mV}$", size=30, labelpad=37)
        #
        ax[1].plot(
            I1_PYd.t / ms,
            I1_PYd.IEPSPs_PY_PY[n_start] / (0.001 * amp * meter**-2),
            linewidth=1.5,
            color="black",
            label="EPSPs from PY",
        )
        ax[1].plot(
            I1_PYd.t / ms,
            I1_PYd.IEPSPs_IN_PY[n_start] / (0.001 * amp * meter**-2),
            linewidth=1.5,
            color="gray",
            label="IPSPs from IN",
        )
        ax[1].spines["top"].set_visible(False)
        ax[1].spines["right"].set_visible(False)
        ax[1].spines["bottom"].set_visible(True)
        ax[1].spines["left"].set_visible(True)
        ax[1].tick_params(axis="both", which="major", labelsize=25, width=2)
        ax[1].set_title("miniPSPs", size=30, loc="left")
        ax[1].set_ylabel(r"$1 \, \mathrm{mA} \cdot \mathrm{m}^{-2}$", size=30, labelpad=25)
        ax[1].legend(fontsize=25, loc="lower right")
        #
        ax[2].plot(
            I1_PYd.t / ms,
            I1_PYd.I_nap[n_start] / (0.001 * amp * meter**-2),
            color="black",
            linewidth=1.5,
            label="10*I_Na(p)",
        )
        ax[2].plot(
            I1_PYd.t / ms,
            I1_PYd.I_na[n_start] / (0.01 * amp * meter**-2),
            color="gray",
            linewidth=1.5,
            label="I_Na",
        )
        ax[2].spines["top"].set_visible(False)
        ax[2].spines["right"].set_visible(False)
        ax[2].spines["bottom"].set_visible(True)
        ax[2].spines["left"].set_visible(True)
        ax[2].tick_params(axis="both", which="major", labelsize=25, width=2)
        ax[2].set_title("Sodium currents", size=30, loc="left")
        ax[2].set_ylabel(r"$10 \, \mathrm{mA} \cdot \mathrm{m}^{-2}$", size=30, labelpad=30)
        ax[2].legend(fontsize=25, loc="lower right")
        #
        ax[3].plot(
            I1_PYd.t / ms,
            I1_PYd.IsynAMPA_PY_PY[n_start] / (0.01 * amp * meter**-2),
            color="black",
            linewidth=1.5,
            label="I_AMPAs",
        )
        ax[3].plot(
            I1_PYd.t / ms,
            I1_PYd.IsynNMDA_PY_PY[n_start] / (0.001 * amp * meter**-2),
            color="gray",
            linewidth=1.5,
            label="10*I_NMDAs",
        )
        ax[3].spines["top"].set_visible(False)
        ax[3].spines["right"].set_visible(False)
        ax[3].spines["bottom"].set_visible(True)
        ax[3].spines["left"].set_visible(True)
        ax[3].tick_params(axis="both", which="major", labelsize=25, width=2)
        ax[3].set_title("Synaptic currents PY-PY", size=30, loc="left")
        ax[3].set_ylabel(r"$10 \, \mathrm{mA} \cdot \mathrm{m}^{-2}$", size=30, labelpad=44)
        ax[3].legend(fontsize=25, loc="lower right")
        #
        ax[4].plot(
            I1_PYd.t / ms, I1_PYd.I_kca[n_start] / (0.001 * amp * meter**-2), color="black", linewidth=1.5
        )
        ax[4].spines["top"].set_visible(False)
        ax[4].spines["right"].set_visible(False)
        ax[4].spines["bottom"].set_visible(True)
        ax[4].spines["left"].set_visible(True)
        ax[4].tick_params(axis="both", which="major", labelsize=25, width=2)
        ax[4].set_title("Calcium-dependent potassium current", size=30, loc="left")
        ax[4].set_ylabel(r"$1 \, \mathrm{mA} \cdot \mathrm{m}^{-2}$", size=30, labelpad=42)
        #
        ax[5].plot(
            I1_PYd.t / ms,
            I1_PYd.IsynGABAA_IN_PY[n_start] / (0.01 * amp * meter**-2),
            color="black",
            linewidth=1.5,
        )
        ax[5].spines["top"].set_visible(False)
        ax[5].spines["right"].set_visible(False)
        ax[5].spines["bottom"].set_visible(True)
        ax[5].spines["left"].set_visible(True)
        ax[5].tick_params(axis="both", which="major", labelsize=25, width=2)
        ax[5].set_title("Synaptic GABAA current IN-PY", size=30, loc="left")
        ax[5].set_ylabel(r"$10 \, \mathrm{mA} \cdot \mathrm{m}^{-2}$", size=30, labelpad=65)
        ax[5].set_ylim([0, 4])
        #
        ax[6].plot(
            S1.t / ms, S1.D[(S_AMPA_PY_PY.j[:] == n_start).nonzero()[0][0]], color="black", linewidth=1.5
        )  # see the first synapse that target the neuron of interest
        ax[6].plot(
            S1.t / ms,
            S1.D[(S_AMPA_PY_PY.j[:] == n_start).nonzero()[0][1]],
            color="black",
            linewidth=0.4,
            alpha=0.6,
        )
        ax[6].plot(
            S1.t / ms,
            S1.D[(S_AMPA_PY_PY.j[:] == n_start).nonzero()[0][2]],
            color="black",
            linewidth=0.4,
            alpha=0.6,
        )
        ax[6].plot(
            S1.t / ms,
            S1.D[(S_AMPA_PY_PY.j[:] == n_start).nonzero()[0][3]],
            color="black",
            linewidth=0.4,
            alpha=0.6,
        )
        ax[6].spines["top"].set_visible(False)
        ax[6].spines["right"].set_visible(False)
        ax[6].spines["bottom"].set_visible(True)
        ax[6].spines["left"].set_visible(True)
        ax[6].tick_params(axis="both", which="major", labelsize=25, width=2)
        ax[6].set_title("AMPA synaptic depression from input synapses", size=30, loc="left")
        ax[6].set_xlabel("ms", size=30, labelpad=30)
        plt.show()

    # Figure 9&11
    elif fig_number in ["9-A1", "9-A2", "9-A3", "9-B1", "9-B2", "9-B3", "11-1", "11-2"]:
        fig, ax = plt.subplots(figsize=(15, 5))
        ax.plot(V2_PYs.t, V2_PYs.v[N // 2] / mV, color="tab:blue", linewidth=1.5)
        ax.set_title(f"PY N={N}", size=35, loc="left")
        ax.set_ylabel("mV", size=30, labelpad=25)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(True)
        ax.spines["left"].set_visible(True)
        ax.yaxis.set_major_locator(MultipleLocator(base=25))
        ax.set_xlim([0, 8])
        ax.tick_params(axis="both", which="major", labelsize=30, width=2)
        plt.show()

    # Figure 10
    elif fig_number in ["10-1", "10-2", "10-3"]:
        fig, ax = plt.subplots(figsize=(15, 5))
        ax.plot(V2_PYs.t, V2_PYs.v[N // 2] / mV, color="tab:blue", linewidth=1.5)
        ax.set_title("PY", size=35, loc="left")
        ax.set_ylabel("mV", size=30, labelpad=25)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(True)
        ax.spines["left"].set_visible(True)
        ax.yaxis.set_major_locator(MultipleLocator(base=25))
        ax.set_xlim([0, 20])
        ax.tick_params(axis="both", which="major", labelsize=30, width=2)
        plt.show()

    # Figure 12
    elif fig_number in ["12-A1", "12-A2"]:
        fig = plt.figure(figsize=(15, 12))
        gs = gridspec.GridSpec(2, 2, width_ratios=[20, 0.5], height_ratios=[1, 1], wspace=0.05, hspace=0.3)
        ax1 = plt.subplot(gs[0, 0])
        im1 = ax1.imshow(
            V2_PYs.v / mV,
            aspect="auto",
            cmap="Greys",
            vmax=-60,
            vmin=-75,
            extent=[0, runtime, N - 1, 0],
            interpolation="bicubic",
        )
        ax1.set_title(r"Conductance PY-IN = 0.02$\mu$S", size=30, loc="left")
        ax1.set_ylabel("Neuron index", size=30, labelpad=30)
        ax1.set_xlabel("Time (s)", size=30, labelpad=30)
        ax1.set_xlim(0.1, 5.1)
        ax1.yaxis.set_major_locator(MultipleLocator(base=25))
        ax1.tick_params(axis="both", which="major", labelsize=25, width=2)
        cbar_ax = plt.subplot(gs[:, 1])
        cbar = fig.colorbar(im1, cax=cbar_ax, orientation="vertical")
        cbar.set_label("Membrane potential (mV)", size=30, labelpad=30)
        cbar.ax.tick_params(labelsize=25, width=2)
        plt.show()

    # Figure 13
    elif fig_number in ["13-A", "13-B"]:
        fig1, ax1 = plt.subplots(figsize=(19, 9))
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)
        ax1.spines["bottom"].set_visible(True)
        ax1.spines["left"].set_visible(True)
        ax1.plot(R2_PYs.t, R2_PYs.i, ".", markersize=2, alpha=0.5, color="tab:blue")
        ax1.set_ylabel("Neuron index", fontsize=30)
        ax1.set_xlabel("Time (s)", size=30)
        ax1.set_xlim([0, 10])
        ax1.tick_params(axis="both", which="major", labelsize=25, width=2)
        ax1.set_title("Raster plot, PY", size=30, loc="left")
        fig1.tight_layout()
        plt.show()

    # Figure 14
    elif fig_number == "14":
        n_choice = 22
        fig, ax = subplots(5, 1, sharex=True, figsize=(19, 29))
        ax[0].plot(V2_PYs.t / ms, V2_PYs.v[n_choice - 4] / mV, linewidth=0.4, color="tab:blue", alpha=0.6)
        ax[0].plot(V2_PYs.t / ms, V2_PYs.v[n_choice - 3] / mV, linewidth=0.4, color="tab:blue", alpha=0.6)
        ax[0].plot(V2_PYs.t / ms, V2_PYs.v[n_choice - 2] / mV, linewidth=0.4, color="tab:blue", alpha=0.6)
        ax[0].plot(V2_PYs.t / ms, V2_PYs.v[n_choice - 1] / mV, linewidth=0.4, color="tab:blue", alpha=0.6)
        ax[0].plot(V2_PYs.t / ms, V2_PYs.v[n_choice] / mV, linewidth=2, color="tab:blue")
        ax[0].plot(V2_PYs.t / ms, V2_PYs.v[n_choice + 1] / mV, linewidth=0.4, color="tab:blue", alpha=0.6)
        ax[0].plot(V2_PYs.t / ms, V2_PYs.v[n_choice + 2] / mV, linewidth=0.4, color="tab:blue", alpha=0.6)
        ax[0].plot(V2_PYs.t / ms, V2_PYs.v[n_choice + 3] / mV, linewidth=0.4, color="tab:blue", alpha=0.6)
        ax[0].plot(V2_PYs.t / ms, V2_PYs.v[n_choice + 4] / mV, linewidth=0.4, color="tab:blue", alpha=0.6)
        ax[0].spines["top"].set_visible(False)
        ax[0].spines["right"].set_visible(False)
        ax[0].spines["bottom"].set_visible(True)
        ax[0].spines["left"].set_visible(True)
        ax[0].set_xlim([0, 10000])
        ax[0].yaxis.set_major_locator(MultipleLocator(base=25))
        ax[0].tick_params(axis="both", which="major", labelsize=25, width=2)
        ax[0].set_title("A. Pyramidal cell", size=30, loc="left")
        ax[0].set_ylabel(r"$\mathrm{mV}$", size=30, labelpad=37)
        #
        ax[1].plot(V2_TC.t / ms, V2_TC.v[n_choice - 4] / mV, linewidth=0.4, color="tab:orange", alpha=0.6)
        ax[1].plot(V2_TC.t / ms, V2_TC.v[n_choice - 3] / mV, linewidth=0.4, color="tab:orange", alpha=0.6)
        ax[1].plot(V2_TC.t / ms, V2_TC.v[n_choice - 2] / mV, linewidth=0.4, color="tab:orange", alpha=0.6)
        ax[1].plot(V2_TC.t / ms, V2_TC.v[n_choice - 1] / mV, linewidth=0.4, color="tab:orange", alpha=0.6)
        ax[1].plot(V2_TC.t / ms, V2_TC.v[n_choice] / mV, linewidth=2, color="tab:orange")
        ax[1].plot(V2_TC.t / ms, V2_TC.v[n_choice + 1] / mV, linewidth=0.4, color="tab:orange", alpha=0.6)
        ax[1].plot(V2_TC.t / ms, V2_TC.v[n_choice + 2] / mV, linewidth=0.4, color="tab:orange", alpha=0.6)
        ax[1].plot(V2_TC.t / ms, V2_TC.v[n_choice + 3] / mV, linewidth=0.4, color="tab:orange", alpha=0.6)
        ax[1].plot(V2_TC.t / ms, V2_TC.v[n_choice + 4] / mV, linewidth=0.4, color="tab:orange", alpha=0.6)
        ax[1].spines["top"].set_visible(False)
        ax[1].spines["right"].set_visible(False)
        ax[1].spines["bottom"].set_visible(True)
        ax[1].spines["left"].set_visible(True)
        ax[1].set_xlim([0, 10000])
        ax[1].yaxis.set_major_locator(MultipleLocator(base=25))
        ax[1].tick_params(axis="both", which="major", labelsize=25, width=2)
        ax[1].set_title("B. Thalamic relay cell", size=30, loc="left")
        ax[1].set_ylabel(r"$\mathrm{mV}$", size=30, labelpad=37)
        #
        ax[2].plot(V1_RE.t / ms, V1_RE.v[n_choice - 4] / mV, linewidth=0.4, color="tab:red", alpha=0.6)
        ax[2].plot(V1_RE.t / ms, V1_RE.v[n_choice - 3] / mV, linewidth=0.4, color="tab:red", alpha=0.6)
        ax[2].plot(V1_RE.t / ms, V1_RE.v[n_choice - 2] / mV, linewidth=0.4, color="tab:red", alpha=0.6)
        ax[2].plot(V1_RE.t / ms, V1_RE.v[n_choice - 1] / mV, linewidth=0.4, color="tab:red", alpha=0.6)
        ax[2].plot(V1_RE.t / ms, V1_RE.v[n_choice] / mV, linewidth=2, color="tab:red")
        ax[2].plot(V1_RE.t / ms, V1_RE.v[n_choice + 1] / mV, linewidth=0.4, color="tab:red", alpha=0.6)
        ax[2].plot(V1_RE.t / ms, V1_RE.v[n_choice + 2] / mV, linewidth=0.4, color="tab:red", alpha=0.6)
        ax[2].plot(V1_RE.t / ms, V1_RE.v[n_choice + 3] / mV, linewidth=0.4, color="tab:red", alpha=0.6)
        ax[2].plot(V1_RE.t / ms, V1_RE.v[n_choice + 4] / mV, linewidth=0.4, color="tab:red", alpha=0.6)
        ax[2].spines["top"].set_visible(False)
        ax[2].spines["right"].set_visible(False)
        ax[2].spines["bottom"].set_visible(True)
        ax[2].spines["left"].set_visible(True)
        ax[2].set_xlim([0, 10000])
        ax[2].yaxis.set_major_locator(MultipleLocator(base=25))
        ax[2].tick_params(axis="both", which="major", labelsize=25, width=2)
        ax[2].set_title("C. Reticular cell", size=30, loc="left")
        ax[2].set_ylabel(r"$\mathrm{mV}$", size=30, labelpad=37)
        #
        ax[3].plot(I2_TC.t / ms, I2_TC.I_t[n_choice] / (0.01 * amp * meter**-2), color="black", linewidth=1.5)
        ax[3].spines["top"].set_visible(False)
        ax[3].spines["right"].set_visible(False)
        ax[3].spines["bottom"].set_visible(True)
        ax[3].spines["left"].set_visible(True)
        ax[3].tick_params(axis="both", which="major", labelsize=25, width=2)
        ax[3].set_title("D. Low-threshold calcium current", size=30, loc="left")
        ax[3].set_ylabel(r"$1 \, \mathrm{mA} \cdot \mathrm{m}^{-2}$", size=30, labelpad=25)
        ax[3].legend(fontsize=25, loc="lower right")
        ax[3].set_xlim([0, 10000])
        fig.tight_layout()
        #
        ax[4].plot(
            I1_RE.t / ms,
            I1_RE.IsynAMPA_PY_RE[n_choice] / (0.01 * amp * meter**-2),
            color="black",
            linewidth=1.5,
        )
        ax[4].spines["top"].set_visible(False)
        ax[4].spines["right"].set_visible(False)
        ax[4].spines["bottom"].set_visible(True)
        ax[4].spines["left"].set_visible(True)
        ax[4].tick_params(axis="both", which="major", labelsize=25, width=2)
        ax[4].set_title("E. Synaptic AMPA current from PY to RE", size=30, loc="left")
        ax[4].set_ylabel(r"$1 \, \mathrm{mA} \cdot \mathrm{m}^{-2}$", size=30, labelpad=25)
        ax[4].legend(fontsize=25, loc="lower right")
        ax[4].set_xlabel("ms", size=30, labelpad=30)
        fig.tight_layout()
        plt.show()

    # Figure 15
    elif fig_number in ["15-A", "15-B", "15-C", "15-D"]:
        fig, ax = subplots(3, 1, sharex=True, figsize=(12, 15))
        ax[0].plot(V2_PYs.t, V2_PYs.v[N // 2] / mV, color="tab:blue", linewidth=1.5)
        ax[0].set_title("PY", size=35, loc="left")
        ax[0].set_ylabel("mV", size=30, labelpad=25)
        ax[0].spines["top"].set_visible(False)
        ax[0].spines["right"].set_visible(False)
        ax[0].spines["bottom"].set_visible(True)
        ax[0].spines["left"].set_visible(True)
        ax[0].yaxis.set_major_locator(MultipleLocator(base=25))
        ax[0].tick_params(axis="both", which="major", labelsize=30, width=2)
        ax[1].plot(V1_RE.t, V1_RE.v[N // 4] / mV, color="tab:red", linewidth=1.5)
        ax[1].set_title("RE", size=35, loc="left")
        ax[1].set_ylabel("mV", size=30, labelpad=25)
        ax[1].spines["top"].set_visible(False)
        ax[1].spines["right"].set_visible(False)
        ax[1].spines["bottom"].set_visible(True)
        ax[1].spines["left"].set_visible(True)
        ax[1].yaxis.set_major_locator(MultipleLocator(base=25))
        ax[1].tick_params(axis="both", which="major", labelsize=30, width=2)
        ax[2].plot(V2_TC.t, V2_TC.v[N // 4] / mV, color="tab:orange", linewidth=1.5)
        ax[2].set_title("TC", size=35, loc="left")
        ax[2].set_ylabel("mV", size=30, labelpad=25)
        ax[2].set_xlabel("Time (s)", size=30, labelpad=25)
        ax[2].spines["top"].set_visible(False)
        ax[2].spines["right"].set_visible(False)
        ax[2].spines["bottom"].set_visible(True)
        ax[2].spines["left"].set_visible(True)
        ax[2].yaxis.set_major_locator(MultipleLocator(base=25))
        ax[2].set_xlim([15, 18])
        ax[2].tick_params(axis="both", which="major", labelsize=30, width=2)
        plt.show()

    # Figure 16 & S4
    elif fig_number in [
        "16-A1",
        "16-A2",
        "16-A3",
        "16-A4",
        "16-B1",
        "16-B2",
        "16-B3",
        "16-B4",
        "16-B5",
        "16-B6",
        "16-B7",
        "S4-A1",
        "S4-B1",
        "S4-A2",
        "S4-B3",
        "S4-A3",
        "S4-B5",
        "S4-A4",
        "S4-B7",
        "S4-B2",
        "S4-B4",
        "S4-B6"
    ]:
        fig, ax = plt.subplots(figsize=(15, 5))
        ax.plot(V2_PYs.t, V2_PYs.v[N // 2] / mV, color="tab:blue", linewidth=1.5)
        ax.set_title("PY, N=20", size=35, loc="left")
        ax.set_ylabel("mV", size=30, labelpad=25)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(True)
        ax.spines["left"].set_visible(True)
        ax.yaxis.set_major_locator(MultipleLocator(base=25))
        ax.set_xlim([0, 30])
        ax.tick_params(axis="both", which="major", labelsize=30, width=2)
        plt.show()

    # Figure 17
    elif fig_number in [
        "17-A1",
        "18-A1",
        "17-A2",
        "18-A2",
        "17-A3",
        "18-A3",
        "17-B1",
        "18-B1",
        "17-B2",
        "18-B2",
        "17-B3",
        "18-B3",
    ]:

        OVERLAP = 50 * ms

        def RHS(raster_i, raster_t, neuron_range, time_window=100 * ms, overlap=OVERLAP):
            # first, remove all spikes that are not emitted by the neurons in neuron_range from raster
            range_raster = where(logical_and(raster_i >= neuron_range[0], raster_i <= neuron_range[-1]))[0]
            small_raster_i = raster_i[range_raster]
            small_raster_t = raster_t[range_raster]
            N_neurons = len(neuron_range)
            # then, compute histogram
            rhs = []
            for t in arange(0 * ms, 15000 * ms, overlap):  # simulation time goes from 0ms to 15000ms
                rhs.append(
                    len(where(logical_and(small_raster_t >= t, small_raster_t < t + time_window))[0])
                    / N_neurons
                )
            return rhs

        rhs_PY = RHS(R2_PYs.i, R2_PYs.t, list(range(38, 63)))
        rhs_TC = RHS(R2_TC.i, R2_TC.t, list(range(19, 32)))
        rhs_input = RHS(monitor_poisson.i, monitor_poisson.t, list(range(12)))
        close("all")
        time_rhs = arange(0 * ms, 15000 * ms, OVERLAP) / 1000

        if fig_number in ["17-A1", "17-A2", "17-A3", "17-B1", "17-B2", "17-B3"]:
            plt.figure(figsize=(20, 20))
            plt.subplots_adjust(hspace=0.5)

            def plot_with_fill(ax, x, y, label, color="black", show_xlabel=False, show_xticks=True):
                ax.plot(x, y, color=color)
                ax.fill_between(x, y, color=color, alpha=0.8)
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
                ax.tick_params(axis="x", labelsize=15 if show_xticks else 0)
                ax.tick_params(axis="y", labelsize=25)
                ax.set_xlim([0.005, 0.015])
                if show_xlabel:
                    ax.set_xlabel("Time", fontsize=30)
                else:
                    ax.set_xlabel("")  # Remove xlabel if not showing

            ax1 = plt.subplot(631)
            plot_with_fill(ax1, time_rhs, rhs_PY, "PY", show_xlabel=False, show_xticks=False)
            ax1.set_ylabel("PY", fontsize=30)
            ax1.set_title("", fontsize=35)
            plt.subplot(634, sharex=ax1)
            plot_with_fill(plt.gca(), time_rhs, rhs_TC, "TC", show_xlabel=False, show_xticks=False)
            plt.gca().set_ylabel("TC", fontsize=30)
            plt.subplot(637, sharex=ax1)
            plot_with_fill(plt.gca(), time_rhs, rhs_input, "Input", show_xlabel=True, show_xticks=True)
            plt.gca().set_ylabel("Input", fontsize=30)
            plt.show()

        if fig_number in ["18-A1", "18-A2", "18-A3", "18-B1", "18-B2", "18-B3"]:
            # For figure 18, we're going to need the crosscorrelation of PY and input RSHs
            # and the power spectra of all RSHs
            OVERLAP = 50
            time_rhs = arange(0, 15000, OVERLAP) / 1000
            crosscorr = signal.correlate(rhs_PY, rhs_input)
            lags = signal.correlation_lags(len(rhs_PY), len(rhs_input))
            # print(lags[0])
            rhs_sampling_freq = 1 / (time_rhs[1])  # sampling frequency of the running spike histogram in Hz
            lags = lags / rhs_sampling_freq  # lags in second
            # print(lags[0])
            freqs, Spectrum_PY = signal.periodogram(rhs_PY, rhs_sampling_freq, scaling="spectrum")
            _, Spectrum_TC = signal.periodogram(rhs_TC, rhs_sampling_freq, scaling="spectrum")
            _, Spectrum_input = signal.periodogram(rhs_input, rhs_sampling_freq, scaling="spectrum")
            plt.figure(figsize=(22, 25))
            plt.subplots_adjust(hspace=1.4, wspace=0.3)
            # Compute range values for cross-correlation data

            print(f"rhs_sampling_freq (dimensionless): {rhs_sampling_freq}")
            rangey = [
                min(
                    crosscorr[
                        int(len(crosscorr) // 2 - 5 * rhs_sampling_freq) : int(
                            len(crosscorr) // 2 + 5 * rhs_sampling_freq
                        )
                    ]
                ),
                max(
                    crosscorr[
                        int(len(crosscorr) // 2 - 5 * rhs_sampling_freq) : int(
                            len(crosscorr) // 2 + 5 * rhs_sampling_freq
                        )
                    ]
                ),
            ]

            def plot_without_fill(
                ax, x, y, xlabel, ylabel, title, show_xlabel=False, x_major_locator_base=None
            ):
                ax.plot(x, y, color="black")
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
                # Configure x-axis ticks if specified
                if x_major_locator_base is not None:
                    ax.xaxis.set_major_locator(MultipleLocator(base=x_major_locator_base))
                ax.tick_params(axis="x", width=2)
                if show_xlabel:
                    ax.tick_params(axis="x", labelsize=25)
                    ax.set_xlabel(xlabel, fontsize=30)
                else:
                    ax.tick_params(axis="x", labelsize=0)
                    ax.set_xlabel("")
                ax.tick_params(axis="y", labelsize=25, width=2)
                ax.set_ylabel(ylabel, fontsize=30)
                ax.set_title(title, fontsize=35)

            # Plot cross-correlation data
            plt.figure(figsize=(15, 15))
            # First plot (cross-correlation)
            plt.subplot(4, 1, 1)
            plot_without_fill(plt.gca(), lags, crosscorr, "Lag (s)", "Crosscorr", "", show_xlabel=True)
            plt.xlim(-5, 5)
            plt.ylim(rangey[0] * 0.9, rangey[1] * 1.1)
            # Second, third, and fourth plots (spectra)
            plt.subplot(4, 1, 2)
            plot_without_fill(
                plt.gca(),
                freqs,
                Spectrum_PY,
                "Frequency (Hz)",
                "PY",
                "",
                show_xlabel=False,
                x_major_locator_base=1,
            )
            plt.ylabel("PY", fontsize=30)
            plt.subplot(4, 1, 3, sharex=plt.gca())
            plot_without_fill(
                plt.gca(),
                freqs,
                Spectrum_TC,
                "Frequency (Hz)",
                "TC",
                "",
                show_xlabel=False,
                x_major_locator_base=1,
            )
            plt.subplot(4, 1, 4, sharex=plt.gca())
            plot_without_fill(
                plt.gca(),
                freqs,
                Spectrum_input,
                "Frequency (Hz)",
                "Input",
                "",
                show_xlabel=True,
                x_major_locator_base=1,
            )
            plt.xlim(0, 4)
            plt.tight_layout()
            plt.show()

    # Figure 19
    elif fig_number in ["19"]:
        fig, ax = subplots(2, 1, sharex=True, figsize=(12, 15))
        ax[0].plot(V2_PYs.t, V2_PYs.v[0] / mV, color="tab:blue", linewidth=1.5)
        ax[0].set_title("Away from stimulus", size=35, loc="left")
        ax[0].set_ylabel("mV", size=30, labelpad=25)
        ax[0].spines["top"].set_visible(False)
        ax[0].spines["right"].set_visible(False)
        ax[0].spines["bottom"].set_visible(True)
        ax[0].spines["left"].set_visible(True)
        ax[0].yaxis.set_major_locator(MultipleLocator(base=25))
        ax[0].tick_params(axis="both", which="major", labelsize=30, width=2)
        ax[1].plot(V2_PYs.t, V2_PYs.v[50] / mV, color="tab:blue", linewidth=1.5)
        ax[1].set_title("Close to stimulus", size=35, loc="left")
        ax[1].set_ylabel("mV", size=30, labelpad=25)
        ax[1].spines["top"].set_visible(False)
        ax[1].spines["right"].set_visible(False)
        ax[1].spines["bottom"].set_visible(True)
        ax[1].spines["left"].set_visible(True)
        ax[1].yaxis.set_major_locator(MultipleLocator(base=25))
        ax[1].tick_params(axis="both", which="major", labelsize=30, width=2)
        ax[1].set_xlim([0, 5])
        plt.show()

    # Supplementary 1
    elif fig_number == "S1":
        # S1-A
        fig, ax = subplots(4, 1, sharex=True, figsize=(19, 24))
        ax[0].plot(V2_PYs.t, V2_PYs.v[50] / mV, color="tab:blue", linewidth=1.5)
        ax[0].set_title("PY", size=30, loc="left")
        ax[0].set_ylabel("mV", size=30, labelpad=30)
        # ax[0].set_xlim([15000,30000])
        ax[0].spines["top"].set_visible(False)
        ax[0].spines["right"].set_visible(False)
        ax[0].spines["bottom"].set_visible(True)
        ax[0].spines["left"].set_visible(True)
        ax[0].yaxis.set_major_locator(MultipleLocator(base=25))
        ax[0].tick_params(axis="both", which="major", labelsize=25, width=2)
        #
        ax[1].plot(V4_INs.t, V4_INs.v[15] / mV, color="tab:green", linewidth=1.5)
        ax[1].set_title("IN", size=30, loc="left")
        ax[1].set_ylabel("mV", size=30, labelpad=30)
        ax[1].spines["top"].set_visible(False)
        ax[1].spines["right"].set_visible(False)
        ax[1].spines["bottom"].set_visible(True)
        ax[1].spines["left"].set_visible(True)
        ax[1].yaxis.set_major_locator(MultipleLocator(base=25))
        ax[1].tick_params(axis="both", which="major", labelsize=25, width=2)
        #
        ax[2].plot(V1_RE.t, V1_RE.v[25] / mV, color="tab:red", linewidth=1.5)
        ax[2].set_title("RE", size=30, loc="left")
        ax[2].set_ylabel("mV", size=30, labelpad=30)
        ax[2].spines["top"].set_visible(False)
        ax[2].spines["right"].set_visible(False)
        ax[2].spines["bottom"].set_visible(True)
        ax[2].spines["left"].set_visible(True)
        ax[2].yaxis.set_major_locator(MultipleLocator(base=25))
        ax[2].tick_params(axis="both", which="major", labelsize=25, width=2)
        #
        ax[3].plot(V2_TC.t, V2_TC.v[25] / mV, color="tab:orange", linewidth=1.5)
        ax[3].set_title("TC", size=30, loc="left")
        ax[3].set_ylabel("mV", size=30, labelpad=30)
        ax[3].set_xlabel("second", size=30)
        ax[3].tick_params(axis="both", labelbottom=True)
        ax[3].spines["top"].set_visible(False)
        ax[3].spines["right"].set_visible(False)
        ax[3].spines["bottom"].set_visible(True)
        ax[3].spines["left"].set_visible(True)
        ax[3].yaxis.set_major_locator(MultipleLocator(base=25))
        ax[3].tick_params(axis="both", which="major", labelsize=25, width=2)
        # plt.savefig('Supplementary1A.png', dpi=300, bbox_inches='tight')
        # S1-B
        fig, ax = subplots(4, 1, sharex=True, figsize=(19, 24))
        ax[0].plot(V2_PYs.t[::50], V2_PYs.v[50][::50] / mV, color="tab:blue", linewidth=1.5)
        ax[0].set_title("PY", size=30, loc="left")
        ax[0].set_ylabel("mV", size=30, labelpad=30)
        # ax[0].set_xlim([15000,30000])
        ax[0].spines["top"].set_visible(False)
        ax[0].spines["right"].set_visible(False)
        ax[0].spines["bottom"].set_visible(True)
        ax[0].spines["left"].set_visible(True)
        ax[0].yaxis.set_major_locator(MultipleLocator(base=25))
        ax[0].tick_params(axis="both", which="major", labelsize=25, width=2)
        #
        ax[1].plot(V4_INs.t[::50], V4_INs.v[15][::50] / mV, color="tab:green", linewidth=1.5)
        ax[1].set_title("IN", size=30, loc="left")
        ax[1].set_ylabel("mV", size=30, labelpad=30)
        ax[1].spines["top"].set_visible(False)
        ax[1].spines["right"].set_visible(False)
        ax[1].spines["bottom"].set_visible(True)
        ax[1].spines["left"].set_visible(True)
        ax[1].yaxis.set_major_locator(MultipleLocator(base=25))
        ax[1].tick_params(axis="both", which="major", labelsize=25, width=2)
        #
        ax[2].plot(V1_RE.t[::50], V1_RE.v[25][::50] / mV, color="tab:red", linewidth=1.5)
        ax[2].set_title("RE", size=30, loc="left")
        ax[2].set_ylabel("mV", size=30, labelpad=30)
        ax[2].spines["top"].set_visible(False)
        ax[2].spines["right"].set_visible(False)
        ax[2].spines["bottom"].set_visible(True)
        ax[2].spines["left"].set_visible(True)
        ax[2].yaxis.set_major_locator(MultipleLocator(base=25))
        ax[2].tick_params(axis="both", which="major", labelsize=25, width=2)
        #
        ax[3].plot(V2_TC.t[::50], V2_TC.v[25][::50] / mV, color="tab:orange", linewidth=1.5)
        ax[3].set_title("TC", size=30, loc="left")
        ax[3].set_ylabel("mV", size=30, labelpad=30)
        ax[3].set_xlabel("second", size=30)
        ax[3].tick_params(axis="both", labelbottom=True)
        ax[3].spines["top"].set_visible(False)
        ax[3].spines["right"].set_visible(False)
        ax[3].spines["bottom"].set_visible(True)
        ax[3].spines["left"].set_visible(True)
        ax[3].yaxis.set_major_locator(MultipleLocator(base=25))
        ax[3].tick_params(axis="both", which="major", labelsize=25, width=2)

    # Supplementary 2
    elif fig_number == "S2":
        fig, ax = subplots(4, 1, sharex=True, figsize=(19, 24))
        ax[0].plot(time_s, VPY_s[50], color="tab:blue", linewidth=1.5)
        ax[0].set_title("PY axosomatic compartment", size=30, loc="left")
        ax[0].set_ylabel("mV", size=30, labelpad=30)
        # ax[0].set_xlim([15000,30000])
        ax[0].spines["top"].set_visible(False)
        ax[0].spines["right"].set_visible(False)
        ax[0].spines["bottom"].set_visible(True)
        ax[0].spines["left"].set_visible(True)
        ax[0].set_ylim(-80, 50)
        ax[0].yaxis.set_major_locator(MultipleLocator(base=25))
        ax[0].tick_params(axis="both", which="major", labelsize=25, width=2)
        #
        ax[1].plot(time_s, VPY_d[50], color="tab:blue", linewidth=1.5)
        ax[1].set_title("PY dendritic compartment", size=30, loc="left")
        ax[1].set_ylabel("mV", size=30, labelpad=30)
        ax[1].spines["top"].set_visible(False)
        ax[1].spines["right"].set_visible(False)
        ax[1].spines["bottom"].set_visible(True)
        ax[1].spines["left"].set_visible(True)
        ax[1].set_ylim(-80, 50)
        ax[1].yaxis.set_major_locator(MultipleLocator(base=25))
        ax[1].tick_params(axis="both", which="major", labelsize=25, width=2)
        #
        ax[2].plot(time_s, VIN[15], color="tab:green", linewidth=1.5)
        ax[2].set_title("IN axosomatic compartment", size=30, loc="left")
        ax[2].set_ylabel("mV", size=30, labelpad=30)
        ax[2].spines["top"].set_visible(False)
        ax[2].spines["right"].set_visible(False)
        ax[2].spines["bottom"].set_visible(True)
        ax[2].spines["left"].set_visible(True)
        ax[2].set_ylim(-100, 50)
        ax[2].yaxis.set_major_locator(MultipleLocator(base=25))
        ax[2].tick_params(axis="both", which="major", labelsize=25, width=2)
        #
        ax[3].plot(time_s, VIN_d[15], color="tab:green", linewidth=1.5)
        ax[3].set_title("IN dendritic compartment", size=30, loc="left")
        ax[3].set_ylabel("mV", size=30, labelpad=30)
        ax[3].set_xlabel("second", size=30)
        ax[3].tick_params(axis="both", labelbottom=True)
        ax[3].spines["top"].set_visible(False)
        ax[3].spines["right"].set_visible(False)
        ax[3].spines["bottom"].set_visible(True)
        ax[3].spines["left"].set_visible(True)
        ax[3].set_ylim(-100, 50)
        ax[3].yaxis.set_major_locator(MultipleLocator(base=25))
        ax[3].tick_params(axis="both", which="major", labelsize=25, width=2)
        fig.tight_layout()

    elif fig_number in ["S3-B", "S3-C", "S3-D", "S3-E"]:
        fig = plt.figure(figsize=(12, 4))
        gs = gridspec.GridSpec(1, 2, width_ratios=[20, 0.5], height_ratios=[1], wspace=0.05)
        ax1 = plt.subplot(gs[0])
        im1 = ax1.imshow(
            V2_PYs.v / mV,
            aspect="auto",
            cmap="Greys",
            vmin=-75,
            vmax=-60,
            extent=[0, runtime, 100, 0],
            interpolation="bicubic",
        )

        ax1.set_title("PY Neurons", size=30, loc="left", pad=30)
        ax1.set_xlabel("Time (s)", size=25, labelpad=10)
        ax1.set_ylabel("Neuron index", size=25, labelpad=30)
        ax1.set_xlim(0, 15)
        ax1.yaxis.set_major_locator(MultipleLocator(base=25))
        ax1.tick_params(axis="both", which="major", labelsize=25, width=2)

        cbar_ax = plt.subplot(gs[1])
        cbar = fig.colorbar(im1, cax=cbar_ax, orientation="vertical")
        cbar.set_label("Membrane potential (mV)", size=25, labelpad=30)
        cbar.ax.tick_params(labelsize=25, width=2)
        plt.show()

    else:
        fig = plt.figure(figsize=(20, 12))
        gs = gridspec.GridSpec(2, 2, width_ratios=[20, 0.5], height_ratios=[2, 1], wspace=0.05, hspace=0.3)
        ax1 = plt.subplot(gs[0, 0])
        im1 = ax1.imshow(
            V2_PYs.v / mV,
            aspect="auto",
            cmap="Greys",
            vmax=-60,
            vmin=-75,
            extent=[0, runtime, N - 1, 0],
            interpolation="bicubic",
        )
        ax1.set_title("PY", size=30, loc="left")
        ax1.set_ylabel("Neuron index", size=30, labelpad=30)
        ax1.set_xlim(0, 30)
        ax1.yaxis.set_major_locator(MultipleLocator(base=25))
        ax1.tick_params(axis="both", which="major", labelsize=25, width=2)
        ax2 = plt.subplot(gs[1, 0])
        im2 = ax2.imshow(
            V2_TC.v / mV,
            aspect="auto",
            cmap="Greys",
            vmax=-60,
            vmin=-75,
            extent=[0, runtime, 50, 0],
            interpolation="bicubic",
        )
        ax2.set_title("TC", size=30, loc="left")
        ax2.set_xlabel("Time (s)", size=30, labelpad=10)
        ax2.set_ylabel("Neuron index", size=30, labelpad=30)
        ax2.set_xlim(0, 30)
        ax2.yaxis.set_major_locator(MultipleLocator(base=25))
        ax2.tick_params(axis="both", which="major", labelsize=25, width=2)

        cbar_ax = plt.subplot(gs[:, 1])
        cbar = fig.colorbar(im1, cax=cbar_ax, orientation="vertical")
        cbar.set_label("Membrane potential (mV)", size=30, labelpad=30)
        cbar.ax.tick_params(labelsize=25, width=2)
        plt.show()
