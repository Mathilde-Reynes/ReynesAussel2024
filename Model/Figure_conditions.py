#!/usr/bin/env python
# coding: utf-8

from brian2 import *
import numpy as np
import matplotlib.pyplot as plt
from brian2.units.constants import *
import matplotlib.gridspec as gridspec


def figure_conditions_pre(fig_number, A_PY_PY, A_PY_IN, N):
    N_PY = N
    N_TC = N_PY // 2
    N_RE = N_PY // 2
    N_IN = N_PY // 4

    # Figure 9
    if fig_number == "9-A1":
        N = 20
    elif fig_number == "9-A2":
        N = 60
    elif fig_number == "9-A3":
        N = 100
    elif fig_number == "9-B1":
        N = 20
        A_PY_PY = 0.00006 * msiemens * 1.2
    elif fig_number == "9-B2":
        N = 60
        A_PY_PY = 0.00006 * msiemens * 1.4
    elif fig_number == "9-B3":
        N = 100
        A_PY_PY = 0.00006 * msiemens * 1.6

    # Figure 10
    elif fig_number == "10-3":
        A_PY_PY = 0.00006 * msiemens * 0.5
        A_PY_IN = 0.000025 * msiemens * 0.5

    elif fig_number == "S4-D":
        A_PY_PY = 0.00006 * msiemens * 0.8

    else:
        None

    return A_PY_PY, A_PY_IN, N


def figure_conditions(
    fig_number,
    all_synapses,
    all_synapses_T,
    all_neurons_T,
    all_neurons,
    g_syn_ampa_tcpy,
    g_syn_ampa_tcin,
    g_syn_ampa_pytc,
    g_syn_ampa_pyre,
):
    RE, TC = all_neurons_T
    PY_dendrite, PY_soma, IN_dendrite, IN_soma = all_neurons

    # Default values for synapses before the conditionals
    syn_PYPY = None
    syn_PYPY_nmda = None
    syn_RETC = None
    syn_TCRE = None
    syn_PYIN = None
    syn_PYIN_nmda = None
    syn_INPY = None
    modulation = None
    base_rate = None
    monitor_poisson = None
    g_syn_ampa_stim = None

    # Dictionnaries

    params_16 = {
        "16-A1": {
            "runtime": 30 * second,
            "g_kl_PY": 0.0025 * msiemens * cm**-2,
            "g_kl_TC": 0.003 * msiemens * cm**-2,
            "syn_PYPY": 0.00015 * msiemens,
            "syn_RETC": 0.0002 * msiemens,
            "syn_TCRE": 0.0004 * msiemens,
        },
        "16-B1": {
            "runtime": 30 * second,
            "g_kl_PY": 0.0025 * msiemens * cm**-2,
            "g_kl_TC": 0.003 * msiemens * cm**-2,
            "syn_PYPY": 0.00015 * msiemens,
            "syn_RETC": 0.0002 * msiemens,
            "syn_TCRE": 0.0004 * msiemens,
        },
        "16-A2": {
            "runtime": 30 * second,
            "g_kl_PY": 0.0025 * 4 / 6 * msiemens * cm**-2,
            "g_kl_TC": 0.002 * msiemens * cm**-2,
            "syn_PYPY": 0.0001266 * msiemens,
            "syn_RETC": 0.0001666 * msiemens,
            "syn_TCRE": 0.000333 * msiemens,
        },
        "16-B3": {
            "runtime": 30 * second,
            "g_kl_PY": 0.0025 * 4 / 6 * msiemens * cm**-2,
            "g_kl_TC": 0.002 * msiemens * cm**-2,
            "syn_PYPY": 0.0001266 * msiemens,
            "syn_RETC": 0.0001666 * msiemens,
            "syn_TCRE": 0.000333 * msiemens,
        },
        "16-A3": {
            "runtime": 30 * second,
            "g_kl_PY": 0.0025 * 2 / 6 * msiemens * cm**-2,
            "g_kl_TC": 0.001 * msiemens * cm**-2,
            "syn_PYPY": 0.0001033 * msiemens,
            "syn_RETC": 0.0001333 * msiemens,
            "syn_TCRE": 0.0002666 * msiemens,
        },
        "16-B5": {
            "runtime": 30 * second,
            "g_kl_PY": 0.0025 * 2 / 6 * msiemens * cm**-2,
            "g_kl_TC": 0.001 * msiemens * cm**-2,
            "syn_PYPY": 0.0001033 * msiemens,
            "syn_RETC": 0.0001333 * msiemens,
            "syn_TCRE": 0.0002666 * msiemens,
        },
        "16-A4": {
            "runtime": 30 * second,
            "g_kl_PY": 0 * msiemens * cm**-2,
            "g_kl_TC": 0 * msiemens * cm**-2,
            "syn_PYPY": 0.00009 * msiemens,
            "syn_RETC": 0.0001 * msiemens,
            "syn_TCRE": 0.0002 * msiemens,
        },
        "16-B7": {
            "runtime": 30 * second,
            "g_kl_PY": 0 * msiemens * cm**-2,
            "g_kl_TC": 0 * msiemens * cm**-2,
            "syn_PYPY": 0.00009 * msiemens,
            "syn_RETC": 0.0001 * msiemens,
            "syn_TCRE": 0.0002 * msiemens,
        },
        "16-B2": {
            "runtime": 30 * second,
            "g_kl_PY": 0.0025 * 5 / 6 * msiemens * cm**-2,
            "g_kl_TC": 0.0025 * msiemens * cm**-2,
            "syn_PYPY": 0.00013833 * msiemens,
            "syn_RETC": 0.00018333 * msiemens,
            "syn_TCRE": 0.0003666 * msiemens,
        },
        "16-B4": {
            "runtime": 30 * second,
            "g_kl_PY": 0.0025 * 3 / 6 * msiemens * cm**-2,
            "g_kl_TC": 0.0015 * msiemens * cm**-2,
            "syn_PYPY": 0.000115 * msiemens,
            "syn_RETC": 0.00015 * msiemens,
            "syn_TCRE": 0.0003 * msiemens,
        },
        "16-B6": {
            "runtime": 30 * second,
            "g_kl_PY": 0.00025 * 1 / 6 * msiemens * cm**-2,
            "g_kl_TC": 0.0005 * msiemens * cm**-2,
            "syn_PYPY": 0.00009166 * msiemens,
            "syn_RETC": 0.0001166 * msiemens,
            "syn_TCRE": 0.0002333 * msiemens,
        },
    }

    stim_params = {
        "17-A1": {
            "runtime": 15 * second,
            "modulation": 0.4 * Hz,
            "g_syn": 0.00015 * msiemens,
            "base_rate": 25 * Hz,
        },
        "18-A1": {
            "runtime": 15 * second,
            "modulation": 0.4 * Hz,
            "g_syn": 0.00015 * msiemens,
            "base_rate": 25 * Hz,
        },
        "17-A2": {
            "runtime": 15 * second,
            "modulation": 1 * Hz,
            "g_syn": 0.00015 * msiemens,
            "base_rate": 25 * Hz,
        },
        "18-A2": {
            "runtime": 15 * second,
            "modulation": 1 * Hz,
            "g_syn": 0.00015 * msiemens,
            "base_rate": 25 * Hz,
        },
        "17-A3": {
            "runtime": 15 * second,
            "modulation": 2.5 * Hz,
            "g_syn": 0.00015 * msiemens,
            "base_rate": 25 * Hz,
        },
        "18-A3": {
            "runtime": 15 * second,
            "modulation": 2.5 * Hz,
            "g_syn": 0.00015 * msiemens,
            "base_rate": 25 * Hz,
        },
        "17-B1": {
            "runtime": 15 * second,
            "modulation": 0.4 * Hz,
            "g_syn": 0.00009 * msiemens,
            "base_rate": 25 * Hz,
        },
        "18-B1": {
            "runtime": 15 * second,
            "modulation": 0.4 * Hz,
            "g_syn": 0.00009 * msiemens,
            "base_rate": 25 * Hz,
        },
        "17-B2": {
            "runtime": 15 * second,
            "modulation": 1 * Hz,
            "g_syn": 0.00009 * msiemens,
            "base_rate": 25 * Hz,
        },
        "18-B2": {
            "runtime": 15 * second,
            "modulation": 1 * Hz,
            "g_syn": 0.00009 * msiemens,
            "base_rate": 25 * Hz,
        },
        "17-B3": {
            "runtime": 15 * second,
            "modulation": 2.5 * Hz,
            "g_syn": 0.00009 * msiemens,
            "base_rate": 25 * Hz,
        },
        "18-B3": {
            "runtime": 15 * second,
            "modulation": 2.5 * Hz,
            "g_syn": 0.00009 * msiemens,
            "base_rate": 25 * Hz,
        },
        "19": {
            "runtime": 15 * second,
            "modulation": 2.5 * Hz,
            "g_syn": 0.00009 * msiemens,
            "base_rate": 25 * Hz,
        },
    }

    stim_params_S3 = {
        "S3-A1": {
            "g_kl": 0.0025 * msiemens * cm**-2,
            "g_syn_PYPY": 0.00015 * msiemens,
            "g_syn_RETC": 0.0002 * msiemens,
            "g_syn_TCRE": 0.0004 * msiemens,
            "g_kl_TC": 0.003 * msiemens * cm**-2,
        },
        "S3-B1": {
            "g_kl": 0.0025 * msiemens * cm**-2,
            "g_syn_PYPY": 0.00015 * msiemens,
            "g_syn_RETC": 0.0002 * msiemens,
            "g_syn_TCRE": 0.0004 * msiemens,
            "g_kl_TC": 0.003 * msiemens * cm**-2,
        },
        "S3-A2": {
            "g_kl": 0.0025 * 4 / 6 * msiemens * cm**-2,
            "g_syn_PYPY": 0.0001266 * msiemens,
            "g_syn_RETC": 0.0002 * msiemens,
            "g_syn_TCRE": 0.0004 * msiemens,
            "g_kl_TC": 0.002 * msiemens * cm**-2,
        },
        "S3-B3": {
            "g_kl": 0.0025 * 4 / 6 * msiemens * cm**-2,
            "g_syn_PYPY": 0.0001266 * msiemens,
            "g_syn_RETC": 0.0002 * msiemens,
            "g_syn_TCRE": 0.0004 * msiemens,
            "g_kl_TC": 0.002 * msiemens * cm**-2,
        },
        "S3-A3": {
            "g_kl": 0.0025 * 2 / 6 * msiemens * cm**-2,
            "g_syn_PYPY": 0.0001033 * msiemens,
            "g_syn_RETC": 0.0002 * msiemens,
            "g_syn_TCRE": 0.0004 * msiemens,
            "g_kl_TC": 0.001 * msiemens * cm**-2,
        },
        "S3-B5": {
            "g_kl": 0.0025 * 2 / 6 * msiemens * cm**-2,
            "g_syn_PYPY": 0.0001033 * msiemens,
            "g_syn_RETC": 0.0002 * msiemens,
            "g_syn_TCRE": 0.0004 * msiemens,
            "g_kl_TC": 0.001 * msiemens * cm**-2,
        },
        "S3-A4": {
            "g_kl": 0 * msiemens * cm**-2,
            "g_syn_PYPY": 0.00008 * msiemens,
            "g_syn_RETC": 0.0002 * msiemens,
            "g_syn_TCRE": 0.0004 * msiemens,
            "g_kl_TC": 0 * msiemens * cm**-2,
        },
        "S3-B7": {
            "g_kl": 0 * msiemens * cm**-2,
            "g_syn_PYPY": 0.00008 * msiemens,
            "g_syn_RETC": 0.0002 * msiemens,
            "g_syn_TCRE": 0.0004 * msiemens,
            "g_kl_TC": 0 * msiemens * cm**-2,
        },
        "S3-B2": {
            "g_kl": 0.0025 * 5 / 6 * msiemens * cm**-2,
            "g_syn_PYPY": 0.00013833 * msiemens,
            "g_syn_RETC": 0.0002 * msiemens,
            "g_syn_TCRE": 0.0004 * msiemens,
            "g_kl_TC": 0.0025 * msiemens * cm**-2,
        },
        "S3-B4": {
            "g_kl": 0.0025 * 3 / 6 * msiemens * cm**-2,
            "g_syn_PYPY": 0.000115 * msiemens,
            "g_syn_RETC": 0.0002 * msiemens,
            "g_syn_TCRE": 0.0004 * msiemens,
            "g_kl_TC": 0.0015 * msiemens * cm**-2,
        },
        "S3-B6": {
            "g_kl": 0.00025 * 1 / 6 * msiemens * cm**-2,
            "g_syn_PYPY": 0.00009166 * msiemens,
            "g_syn_RETC": 0.0002 * msiemens,
            "g_syn_TCRE": 0.0004 * msiemens,
            "g_kl_TC": 0.0005 * msiemens * cm**-2,
        },
    }

    # For figure 2 to 4, please refer to the corresponding .py files as values are computed using Bazhenov et al. (2002) original results

    if fig_number == "5":
        runtime = 15 * second  # 30
    elif fig_number == "7":
        runtime = 30 * second
    elif fig_number == "8":
        runtime = 10 * second
    elif fig_number in ["9-A1", "9-A2", "9-A3", "9-B1", "9-B2", "9-B3"]:
        runtime = 8 * second
        g_syn_ampa_tcpy = 0 * msiemens
        g_syn_ampa_tcin = 0 * msiemens
        g_syn_ampa_pytc = 0 * msiemens
        g_syn_ampa_pyre = 0 * msiemens
    elif fig_number in ["10-1", "10-3"]:
        runtime = 20 * second
        g_syn_ampa_tcpy = 0 * msiemens
        g_syn_ampa_tcin = 0 * msiemens
        g_syn_ampa_pytc = 0 * msiemens
        g_syn_ampa_pyre = 0 * msiemens
    elif fig_number == "10-2":
        runtime = 20 * second
        g_syn_ampa_tcpy = 0 * msiemens
        g_syn_ampa_tcin = 0 * msiemens
        g_syn_ampa_pytc = 0 * msiemens
        g_syn_ampa_pyre = 0 * msiemens
        syn_PYPY = all_synapses[0]
        syn_PYPY.g_syn = 0.00015 * msiemens * 0.5
        syn_PYPY_nmda = all_synapses[2]
        syn_PYPY_nmda.g_syn = 0.00001 * msiemens * 0.5
        syn_PYIN = all_synapses[1]
        syn_PYIN.g_syn = 0.00005 * msiemens * 0.5
        syn_PYIN_nmda = all_synapses[3]
        syn_PYIN_nmda.g_syn = 0.000008 * msiemens * 0.5
        syn_INPY = all_synapses[4]
        syn_INPY.g_syn = 0.00005 * msiemens * 0.5
    elif fig_number == "11-1":
        runtime = 8 * second
    elif fig_number == "11-2":
        runtime = 8 * second
    elif fig_number == "12-A1":
        runtime = 8 * second
        g_syn_ampa_tcpy = 0 * msiemens
        g_syn_ampa_tcin = 0 * msiemens
        g_syn_ampa_pytc = 0 * msiemens
        g_syn_ampa_pyre = 0 * msiemens
        syn_PYIN = all_synapses[1]
        syn_PYIN.g_syn = 0.00002 * msiemens
    elif fig_number == "12-A2":
        runtime = 8 * second
        g_syn_ampa_tcpy = 0 * msiemens
        g_syn_ampa_tcin = 0 * msiemens
        g_syn_ampa_pytc = 0 * msiemens
        g_syn_ampa_pyre = 0 * msiemens
        syn_PYIN = all_synapses[1]
        syn_PYIN.g_syn = 0.00007 * msiemens
    elif fig_number == "12-B1":
        runtime = 20 * second
        g_syn_ampa_tcpy = 0 * msiemens
        g_syn_ampa_tcin = 0 * msiemens
        g_syn_ampa_pytc = 0 * msiemens
        g_syn_ampa_pyre = 0 * msiemens
        syn_PYPY = all_synapses[0]
        syn_PYPY.g_syn = (
            0.00012 * msiemens
        )  # 0.00012*msiemens ; 0.00013*msiemens ; 0.00014*msiemens ; 0.00015*msiemens ; 0.00016*msiemens
    elif fig_number == "12-B2":
        runtime = 20 * second
        g_syn_ampa_tcpy = 0 * msiemens
        g_syn_ampa_tcin = 0 * msiemens
        g_syn_ampa_pytc = 0 * msiemens
        g_syn_ampa_pyre = 0 * msiemens
        syn_PYIN = all_synapses[1]
        syn_PYIN.g_syn = (
            0.00002 * msiemens
        )  # 0.00002*msiemens ; 0.00003*msiemens ; 0.00004*msiemens ; 0.00005*msiemens ; 0.00006*msiemens ; 0.00007*msiemens ; 0.00008*msiemens
    elif fig_number == "13-A":
        runtime = 10 * second
    elif fig_number == "13-B":
        runtime = 10 * second
        g_syn_ampa_tcpy = 0 * msiemens
        g_syn_ampa_tcin = 0 * msiemens
        g_syn_ampa_pytc = 0 * msiemens
        g_syn_ampa_pyre = 0 * msiemens
    elif fig_number == "14":
        runtime = 10 * second
    elif fig_number == "15-A":
        runtime = 20 * second
        PY_dendrite.g_kl = 0 * msiemens * cm**-2
        TC.g_kl_TC = 0 * msiemens * cm**-2
        RE.g_kl_RE = 0 * msiemens * cm**-2
        syn_PYPY = all_synapses[0]
        syn_PYPY.g_syn = 0.00015 * msiemens
        syn_RETC = all_synapses_T[0]
        syn_RETC.g_syn = 0.0002 * msiemens
        syn_TCRE = all_synapses_T[-1]
        syn_TCRE.g_syn = 0.0004 * msiemens
    elif fig_number == "15-B":
        runtime = 20 * second
        PY_dendrite.g_kl = 0 * msiemens * cm**-2
        TC.g_kl_TC = 0 * msiemens * cm**-2
        RE.g_kl_RE = 0 * msiemens * cm**-2
        syn_PYPY = all_synapses[0]
        syn_PYPY.g_syn = 0.00009 * msiemens
        syn_RETC = all_synapses_T[0]
        syn_RETC.g_syn = 0.0001 * msiemens
        syn_TCRE = all_synapses_T[-1]
    elif fig_number == "15-C":
        runtime = 20 * second
        PY_dendrite.g_kl = 0 * msiemens * cm**-2
        TC.g_kl_TC = 0 * msiemens * cm**-2
        RE.g_kl_RE = 0 * msiemens * cm**-2
        syn_PYPY = all_synapses[0]
        syn_PYPY.g_syn = 0.00009 * msiemens
        syn_RETC = all_synapses_T[0]
        syn_RETC.g_syn = 0.0002 * msiemens
        syn_TCRE = all_synapses_T[-1]
        syn_TCRE.g_syn = 0.0004 * msiemens
    elif fig_number == "15-D":
        runtime = 20 * second
        syn_PYPY = all_synapses[0]
        syn_PYPY.g_syn = (
            0.0001 * msiemens
        )  # 0.00008*msiemens ; 0.0001*msiemens ; 0.00012*msiemens ; 0.00014*msiemens ; 0.00016*msiemens
    elif fig_number in params_16:
        params = params_16[fig_number]
        runtime = params["runtime"]
        PY_dendrite.g_kl = params["g_kl_PY"]
        TC.g_kl_TC = params["g_kl_TC"]
        syn_PYPY = all_synapses[0]
        syn_PYPY.g_syn = params["syn_PYPY"]
        syn_RETC = all_synapses_T[0]
        syn_RETC.g_syn = params["syn_RETC"]
        syn_TCRE = all_synapses_T[-1]
        syn_TCRE.g_syn = params["syn_TCRE"]
    elif fig_number in stim_params:
        params = stim_params[fig_number]
        runtime = params["runtime"]
        modulation = params["modulation"]
        g_syn = params["g_syn"]
        base_rate = params["base_rate"]
        PY_dendrite.g_kl = 0.003 * msiemens * cm**-2
        TC.g_kl_TC = 0.003 * msiemens * cm**-2
        syn_PYPY = all_synapses[0]
        syn_PYPY.g_syn = g_syn
        syn_RETC = all_synapses_T[0]
        syn_RETC.g_syn = 0.0002 * msiemens
        syn_TCRE = all_synapses_T[-1]
        syn_TCRE.g_syn = 0.0004 * msiemens
        g_syn_ampa_stim = 0.0004 * msiemens
    elif fig_number == "S1" or fig_number == "S2":
        runtime = 10 * second
    elif fig_number in stim_params_S3:
        params = stim_params_S3[fig_number]
        PY_dendrite.g_kl = params["g_kl"]
        TC.g_kl_TC = params["g_kl_TC"]
        syn_PYPY = all_synapses[0]
        syn_PYPY.g_syn = params["g_syn_PYPY"]
        syn_RETC = all_synapses_T[0]
        syn_RETC.g_syn = params["g_syn_RETC"]
        syn_TCRE = all_synapses_T[-1]
        syn_TCRE.g_syn = params["g_syn_TCRE"]
    elif fig_number == "S4-B":
        runtime = 15 * second
    elif fig_number == "S4-C":
        g_syn_ampa_tcpy = 0.00001 * msiemens
        g_syn_ampa_tcin = 0.00001 * msiemens
        runtime = 15 * second
    elif fig_number == "S4-D":
        runtime = 15 * second
    elif fig_number == "S4-E":
        PY_dendrite.g_kca = 0.06 * msiemens * cm**-2
        runtime = 15 * second
    else:
        runtime = 10 * second
        print(
            "As the fig_number provided was not recognized, a 10-seconds simulation with standard parameters will run (as per done for figure 5). For reference, possible fig_number are: 5, 7, 8, 9-A1, 9-A2, 9-A3, 9-B1, 9-B2, 9-B3, 10-1, 10-2, 10-3, 11-1, 11-2, 12-A1, 12-A2, 12-B1, 12-B2, 13-A, 13-B, 14, 15-A, 15-B, 15-C, 15-D, 16-A1, 16-A2, 16-A3, 16-A4, 16-B1, 16-B2, 16-B3, 16-B5, 16-B4, 16-B6, 16-B7, 17-A1, 17-A2, 17-A3, 18-A1, 18-A2, 18-A3, 17-B1, 17-B2, 17-B3, 18-B1, 18-B2, 18-B3, 19, S1, S2, S3-A1, S3-A2, S3-A3, S3-A4, S3-B1, S3-B2, S3-B3, S3-B4, S3-B5, S3-B6, S3-B7, S4-B, S4-C, S4-D, S4-E"
        )

    return (
        runtime,
        g_syn_ampa_tcpy,
        g_syn_ampa_tcin,
        g_syn_ampa_pytc,
        g_syn_ampa_pyre,
        modulation,
        base_rate,
        g_syn_ampa_stim,
        monitor_poisson,
    )
