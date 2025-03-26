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

    # For figure 2 to 4, please refer to the corresponding .py files as values are computed using Bazhenov et al. (2002) original results

    # Figure 5
    if fig_number == "5":
        runtime = 5 * second  # 30
    else:
        runtime = 10 * second

    # Figure 7
    if fig_number == "7":
        runtime = 30 * second

    # Figure 8
    if fig_number == "8":
        runtime = 10 * second

    # Figure 9
    if fig_number in ["9-A1", "9-A2", "9-A3", "9-B1", "9-B2", "9-B3"]:
        runtime = 8 * second
        g_syn_ampa_tcpy = 0 * msiemens
        g_syn_ampa_tcin = 0 * msiemens
        g_syn_ampa_pytc = 0 * msiemens
        g_syn_ampa_pyre = 0 * msiemens

    # Figure 10
    if fig_number in ["10-1", "10-3"]:
        runtime = 20 * second
        g_syn_ampa_tcpy = 0 * msiemens
        g_syn_ampa_tcin = 0 * msiemens
        g_syn_ampa_pytc = 0 * msiemens
        g_syn_ampa_pyre = 0 * msiemens
    if fig_number == "10-2":
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

    # Figure 11
    if fig_number == "11-1":
        runtime = 8 * second
    if fig_number == "11-2":
        runtime = 8 * second

    # Figure 12
    if fig_number == "12-A1":
        runtime = 8 * second
        g_syn_ampa_tcpy = 0 * msiemens
        g_syn_ampa_tcin = 0 * msiemens
        g_syn_ampa_pytc = 0 * msiemens
        g_syn_ampa_pyre = 0 * msiemens
        g_syn_ampa_pyin = 0.00002 * msiemens
    if fig_number == "12-A2":
        runtime = 8 * second
        g_syn_ampa_tcpy = 0 * msiemens
        g_syn_ampa_tcin = 0 * msiemens
        g_syn_ampa_pytc = 0 * msiemens
        g_syn_ampa_pyre = 0 * msiemens
        g_syn_ampa_pyin = 0.00007 * msiemens
    if fig_number == "12-B1":
        runtime = 20 * second
        g_syn_ampa_tcpy = 0 * msiemens
        g_syn_ampa_tcin = 0 * msiemens
        g_syn_ampa_pytc = 0 * msiemens
        g_syn_ampa_pyre = 0 * msiemens
        g_syn_ampa_pypy = (
            0.00012 * msiemens
        )  # 0.00012*msiemens ; 0.00013*msiemens ; 0.00014*msiemens ; 0.00015*msiemens ; 0.00016*msiemens
    if fig_number == "12-B2":
        runtime = 20 * second
        g_syn_ampa_tcpy = 0 * msiemens
        g_syn_ampa_tcin = 0 * msiemens
        g_syn_ampa_pytc = 0 * msiemens
        g_syn_ampa_pyre = 0 * msiemens
        g_syn_ampa_pyin = (
            0.00002 * msiemens
        )  # 0.00002*msiemens ; 0.00003*msiemens ; 0.00004*msiemens ; 0.00005*msiemens ; 0.00006*msiemens ; 0.00007*msiemens ; 0.00008*msiemens

    # Figure 13
    if fig_number == "13-A":
        runtime = 10 * second
    if fig_number == "13-B":
        runtime = 10 * second
        g_syn_ampa_tcpy = 0 * msiemens
        g_syn_ampa_tcin = 0 * msiemens
        g_syn_ampa_pytc = 0 * msiemens
        g_syn_ampa_pyre = 0 * msiemens

    # Figure 14
    if fig_number == "14":
        runtime = 10 * second

    # Figure 15
    if fig_number == "15-A":
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
    if fig_number == "15-B":
        runtime = 20 * second
        PY_dendrite.g_kl = 0 * msiemens * cm**-2
        TC.g_kl_TC = 0 * msiemens * cm**-2
        RE.g_kl_RE = 0 * msiemens * cm**-2
        syn_PYPY = all_synapses[0]
        syn_PYPY.g_syn = 0.00009 * msiemens
        syn_RETC = all_synapses_T[0]
        syn_RETC.g_syn = 0.0001 * msiemens
        syn_TCRE = all_synapses_T[-1]
    if fig_number == "15-C":
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
    if fig_number == "15-D":
        runtime = 20 * second
        syn_PYPY = all_synapses[0]
        syn_PYPY.g_syn = (
            0.0001 * msiemens
        )  # 0.00008*msiemens ; 0.0001*msiemens ; 0.00012*msiemens ; 0.00014*msiemens ; 0.00016*msiemens

    # Figure 16
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
    if fig_number in params_16:
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

    # Figure 17/18/19
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
    if fig_number in stim_params:
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

    # Supplementary 1 & 2
    if fig_number == "S1" or fig_number == "S2":
        runtime = 10 * second

    # Supplementary 3
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
    if fig_number in stim_params_S3:
        params = stim_params_S3[fig_number]
        PY_dendrite.g_kl = params["g_kl"]
        TC.g_kl_TC = params["g_kl_TC"]
        syn_PYPY = all_synapses[0]
        syn_PYPY.g_syn = params["g_syn_PYPY"]
        syn_RETC = all_synapses_T[0]
        syn_RETC.g_syn = params["g_syn_RETC"]
        syn_TCRE = all_synapses_T[-1]
        syn_TCRE.g_syn = params["g_syn_TCRE"]

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
