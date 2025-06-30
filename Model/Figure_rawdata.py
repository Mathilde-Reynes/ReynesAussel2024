#!/usr/bin/env python
# coding: utf-8

from brian2 import *
import numpy as np
import matplotlib.pyplot as plt
from brian2.units.constants import *
import matplotlib.gridspec as gridspec

import os

print(os.getcwd())


def figure_rawdata(fig_number, all_monitors, all_monitors_T, monitor_poisson):
    V1_PYd, V2_PYs, V3_INd, V4_INs, R2_PYs, R4_INs, I1_PYd, I2_INd, S1, S2, M0, M1 = all_monitors
    V1_RE, V2_TC, R1_RE, R2_TC, I1_RE, I2_TC = all_monitors_T

    # Figure 5
    if fig_number == "5":
        folder_name = f"Figure {fig_number}"
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
        try:
            np.savetxt(os.path.join(folder_name, "PYs_v.txt"), V2_PYs.v / mV)
            np.savetxt(os.path.join(folder_name, "TC_v.txt"), V2_TC.v / mV)
            print(f"Raw data saved to {os.path.abspath(folder_name)}")
        except Exception as e:
            print(f"Error saving raw data: {e}")

    # Figure 7
    elif fig_number == "7":
        folder_name = f"Figure {fig_number}"
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
        savetxt(os.path.join(folder_name, "PYs_v.txt"), V2_PYs.v / mV)
        savetxt(os.path.join(folder_name, "INs_v.txt"), V4_INs.v / mV)
        savetxt(os.path.join(folder_name, "PYd_v.txt"), V1_PYd.v / mV)
        savetxt(os.path.join(folder_name, "INd_v.txt"), V3_INd.v / mV)
        savetxt(os.path.join(folder_name, "TC_v.txt"), V2_TC.v / mV)
        savetxt(os.path.join(folder_name, "RE_v.txt"), V1_RE.v / mV)
        print(f"Raw data saved to {os.path.abspath(folder_name)}")

    # Figure 8
    elif fig_number == "8":
        folder_name = f"Figure {fig_number}"
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
        savetxt(os.path.join(folder_name, "PYs_v.txt"), V2_PYs.v / mV)
        savetxt(os.path.join(folder_name, "IEPSPs_PY_PY.txt"), I1_PYd.IEPSPs_PY_PY)
        savetxt(os.path.join(folder_name, "IEPSPs_PY_IN.txt"), I1_PYd.IEPSPs_PY_IN)
        savetxt(os.path.join(folder_name, "PYd_Inap.txt"), I1_PYd.I_nap)
        savetxt(os.path.join(folder_name, "PYd_Ina.txt"), I1_PYd.I_na)
        savetxt(os.path.join(folder_name, "AMPA_PY_PY.txt"), I1_PYd.IsynAMPA_PY_PY)
        savetxt(os.path.join(folder_name, "NMDA_PY_PY.txt"), I1_PYd.IsynNMDA_PY_PY)
        savetxt(os.path.join(folder_name, "PYd_Ikca.txt"), I1_PYd.I_kca)
        savetxt(os.path.join(folder_name, "GABAA_IN_PY.txt"), I1_PYd.IsynGABAA_IN_PY)
        print(f"Raw data saved to {os.path.abspath(folder_name)}")

    # Figure 9/10/11/19
    elif fig_number in [
        "9-A1",
        "9-A2",
        "9-A3",
        "9-B1",
        "9-B2",
        "9-B3",
        "10-1",
        "10-2",
        "10-3",
        "11-1",
        "11-2",
        "19",
    ]:
        folder_name = f"Figure {fig_number}"
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
        savetxt(os.path.join(folder_name, "PYs_v.txt"), V2_PYs.v / mV)
        print(f"Raw data saved to {os.path.abspath(folder_name)}")

    # Figure 12/13
    elif fig_number in ["12-A1", "12-A2", "13-A", "13-B"]:
        folder_name = f"Figure {fig_number}"
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
        savetxt(os.path.join(folder_name, "PY_raster.txt"), R2_PYs.i)
        savetxt(os.path.join(folder_name, "PY_raster.txt"), R2_PYs.t)
        print(f"Raw data saved to {os.path.abspath(folder_name)}")

    # Figure 14
    elif fig_number == "14":
        folder_name = f"Figure {fig_number}"
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
        savetxt(os.path.join(folder_name, "PYs_v.txt"), V2_PYs.v / mV)
        savetxt(os.path.join(folder_name, "TC_v.txt"), V2_TC.v / mV)
        savetxt(os.path.join(folder_name, "RE_v.txt"), V1_RE.v / mV)
        savetxt(os.path.join(folder_name, "It_TC.txt"), I2_TC.I_t)
        savetxt(os.path.join(folder_name, "IsynAMPA_PY_RE.txt"), I1_RE.IsynAMPA_PY_RE)
        print(f"Raw data saved to {os.path.abspath(folder_name)}")

    # Figure 15
    elif fig_number in ["15-A", "15-B", "15-C", "15-D"]:
        folder_name = f"Figure {fig_number}"
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
        savetxt(os.path.join(folder_name, "PYs_v.txt"), V2_PYs.v / mV)
        savetxt(os.path.join(folder_name, "TC_v.txt"), V2_TC.v / mV)
        savetxt(os.path.join(folder_name, "RE_v.txt"), V1_RE.v / mV)
        print(f"Raw data saved to {os.path.abspath(folder_name)}")

    # Figure 16/S5
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
        "S5-A1",
        "S5-B1",
        "S5-A2",
        "S5-B3",
        "S5-A3",
        "S5-B5",
        "S5-A4",
        "S5-B7",
        "S5-B2",
        "S5-B4",
        "S5-B6"
    ]:
        folder_name = f"Figure {fig_number}"
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
        savetxt(os.path.join(folder_name, "PYs_v.txt"), V2_PYs.v / mV)
        print(f"Raw data saved to {os.path.abspath(folder_name)}")

    # Figure 17/18
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
        folder_name = f"Figure {fig_number}"
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
        savetxt(os.path.join(folder_name, "PY_raster_i.txt"), R2_PYs.i)
        savetxt(os.path.join(folder_name, "PY_raster_t.txt"), R2_PYs.t)
        savetxt(os.path.join(folder_name, "TC_raster_i.txt"), R2_TC.i)
        savetxt(os.path.join(folder_name, "TC_raster_t.txt"), R2_TC.t)
        savetxt(os.path.join(folder_name, "input_i.txt"), monitor_poisson.i)
        savetxt(os.path.join(folder_name, "input_t.txt"), monitor_poisson.t)
        print(f"Raw data saved to {os.path.abspath(folder_name)}")

    # S1
    elif fig_number == "S1":
        folder_name = f"Figure {fig_number}"
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
        savetxt(os.path.join(folder_name, "PYs_v.txt"), V2_PYs.v / mV)
        savetxt(os.path.join(folder_name, "INs_v.txt"), V4_INs.v / mV)
        savetxt(os.path.join(folder_name, "TC_v.txt"), V2_TC.v / mV)
        savetxt(os.path.join(folder_name, "RE_v.txt"), V1_RE.v / mV)
        print(f"Raw data saved to {os.path.abspath(folder_name)}")

    # S2
    elif fig_number == "S2":
        folder_name = f"Figure {fig_number}"
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
        savetxt(os.path.join(folder_name, "PYs_v.txt"), V2_PYs.v / mV)
        savetxt(os.path.join(folder_name, "INs_v.txt"), V4_INs.v / mV)
        savetxt(os.path.join(folder_name, "PYd_v.txt"), V1_PYd.v / mV)
        savetxt(os.path.join(folder_name, "INd_v.txt"), V3_INd.v / mV)
        print(f"Raw data saved to {os.path.abspath(folder_name)}")

    # S3
    elif fig_number in ["S3-B", "S3-C", "S3-D", "S3-E"]:
        folder_name = f"Figure {fig_number}"
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
        try:
            np.savetxt(os.path.join(folder_name, "PYs_v.txt"), V2_PYs.v / mV)
            np.savetxt(os.path.join(folder_name, "TC_v.txt"), V2_TC.v / mV)
            print(f"Raw data saved to {os.path.abspath(folder_name)}")
        except Exception as e:
            print(f"Error saving raw data: {e}")

    else:
        None
