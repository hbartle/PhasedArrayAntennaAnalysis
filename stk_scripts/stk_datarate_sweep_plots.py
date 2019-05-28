#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Hannes Bartle <mail@hbartle.de>
#
# Distributed under terms of the MIT license.

import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.dates as mdates
#plt.style.use('classic')
#plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rcParams.update({'font.size': 14})
import matplotlib.gridspec as gridspec
import pandas
import datetime

def dateparse (time_in_secs):
    return datetime.datetime.fromtimestamp(float(time_in_secs))


"""
    Plot the results of the STK datarate sweep
"""


if __name__ == "__main__":

    directory = "C:/Users/hanne/Nextcloud/studies/au/courses/master_thesis/stk_analysis_results/sat2gnd/varying_datarate/"


    filenames = ['Satellite-DelphiniX-Transmitter-PAA2x2-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_0.5Mbps.csv',
                 'Satellite-DelphiniX-Transmitter-PAA2x2-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_1Mbps.csv',
                 'Satellite-DelphiniX-Transmitter-PAA2x2-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_2Mbps.csv',
                 'Satellite-DelphiniX-Transmitter-PAA2x2-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_5Mbps.csv',
                 'Satellite-DelphiniX-Transmitter-PAA2x2-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_10Mbps.csv',
                 'Satellite-DelphiniX-Transmitter-PAA2x2-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_20Mbps.csv',
                 'Satellite-DelphiniX-Transmitter-PAA2x2-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_50Mbps.csv',
                 'Satellite-DelphiniX-Transmitter-PAA2x2-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_100Mbps.csv',
                'Satellite-DelphiniX-Transmitter-PAA2x2-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_150Mbps.csv',
                'Satellite-DelphiniX-Transmitter-PAA2x2-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_200Mbps.csv']



    data = []

    for filename in filenames:
        d = pandas.read_csv(directory + filename,
                            delimiter=',',
                            parse_dates=["Time (UTCG)"],
                            date_parser=lambda dates: [pandas.datetime.strptime(d, '%d %b %Y %H:%M:%S.%f') for d in
                                                       dates])
        data.append(d)

    SNR_avg = []
    SNR_std = []
    for i, d in enumerate(data):
        accesses = np.unique(d["Access Number"].values)
        access_avg = d.groupby(["Access Number"]).mean()
        access_max = d.groupby(["Access Number"]).max()
        access_min = d.groupby(["Access Number"]).min()
        access_std = d.groupby(["Access Number"]).std()
        SNR_avg.append(access_avg["C/N (dB)"].values[0:20].mean())
        SNR_std.append(access_std["C/N (dB)"].values[0:20].mean())

    SNR_avg_paa2 = np.asarray(SNR_avg)
    SNR_std_paa2 = np.asarray(SNR_std)



    filenames = ['Satellite-DelphiniX-Transmitter-PAA4x4-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_0.5Mbps.csv',
                 'Satellite-DelphiniX-Transmitter-PAA4x4-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_1Mbps.csv',
                 'Satellite-DelphiniX-Transmitter-PAA4x4-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_2Mbps.csv',
                 'Satellite-DelphiniX-Transmitter-PAA4x4-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_5Mbps.csv',
                 'Satellite-DelphiniX-Transmitter-PAA4x4-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_10Mbps.csv',
                 'Satellite-DelphiniX-Transmitter-PAA4x4-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_20Mbps.csv',
                 'Satellite-DelphiniX-Transmitter-PAA4x4-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_50Mbps.csv',
                 'Satellite-DelphiniX-Transmitter-PAA4x4-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_100Mbps.csv',
                'Satellite-DelphiniX-Transmitter-PAA4x4-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_150Mbps.csv',
                'Satellite-DelphiniX-Transmitter-PAA4x4-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_200Mbps.csv']

    data = []

    for filename in filenames:
        d = pandas.read_csv(directory + filename,
                            delimiter=',',
                            parse_dates=["Time (UTCG)"],
                            date_parser=lambda dates: [pandas.datetime.strptime(d, '%d %b %Y %H:%M:%S.%f') for d in
                                                       dates])
        data.append(d)

    SNR_avg = []
    SNR_std = []
    for i, d in enumerate(data):
        accesses = np.unique(d["Access Number"].values)
        access_avg = d.groupby(["Access Number"]).mean()
        access_max = d.groupby(["Access Number"]).max()
        access_min = d.groupby(["Access Number"]).min()
        access_std = d.groupby(["Access Number"]).std()
        SNR_avg.append(access_avg["C/N (dB)"].values[0:20].mean())
        SNR_std.append(access_std["C/N (dB)"].values[0:20].mean())

    SNR_avg_paa4 = np.asarray(SNR_avg)
    SNR_std_paa4 = np.asarray(SNR_std)

    #### Patch Antenna #####


    filenames = ['Satellite-DelphiniX-Transmitter-PatchTransmitter-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_0.5Mbps.csv',
                 'Satellite-DelphiniX-Transmitter-PatchTransmitter-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_1Mbps.csv',
                 'Satellite-DelphiniX-Transmitter-PatchTransmitter-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_2Mbps.csv',
                 'Satellite-DelphiniX-Transmitter-PatchTransmitter-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_5Mbps.csv',
                 'Satellite-DelphiniX-Transmitter-PatchTransmitter-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_10Mbps.csv',
                 'Satellite-DelphiniX-Transmitter-PatchTransmitter-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_20Mbps.csv',
                 'Satellite-DelphiniX-Transmitter-PatchTransmitter-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_50Mbps.csv',
                 'Satellite-DelphiniX-Transmitter-PatchTransmitter-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_100Mbps.csv',
                'Satellite-DelphiniX-Transmitter-PatchTransmitter-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_150Mbps.csv',
                'Satellite-DelphiniX-Transmitter-PatchTransmitter-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison_200Mbps.csv']

    data = []

    for filename in filenames:
        d = pandas.read_csv(directory + filename,
                        delimiter=',',
                        parse_dates=["Time (UTCG)"],
                        date_parser=lambda dates: [pandas.datetime.strptime(d, '%d %b %Y %H:%M:%S.%f') for d in dates])
        data.append(d)


    SNR_avg = []
    SNR_std = []
    for i,d in enumerate(data):
        accesses = np.unique(d["Access Number"].values)
        access_avg = d.groupby(["Access Number"]).mean()
        access_max = d.groupby(["Access Number"]).max()
        access_min = d.groupby(["Access Number"]).min()
        access_std = d.groupby(["Access Number"]).std()
        SNR_avg.append(access_avg["C/N (dB)"].values[0:20].mean())
        SNR_std.append(access_std["C/N (dB)"].values[0:20].mean())

    SNR_avg_patch = np.asarray(SNR_avg)
    SNR_std_patch = np.asarray(SNR_std)



    datarates = [0.5,1,2,5,10,20,50,100,150,200]

    gs = gridspec.GridSpec(1, 1)
    plt.figure()
    ax0 = plt.subplot(gs[0])

    ax0.plot(datarates,SNR_avg_patch,".--",label="Patch")
    ax0.fill_between(datarates,
                     SNR_avg_patch - SNR_std_patch,
                     SNR_avg_patch + SNR_std_patch,
                     alpha=0.3)
    ax0.plot(datarates,SNR_avg_paa2,".--",label="PAA2x2")
    ax0.fill_between(datarates,
                     SNR_avg_paa2 - SNR_std_paa2,
                     SNR_avg_paa2 + SNR_std_paa2,
                     alpha=0.3)
    ax0.plot(datarates,SNR_avg_paa4,".--",label="PAA4x4")
    ax0.fill_between(datarates,
                     SNR_avg_paa4 - SNR_std_paa4,
                     SNR_avg_paa4 + SNR_std_paa4,
                     alpha=0.3)
    ax0.set_ylabel("SNR [dB]")
    ax0.set_xlabel("Datarate [Mbps]")
    ax0.grid(True)
    plt.legend(loc='upper right')
    ax0.set_ylim([0,30])

    dic = {"Data Rate":datarates,
           "SNR Patch": SNR_avg_patch,
           "SNR Patch Std": SNR_std_patch,
           "SNR PAA2x2": SNR_avg_paa2,
           "SNR PAA2x2 Std": SNR_std_paa2,
           "SNR PAA4x4": SNR_avg_paa4,
           "SNR PAA4x4 Std": SNR_std_paa4}

    df =  pandas.DataFrame(data=dic)
    df.to_csv(directory + "snr_datarate_sweep.txt", index=None, sep=",", float_format='%.2f')

    plt.show()
