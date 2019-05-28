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
#plt.style.use('seaborn')
#plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rcParams.update({'font.size': 14})
import matplotlib.gridspec as gridspec
import pandas
import datetime

def dateparse (time_in_secs):
    return datetime.datetime.fromtimestamp(float(time_in_secs))


"""
    Plot the average link budget over all accesses 
"""

if __name__ == "__main__":

    # Load Link Budget Comparison
    directory= "C:/Users/hanne/Nextcloud/studies/au/courses/master_thesis/stk_analysis_results/sat2gnd/linkbudgetcomparison_10Mbps/"
    filenames = ["Satellite-DelphiniX-Transmitter-PAA2x2-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison.csv",
                "Satellite-DelphiniX-Transmitter-PAA3x3-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison.csv",
                "Satellite-DelphiniX-Transmitter-PAA4x4-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison.csv",
                "Satellite-DelphiniX-Transmitter-PAA5x5-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison.csv",
                "Satellite-DelphiniX-Transmitter-PatchTransmitter-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison.csv"]

    labels = ["2x2","3x3","4x4","5x5","Patch"]



    data_output_directory = directory + "access_avg/"

    os.makedirs(data_output_directory, exist_ok=True)

    data = []

    for filename in filenames:
        d = pandas.read_csv(directory + filename,
                        delimiter=',',
                        parse_dates=["Time (UTCG)"],
                        date_parser=lambda dates: [pandas.datetime.strptime(d, '%d %b %Y %H:%M:%S.%f') for d in dates])
        data.append(d)

    gs = gridspec.GridSpec(2, 1)
    plt.figure()
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1], sharex=ax0)

    for i,d in enumerate(data):
        accesses = np.unique(d["Access Number"].values)
        access_avg = d.groupby(["Access Number"]).mean()
        access_max = d.groupby(["Access Number"]).max()
        access_min = d.groupby(["Access Number"]).min()
        access_std = d.groupby(["Access Number"]).std()


        ax0.plot(accesses[0:20],access_avg["EIRP (dBW)"].values[0:20],".--",label=labels[i])
        ax0.fill_between(accesses[0:20],
                         access_avg["EIRP (dBW)"].values[0:20] - access_std["EIRP (dBW)"].values[0:20],
                         access_avg["EIRP (dBW)"].values[0:20] + access_std["EIRP (dBW)"].values[0:20],
                         alpha=0.3)
        ax0.set_ylabel("EIRP", )
        ax0.grid(True)

        dic = {"Access Number": accesses[0:20],
               "EIRP Avg": access_avg["EIRP (dBW)"].values[0:20],
               "EIRP Std": access_std["EIRP (dBW)"].values[0:20]}

        df =  pandas.DataFrame(data=dic)
        df.to_csv(data_output_directory + "eirp_access" + labels[i] + ".txt", index=None, sep=",",float_format='%.2f')


        ax1.plot(accesses[0:20],access_avg["C/N (dB)"].values[0:20],".--",label=labels[i])
        ax1.fill_between(accesses[0:20],
                         access_avg["C/N (dB)"].values[0:20] - access_std["C/N (dB)"].values[0:20],
                         access_avg["C/N (dB)"].values[0:20] + access_std["C/N (dB)"].values[0:20],
                         alpha=0.3)
        ax1.set_ylabel("C/N", )
        ax1.grid(True)

        dic = {"Access Number": accesses[0:20],
               "SNR Avg": access_avg["C/N (dB)"].values[0:20],
               "SNR Std": access_std["C/N (dB)"].values[0:20]}

        df =  pandas.DataFrame(data=dic)
        df.to_csv(data_output_directory + "snr_access" + labels[i] + ".txt", index=None, sep=",",float_format='%.2f')

    plt.legend(loc='lower left', bbox_to_anchor=(1, 0.5))
    plt.show()
