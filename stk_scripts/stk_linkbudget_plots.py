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
    Plot several link budget simulation results
"""


if __name__ == "__main__":

    # Load Data
    directory= "C:/Users/hanne/Nextcloud/studies/au/courses/master_thesis/stk_analysis_results/sat2gnd/linkbudgetcomparison_10Mbps/"
    filenames = ["Satellite-DelphiniX-Transmitter-PAA2x2-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison.csv",
                "Satellite-DelphiniX-Transmitter-PAA3x3-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison.csv",
                "Satellite-DelphiniX-Transmitter-PAA4x4-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison.csv",
                "Satellite-DelphiniX-Transmitter-PAA5x5-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison.csv",
                "Satellite-DelphiniX-Transmitter-PatchTransmitter-To-Place-Aarhus-Receiver-Receiver2 LinkBudgetCustomComparison.csv"]

    labels = ["2x2","3x3","4x4","5x5","Patch"]


    data = []

    for filename in filenames:
        d = pandas.read_csv(directory + filename,
                        delimiter=',',
                        parse_dates=["Time (UTCG)"],
                        date_parser=lambda dates: [pandas.datetime.strptime(d, '%d %b %Y %H:%M:%S.%f') for d in dates])
        data.append(d)



    # Select which access to plot

    access=18 # Good Pass
    data_output_directory = directory + "goodpass/"
    #access=16  #Bad Pass
    #data_output_directory = directory + "badpass/"

    os.makedirs(data_output_directory, exist_ok=True)

    if not access == []:
        a = data[0].loc[d["Access Number"] == access]
    else:
        a = d


    r = a["Range (km)"]
    el = a["Elevation (deg)"]
    t = a["Time (UTCG)"]

    cn = []
    eirp = []
    for d in data:
        if not access == []:
            a = d.loc[d["Access Number"] == access]
        else:
            a = d
        cn.append(a["C/N (dB)"])
        eirp.append(a["EIRP (dBW)"])



    # Plot
    gs = gridspec.GridSpec(4, 1)
    plt.figure()
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1],sharex=ax0)
    ax2 = plt.subplot(gs[2],sharex=ax0)
    ax3 = plt.subplot(gs[3],sharex=ax0)

    ax0.plot(t, r, "r")
    ax0.set_ylabel("Range [km]", )
    ax0.grid(True)
    df = pandas.concat([t, r], axis=1)
    df.to_csv(data_output_directory + "range.txt", index=None, sep=",", float_format='%.2f')

    ax1.plot(t, el, "r")
    ax1.set_ylabel("Elevation [deg]", )
    ax1.grid(True)
    df = pandas.concat([t, el], axis=1)
    df.to_csv(data_output_directory + "elevation.txt", index=None, sep=",", float_format='%.2f')


    for i,x in enumerate(eirp):
        ax2.plot(t,x,label=labels[i])
        df =  pandas.concat([t, x],axis=1)
        df.to_csv(data_output_directory + "eirp" + labels[i] + ".txt", index=None, sep="," ,float_format='%.2f')

    ax2.set_ylabel("EIRP [dBW]")
    ax2.set_xlabel("Time")
    #ax1.set_ylim([0,20])
    ax2.grid(True)


    for i,x in enumerate(cn):
        ax3.plot(t,x,label=labels[i])
        df =  pandas.concat([t, x],axis=1)
        df.to_csv(data_output_directory + "snr" + labels[i] + ".txt", index=None, sep="," ,float_format='%.2f')

    ax3.set_ylabel("C/N")
    ax3.set_xlabel("Time")
    #ax2.set_ylim([0,60])
    ax3.grid(True)


    ax3.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax3.xaxis.set_major_locator(mdates.MinuteLocator())
    plt.gcf().autofmt_xdate()

    plt.legend(loc='lower left', bbox_to_anchor=(1, 0.5))



    plt.show()

