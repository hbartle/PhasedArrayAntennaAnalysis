#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Hannes Bartle <mail@hbartle.de>
#
# Distributed under terms of the MIT license.

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.dates as mdates
#plt.style.use('classic')
#plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rcParams.update({'font.size': 14})
import matplotlib.gridspec as gridspec
cmap = plt.get_cmap('jet')
import pandas
import datetime



"""
    Plots for Sat2Sat scenario
"""


# Load Data
directory = "C:/Users/hanne/Nextcloud/studies/au/courses/master_thesis/stk_analysis_results/sat2sat/"
filename = "Satellite-LeadingSat-Transmitter-Fixed2x2-To-Satellite-TrailingSat-Receiver-Fixed2x2 AER.csv"

d = pandas.read_csv(directory + filename,
                    delimiter=',',
                    parse_dates=["Time (UTCG)"],
                    date_parser=lambda dates: [pandas.datetime.strptime(d, '%d %b %Y %H:%M:%S.%f') for d in dates])
r = d["Range (km)"]
if __name__ == "__main__":

    directory = "C:/Users/hanne/Nextcloud/studies/au/courses/master_thesis/stk_analysis_results/sat2sat/comparison_fixed_paa/"


    filenames_fixed = ['Satellite-LeadingSat-Transmitter-Fixed2x2-To-Satellite-TrailingSat-Receiver-Fixed2x2 Link Budget - Detailed.csv',
'Satellite-LeadingSat-Transmitter-Fixed3x3-To-Satellite-TrailingSat-Receiver-Fixed3x3 Link Budget - Detailed.csv',
'Satellite-LeadingSat-Transmitter-Fixed4x4-To-Satellite-TrailingSat-Receiver-Fixed4x4 Link Budget - Detailed.csv',
'Satellite-LeadingSat-Transmitter-Fixed5x5-To-Satellite-TrailingSat-Receiver-Fixed5x5 Link Budget - Detailed.csv']
    filenames_paa = ['Satellite-LeadingSat-Transmitter-PAA2x2-To-Satellite-TrailingSat-Receiver-PAA2x2 Link Budget - Detailed.csv',
'Satellite-LeadingSat-Transmitter-PAA3x3-To-Satellite-TrailingSat-Receiver-PAA3x3 Link Budget - Detailed.csv',
'Satellite-LeadingSat-Transmitter-PAA4x4-To-Satellite-TrailingSat-Receiver-PAA4x4 Link Budget - Detailed.csv',
'Satellite-LeadingSat-Transmitter-PAA5x5-To-Satellite-TrailingSat-Receiver-PAA5x5 Link Budget - Detailed.csv']

    labels = ["2x2","3x3","4x4","5x5"]


    gs = gridspec.GridSpec(2, 1)
    plt.figure()
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1],sharex=ax0)


    t=[]
    snr_fixed = []
    eirp_fixed = []
    snr_paa = []
    eirp_paa = []
    for i,label in enumerate(labels):
        color = cmap(float(i) / 4)

        d_fixed = pandas.read_csv(directory + filenames_fixed[i],
                            delimiter=',',
                            parse_dates=["Time (UTCG)"],
                            date_parser=lambda dates: [pandas.datetime.strptime(d, '%d %b %Y %H:%M:%S.%f') for d in
                                                       dates])

        t = d_fixed["Time (UTCG)"]
        snr_fixed.append(d_fixed["C/N (dB)"])
        eirp_fixed.append(d_fixed["EIRP (dBW)"])

        ax0.plot(r,eirp_fixed[i],label=label,c=color)
        ax1.plot(r,snr_fixed[i],label=label,c=color)



        d_paa = pandas.read_csv(directory + filenames_paa[i],
                                  delimiter=',',
                                  parse_dates=["Time (UTCG)"],
                                  date_parser=lambda dates: [pandas.datetime.strptime(d, '%d %b %Y %H:%M:%S.%f') for d
                                                             in
                                                             dates])

        t = d_paa["Time (UTCG)"]
        snr_paa.append(d_paa["C/N (dB)"])
        eirp_paa.append(d_paa["EIRP (dBW)"])

        ax0.plot(r,eirp_paa[i],"--",label=label,c=color)
        ax1.plot(r,snr_paa[i],"--",label=label,c=color)



    ax1.set_ylabel("SNR [dB]")
    ax0.set_ylabel("EIRP [dBW]")

    ax0.grid(True)
    ax1.grid(True)
    plt.legend(loc='upper right')


    dic = {"Time": t,
           "Range": r,
           "SNR Fixed2x2": snr_fixed[0],
           "SNR PAA2x2": snr_paa[0],
           "EIRP Fixed2x2": eirp_fixed[0],
           "EIRP PAA2x2": eirp_paa[0],
           "SNR Fixed3x3": snr_fixed[1],
           "SNR PAA3x3": snr_paa[1],
           "EIRP Fixed3x3": eirp_fixed[1],
           "EIRP PAA3x3": eirp_paa[1],
           "SNR Fixed4x4": snr_fixed[2],
           "SNR PAA4x4": snr_paa[2],
           "EIRP Fixed4x4": eirp_fixed[2],
           "EIRP PAA4x4": eirp_paa[2],
           "SNR Fixed5x5": snr_fixed[3],
           "SNR PAA5x5": snr_paa[3],
           "EIRP Fixed5x5": eirp_fixed[3],
           "EIRP PAA5x5": eirp_paa[3]}

    df =  pandas.DataFrame(data=dic)
    df.to_csv(directory + "comparison_fixed_paa.txt", index=None, sep=",", float_format='%.2f')

    plt.show()
