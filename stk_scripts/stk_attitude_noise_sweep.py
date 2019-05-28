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
from util import *


"""
    Parse and plot STK results, Attitude Noise Sweep
"""

if __name__ == "__main__":

    output_directory = "C:/Users/hanne/Nextcloud/studies/au/courses/master_thesis/stk_analysis_results/sat2sat_noise/"
    directories = ["C:/Users/hanne/Nextcloud/studies/au/courses/master_thesis/stk_analysis_results/sat2sat_noise/sigma1/",
                   "C:/Users/hanne/Nextcloud/studies/au/courses/master_thesis/stk_analysis_results/sat2sat_noise/sigma2/",
                   "C:/Users/hanne/Nextcloud/studies/au/courses/master_thesis/stk_analysis_results/sat2sat_noise/sigma3/",
                   "C:/Users/hanne/Nextcloud/studies/au/courses/master_thesis/stk_analysis_results/sat2sat_noise/sigma4/",
                   "C:/Users/hanne/Nextcloud/studies/au/courses/master_thesis/stk_analysis_results/sat2sat_noise/sigma5/",
                   "C:/Users/hanne/Nextcloud/studies/au/courses/master_thesis/stk_analysis_results/sat2sat_noise/sigma6/",
                   "C:/Users/hanne/Nextcloud/studies/au/courses/master_thesis/stk_analysis_results/sat2sat_noise/sigma7/",
                   "C:/Users/hanne/Nextcloud/studies/au/courses/master_thesis/stk_analysis_results/sat2sat_noise/sigma8/",
                   "C:/Users/hanne/Nextcloud/studies/au/courses/master_thesis/stk_analysis_results/sat2sat_noise/sigma9/",
                   "C:/Users/hanne/Nextcloud/studies/au/courses/master_thesis/stk_analysis_results/sat2sat_noise/sigma10/"]



    filenames_fixed = ['Satellite-LeadingSat-Transmitter-Fixed2x2-To-Satellite-TrailingSat-Receiver-Fixed2x2 NoisyLinkBudgetAnalysis.csv',
                'Satellite-LeadingSat-Transmitter-Fixed3x3-To-Satellite-TrailingSat-Receiver-Fixed3x3 NoisyLinkBudgetAnalysis.csv',
                'Satellite-LeadingSat-Transmitter-Fixed4x4-To-Satellite-TrailingSat-Receiver-Fixed4x4 NoisyLinkBudgetAnalysis.csv',
                'Satellite-LeadingSat-Transmitter-Fixed5x5-To-Satellite-TrailingSat-Receiver-Fixed5x5 NoisyLinkBudgetAnalysis.csv']

    filenames_paa=[
                'Satellite-LeadingSat-Transmitter-PAA2x2-To-Satellite-TrailingSat-Receiver-PAA2x2 NoisyLinkBudgetAnalysis.csv',
                'Satellite-LeadingSat-Transmitter-PAA3x3-To-Satellite-TrailingSat-Receiver-PAA3x3 NoisyLinkBudgetAnalysis.csv',
                'Satellite-LeadingSat-Transmitter-PAA4x4-To-Satellite-TrailingSat-Receiver-PAA4x4 NoisyLinkBudgetAnalysis.csv',
                'Satellite-LeadingSat-Transmitter-PAA5x5-To-Satellite-TrailingSat-Receiver-PAA5x5 NoisyLinkBudgetAnalysis.csv']

    sigma = [1,2,3,4,5,6,7,8,9,10]
    labels = ["2x2","3x3","4x4","5x5"]
    quantity = "Xmtr Gain (dB)"
    #quantity = "Rcvr Gain (dB)"
    #quantity = "C/N (dB)"
    #quantity = "log(BER)"

    d22_fixed = []
    d22_paa = []
    d33_fixed = []
    d33_paa = []
    d44_fixed = []
    d44_paa = []
    d55_fixed = []
    d55_paa = []

    for directory in directories:

        d22_fixed.append(pandas.read_csv(directory + filenames_fixed[0],
                                delimiter=',',
                                parse_dates=["Time (UTCG)"],
                                date_parser=lambda dates: [pandas.datetime.strptime(d, '%d %b %Y %H:%M:%S.%f')
                                                           for d in dates]))

        d22_paa.append(pandas.read_csv(directory + filenames_paa[0],
                                 delimiter=',',
                                 parse_dates=["Time (UTCG)"],
                                 date_parser=lambda dates: [pandas.datetime.strptime(d, '%d %b %Y %H:%M:%S.%f')
                                                            for d in dates]))
        d33_fixed.append(pandas.read_csv(directory + filenames_fixed[1],
                                         delimiter=',',
                                         parse_dates=["Time (UTCG)"],
                                         date_parser=lambda dates: [pandas.datetime.strptime(d, '%d %b %Y %H:%M:%S.%f')
                                                                    for d in dates]))

        d33_paa.append(pandas.read_csv(directory + filenames_paa[1],
                                       delimiter=',',
                                       parse_dates=["Time (UTCG)"],
                                       date_parser=lambda dates: [pandas.datetime.strptime(d, '%d %b %Y %H:%M:%S.%f')
                                                                  for d in dates]))
        d44_fixed.append(pandas.read_csv(directory + filenames_fixed[2],
                                         delimiter=',',
                                         parse_dates=["Time (UTCG)"],
                                         date_parser=lambda dates: [pandas.datetime.strptime(d, '%d %b %Y %H:%M:%S.%f')
                                                                    for d in dates]))

        d44_paa.append(pandas.read_csv(directory + filenames_paa[2],
                                       delimiter=',',
                                       parse_dates=["Time (UTCG)"],
                                       date_parser=lambda dates: [pandas.datetime.strptime(d, '%d %b %Y %H:%M:%S.%f')
                                                                  for d in dates]))
        d55_fixed.append(pandas.read_csv(directory + filenames_fixed[3],
                                         delimiter=',',
                                         parse_dates=["Time (UTCG)"],
                                         date_parser=lambda dates: [pandas.datetime.strptime(d, '%d %b %Y %H:%M:%S.%f')
                                                                    for d in dates]))

        d55_paa.append(pandas.read_csv(directory + filenames_paa[3],
                                       delimiter=',',
                                       parse_dates=["Time (UTCG)"],
                                       date_parser=lambda dates: [pandas.datetime.strptime(d, '%d %b %Y %H:%M:%S.%f')
                                                                  for d in dates]))

    ## Plot Dedicated Array Over Time
    for idx,s in enumerate(sigma):
        fig = plt.figure()
        fig.suptitle("Sigma=" + str(s))
        ax0 = plt.subplot()


        end = 400
        ax0.plot(d55_fixed[idx]["Time (UTCG)"][:end],d55_fixed[idx][quantity][:end],label=labels[3])
        ax0.plot(d44_fixed[idx]["Time (UTCG)"][:end],d44_fixed[idx][quantity][:end],label=labels[2])
        ax0.plot(d33_fixed[idx]["Time (UTCG)"][:end],d33_fixed[idx][quantity][:end],label=labels[1])
        ax0.plot(d22_fixed[idx]["Time (UTCG)"][:end],d22_fixed[idx][quantity][:end],label=labels[0])


        ax0.set_ylim([0,20])
        ax0.set_ylabel(quantity)
        ax0.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
        ax0.xaxis.set_major_locator(mdates.AutoDateLocator())
        plt.gcf().autofmt_xdate()

        dic = {"Time": d55_fixed[idx]["Time (UTCG)"][:end],
               "Xmtr Gain 5x5 fixed": d55_fixed[idx][quantity][:end],
               "Xmtr Gain 4x4 fixed": d44_fixed[idx][quantity][:end],
               "Xmtr Gain 3x3 fixed": d33_fixed[idx][quantity][:end],
               "Xmtr Gain 2x2 fixed": d22_fixed[idx][quantity][:end],
               "Xmtr Gain 5x5 PAA": d55_paa[idx][quantity][:end],
               "Xmtr Gain 4x4 PAA": d44_paa[idx][quantity][:end],
               "Xmtr Gain 3x3 PAA": d33_paa[idx][quantity][:end],
               "Xmtr Gain 2x2 PAA": d22_paa[idx][quantity][:end]}

        df = pandas.DataFrame(data=dic)
        df.to_csv(output_directory + "xmtr_gain_sigma=" + str(s) + ".txt", index=None, sep=",", float_format='%.2f')



    ## Compare Array for fixed Sigma
    end = 400
    idx = 3

    fig = plt.figure()
    fig.suptitle("Comparison Fixed Array vs PAA 5x5")
    ax0 = plt.subplot()
    ax0.plot(d55_fixed[idx]["Time (UTCG)"][:end],d55_fixed[idx][quantity][:end],label=labels[3])
    ax0.plot(d55_paa[idx]["Time (UTCG)"][:end],d55_paa[idx][quantity][:end],label=labels[3])
    ax0.set_ylim([0,20])
    ax0.set_ylabel(quantity)
    ax0.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax0.xaxis.set_major_locator(mdates.AutoDateLocator())
    plt.gcf().autofmt_xdate()

    fig = plt.figure()
    fig.suptitle("Comparison Fixed Array vs PAA 4x4")
    ax0 = plt.subplot()
    ax0.plot(d44_fixed[idx]["Time (UTCG)"][:end],d44_fixed[idx][quantity][:end],label=labels[2])
    ax0.plot(d44_paa[idx]["Time (UTCG)"][:end],d44_paa[idx][quantity][:end],label=labels[2])
    ax0.set_ylim([0,20])
    ax0.set_ylabel(quantity)
    ax0.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax0.xaxis.set_major_locator(mdates.AutoDateLocator())
    plt.gcf().autofmt_xdate()

    fig = plt.figure()
    fig.suptitle("Comparison Fixed Array vs PAA 3x3")
    ax0 = plt.subplot()
    ax0.plot(d33_fixed[idx]["Time (UTCG)"][:end],d33_fixed[idx][quantity][:end],label=labels[1])
    ax0.plot(d33_paa[idx]["Time (UTCG)"][:end],d33_paa[idx][quantity][:end],label=labels[1])
    ax0.set_ylim([0,20])
    ax0.set_ylabel(quantity)
    ax0.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax0.xaxis.set_major_locator(mdates.AutoDateLocator())
    plt.gcf().autofmt_xdate()

    fig = plt.figure()
    fig.suptitle("Comparison Fixed Array vs PAA 2x2")
    ax0 = plt.subplot()
    ax0.plot(d22_fixed[idx]["Time (UTCG)"][:end],d22_fixed[idx][quantity][:end],label=labels[0])
    ax0.plot(d22_paa[idx]["Time (UTCG)"][:end],d22_paa[idx][quantity][:end],label=labels[0])
    ax0.set_ylim([0,20])
    ax0.set_ylabel(quantity)
    ax0.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax0.xaxis.set_major_locator(mdates.AutoDateLocator())
    plt.gcf().autofmt_xdate()



    ### Plot Averages
    d22_fixed_avg = []
    d22_paa_avg = []
    d33_fixed_avg = []
    d33_paa_avg = []
    d44_fixed_avg = []
    d44_paa_avg = []
    d55_fixed_avg = []
    d55_paa_avg = []
    for i,s in enumerate(sigma):
        d22_fixed_avg.append(pow2db(db2pow(d22_fixed[i][quantity]).mean()))
        d22_paa_avg.append(pow2db(db2pow(d22_paa[i][quantity]).mean()))

        d33_fixed_avg.append(pow2db(db2pow(d33_fixed[i][quantity]).mean()))
        d33_paa_avg.append(pow2db(db2pow(d33_paa[i][quantity]).mean()))

        d44_fixed_avg.append(pow2db(db2pow(d44_fixed[i][quantity]).mean()))
        d44_paa_avg.append(pow2db(db2pow(d44_paa[i][quantity]).mean()))

        d55_fixed_avg.append(pow2db(db2pow(d55_fixed[i][quantity]).mean()))
        d55_paa_avg.append(pow2db(db2pow(d55_paa[i][quantity]).mean()))

    fig = plt.figure()
    fig.suptitle("Average " + quantity)
    ax0 = plt.subplot()
    ax0.plot(sigma,d22_fixed_avg,label=labels[0] + ",fixed")
    ax0.plot(sigma,d22_paa_avg,label=labels[0] + ",PAA")
    ax0.plot(sigma,d33_fixed_avg,label=labels[1] + ",fixed")
    ax0.plot(sigma,d33_paa_avg,label=labels[1] + ",PAA")
    ax0.plot(sigma,d44_fixed_avg,label=labels[2] + ",fixed")
    ax0.plot(sigma,d44_paa_avg,label=labels[2] + ",PAA")
    ax0.plot(sigma,d55_fixed_avg,label=labels[3] + ",fixed")
    ax0.plot(sigma,d55_paa_avg,label=labels[3] + ",PAA")
    ax0.set_ylabel(quantity)
    ax0.grid(True)

    dic = {"Std Deviation": sigma,
           "Avg Xmtr Gain 5x5 fixed": d55_fixed_avg,
           "Avg Xmtr Gain 5x5 PAA": d55_paa_avg,
           "Avg Xmtr Gain 4x4 fixed": d44_fixed_avg,
           "Avg Xmtr Gain 4x4 PAA": d44_paa_avg,
           "Avg Xmtr Gain 3x3 fixed": d33_fixed_avg,
           "Avg Xmtr Gain 3x3 PAA": d33_paa_avg,
           "Avg Xmtr Gain 2x2 fixed": d22_fixed_avg,
           "Avg Xmtr Gain 2x2 PAA": d22_paa_avg,
           }

    df = pandas.DataFrame(data=dic)
    df.to_csv(output_directory + "avg_xmtr_gain.txt", index=None, sep=",", float_format='%.2f')

    plt.legend()

    plt.show()
