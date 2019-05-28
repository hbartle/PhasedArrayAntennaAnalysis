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
    Plot Azimuth, Elevation and Range for the Sat2Sat scenario
"""


if __name__ == "__main__":

    # Load Data
    directory= "C:/Users/hanne/Nextcloud/studies/au/courses/master_thesis/stk_analysis_results/sat2sat/"
    filename = "Satellite-LeadingSat-Transmitter-Fixed2x2-To-Satellite-TrailingSat-Receiver-Fixed2x2 AER.csv"

    d = pandas.read_csv(directory + filename,
                        delimiter=',',
                        parse_dates=["Time (UTCG)"],
                        date_parser=lambda dates: [pandas.datetime.strptime(d, '%d %b %Y %H:%M:%S.%f') for d in dates])


    start = 0
    t = d["Time (UTCG)"][start:]
    r = d["Range (km)"][start:]
    el = d["Elevation (deg)"][start:]

    # Plot
    gs = gridspec.GridSpec(2, 1)
    plt.figure()
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1],sharex=ax0)

    ax0.plot(t, r)
    ax0.set_ylabel("Range [km]")
    ax0.grid(True)

    ax1.plot(t, el+90)
    ax1.set_ylabel("Elevation [deg]")
    ax1.grid(True)

    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%d %b %Y '))
    ax1.xaxis.set_major_locator(mdates.DayLocator())
    plt.gcf().autofmt_xdate()

    d["Elevation (deg)"] = d["Elevation (deg)"] + 90
    d.to_csv(directory + "AER_sat2sat.txt", index=None, sep=",", float_format='%.2f')


    plt.show()

