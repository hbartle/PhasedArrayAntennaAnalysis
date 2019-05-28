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

from util import *

directories = ["C:/Users/hanne/Nextcloud/studies/au/courses/master_thesis/sim_data/quantization_analysis/5x5",
            "C:/Users/hanne/Nextcloud/studies/au/courses/master_thesis/sim_data/quantization_analysis/4x4",
            "C:/Users/hanne/Nextcloud/studies/au/courses/master_thesis/sim_data/quantization_analysis/3x3",
            "C:/Users/hanne/Nextcloud/studies/au/courses/master_thesis/sim_data/quantization_analysis/2x2"]
output_directory =  "C:/Users/hanne/Nextcloud/studies/au/courses/master_thesis/sim_data/quantization_analysis/"

cuts_file = 'Directivity,Phi=0.txt'
hpbw_file = 'Realized Gain,Phi=PAA_FA_SCANPHI, HPBW.txt'
mainlobe_file = 'Realized Gain,Phi=PAA_FA_SCANPHI, MainLobeDirection.txt'
maxvalue_file = 'Realized Gain,Phi=PAA_FA_SCANPHI, MaxValue.txt'
sll_file = 'Realized Gain,Phi=PAA_FA_SCANPHI, SLL.txt'
value_file = 'Realized Gain,Theta=PAA_FA_SCANTHETA,Phi=PAA_FA_SCANPHI.txt'
resnav_file =  "result_navigator.csv"

arraysize = ["5x5","4x4","3x3","2x2"]

for k,directory in enumerate(directories):
    ### Cuts ###
    with open(directory + "_analog/" + cuts_file, newline='') as f:
        res = []
        start = 0
        for line in f:

            if "----------------------------------------------------------------------" in line:

                start = 1
            elif line in ['\n', '\r\n']:
                start = 0
            elif start:
                res.append(list(map(float, line.split())))
    results_analog = np.asarray(res)
    results_analog = np.split(results_analog, 11)



    with open(directory + "_digital/" + cuts_file, newline='') as f:
        res = []
        start = 0
        for line in f:

            if "----------------------------------------------------------------------" in line:

                start = 1
            elif line in ['\n', '\r\n']:
                start = 0
            elif start:
                res.append(list(map(float, line.split())))
    results = np.asarray(res)
    results = np.split(results, 66)


    result_navigator = np.genfromtxt(directory + "_digital/" + resnav_file,skip_header=1,delimiter=",")


    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(results_analog[0][:, 0], results_analog[0][:, 1],label="Theta=45" + ",Analog")

    indices = np.arange(0,6,1)
    for idx in indices:
        ax.plot(results[idx][:,0],results[idx][:,1],label="Theta=" + str(result_navigator[idx,1])+ ",Bits=" + str(result_navigator[idx,2]))

    plt.legend()

    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(results_analog[10][:, 0], results_analog[10][:, 1],label="Theta=90" + ",Analog")

    indices = np.arange(15,66,10)
    for idx in indices:
        ax.plot(results[idx][:,0],results[idx][:,1],label="Theta=" + str(result_navigator[idx,1])+ ",Bits=" + str(result_navigator[idx,2]))


    plt.legend()

    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(results_analog[3][:, 0], results_analog[3][:, 1],label="Theta=20" + ",Analog")

    indices = np.arange(8,66,10)
    for idx in indices:
        ax.plot(results[idx][:,0],results[idx][:,1],label="Theta=" + str(result_navigator[idx,1])+ ",Bits=" + str(result_navigator[idx,2]))


    plt.legend()
    #plt.show()

    dic = {"Theta": results[8][:,0],
           "Analog": results_analog[3][:, 1],
           "Bit1": results[8][:,1],
           "Bit2": results[18][:, 1],
           "Bit3": results[28][:, 1],
           "Bit4": results[38][:, 1],
           "Bit5": results[48][:, 1],
           "Bit6": results[58][:, 1]}

    df = pandas.DataFrame(data=dic)
    df.to_csv(output_directory + "quantization_cuts_theta=20," + arraysize[k] + ".txt",
              index=None, sep=",", float_format='%.2f')

    dic = {"Theta": results[11][:,0],
           "Analog": results_analog[6][:, 1],
           "Bit1": results[11][:,1],
           "Bit2": results[21][:, 1],
           "Bit3": results[31][:, 1],
           "Bit4": results[41][:, 1],
           "Bit5": results[51][:, 1],
           "Bit6": results[61][:, 1]}

    df = pandas.DataFrame(data=dic)
    df.to_csv(output_directory + "quantization_cuts_theta=50," + arraysize[k] + ".txt",
              index=None, sep=",", float_format='%.2f')

    dic = {"Theta": results[15][:,0],
           "Analog": results_analog[10][:, 1],
           "Bit1": results[15][:,1],
           "Bit2": results[25][:, 1],
           "Bit3": results[35][:, 1],
           "Bit4": results[45][:, 1],
           "Bit5": results[55][:, 1],
           "Bit6": results[65][:, 1]}

    df = pandas.DataFrame(data=dic)
    df.to_csv(output_directory + "quantization_cuts_theta=90," + arraysize[k] + ".txt",
              index=None, sep=",", float_format='%.2f')

    ### Single Value Plots ###
    file = directory + "_digital/" + value_file
    value = np.genfromtxt(file,skip_header=2)
    value_bits = np.split(value[6:,1],6)

    file_analog = directory + "_analog/" + value_file
    value_analog = np.genfromtxt(file_analog,skip_header=2)
    mask = np.ones(value_analog.shape[0],dtype=bool)
    mask[5] = 0
    value_analog = value_analog[mask]


    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(value_analog[:,1],label="Analog")
    for i in range(0,6):
        ax.plot(value_bits[i],label="Bits=" + str(i+1))
    plt.legend()


    file = directory + "_digital/" + mainlobe_file
    mainlobe_direction = np.genfromtxt(file,skip_header=2)
    mainlobe_direction_bits = np.split(mainlobe_direction[6:,1],6)

    file_analog = directory + "_analog/" + mainlobe_file
    mainlobe_analog = np.genfromtxt(file_analog,skip_header=2)
    mask = np.ones(mainlobe_analog.shape[0],dtype=bool)
    mask[5] = 0
    mainlobe_analog = mainlobe_analog[mask]



    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(mainlobe_analog[:,1] - np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90]), label="Analog")

    for i in range(0,6):
        ax.plot(mainlobe_direction_bits[i] - np.array([0, 10 ,20 ,30 ,40 ,50 ,60 ,70 ,80,90]),label="Bits=" + str(i+1))
    plt.legend()

    dic = {"Theta": np.arange(0,100,10),
           "Analog": value_analog[:, 1],
           "Bit1": value_bits[0],
           "Bit2": value_bits[1],
           "Bit3": value_bits[2],
           "Bit4": value_bits[3],
           "Bit5": value_bits[4],
           "Bit6": value_bits[5]}

    df = pandas.DataFrame(data=dic)
    df.to_csv(output_directory + "realized_gain_in_steering_direction," + arraysize[k] + ".txt",
              index=None, sep=",", float_format='%.2f')


    theta = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90])
    dic = {"Theta": theta,
           "Analog": mainlobe_analog[:,1] - theta,
           "Bit1": mainlobe_direction_bits[0]  - theta,
           "Bit2": mainlobe_direction_bits[1]  - theta,
           "Bit3": mainlobe_direction_bits[2] - theta,
           "Bit4": mainlobe_direction_bits[3] - theta,
           "Bit5": mainlobe_direction_bits[4] - theta,
           "Bit6": mainlobe_direction_bits[5] - theta}

    df = pandas.DataFrame(data=dic)
    df.to_csv(output_directory + "steering_offset," + arraysize[k] + ".txt",
              index=None, sep=",", float_format='%.2f')




plt.show()

