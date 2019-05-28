#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Hannes Bartle <mail@hbartle.de>
#
# Distributed under terms of the MIT license.



import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

plt.style.use('classic')
#plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 16})

import matplotlib.gridspec as gridspec
from scipy.interpolate import interp2d
import scipy.constants as const
from util import *


if __name__=="__main__":



    number_of_elements_max = np.arange(2,20,1)
    array_length_max = 0.1 # See Gomspace ANT2000 (S-band Patch) for reference

    max_directivity = lin2db(number_of_elements_max)

    freq_max = const.speed_of_light/2 * number_of_elements_max / array_length_max


    # Plot
    ax0 = plt.gca()
    ax1 = ax0.twinx()
    ax0.plot(number_of_elements_max,freq_max/1e9,"ko--")
    ax0.set_ylabel("Frequency [GHz]",color="k")
    ax0.set_xlabel("Array Size")
    ax0.set_xlim([number_of_elements_max.min(),number_of_elements_max.max() ])
    ax1.plot(number_of_elements_max,max_directivity,"mo--")
    ax1.set_ylabel("Approx. Increase in Directivity [dB]",color="m")

    labels = [str(n) + "x" + str(n) for n in number_of_elements_max]
    ax0.set_xticks(number_of_elements_max)
    ax0.set_xticklabels(labels)
    #ax0.yaxis.set_major_locator(ticker.MultipleLocator(3))
    labels = np.arange(3,36,3)
    ax0.set_yticks(labels[:-1])
    ax0.set_yticklabels(labels[1:])

    np.set_printoptions(precision=2)
    for i,n in enumerate(number_of_elements_max):
        print("(" + str(n)+"," + str(max_directivity[i]) + ")")

    labels = np.arange(6,33,3)
    ax1.set_yticks(labels[:-1])
    ax1.set_yticklabels(labels[1:])

    ax0.grid(True)
    ax1.grid(0)

    plt.show()
