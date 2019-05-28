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
plt.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 16})
plt.rc('grid', linestyle="--")

import matplotlib.gridspec as gridspec
from scipy.interpolate import interp2d
import scipy.constants as const
from util import *


if __name__=="__main__":

    directory = "../sim_data/element_tuning_5x5/"
    files = [directory + "s11_parameters_5x5_patch_untuned.txt",
             directory + "s_parameter_2x2_array_tuned.txt",
             directory + "s_parameter_3x3_array_tuned.txt",
             directory + "s_parameter_4x4_array_tuned.txt",
             directory + "s_parameter_5x5_array_tuned.txt"]

    number_of_elements = [25, 4, 9, 16, 25]

    for i,file in enumerate(files):
        plot_name = file.split("/")[3]
        plot_name = plot_name[:-4]
        with open(file, newline='') as f:
            s_param = []
            start = 0
            for line in f:
                if "----------------------------------------------------------------------" in line:
                    start = 1
                elif line in ['\n', '\r\n']:
                    start = 0
                elif start:
                    s_param.append(list(map(float, line.split())))
        s_param= np.asarray(s_param)
        s_param = np.split(s_param, number_of_elements[i])
        freq = s_param[0][:,0]

        # Plot
        fig = plt.figure()
        fig.canvas.set_window_title(plot_name)
        ax0 = plt.subplot()
        for i,s in enumerate(s_param):
            ax0.plot(freq,s[:,1], "-",linewidth=0.5)

        ax0.set_xlim([6,10])
        ax0.set_ylim([-35,0])
        ax0.set_yticks(np.arange(-30,5,5))
        ax0.set_xlabel("Frequency [GHz]")
        ax0.set_ylabel("S11 [dB]")

        plt.grid(True)
        plt.savefig(directory + plot_name + ".png")

    plt.show()
