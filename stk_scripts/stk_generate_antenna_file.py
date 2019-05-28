#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Hannes Bartle <mail@hbartle.de>
#
# Distributed under terms of the MIT license.


import sys
import math
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
plt.style.use('seaborn-whitegrid')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 14})
from scipy import interpolate
from scipy.signal import argrelextrema

from CSTResultLoader import *
from util import *


"""
    Convert the CST simulation results to a STK antenna file
"""

if __name__ == "__main__":
    result_loader = CSTResultLoader()

    directory = "../../sim_data/reference_antennas/"
    filenames = ["2x2_uniform_broadside_realized_gain.txt",
                 "3x3_uniform_broadside_realized_gain.txt",
                 "4x4_uniform_broadside_realized_gain.txt",
                 "5x5_uniform_broadside_realized_gain.txt"]

    # Create STK antenna file for reference antennas
    for filename in filenames:
        result_loader.convertASCII2STKAntennaFile(directory,filename)


