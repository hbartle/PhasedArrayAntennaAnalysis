#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Hannes Bartle <mail@hbartle.de>
#
# Distributed under terms of the MIT license.

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
plt.style.use('seaborn-whitegrid')
# plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 14})

# Own Modules
from antenna_toolbox.Antenna import Antenna
from antenna_toolbox.RectangularPatch import RectangularPatch
from util import *
import pandas

if __name__ == "__main__":

    f = 8e9
    e_r = 3.48
    h = 1.55e-3

    ### Analytical Patch ###
    patch = RectangularPatch(e_r, h, f, rolloff=0.1)
    patch.getField(0, 90, 1, 0, 360, 1)
    patch.calcDirectivity()


    gs = gridspec.GridSpec(2, 2)
    #ax1 = plt.subplot(gs[0],polar=True)
    #ax2 = plt.subplot(gs[1],projection='3d')
    #ax3 = plt.subplot(gs[2],polar=True)
    #ax4 = plt.subplot(gs[3],projection='3d')

    fig = plt.figure()
    ax1 = fig.add_subplot(111, polar=True)
    fig.suptitle('Directivity, Phi= 0°')


    fig = plt.figure()
    ax2 = fig.add_subplot(111, projection='3d')
    fig.suptitle('Normalized Pattern, Analytical')

    fig = plt.figure()
    ax3 = fig.add_subplot(111, polar=True)
    fig.suptitle('Directivity, Phi= 90°')

    fig = plt.figure()
    ax4 = fig.add_subplot(111, projection='3d')
    fig.suptitle('Normalized Pattern, Simulated')

    #patch.plotPatternCut("phi", [0,90],ax=ax1,labels=["0Deg","90Deg"])
    patch.plotPattern3D(ax=ax2,logarithmic=True)


    ### Simulated Patch ###

    patch_sim = Antenna()
    patch_sim.loadPatternFromFile("../sim/patch_8Ghz_efield_ro5880_0.254.txt")
    patch_sim.calcDirectivity()


    #patch_sim.plotPatternCut("phi", [0,90],ax=ax3,labels=["0Deg","90Deg"])
    patch_sim.plotPattern3D(ax=ax4,logarithmic=True)

    _,dir_phi0,theta_phi0 = patch.plotPatternCut("phi", [0],ax=ax1,labels=["Analytical"])
    _, dir_sim_phi0, theta_sim_phi0 = patch_sim.plotPatternCut("phi", [0],ax=ax1,labels=["Simulated"])
    _, dir_phi90, theta_phi90 = patch.plotPatternCut("phi", [90],ax=ax3,labels=["Analytical"])
    _, dir_sim_phi90, theta_phi90 = patch_sim.plotPatternCut("phi", [90],ax=ax3,labels=["Simulated"])


    dic = {"Theta": rad2deg(theta_phi0),
           "Dir Phi0": lin2db(dir_phi0).clip(min=-40),
           "Dir Phi90": lin2db(dir_phi90).clip(min=-40)}
    dic_sim = {"Theta":rad2deg(theta_sim_phi0),
                "Dir Sim Phi0": lin2db(dir_sim_phi0).clip(min=-40),
                "Dir Sim Phi90": lin2db(dir_sim_phi90).clip(min=-40)}

    df = pandas.DataFrame(data=dic)
    df.to_csv("../sim/directivity_patch_analytical.txt", index=None, sep=",", float_format='%.2f')
    df_sim = pandas.DataFrame(data=dic_sim)
    df_sim.to_csv("../sim/directivity_patch_sim.txt", index=None, sep=",", float_format='%.2f')

    #plt.show()

