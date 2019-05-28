#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Hannes Bartle <mail@hbartle.de>
#
# Distributed under terms of the MIT license.


import matplotlib.pyplot as plt

# Own Modules
from antenna_toolbox.Antenna import Antenna
from antenna_toolbox.RectangularArray import RectangularArray

if __name__ == "__main__":



    f = 8e9
    patch = Antenna()
    patch.loadPatternFromFile("sim/patch_8Ghz_efield_ro5880_0.254.txt")

    paa1 = RectangularArray(4,4,0.5,0.5,f,arraytype="uniform",sidelobe_level=30)
    paa1.setElement(patch)
    paa1.getArrayFactor(0,180,1,0,360,1,0,0)
    paa1.getField(0,0)
    paa1.calcDirectivity()

    paa1_sim = Antenna()
    paa1_sim.loadPatternFromFile("sim/uniform/farfield_array_4x4_0-0.txt")
    paa1_sim.calcDirectivity()


    paa1.plotDirectivity3D(logarithmic=True,log_range=-40)
    paa1_sim.plotDirectivity3D(logarithmic=True,log_range=-40)

    plt.figure()
    plr = True
    ax = plt.subplot(111,polar=plr)
    paa1.plotDirectivityCut("phi",[0],ax=ax,polar=plr,labels=["Sim + Analytical"])
    paa1_sim.plotDirectivityCut("phi",[0],ax=ax,polar=plr,labels=["Sim"])

 #   f = 8e9
 #   e_r = 3.66
 #   h = 0.101e-3
 #   patch = RectangularPatch(e_r,h,f)
 #   patch.getField(0,90,1,0,360,1)
 #
 #   # Generate Arrays
 #   dims = [2,3,4,5]
 #   arrays = []
 #   for dim in dims:
 #       paa = RectangularArray(dim,dim,0.5,0.5,f,arraytype="uniform")
 #       paa.setElement(patch)
 #       paa.getField(0,0)
 #       paa.calcDirectivity()
 #       arrays.append(paa)
 #
 #
 #   ## Plotting
 #   plr = True
 #   ax = plt.subplot(111,polar=plr)
 #
 #   for array in arrays:
 #       array.plotDirectivityCut("phi", [0],ax=ax,polar=plr,labels=["Uniform, "+str(array.Nelements_x) + "x" + str(array.Nelements_y)])
 #       #paa.plotPattern3D(logarithmic=False,log_range=-40)
 #   arrays[0].plotDirectivity3D()
    plt.show()


