#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Hannes Bartle <mail@hbartle.de>
#
# Distributed under terms of the MIT license.

import matplotlib.pyplot as plt

# Own Modules
from antenna_toolbox.RectangularArray import RectangularArray

if __name__ == "__main__":

    f = 8e9
    sdlb = 30

#    dims = [3,3,4,5]
#    arrays = []
#    for dim in dims:
#        paa = RectangularArray(dim,dim,0.5,0.5,f,arraytype="uniform")
#        paa.getArrayFactor(0,180,1,0,360,1,0,0)
#        arrays.append(paa)
#
#
#    plr = True
#    ax = plt.subplot(111,polar=plr)
#    for array in arrays:
#        array.plotArrayFactorCut("phi", [0],ax=ax,polar=plr,label="Uniform, " + str(array.Nelements_x) + "x" + str(array.Nelements_y))
#
    angles = [0, 22.5, 45, 67.5, 90]
    arrays = []
    for angle in angles:
        paa = RectangularArray(5,5,0.5,0.5,f,arraytype="dt",sidelobe_level=30)
        paa.getArrayFactor(0,180,1,0,360,1,angle,0)
        arrays.append(paa)


    plr = True
    ax = plt.subplot(111,polar=plr)
    for i,array in enumerate(arrays):
        array.plotArrayFactorCut("phi", [0],ax=ax,polar=plr,label="Uniform, " + str(angles[i]) + "Deg")

    arrays[0].plotArrayFactor3D(logarithmic=False)


    plt.show()


