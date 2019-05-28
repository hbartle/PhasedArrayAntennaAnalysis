#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Hannes Bartle <mail@hbartle.de>
#
# Distributed under terms of the MIT license.

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

plt.style.use('seaborn-whitegrid')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 14})

from antenna_toolbox.RectangularArray import RectangularArray
from util import *

"""
    Plot the phase shift for each array element given a certain steering angle
"""

if __name__=="__main__":

    theta_angles = np.arange(0,90,1)
    phi_angles = np.arange(0,360,1)
    theta_angles = np.asarray([10,20,30])
    phi_angles = np.asarray([45])




    array = RectangularArray(10, 10, 0.5, 0.5, 8e9)
    for theta in theta_angles:
        for phi in phi_angles:
            array.getArrayFactor(0,90,10,0,360,20,theta,phi)
            phase = rad2deg(array.element_phase_shift)
            phase = np.mod(phase, 360)



            patches = []
            color_list = phase.flatten()

            # Add rectangles
            width = 0.01*1000
            height = 0.01*1000
            for a_x, a_y in zip(array.elements_x.flatten()*1000,array.elements_y.flatten()*1000):
                patch = Rectangle(xy=(a_x - width / 2, a_y - height / 2), width=width, height=height)
                patches.append(patch)


            collection = PatchCollection(patches,cmap=plt.cm.hot)
            collection.set_array(color_list)

            # Figure
            fig = plt.figure()
            fig.canvas.set_window_title("Element Phase Shift" + " Theta: " + str(theta) + "°, Phi: " + str(phi) + "°")
            ax = plt.subplot()

            ax.add_collection(collection)


            ax.set_title("Element Phase Shift\n" + "Theta: " + str(theta) + "°, Phi: " + str(phi) + "°")
            ax.set_xlim(array.elements_x.min()*1000*1.2,array.elements_x.max()*1000*1.2)
            ax.set_ylim(array.elements_y.min()*1000*1.2,array.elements_y.max()*1000*1.2)
            ax.set_xticks(np.unique(array.elements_x)*1000)
            ax.set_yticks(np.unique(array.elements_y)*1000)
            ax.set_facecolor('grey')
            plt.xticks(rotation=45)
            plt.colorbar(collection)

    plt.show()

