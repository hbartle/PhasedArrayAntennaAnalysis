#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Hannes Bartle <mail@hbartle.de>
#
# Distributed under terms of the MIT license.

import os
import matplotlib.pyplot as plt

from scipy.interpolate import interp2d

# Own Modules
from antenna_toolbox.Antenna import Antenna
from antenna_toolbox.RectangularArray import RectangularArray
from util import *

if __name__ == "__main__":
    # Load Unit Cell
    f = 8e9
    patch = Antenna()
    patch.loadPatternFromFile("sim/patch_8Ghz_efield_ro5880_0.254.txt")

    # Array Dimensions
    dims = [2,3,4,5]
    #dims = [3,4]

    # Steering angles
    theta_angles = np.linspace(0,90,19)
    phi_angles = np.linspace(0,360,73)
    #theta_angles = np.linspace(0,90,4)
    #phi_angles = np.linspace(0,90,4)
    theta_mesh,phi_mesh = np.meshgrid(theta_angles,phi_angles)

    dir_heatmap = np.zeros((theta_angles.size,phi_angles.size))
    arrays = []
    dir_heatmaps = []


    for dim in dims:
        result_path = "results/steering_sweep/uniform/array"+ str(dim) + "x" + str(dim) + "/"
        if not os.path.exists(result_path):
            os.makedirs(result_path)

        max_dir_filename = "max_directivity_" + str(dim) + "x" + str(dim) + "_array.txt"

        if not os.path.isfile(result_path + max_dir_filename):
            print("Calculating " + str(dim) + "x" + str(dim) + "-Array")
            # Create Array
            paa = RectangularArray(dim,dim,0.5,0.5,f,arraytype="uniform")
            paa.setElement(patch)


            # Loop through steering angles
            for t,theta in enumerate(theta_angles):
                for p,phi in enumerate(phi_angles):
                    paa.getField(theta,phi)
                    paa.calcDirectivity()

                    dir_filename = "directivity_" + str(dim) + "x" + str(dim) + "_array_" + str(theta) + "-" + str(phi) +".txt"
                    # Store directivity pattern
                    np.savetxt(result_path + dir_filename,np.array((paa.theta_mesh.flatten(),paa.phi_mesh.flatten(),paa.directivity.flatten())).T)

                    dir_heatmap[t,p]= paa.max_directivity
                    print("Theta= " +str(theta)+", Phi= " + str(phi) + ": Max. Dir= " + str(pow2db(paa.max_directivity)))

            #Store max directivity map
            np.savetxt(result_path + max_dir_filename,np.array((theta_mesh.flatten(),phi_mesh.flatten(),dir_heatmap.flatten())).T)
            arrays.append(paa)

        else:
            print("Loading stored results...")
            theta_mesh, phi_mesh, dir_heatmap = np.loadtxt(result_path + max_dir_filename, unpack=True)
            ntheta = np.unique(theta_mesh).size
            nphi = np.unique(phi_mesh).size
            theta_mesh = np.reshape(theta_mesh,(nphi,ntheta))
            phi_mesh = np.reshape(phi_mesh,(nphi,ntheta))
            dir_heatmap = np.reshape(dir_heatmap,(nphi,ntheta))


        dir_heatmaps.append(dir_heatmap)


    x = np.arange(0, 90, 0.5)
    y = np.arange(0, 90, 0.5)
    X, Y = np.meshgrid(x, y)

    for heatmap in dir_heatmaps:

        h = pow2db(heatmap)
        f = interp2d(theta_mesh, phi_mesh, h, kind='cubic')

        fig,ax = plt.subplots()
        c = ax.pcolormesh(theta_mesh,phi_mesh,h, cmap='jet',vmin=h.min(),vmax=h.max())
        ax.set_title("Maximum Directivity")
        fig.colorbar(c, ax=ax)

    plt.show()
