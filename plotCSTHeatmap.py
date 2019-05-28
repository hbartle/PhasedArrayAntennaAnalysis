#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Hannes Bartle <mail@hbartle.de>
#
# Distributed under terms of the MIT license.



import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
plt.style.use('seaborn-whitegrid')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 14})
import matplotlib.gridspec as gridspec


"""
    Plot Heatmaps (Theta/Phi) of the different antenna performance indicators
"""


if __name__=="__main__":

    plot_directory = "../plots/steering_sweep_results/"

    directories = ["../sim_data/5x5/chebychev/",
                   "../sim_data/4x4/chebychev/",
                   "../sim_data/3x3/chebychev/",
                   "../sim_data/2x2/chebychev/",
                   "../sim_data/5x5/uniform/",
                   "../sim_data/4x4/uniform/",
                   "../sim_data/3x3/uniform/",
                   "../sim_data/2x2/uniform/",
                   "../sim_data/5x5/binomial/",
                   "../sim_data/4x4/binomial/",
                   "../sim_data/3x3/binomial/",
                   "../sim_data/2x2/binomial/"]

    directories = ["../sim_data/2x2/uniform/",
                   "../sim_data/3x3/uniform/",
                   "../sim_data/4x4/uniform/",
                   "../sim_data/5x5/uniform/",
                   "../sim_data/2x2/binomial/",
                   "../sim_data/3x3/binomial/",
                   "../sim_data/4x4/binomial/",
                   "../sim_data/5x5/binomial/",
                   "../sim_data/2x2/chebychev/",
                   "../sim_data/3x3/chebychev/",
                   "../sim_data/4x4/chebychev/",
                   "../sim_data/5x5/chebychev/"]

    min_sll = []
    max_sll = []
    rms_sll = []

    min_hpbw =[]
    max_hpbw =[]

    max_direction_offset =[]


    for directory in directories:



        # Max Directivity
        #filename = "Directivity,Theta=PAA_FA_SCANTHETA,Phi=PAA_FA_SCANPHI.txt"
        #filename = "Directivity,Phi=PAA_FA_SCANPHI,Max. Value.txt"
        #phi_scan, theta_scan, max_dir  = np.loadtxt(directory + filename, skiprows=2, unpack=True)
        #phi_mesh = np.reshape(phi_scan,(-1,np.unique(phi_scan).size))
        #theta_mesh = np.reshape(theta_scan,(-1,np.unique(phi_scan).size))
        #max_dir_mesh = np.reshape(max_dir,(-1,np.unique(phi_scan).size))



        filename = "result_navigator.csv"
        run_id, phi_scan, theta_scan  = np.loadtxt('../sim_data/' + filename,delimiter=',', skiprows=1, unpack=True, )
        phi_mesh = np.reshape(phi_scan, (-1, np.unique(phi_scan).size))
        theta_mesh = np.reshape(theta_scan, (-1, np.unique(phi_scan).size))



        # Side Lobe Level
        filename = "Directivity,Phi=PAA_FA_SCANPHI,Side Lobe Level.txt"

        run_id, sll_phi  = np.loadtxt(directory + filename, skiprows=2, unpack=True)
        sll_phi_mesh = np.reshape(sll_phi, (-1, np.unique(phi_scan).size))

        # HPBW
        filename = "Directivity,Phi=PAA_FA_SCANPHI,HPBW.txt"
        run_id, hpbw  = np.loadtxt(directory + filename, skiprows=2, unpack=True)
        hpbw_mesh = np.reshape(hpbw, (-1, np.unique(phi_scan).size))

        # Direction Offset
        filename = "Directivity,Phi=PAA_FA_SCANPHI,Main Lobe Direction.txt"
        run_id, main_lobe_direction  = np.loadtxt(directory + filename, skiprows=2, unpack=True)
        main_lobe_direction = np.absolute(main_lobe_direction)
        direction_offset = np.reshape(main_lobe_direction, (-1, np.unique(phi_scan).size))  - theta_mesh

        #
        # fig2 = plt.figure()
        # fig2.canvas.set_window_title(directory + "- Maximum Side Lobe Level")
        # ax2 = plt.subplot()
        # fig3 = plt.figure()
        # fig3.canvas.set_window_title(directory + "- HPBW")
        # ax3 = plt.subplot()
        # fig4 = plt.figure()
        # fig4.canvas.set_window_title(directory + "- Steering Direction Offset")
        # ax4 = plt.subplot()
        #
        #
        #
        # c = ax2.pcolormesh(phi_mesh, theta_mesh, sll_phi_mesh, cmap='jet', vmin=sll_phi.min(), vmax=sll_phi.max())
        # ax2.set_title("Maximum Side Lobe Level")
        # ax2.set_xlabel("Phi Scan [°]")
        # ax2.set_ylabel("Theta Scan [°]")
        # ax2.set_xticks(np.arange(0,200,20))
        # ax2.set_yticks(np.arange(0,100,10))
        # cbar = fig2.colorbar(c, ax=ax2)
        # cbar.set_label('dB', rotation=270,labelpad=20)
        # fig2.savefig(plot_directory + directory.split('/')[2] + "-" + directory.split('/')[3] + "- Maximum Side Lobe Level.png")
        #
        # c = ax3.pcolormesh(phi_mesh, theta_mesh, hpbw_mesh, cmap='jet', vmin=hpbw.min(), vmax=hpbw.max())
        # ax3.set_title("HPBW")
        # ax3.set_xlabel("Phi Scan [°]")
        # ax3.set_ylabel("Theta Scan [°]")
        # ax3.set_xticks(np.arange(0,200,20))
        # ax3.set_yticks(np.arange(0,100,10))
        # cbar = fig3.colorbar(c, ax=ax3)
        # cbar.set_label('Degrees', rotation=270,labelpad=20)
        # fig3.savefig(plot_directory + directory.split('/')[2] + "-" + directory.split('/')[3] + "- HPBW.png")
        #
        #
        # c = ax4.pcolormesh(phi_mesh, theta_mesh, direction_offset, cmap='jet', vmin=direction_offset.min(), vmax=direction_offset.max())
        # ax4.set_title("Steering Direction Offset")
        # ax4.set_xlabel("Phi Scan [°]")
        # ax4.set_ylabel("Theta Scan [°]")
        # ax4.set_xticks(np.arange(0,200,20))
        # ax4.set_yticks(np.arange(0,100,10))
        # cbar = fig4.colorbar(c, ax=ax4)
        # cbar.set_label('Degrees', rotation=270,labelpad=20)
        # fig4.savefig(plot_directory + directory.split('/')[2] + "-" + directory.split('/')[3] + "- Steering Direction Offset.png")


        min_sll.append(sll_phi.min())
        max_sll.append(sll_phi.max())
        rms_sll.append(np.sqrt(np.mean(sll_phi**2)))
        #rms_sll.append(np.mean(sll_phi))

        min_hpbw.append(hpbw.min())
        max_hpbw.append(hpbw.max())
        max_direction_offset.append(direction_offset.min())

    for i, n in enumerate(rms_sll):
        print("(" + str(i % 4) + "," + str(rms_sll[i]) + ")")
    #for i,n in enumerate(min_sll):
     #   print("(" + str(i%4)+"," + str(min_sll[i]) + ")")
    #for i,n in enumerate(max_sll):
     #   print("(" + str(i%4)+"," + str(max_sll[i]) + ")")
    #for i,n in enumerate(min_hpbw):
     #   print("(" + str(i%4)+"," + str(min_hpbw[i]) + ")")
    #for i,n in enumerate(max_hpbw):
     #   print("(" + str(i%4)+"," + str(max_hpbw[i]) + ")")
    #for i,n in enumerate(max_direction_offset):
     #   print("(" + str(i%4)+"," + str(max_direction_offset[i]) + ")")

    #plt.show()


