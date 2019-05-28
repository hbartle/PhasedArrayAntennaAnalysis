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

from util import *


class CSTResult:
    def __init__(self,run_id,theta_scan,phi_scan,field_cut):
        self.run_id = run_id
        self.theta_scan = theta_scan
        self.phi_scan = phi_scan
        self.field_cut = field_cut

    def calcMaxDir(self):
        """ Calculate maximum directivity and direction thereof """
        ind = np.unravel_index(np.argmax(self.field_cut[:,1]), self.field_cut[:,1].shape)
        max_dir = self.field_cut[ind,1]
        max_dir_theta =  self.field_cut[ind,0]

        return max_dir[0],max_dir_theta[0]

    def calcHPBW(self):
        """ Calculate Half-Power Beamwidth """
        max_dir,max_dir_theta = self.calcMaxDir()
        half_power_level = max_dir - 3
        index = (np.abs(self.field_cut[:,1] - half_power_level)).argmin()
        hpbw = 2*np.abs(max_dir_theta-self.field_cut[index,0])

        return hpbw


    def calcSLL(self):
        """ Calculate Side-Lobe Level """
        max_dir,max_dir_theta = self.calcMaxDir()
        ind = argrelextrema(self.field_cut[:, 1], np.greater)
        maxima = self.field_cut[ind, 1].flatten()
        maxima_sorted = maxima[maxima.argsort()]
        max_sidelobe = maxima_sorted[-2:-1]
        sll = -(max_dir -max_sidelobe)
        return sll[0],max_sidelobe[0]

    def calcSteeringOffset(self):
        """ Calculate the difference between steering angle and actual main beam direction """
        ind = np.unravel_index(np.argmax(self.field_cut[:,1]), self.field_cut[:,1].shape)
        max_dir_theta =  self.field_cut[ind,0]

        steering_offset = np.abs(max_dir_theta-self.theta_scan)
        return steering_offset[0]

class CSTSweepResult:
    """ Class for results obtain by a CST scan sweep"""

    def __init__(self,theta=[],phi=[],directivity=[],gain_ieee=[],gain_realized=[],array_size=[],array_type=[]):
        self.theta_scan = theta
        self.phi_scan = phi
        self.directivity = directivity
        self.gain_ieee = gain_ieee
        self.gain_realized = gain_realized
        self.array_size=array_size
        self.array_type = array_type


    def constructSTKAntennaFile(self,output_file_dir):
        """ Generate the External Antenna File for STK simulations """
        print("Construct STK Antenna file from scanning sweep data...")
        print("File: " + output_file_dir)
        filename = "stk_antenna_file_" + self.array_size + "_" + self.array_type + "_gain_ieee.txt"
        with open(output_file_dir + filename, "w") as file:
            file.write("stk.v.6.0\n" +
                       "ThetaPhiPattern\n" +
                       "AngleUnits Degrees\n" +
                       "NumberOfPoints " + str(self.theta_scan.size * 2) + "\n" +
                       "PatternData\n")
            np.savetxt(file, np.array([self.theta_scan.astype(np.float),
                                       self.phi_scan.astype(np.float),
                                       self.gain_ieee]).T, fmt="%.2f")
            np.savetxt(file, np.array([self.theta_scan.astype(np.float),
                                       self.phi_scan.astype(np.float)+180,
                                       self.gain_ieee]).T, fmt="%.2f")

        filename = "stk_antenna_file_" + self.array_size + "_" + self.array_type + "_gain_realized.txt"
        with open(output_file_dir + filename, "w") as file:
            file.write("stk.v.6.0\n" +
                       "ThetaPhiPattern\n" +
                       "AngleUnits Degrees\n" +
                       "NumberOfPoints " + str(self.theta_scan.size * 2) + "\n" +
                       "PatternData\n")

            np.savetxt(file, np.array([self.theta_scan.astype(np.float),
                                       self.phi_scan.astype(np.float),
                                       self.gain_realized]).T, fmt="%.2f")
            np.savetxt(file, np.array([self.theta_scan.astype(np.float),
                                       self.phi_scan.astype(np.float)+180,
                                       self.gain_realized]).T, fmt="%.2f")
        print("Done.")
class CSTResultLoader:

    def loadSweepResultFile(self, results_file, run_id_file):
        """
        :param results_file: Path to CST results file
        :param run_id_file: Path to Result navigator file
        :return: Result object containing Directivity, Gain (IEEE/Realized) of scanning sweep
        """
        with open(results_file,newline='') as f:
            res = []
            start=0
            for line in f:
                if "----------------------------------------------------------------------" in line:
                    start=1
                elif line in ['\n', '\r\n']:
                    start=0
                elif start:
                    res.append(list(map(float, line.split())))
        results = np.asarray(res)
        sim_results = np.split(results,3)
        directivity = sim_results[0]
        gain_ieee = sim_results[1]
        gain_realized = sim_results[2]

        theta = []
        phi = []
        with open(run_id_file,'r') as f:
            reader = csv.reader(f, delimiter=",")
            next(reader)
            for line in reader:
                if not line[0] == '0':
                    phi.append(line[1])
                    theta.append(line[2])
        theta = np.asarray(theta)
        phi = np.asarray(phi)




        array_type = directory.split("/")
        result_object = CSTSweepResult(theta=theta, phi=phi,
                                       directivity=directivity[:,1],
                                       gain_ieee=gain_ieee[:,1],
                                       gain_realized=gain_realized[:,1],
                                       array_size=array_type[2],
                                       array_type=array_type[3])

        return result_object

    def loadCutsResultFile(self, results_file, run_id_file):
        with open(results_file,newline='') as f:
            res = []
            start=0
            for line in f:

                if "----------------------------------------------------------------------" in line:


                    start=1
                elif line in ['\n', '\r\n']:
                    start=0
                elif start:
                    res.append(list(map(float, line.split())))
        results = np.asarray(res)
        results = np.split(results,361)
        result_objects = []
        with open(run_id_file,newline='') as f:
            csvReader = csv.reader(f)
            next(csvReader,None)
            for i,row in enumerate(csvReader):
                if not row[0] == '0':
                    row = list(map(int, row))
                    obj = CSTResult(row[0],row[1],row[2],results[i])
                    result_objects.append(obj)
        return result_objects

    def convertASCII2STKAntennaFile(self, directory, results_file):
        print("Converting far-field ASCII file to STK antenna file...")
        print("File: "+ directory + results_file)
        with open(directory + results_file,'r') as f:
            data = np.genfromtxt(f,skip_header=2)
            theta = data[:,0]
            phi = data[:,1]
            gain = data[:,2]

        filename = "stk_antenna_file_" + results_file
        with open(directory + filename, "w") as file:
            file.write("stk.v.6.0\n" +
                       "ThetaPhiPattern\n" +
                       "AngleUnits Degrees\n" +
                       "NumberOfPoints " + str(gain.size) + "\n" +
                       "PatternData\n")

            np.savetxt(file, np.array([theta,phi,gain]).T, fmt="%.2f")

        print("Done.")

class CSTResultAnalyzer:
    """ Analyze CSTResult Objects """

    def __init__(self,result_objects):
        self.result_objects = result_objects



    def analyzeCut(self,theta_scan, phi_scan):
        """ Analyze a single cut """
        fig = plt.figure()
        ax = fig.add_subplot(111, polar=True)
        for obj in self.result_objects:
            if obj.phi_scan == phi_scan and obj.theta_scan == theta_scan:

                # Plot Cut
                ax.plot(obj.field_cut[:, 0] * math.pi / 180, obj.field_cut[:, 1],
                        label="SCANTHETA=" + str(obj.theta_scan) + ",SCANPHI=" + str(obj.phi_scan))

                plot_min = obj.field_cut[:,1].min()
                plot_max = obj.field_cut[:,1].max()+obj.field_cut[:,1].max()/3
                plt.ylim(plot_min, plot_max)

                # Plot Max Directivity
                max_dir,max_dir_direction = obj.calcMaxDir()
                ax.plot((0, max_dir_direction*math.pi/180),(plot_min,max_dir),"g--",linewidth=0.8)

                # Plot HPBW
                hpbw = obj.calcHPBW()
                ax.plot((0,(max_dir_direction+hpbw/2)*math.pi/180),(plot_min,plot_max), "r--",linewidth=0.8)
                ax.plot((0,(max_dir_direction-hpbw/2)*math.pi/180),(plot_min,plot_max), "r--",linewidth=0.8)

                # Plot SLL
                sll, max_sidelobe = obj.calcSLL()
                ax.plot(np.linspace(0, 2 * np.pi, 100), np.ones(100) * max_sidelobe, "m--",linewidth=0.8)

                # Steering Offset
                steering_offset = obj.calcSteeringOffset()
                ax.plot((0,(obj.theta_scan)*math.pi/180),(plot_min,plot_max), "k--",linewidth=0.8)



        ax.legend()
        print(max_dir)
        ax.annotate('Max. Directivity: ' + str(np.round(max_dir,2)) + "dBi\n" +
                    "HPBW: " + str(np.round(hpbw)) +"°\n" +
                    "Sidelobe-Level: " + str(np.round(sll,2)) + "dB\n" +
                    "Steering Offset: " + str(steering_offset) + "°",
                    xy=(0.1, 0.9),
                    xycoords=('figure fraction', 'figure fraction'),
                    size=10, ha='left', va='center')

        ax.set_theta_zero_location("N")
        plt.grid(True)
        plt.show()

    def plotHeatmap(self,title = "",value="gain_realized",interp=False):
        res = self.result_objects
        if value == "gain_realized":
            v = res.gain_realized
            unit = "dBi"
        elif value == "gain_ieee":
            v = res.gain_ieee
            unit = "dBi"
        elif value == "directivity":
            v = res.directivity
            unit = "dBi"
        elif value == "full_efficiency":
            v = 100*db2pow(res.gain_realized - res.directivity)
            unit = "Percent"
        elif value == "radiation_efficiency":
            v = 100*db2pow(res.gain_ieee - res.directivity)
            unit = "Percent"
        else:
            print("Specify a correct quantity: directivity,gain_ieee,gain_realized")

        phi_mesh = np.reshape(res.phi_scan, (-1, np.unique(res.phi_scan).size))
        theta_mesh = np.reshape(res.theta_scan, (-1, np.unique(res.phi_scan).size))
        v_mesh = np.reshape(v,(-1, np.unique(res.phi_scan).size))

        if interp == True:
            theta_interp = np.arange(0, 90, 1)
            phi_interp = np.arange(0,180, 1)
            f = interpolate.interp2d(theta_mesh.astype(np.float), phi_mesh.astype(np.float), v_mesh, kind='linear')
            theta_mesh, phi_mesh = np.meshgrid(theta_interp,phi_interp)
            v_mesh = f(theta_interp,phi_interp)



        fig = plt.figure()
        fig.canvas.set_window_title(title)
        ax = plt.subplot()

        c = ax.pcolormesh(phi_mesh, theta_mesh, v_mesh,
                          cmap='jet',
                          vmin=v_mesh.min(), vmax=v_mesh.max())
        ax.set_title(title.split("-")[-1])
        ax.set_xlabel("Phi Scan [°]")
        ax.set_ylabel("Theta Scan [°]")
        ax.xaxis.set_major_locator(ticker.MultipleLocator(2))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
        colorbar_ticks = np.linspace(v_mesh.min(),v_mesh.max(),8)
        cbar = fig.colorbar(c, ax=ax,ticks=colorbar_ticks)
        #cbar = fig.colorbar(c, ax=ax)
        cbar.set_label(unit, rotation=270, labelpad=20)

        return fig

        # Find and plot line of constant value
        #idx = np.where(v>=(v.max()-3))[0]
        #const = v[idx]
        #print(const)

    def calcSteeringBandwidth(self):
        gain_realized = self.result_objects.gain_realized
        max_gain_realized = gain_realized.max()

        bandwidth_gain = max_gain_realized - gain_realized
        index = np.where(bandwidth_gain>=3)
        theta = np.asarray(list(map(int,self.result_objects.theta_scan[index])))
        return theta.min()


if __name__ == "__main__":
    result_loader = CSTResultLoader()

    plot_directory = "../plots/steering_sweep_results/"

    directories = ["../sim_data/5x5/uniform/",
                 "../sim_data/5x5/binomial/",
                 "../sim_data/5x5/chebychev/",
                 "../sim_data/4x4/uniform/",
                 "../sim_data/4x4/binomial/",
                 "../sim_data/4x4/chebychev/",
                 "../sim_data/3x3/uniform/",
                 "../sim_data/3x3/binomial/",
                 "../sim_data/3x3/chebychev/",
                 "../sim_data/2x2/uniform/",
                 "../sim_data/2x2/binomial/",
                 "../sim_data/2x2/chebychev/"]

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

    max_gain_realized = []
    min_sll = []
    min_hpbw =[]
    max_steering_offset = []
    steering_bandwidth = []


    for directory in directories:
        result_objects_cut = result_loader.loadCutsResultFile(directory + "Directivity,Phi=0.txt",
                                                "../sim_data/result_navigator.csv")
        res_analyzer = CSTResultAnalyzer(result_objects_cut)
        #res_analyzer.analyzeCut(0,0)

        # Load Scanning Sweep Result and plot Heatmaps
        result_object_sweep = result_loader.loadSweepResultFile(directory + "DGG,Value,Theta=PAA_FA_SCANTHETA,Phi=PAA_FA_SCANPHI.txt",
                                                      "../sim_data/result_navigator.csv")

        analyzer = CSTResultAnalyzer(result_object_sweep)
        # Realized Gain
        title = directory.split('/')[2] + '-' + directory.split('/')[3] + " - Scanning Direction, Realized Gain"
        #fig = analyzer.plotHeatmap(title=title, value="gain_realized")
        #fig.savefig(plot_directory + title + ".png")

        # Gain IEEE
        title = directory.split('/')[2] + '-' + directory.split('/')[3] +" - Scanning Direction, Gain IEEE"
        #fig = analyzer.plotHeatmap(title=title, value="gain_ieee")
        #fig.savefig(plot_directory + title + ".png")

        # Save Max Values for Comparison plots
        max_gain_realized.append(result_object_sweep.gain_realized.max())
        steering_bandwidth.append(analyzer.calcSteeringBandwidth())


        result_object_sweep = result_loader.loadSweepResultFile(directory + "DGG,Max. Value,Phi=PAA_FA_SCANPHI.txt",
                                                      "../sim_data/result_navigator.csv")
        analyzer = CSTResultAnalyzer(result_object_sweep)

        # Realized Gain Max Value
        title = directory.split('/')[2] + '-' + directory.split('/')[3] + " - Max. Value, Realized Gain"
        #fig = analyzer.plotHeatmap(title=title,value="gain_realized")
        #fig.savefig(plot_directory + title + ".png")

        # Gain IEEE Max Value
        title = directory.split('/')[2] + '-' + directory.split('/')[3] + " - Max. Value, Gain IEEE"
        #fig = analyzer.plotHeatmap(title=title,value="gain_ieee")
        #fig.savefig(plot_directory + title + ".png")

        # Efficiency
        title = directory.split('/')[2] + '-' + directory.split('/')[3] + " - Full Efficiency"
        #fig = analyzer.plotHeatmap(title=title, value="full_efficiency")
        #fig.savefig(plot_directory + title + ".png")

        title = directory.split('/')[2] + '-' + directory.split('/')[3] + " - Radiation Efficiency"
        #fig = analyzer.plotHeatmap(title=title, value="radiation_efficiency")
        #fig.savefig(plot_directory + title + ".png")

        #result_object_sweep.constructSTKAntennaFile(directory)



    # Create STK antenna file for reference antennas
    # result_loader.convertASCII2STKAntennaFile("../sim_data/reference_antennas/","8ghz_patch_realized_gain.txt")
    # result_loader.convertASCII2STKAntennaFile("../sim_data/reference_antennas/",
    #                                           "4x4_chebychev_broadside_array_realized_gain.txt")

    #plt.show()
    for i,n in enumerate(max_gain_realized):
        print("(" + str(i%4)+"," + str(max_gain_realized[i]) + ")")
