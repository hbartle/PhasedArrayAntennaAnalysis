#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Hannes Bartle <mail@hbartle.de>
#
# Distributed under terms of the MIT license.
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from matplotlib import cm
from math import cos, sin, sqrt
from mpl_toolkits.mplot3d import Axes3D

from util import *

class Antenna:
    """ Antenna Parent Class

        Implements general functions needed for all types of antennas, e.g. plotting.
    """

    def loadPatternFromFile(self,patternfile):
        """ Loads a .txt file containing the simulated far-field pattern

            Format has to be CST Studio FarField Export.

            Args:
                patternfile: Path to file containing the pattern

            Return:
                pattern: Antenna Pattern

        """

        # Load txt file
        theta,phi,e_abs,e_theta,phase_theta,e_phi,phase_phi,ax_ratio= np.loadtxt(patternfile,skiprows=2,unpack=True)


        self.Ntheta = np.unique(theta).size
        self.Nphi = np.unique(phi).size
        self.theta_step = theta[1]-theta[0]
        self.phi_step = phi[1+self.Nphi] - phi[0]


        # Reshape into nd Array and copy value at phi=0/360 degree
        self.theta_mesh = deg2rad(np.reshape(theta,(self.Nphi,self.Ntheta))).T
        self.theta_mesh = np.concatenate([self.theta_mesh,self.theta_mesh[:,0].reshape((self.Ntheta,1))],axis=1)

        self.phi_mesh = deg2rad(np.reshape(phi,(self.Nphi,self.Ntheta))).T
        self.phi_mesh = np.concatenate([self.phi_mesh,self.phi_mesh[:,0].reshape((self.Ntheta,1))+2*math.pi],axis=1)

        pattern = np.reshape(e_abs,(self.Nphi,self.Ntheta)).T
        pattern = np.concatenate([pattern,pattern[:,0].reshape((self.Ntheta,1))],axis=1)

        self.pattern = normalize(pattern)

        # Add one sample at phi=360 degree
        self.Nphi=self.Nphi+1
        self.theta_start = rad2deg(self.theta_mesh.min())
        self.theta_stop = rad2deg(self.theta_mesh.max())
        self.phi_start = rad2deg(self.phi_mesh.min())
        self.phi_stop = rad2deg(self.phi_mesh.max())

    def calcDirectivity(self):
        """ Calculate Maximum directiviy of antenna

            Returns:
                directivity: Directivity in dB
                max_directivity: Maximum Directivity in specific direction
                theta_max: Theta angle of direction with maximum directivity
                phi_max: Phi angle of direction with maximum directivity

        """
        pattern = normalize(self.pattern)

        rad_intensity = pattern*pattern

        rad_intensity_iso = 0
        dtheta = deg2rad(self.theta_step)
        dphi = deg2rad(self.phi_step)

        #for t,theta in enumerate(self.theta_mesh):
        #    for p,phi in enumerate(self.phi_mesh):
        #        rad_intensity_iso = rad_intensity_iso + rad_intensity[t,p]*np.sin(theta)*dtheta*dphi

        rad_intensity_iso = np.sum(rad_intensity*np.sin(self.theta_mesh)*dtheta*dphi)

        rad_intensity_iso = 1/(4 * math.pi) * rad_intensity_iso

        self.directivity = rad_intensity/rad_intensity_iso

        self.max_directivity = self.directivity.max()

        return self.directivity,self.max_directivity


    def calcHPBW(self):
        """ Calculate Half-Power Beamwidth of Antenna

            Returns:
                hpbw: Half-Power Beamwidth in Degrees
        """
        pass



    def getPatternCut(self,dim,angle):
        if dim =="theta":
            t_idx=round(angle/self.theta_step)

            return self.pattern[t_idx,:], self.phi_mesh[t_idx,:]


        elif dim =="phi":
            p_idx= int(round(angle/self.phi_step))
            offset = int(round(180/self.phi_step))
            pattern_cut =  np.concatenate([np.flip(self.pattern[:,p_idx + offset]),self.pattern[:,p_idx]])
            cut_angles = np.concatenate([-np.flip(self.theta_mesh[:,p_idx + offset]),self.theta_mesh[:,p_idx]])

            return pattern_cut,cut_angles

    def plotPatternCut(self,dim,angles,ax=None,logarithmic=True,polar=True,labels=["Cut"]):
        """ Plot Pattern Cuts in given dimension

            Args:
                dim: Dimension of cut (Theta/Phi).
                angles: List of cutting angles.
                ax: Axis object to plot into
                logarithmic: Plot in dB
                polar: Polar plot

            Return:
                ax: Axes object
        """

        if polar == True:
            if ax==None:
                fig = plt.figure()
                ax = fig.add_subplot(111,polar=True)

            ax.set_theta_zero_location('N')
            ax.grid(True)


            for i,angle in enumerate(angles):

                pattern_cut,cut_angles = self.getPatternCut(dim,angle)

                if dim =="phi":
                    ax.set_xticks(np.pi/180. * np.linspace(180,  -180, 8, endpoint=False))
                    ax.set_thetalim(-np.pi, np.pi)

                if logarithmic == True:
                    ax.plot(cut_angles,lin2db(pattern_cut),label=labels[i])
                    ax.set_yticks(range(-50,0 , 10))
                    ax.set_ylim(-50, 0)
                else:
                    ax.plot(cut_angles,pattern_cut,label=labels[i])
                    ax.set_yticks(range(0,1 , 0.2))
                    ax.set_ylim(0, 1)

        elif polar == False:

            if ax==None:
                fig = plt.figure()
                ax = fig.add_subplot(111,polar=False)
            ax.grid(True)

            for i,angle in enumerate(angles):

                pattern_cut,cut_angles = self.getPatternCut(dim,angle)

                if logarithmic == True:
                    ax.plot(cut_angles,lin2db(pattern_cut),label=labels[i])
                    ax.set_yticks(range(-50,0 , 10))
                    ax.set_ylim(-50, 0)
                else:
                    ax.plot(cut_angles,pattern_cut,label=labels[i])
                    ax.set_yticks(range(0,1 , 0.2))
                    ax.set_ylim(0, 1)

        ax.legend(loc="lower right",bbox_to_anchor=(1.4, 0))

        return ax,pattern_cut,cut_angles



    def plotPattern3D(self,ax=None, logarithmic=True,log_range=-20):
        """ Plot 3D Pattern

            Args:
                ax: Axis object to plot into
                logarithmic: Plot in dB
                log_range: Minimum Value of plotting range in dB

            Return:
                ax: Axes object
        """

        if ax==None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')


        if logarithmic ==True:
            pattern_db = lin2db(self.pattern)
            pattern_min =  np.min(pattern_db)
            pattern_max =  np.max(pattern_db)
            e = pattern_db - (log_range + pattern_max)
            e[e<0.0] = 0.0
        else:
            e = self.pattern
            pattern_min =  np.min(self.pattern)
            pattern_max =  np.max(self.pattern)

        x,y,z = sph2cart(self.phi_mesh,self.theta_mesh,e)

        ax.plot_surface(x, y, z,facecolors=cm.jet(self.pattern))
        plt.xlabel('X')
        plt.ylabel('Y')
        lim_left,lim_right = ax.get_zlim()
        ax.set_xlim(lim_left-lim_right/2,lim_right/2)
        ax.set_ylim(lim_left-lim_right/2,lim_right/2)

        m = cm.ScalarMappable(cmap=plt.cm.jet)
        if logarithmic==True:
            m.set_array((e+(log_range+pattern_max)))
        else:
            m.set_array(e)
        cbar = plt.colorbar(m,ax=ax)
        cbar.set_label('dB', rotation=270,labelpad=12)

        return ax


    def getDirectivityCut(self,dim,angle):
        if dim =="theta":
            t_idx=round(angle/self.theta_step)

            return self.directivity[t_idx,:], self.phi_mesh[t_idx,:]


        elif dim =="phi":
            p_idx= int(round(angle/self.phi_step))
            offset = int(round(180/self.phi_step))
            directivity_cut =  np.concatenate([np.flip(self.directivity[:,p_idx + offset]),self.directivity[:,p_idx]])
            cut_angles = np.concatenate([-np.flip(self.theta_mesh[:,p_idx + offset]),self.theta_mesh[:,p_idx]])

            return directivity_cut,cut_angles

    def plotDirectivityCut(self,dim,angles,ax=None,logarithmic=True,polar=True,labels=["Cut"]):
        """ Plot Directivity Cuts in given dimension

            Args:
                dim: Dimension of cut (Theta/Phi).
                angles: List of cutting angles.
                ax: Axis object to plot into
                logarithmic: Plot in dB
                polar: Polar plot

            Return:
                ax: Axes object
        """

        if polar == True:
            if ax==None:
                fig = plt.figure()
                ax = fig.add_subplot(111,polar=True)

            ax.set_theta_zero_location('N')
            ax.grid(True)


            for i,angle in enumerate(angles):

                directivity_cut,cut_angles = self.getDirectivityCut(dim,angle)
                if dim =="phi":
                    ax.set_xticks(np.pi/180. * np.linspace(180,  -180, 8, endpoint=False))
                    ax.set_thetalim(-np.pi, np.pi)

                if logarithmic == True:
                    ax.plot(cut_angles,pow2db(directivity_cut),label=labels[i])
                    ax.set_yticks(range(-50,int(directivity_cut.max()) , 10))
                    ax.set_ylim(-50, int(directivity_cut.max()))
                else:
                    ax.plot(cut_angles,directivity_cut,label=labels[i])
                    #ax.set_yticks(range(0,1 , 0.2))
                    #ax.set_ylim(0, 1)

        elif polar == False:

            if ax==None:
                fig = plt.figure()
                ax = fig.add_subplot(111,polar=False)
            ax.grid(True)

            for i,angle in enumerate(angles):

                directivity_cut,cut_angles = self.getDirectivityCut(dim,angle)

                if logarithmic == True:
                    ax.plot(cut_angles,pow2db(directivity_cut),label=labels[i])
                    ax.set_yticks(range(-50,int(directivity_cut.max()) , 10))
                    ax.set_ylim(-50, int(directivity_cut.max()))
                else:
                    ax.plot(cut_angles,directivity_cut,label=labels[i])
                    #ax.set_yticks(range(0,1 , 0.2))
                    #ax.set_ylim(0, 1)

        ax.legend()

        return ax

    def plotDirectivity3D(self,ax=None, logarithmic=True,log_range=-20):
        """ Plot 3D Directivity

            Args:
                ax: Axis object to plot into
                logarithmic: Plot in dB
                log_range: Minimum Value of plotting range in dB

            Return:
                ax: Axes object
        """

        if ax==None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')


        if logarithmic ==True:
            directivity_db = pow2db(self.directivity)
            directivity_min =  np.min(directivity_db)
            directivity_max =  np.max(directivity_db)
            e = directivity_db - (log_range + directivity_max)
            e[e<0.0] = 0.0
        else:
            e = self.directivity
            directivity_min =  np.min(self.directivity)
            directivitye_max =  np.max(self.directivity)

        x,y,z = sph2cart(self.phi_mesh,self.theta_mesh,e)

        ax.plot_surface(x, y, z,facecolors=cm.jet(normalize(self.directivity)))
        plt.xlabel('X')
        plt.ylabel('Y')
        lim_left,lim_right = ax.get_zlim()
        ax.set_xlim(lim_left-lim_right/2,lim_right/2)
        ax.set_ylim(lim_left-lim_right/2,lim_right/2)

        m = cm.ScalarMappable(cmap=plt.cm.jet)
        if logarithmic==True:
            m.set_array((e+(log_range+directivity_max)))
        else:
            m.set_array(e)
        plt.colorbar(m)

        return ax
