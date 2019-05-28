#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Hannes Bartle <mail@hbartle.de>
#
# Distributed under terms of the MIT license.
import matplotlib.pyplot as plt

from matplotlib import cm

from util import *
from antenna_toolbox.Antenna import Antenna


class RectangularArray(Antenna):
    """ Rectangular Phased Array"""

    def __init__(self,Nelements_x,Nelements_y,d_x,d_y,freq,arraytype="uniform",sidelobe_level=20):

        self.Nelements_x = Nelements_x
        self.Nelements_y = Nelements_y
        self.element_distance_x = d_x
        self.element_distance_y = d_y
        self.frequency = freq


        wavelength = 3e8/freq

        x = d_x * wavelength * np.linspace(-Nelements_x/2,Nelements_x/2,Nelements_x)
        y = d_y * wavelength * np.linspace(-Nelements_y/2,Nelements_y/2,Nelements_y)

        self.elements_y, self.elements_x = np.meshgrid(y,x)


        self.arraytype = arraytype

        if arraytype == "uniform":
            self.uniformWeights()
        elif arraytype == "binomial":
            self.binomialWeights()
        elif arraytype == "dt":
            self.dtWeights(sidelobe_level)

    def setElement(self,element):
        self.element = element

    def uniformWeights(self):
        """ Uniform Amplitude Weights """
        self.weights = np.ones((self.Nelements_x,self.Nelements_y))
        return self.weights

    def binomialWeights(self):
        n = self.Nelements_x-1
        m = self.Nelements_y-1

        w_x = [binomial(n,i) for i in range(n+1)]
        w_y = [binomial(m,i) for i in range(m+1)]

        w_y_grid,w_x_grid = np.meshgrid(w_y,w_x)
        self.weights = w_y_grid*w_x_grid
        return self.weights


    def dtWeights(self,sidelobe_level):
        """ Calculate Dolph-Tchebychev weights

            See Balanis, Antenna Theory,4ed, Chapter 6.8.4, p338

        """
        R = db2lin(sidelobe_level)
        w_x = dolphTschebychev(self.Nelements_x,R)
        w_y = dolphTschebychev(self.Nelements_y,R)

        w_y_grid,w_x_grid = np.meshgrid(w_y,w_x)
        self.weights = w_y_grid*w_x_grid
        return self.weights


    def getField(self,theta_steering=0,phi_steering=0):
        """ Calculate Resulting Field (Element*AF)

            Uses the same step size and range for far-field calculation as the element.

            Args:
                theta_steering: Theta angle of beam steering direction
                phi_steering: Phi angle of beam steering direction

            Return:
                pattern: Combined Field Pattern
        """
        e = self.element
        self.getArrayFactor(e.theta_start,e.theta_stop,e.theta_step,e.phi_start,e.phi_stop,e.phi_step,theta_steering,phi_steering)
        Ntheta = (e.theta_stop-e.theta_start)/e.theta_step + 1
        Nphi = (e.phi_stop-e.phi_start)/e.phi_step + 1

        theta_range = deg2rad(np.linspace(e.theta_start,e.theta_stop,Ntheta))
        phi_range = deg2rad(np.linspace(e.phi_start,e.phi_stop,Nphi))
        self.pattern = np.ones((theta_range.size,phi_range.size))
        self.theta_mesh_pattern = np.ones((theta_range.size,phi_range.size))
        self.phi_mesh_pattern = np.ones((theta_range.size,phi_range.size))

        # Pattern Multiplication
        self.pattern = e.pattern*self.array_factor

        self.pattern=normalize(self.pattern)
        return self.pattern


    def getArrayFactor(self,theta_start,theta_stop,theta_step,phi_start,phi_stop,phi_step,theta_steering=0,phi_steering=0,pre_comp_phase=False):
        """ Calculate array factor

            Args:
                theta_step: Step size in Theta dimension
                phi_step: Step size in Phi dimension
                theta_steering: Direction of steered beam
                phi_steering: Direction of steered beam

            Return:
                array_factor: Array factor as ndarray in theta/phi-domain

        """

        # Wavenumber
        k = 2*math.pi*self.frequency/3e8

        # Vectors to each element
        rx = self.elements_x.flatten()
        ry = self.elements_y.flatten()
        rz = np.zeros(self.elements_x.size)
        r = np.vstack((rx,ry,rz))

        if pre_comp_phase == False:
            # Get element phase shifts from steering direction
            X0, Y0,Z0 = sph2cart(deg2rad(phi_steering),deg2rad(theta_steering),999)
            R0 = np.array((X0,Y0,Z0))
            self.element_phase_shift = -k * R0 @ r
            self.element_phase_shift = np.reshape(self.element_phase_shift,self.elements_x.shape)

        self.Ntheta = (theta_stop-theta_start)/theta_step + 1
        self.Nphi = (phi_stop-phi_start)/phi_step + 1
        self.theta_step = theta_step
        self.phi_step = phi_step

        theta_range = deg2rad(np.linspace(theta_start,theta_stop,self.Ntheta))
        phi_range = deg2rad(np.linspace(phi_start,phi_stop,self.Nphi))
        self.array_factor = np.ones((theta_range.size,phi_range.size))
        self.theta_mesh = np.ones((theta_range.size,phi_range.size))
        self.phi_mesh = np.ones((theta_range.size,phi_range.size))


        for t,theta in enumerate(theta_range):
            for p,phi in enumerate(phi_range):

                # Far field direction unit vector
                X,Y,Z= sph2cart(phi,theta,1)
                R = np.array((X,Y,Z))

                psi =  k * R @ r + self.element_phase_shift.flatten()
                self.array_factor[t,p] =   np.real(np.dot(self.weights.flatten(),np.exp(1j*psi)))
                self.theta_mesh[t,p] = theta
                self.phi_mesh[t,p] = phi


        self.array_factor = normalize(self.array_factor)

        return self.array_factor


    def getArrayFactorCut(self,dim,angle):
        """ Calculate Cutting plane in Array Factor

            Args:
                dim: Cutting Dimension (theta/phi)
                angle: Angle of Cutting Plane

            Return:
                array_factor_cut: Array of AF values in cutting plane
                cut_angles: Angles in cutting plane

        """
        if dim =="theta":
            t_idx=round(angle/self.theta_step)

            return self.array_factor[t_idx,:], self.phi_mesh[t_idx,:]


        elif dim =="phi":
            p_idx= round(angle/self.phi_step)
            offset = round(180/self.phi_step)
            array_factor_cut =  np.concatenate([np.flip(self.array_factor[:,p_idx + offset]),self.array_factor[:,p_idx]])
            cut_angles = np.concatenate([-np.flip(self.theta_mesh[:,p_idx + offset]),self.theta_mesh[:,p_idx]])

            return array_factor_cut,cut_angles


    def plotArrayFactorCut(self,dim,angles,ax=None,logarithmic=True,polar=True,label="Cut"):
        """ Plot AF Cuts in given dimension

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


            for angle in angles:

                array_factor_cut,cut_angles = self.getArrayFactorCut(dim,angle)

                if dim =="phi":
                    ax.set_xticks(np.pi/180. * np.linspace(180,  -180, 8, endpoint=False))
                    ax.set_thetalim(-np.pi, np.pi)

                if logarithmic == True:
                    ax.plot(cut_angles,lin2db(array_factor_cut),label=label)
                    ax.set_yticks(range(-60,0 , 10))
                    ax.set_ylim(-60, 0)
                else:
                    ax.plot(cut_angles,array_factor_cut,label=label)
                    ax.set_yticks(range(0,1 , 0.2))
                    ax.set_ylim(0, 1)

        elif polar == False:

            if ax==None:
                fig = plt.figure()
                ax = fig.add_subplot(111,polar=False)
            ax.grid(True)

            for angle in angles:

                array_factor_cut,cut_angles = self.getArrayFactorCut(dim,angle)

                if logarithmic == True:
                    ax.plot(cut_angles,lin2db(array_factor_cut),label=label)
                    ax.set_yticks(range(-60,0 , 10))
                    ax.set_ylim(-60, 0)
                else:
                    ax.plot(cut_angles,array_factor_cut,label=label)
                    ax.set_yticks(range(0,1 , 0.2))
                    ax.set_ylim(0, 1)

        plt.legend()

        return ax



    def plotArrayFactor3D(self,ax=None, logarithmic=True,log_range=-40):
        """ Plot 3D Array Factor

            Args:
                ax: Axis object to plot into
                logarithmic: Plot in dB

            Return:
                ax: Axes object
        """


        if ax==None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')


        if logarithmic ==True:
            array_factor_db = lin2db(self.array_factor)
            array_factor_min =  np.min(array_factor_db)
            array_factor_max =  np.max(array_factor_db)
            e = array_factor_db - (log_range + array_factor_max)
            e[e<0.0] = 0.0
        else:
            e = self.array_factor
            array_factor_min =  np.min(self.array_factor)
            array_factor_max =  np.max(self.array_factor)

        x,y,z = sph2cart(self.phi_mesh,self.theta_mesh,e)

        ax.plot_surface(x, y, z,facecolors=cm.jet(self.array_factor))
        plt.xlabel('X')
        plt.ylabel('Y')
        lim_left,lim_right = ax.get_zlim()
        ax.set_xlim(lim_left,lim_right)
        ax.set_ylim(lim_left,lim_right)

        m = cm.ScalarMappable(cmap=plt.cm.jet)
        if logarithmic==True:
            m.set_array(e+(log_range+array_factor_max))
        else:
            m.set_array(e)
        plt.colorbar(m)

        return ax
