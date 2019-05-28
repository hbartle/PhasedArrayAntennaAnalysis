#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Hannes Bartle <mail@hbartle.de>
#
# Distributed under terms of the MIT license.

from math import cos, sin, sqrt

from util import *
from antenna_toolbox.Antenna import Antenna

class RectangularPatch(Antenna):
    """ Rectangular Patch Antenna """

    def __init__(self,epsilon_r,substrate_height,frequency,rolloff=1):
        self.epsilon_r = epsilon_r
        self.substrate_height = substrate_height
        self.frequency = frequency

        self.getDimensions()



        self.rolloff_factor = rolloff

    def getDimensions(self):
        """ Calculate Standard Rectangular Patch Dimensions (Balanis,4th Ed.,14.2"""

        e_r = self.epsilon_r
        h = self.substrate_height
        f = self.frequency

        # Permittivity of free space
        epsilon0 = 8.854185e-12

        wavelength = 3e8 / f
        wavelength_adapt = wavelength / sqrt(e_r)

        # Patch Witdh
        W = (3e8 / (2 * f)) * sqrt(2 / (e_r + 1))

        e_r_eff = ((e_r + 1) / 2) + ((e_r - 1) / 2) * (1 + 12 * (h / W)) ** -0.5

        F1 = (e_r_eff + 0.3) * (W / h + 0.264)
        F2 = (e_r_eff - 0.258) * (W / h + 0.8)
        dL = h * 0.412 * (F1 / F2)

        wavelength_adapt = wavelength / sqrt(e_r_eff)
        self.L = (wavelength_adapt/ 2) - 2 * dL
        self.substrate_height_eff = h * sqrt(e_r)
        self.W = W
        self.epsilon_r_eff = e_r_eff

        return self.W,self.L,self.epsilon_r_eff,self.substrate_height_eff


    def calcFieldValue(self,theta_in,phi_in):
        """ Calculate far-field normalized pattern value """

        # Wavenumber
        wavelength = 3e8 / self.frequency
        k = 2 * math.pi / wavelength

        # Calculate Far-Field Points and rotate coordinate system 90 degrees
        xff, yff, zff = sph2cart(phi_in, theta_in, 999)
    #    #xffd = zff
    #    #yffd = yff
    #    #zffd = xff
    #    xffd = xff
    #    yffd = yff
    #    zffd = zff
    #    phi, theta, r = cart2sph(xffd, yffd, zffd)
    #
    #    if theta == 0:
    #        theta = 1e-9
    #    if phi == 0:
    #        phi = 1e-9
    #
    #    # Calculate Normalized Far-Field Pattern
    #    X = k*self.W/2 *sin(theta)*sin(phi)
    #    Fphi = -sin(X)/X*cos(k*self.L*sin(theta)*cos(phi))*cos(theta)*sin(phi)
    #    Ftheta = sin(X)/X*cos(k*self.L*sin(theta)*cos(phi))*cos(phi)
    #
    #
    #    # Smooth pattern using roll-off
    #    theta_deg = theta_in * 180 / math.pi
    #    F1 = 1 / (((self.rolloff_factor * (abs(theta_deg) - 90)) ** 2) + 0.001)
    #    PatEdgeSF = 1 / (F1 + 1)
    #
    #    UNF = 1.0006


        # Balanis Method
        xffd = zff
        yffd = yff
        zffd = xff
        phi, theta, r = cart2sph(xffd, yffd, zffd)

        if theta == 0:
            theta = 1e-9
        if phi == 0:
            phi = 1e-9



        X = k * self.substrate_height_eff / 2 * sin(theta)*cos(phi)
        Z = k * self.W / 2 * cos(theta)

        if theta_in <= math.pi / 2:
            Ftotal = sin(theta)*sin(X)/X * sin(Z)/Z * cos(k * self.L / 2 * sin(theta) *sin(phi))
        else:
            Ftotal = 1e-9



        return Ftotal


    def getField(self,theta_start,theta_stop,theta_step,phi_start,phi_stop,phi_step):
        """ Calculate Normalized Pattern for angular range """
        self.Ntheta = (theta_stop-theta_start)/theta_step + 1
        self.Nphi = (phi_stop-phi_start)/phi_step + 1
        self.theta_step = theta_step
        self.phi_step = phi_step
        self.theta_start = theta_start
        self.theta_stop = theta_stop
        self.phi_start = phi_start
        self.phi_stop = phi_stop

        theta_range = deg2rad(np.linspace(theta_start,theta_stop,self.Ntheta))
        phi_range = deg2rad(np.linspace(phi_start,phi_stop,self.Nphi))
        self.pattern = np.ones((theta_range.size,phi_range.size))
        self.theta_mesh = np.ones((theta_range.size,phi_range.size))
        self.phi_mesh = np.ones((theta_range.size,phi_range.size))



        for t,theta in enumerate(theta_range):
            for p,phi in enumerate(phi_range):
                self.pattern[t,p] = self.calcFieldValue(theta,phi)
                self.theta_mesh[t,p] = theta
                self.phi_mesh[t,p] = phi

        return self.pattern



