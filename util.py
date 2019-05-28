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

from scipy.misc import factorial

""" Conversion """
def deg2rad(angle):
    return angle*math.pi/180
def rad2deg(angle):
    return angle*180/math.pi

def lin2db(value_lin):
    return 20*np.log10(value_lin)
def db2lin(value_db):
    return 10**(value_db/20)
def pow2db(value_pow):
    return 10*np.log10(value_pow)
def db2pow(value_db):
    return 10**(value_db/10)

def cart2sph(x,y,z):
    r = np.sqrt(x**2 + y**2 + z**2)
    phi = np.arctan2(y,x)
    theta = np.arccos(z/r)
    return phi, theta, r

def sph2cart(phi,theta,r):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z

""" Data Treatment """
def normalize(x):
    return (x - x.min()) / (np.ptp(x))


""" Polynomials """
def binomial(x, y):
    """ Calculate Binomial Coefficient """
    try:
        binom = math.factorial(x) // math.factorial(y) // math.factorial(x - y)
    except ValueError:
        binom = 0
    return binom

def dolphTschebychev(n_elements,R):
    """ Calculate Weights based on Tschebychev Polynomial

        Args:
            n_elements: Number of elements in array
            R: Sidelobe-level (voltage ratio)

        Return:
            w: polynomial weights
    """
    if R <= 0:
        sys.exit('Sidelobe Level has to be positive!')

    N = math.floor(n_elements/2)

    beta = (np.cosh(np.arccosh(R)/(n_elements-1)))**2
    alpha = 1 - 1/beta

    # Even Case
    if n_elements == 2*N:
        nend = N - 1
        w = np.zeros(N)
    else:
        nend = N
        w = np.zeros(N+1)

    w[0] = 1
    for n in range(1,nend+1):
        p = 1
        for m in range(1,n):
            f_m = m * (n_elements-1-2*n + m) / ((n-m) * (n+1-m))
            p = p * alpha * f_m + 1

        w[n] = (n_elements-1)*alpha * p

    if n_elements == 2*N:
        w = np.concatenate((w,np.flipud(w)))
    else:
        w = np.concatenate((w[:-1],np.flipud(w)))


    return w

