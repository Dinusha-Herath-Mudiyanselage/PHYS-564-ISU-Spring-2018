# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 08:29:31 2017

@author: trvsst

class defining the trajectories of an inverse power law potential
"""

import numpy as np
from scipy.integrate import quad
from scipy import optimize

class OrbitGravity(object):
    """
    Trajectories of an inverse power law potential
    """
    
    def __init__(self, exc, theta_0=0, eps=1e-5):
        """
        The constructor
        
        :param exc: excentricity
        :param theta_0: angle at which the turning point occurs
        :param eps: small value to avoid singular values for open orbits
        """
        
        self.e = exc
        self.theta_0 = theta_0
        self.eps = eps
        
        if exc == 0.0:
            self.theta_max = np.pi
            self.theta_min = - np.pi
            self.trajectory = 'Circle'
        
        elif exc < 1:
            self.theta_max = np.pi
            self.theta_min = - np.pi
            self.trajectory = 'Ellipse'
        elif exc==1:
            self.theta_max = np.pi - eps
            self.theta_min = - np.pi + eps
            self.trajectory = 'Parabola'
        else:
            self.theta_max = self.theta_0 + np.arccos(-1/self.e) - eps
            self.theta_min = self.theta_0 - np.arccos(-1/self.e) -eps
            self.trajectory = 'Hyperbole'
        
        if self.e < 1:
            self.period = 2*np.pi/np.sqrt(1-self.e**2)**3
        else:
            self.period = np.inf
    
    def radial(self, theta):
        """Returns the radial position of the trajectory
        
        :param theta: angles to compute the trajectory
        :return: trajectory values
        :rtype: ndarray
        """
        
        theta = np.asarray(theta)
        
        return 1.0/(1+self.e*np.cos(theta-self.theta_0))

    def coordinate(self, time):
        """Returns the coordinates as a function of time
        
        :param time: numpy array with the times 
        :return : positons at the times ginve (in polar coordinates)
        :rtype: ndarray
        """
        
        def intg(x):
            return 1/(1+self.e*np.cos(x))**2
            
        def theta_f(th, val):
            res = quad(intg, 0.0, th)
            # note that the trajectory starts at the perihelion
            return res[0]-val
        
        theta = np.zeros_like(time)
        num_turns = np.zeros_like(time)
        for ind, x in enumerate(time):
            z = x
            a_max = self.theta_max
            if self.e < 1:
                z = np.remainder(z, self.period)
                num_turns[ind] = np.floor(x/self.period)
                a_max = 2*np.pi
            theta[ind] = optimize.brenth(theta_f, 0.0, a_max, args=(z,))
        
        return np.array([theta, self.radial(theta), num_turns])
