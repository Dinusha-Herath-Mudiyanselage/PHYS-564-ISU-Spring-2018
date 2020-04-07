# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 08:29:31 2017

@author: trvsst

class defining the trajectories of an inverse power law potential
"""

import numpy as np

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
            self.theta_min = - np.pi - eps
            self.trajectory = 'Parabola'
        else:
            self.theta_max = self.theta_0 + np.arccos(-1/self.e) - eps
            self.theta_min = self.theta_0 - np.arccos(-1/self.e) -eps
            self.trajectory = 'Hyperbole'
    
    def radial(self, theta):
        """Returns the radial position of the trajectory
        
        :param theta: angles to compute the trajectory
        :return: trajectory values
        :rtype: ndarray
        """
        
        theta = np.asarray(theta)
        
        return 1.0/(1+self.e*np.cos(theta-self.theta_0))
