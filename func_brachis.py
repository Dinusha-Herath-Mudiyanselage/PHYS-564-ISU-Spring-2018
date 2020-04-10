# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 13:56:53 2017

@author: Dinusha Herath Mudiyanselage
"""

import numpy as np
from scipy import optimize

class BrachisFuncs(object):
    """
    Class function to get the Brachistochrone functions
    """
    
    def __init__(self, x_f, y_f):
        """
        The constructor
        
        We assume that the curve starts at the origin (x_i=0, y_i=0)
        
        :param x_f: final x coordinate
        :param y_f: final y coordinate
        
        """
        # obtain the initial conditions

        self.x_f = x_f
        self.y_f = y_f
        
        v = self.x_var(x_f, y_f)
        
        self.a_val = v[0]
        self.theta_f = v[1]

        self.intval = np.array([0, x_f])

    def f(self, x):
        """Function
        
        :param x: value at x
        :return : value of the function
        :rtype : float
        """
    
        x = np.asarray(x)
        
        # a float is read as 0-dimensional array, need to resize it to 1-dim
        if x.ndim == 0:
            x = x.reshape(1)
        
        vals = np.zeros_like(x)

        for ind, v in enumerate(x):
            theta = self.theta(v)
            vals[ind] = -self.y_cycloid(theta)
        
        return vals
        
    def f_prime(self, x):
        """Derivative
        
        :param x: value at x
        :return : value of the function derivative
        :rtype : float
        """
        
        x = np.asarray(x)
        
        # a float is read as 0-dimensional array, need to resize it to 1-dim
        if x.ndim == 0:
            x = x.reshape(1)    
        
        der = np.zeros_like(x)
        
        for ind, v in enumerate(x):
            theta = self.theta(v)
            der[ind] = -np.sin(theta)/(1-np.cos(theta))
        
        return der

    def is_in_interval(self, x):
        """Checks whether a point x is in the domain of the function
        
        param x: value of x
        :return : True if the point x is in the interval intv
        :rtype: bool
        """
        
        bool_intval = False

        if np.logical_and(x >= self.intval[0], x <= self.intval[1]):
            bool_intval = True
        else:
            print(x, self.intval)
            
        return bool_intval
    
    def x_cycloid(self, theta):
        return 0.5*self.a_val*(theta-np.sin(theta))
    
    def y_cycloid(self, theta):
        return 0.5*self.a_val*(1-np.cos(theta))
    
    def theta(self, x, eps=1e-14):
        """returns the value of theta for a given x
        
        :param x: x coordinate
        :return theta: theta value
        :param eps: precision value
        :rtype: float
        """
    
        if x < eps:
            return 0.0
        
        def fun(z):
            return self.x_cycloid(z)-x
        
        return optimize.brenth(fun, eps, 2*np.pi-eps)
    
    @staticmethod   
    def x_var(x, y, eps=1e-4):
        """ 
        Returns the value of a, theta for arbitrary boundary conditions
        
        :param x: x coordinate
        :param y: y coordinate
        :param eps: precision value
        :return : value of (a, theta)
        :rtype: ndarray
        """   
        
        r = y/x
        
        if x < eps:
            theta = 2*x/y
        else:
            def ini_cond(z):
                return r-(1-np.cos(z))/(z - np.sin(z))
    
            theta = optimize.brenth(ini_cond, eps, 2*np.pi-eps)
        
        a_val =  2*y/(1-np.cos(theta))
        
        return np.array([a_val, theta])
