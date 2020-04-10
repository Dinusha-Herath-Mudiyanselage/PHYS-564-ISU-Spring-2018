# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 11:56:43 2017

@author: Dinusha Herath Mudiyanselage

This module defines a function class that puts together the function and its 
derivatives. 

"""
import numpy as np
from scipy.misc import derivative

class FuncDevs(object):
    """
    Class function
    """
    
    def __init__(self, func, nm, intv = (-np.inf, np.inf), xtest = False):
        """
        The constructor
        
        :param func: function list
        :param intv: interval where the function and derivatives are defined 
        
        """
        
        self.num_ders = len(func)
        self.intval = intv
        self.name = nm
        
        if isinstance(xtest, np.ndarray):
            self.xtest = xtest
        else:
            self.xtest = np.array([1.0, 2.0, 3.0, 4.0])
            # define some default values
        
        self.fun = func
        
        # check that derivatives are computed correctly
        dx = 1e-6
        tol = 1e-5
        for val in self.xtest:
            if np.logical_and(val > intv[0], val < intv[1]):
                der = derivative(self.fun[0], val, dx)
                if np.logical_not(np.allclose(der, self.fun[1](val), atol = tol)):
                    print(val, ' ', self.fun[1](val), ' ',  der)
                    raise ValueError('Derivative is not well calculated')
                    
    def f(self, x):
        """Function
        
        :param x: value at x
        :return : value of the function
        :rtype : float
        """
        
        return self.fun[0](x)
        
    def f_prime(self, x):
        """Derivative
        
        :param x: value at x
        :return : value of the function derivative
        :rtype : float
        """
        
        return self.fun[1](x)

    def is_in_interval(self, x):
        """Checks whether a point x is in the domain of the function
        
        param x: value of x
        :return : True if the point x is in the interval intv
        :rtype: bool
        """
        
        bool_intval = False

        if np.logical_and(x > self.intval[0], x < self.intval[1]):
            bool_intval = True
            
        return bool_intval
                
        
