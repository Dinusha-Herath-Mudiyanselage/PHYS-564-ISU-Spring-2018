# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 12:22:44 2017

@author: Dinusha Herath Mudiyanselage
This module contains functions necessary to evaluate the times it takes for an
object initially at rest to reach a certain point. This is called the 
brachistochrone problem. See textbook, section 2.2. This uses the class func.
"""

import numpy as np
from scipy import integrate

def time_end(w, x_f):
    """function that calculates the time it takes, starting from (0,0) to
       reach the point (x_f, f(x_f))
    
    :param w: curve along the trajectory defined from a FuncDevs object 
    :param x_f: value of x_f at the end
    :return: time it takes (units of ..math :: \\frac{1}{\\sqrt{g}} ) and error
    :rtype: tuple
    """
    
    if np.logical_not(w.is_in_interval(x_f)):
        raise ValueError('Proposed x value is not in the function domain')
        
    def ing(x):
        """returns the value of the integrand
        
        :param x: x value
        :return: value of the integrand
        :rtype: float
        """
        return np.sqrt((1+(w.f_prime(x))**2)/(-2*w.f(x)))
        
    val, err = integrate.quad(ing, 0, x_f)
    # compute value and error
    
    return (val, err)
