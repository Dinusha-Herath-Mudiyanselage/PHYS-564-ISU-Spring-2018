# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 12:41:30 2017

@author: trvsst

This module contains functions necessary to evaluate the times it takes for an
object initially at rest to reach a certain point. This is called the 
brachistochrone problem. See textbook, section 2.2
"""

import numpy as np
from scipy import integrate

def time_end(f,g, x_f):
    """function that calculates the time it takes, starting from (0,0) to
       reach the point (x_f, f(x_f))
    
    :param f: curve along the trajectory
    :param g: derivative of the trajectory
    :param x_f: value of x_f at the end
    :return: time it takes (units of ..math :: \\frac{1}{\\sqrt{g}} ) and error
    :rtype: tuple
    """
    
    def ing(x):
        """returns the value of the integrand
        
        :param x: x value
        :return: value of the integrand
        :rtype: float
        """
        return np.sqrt((1+(g(x))**2)/(-2*f(x)))
        
    val, err = integrate.quad(ing, 0, x_f)
    # compute value and error
    
    return (val, err)
