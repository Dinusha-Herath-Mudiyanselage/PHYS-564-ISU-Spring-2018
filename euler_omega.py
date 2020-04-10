# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 13:02:44 2018

@author: Dinusha Herath Mudiyanselage
"""
import numpy as np

def omega_body(a_e, v_e):
    """
    Computes the omega vector in the body frame
    
    :param e_euler: Euler angles :math:`(\\phi, \\theta, \\psi)`
    :param v_euler: Derivative of Euler angles :math:`(\\phi, \\theta, \\psi)`
    :return: omega in the body frame
    :rtype: ndarrray
    """
    
    omega = np.zeros_like(a_e)
    
    omega[0] = np.sin(a_e[2])*np.sin(a_e[1])*v_e[0] + np.cos(a_e[2])*v_e[1]
    omega[1] = np.cos(a_e[2])*np.sin(a_e[1])*v_e[0] - np.sin(a_e[2])*v_e[1]
    omega[2] = v_e[2]+ np.cos(a_e[1])*v_e[0]
    
    return omega

def omega_reference(a_e, v_e):
    """
    Computes the omega vector in the reference frame
    
    :param e_euler: Euler angles :math:`(\\phi, \\theta, \\psi)`
    :param v_euler: Derivative of Euler angles :math:`(\\phi, \\theta, \\psi)`
    :return: omega in the body frame
    :rtype: ndarrray
    """
    
    omega = np.zeros_like(a_e)
    
    omega[0] = np.sin(a_e[0])*np.sin(a_e[1])*v_e[2] + np.cos(a_e[0])*v_e[1]
    omega[1] = -np.cos(a_e[0])*np.sin(a_e[1])*v_e[2] + np.sin(a_e[0])*v_e[1]
    omega[2] = v_e[0]+ np.cos(a_e[1])*v_e[2]
    
    return omega
