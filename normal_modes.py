# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 12:01:18 2018

@author: Dinusha Herath Mudiyanselage
"""

import numpy as np
from numpy import linalg as la

class NormalModes(object):
    """
    Class computing normal modes
    """
    
    def __init__(self, mat_kin, mat_vib):
        """The constructor
        
        :param mat_kin: kinetic energy matrix
        :param mat_vib: vibration matrix
        """
        
        self.kin = mat_kin
        self.vib = mat_vib
        
        dig_t = np.sqrt(mat_kin)
        
        diag = la.inv(dig_t)
        
        tot = np.dot(diag, np.dot(mat_vib, diag))
        
        eig_val, eig_vec = la.eigh(tot)
        
        # omega 
        self.omega = np.sqrt(eig_val)
        
        # make sure the eigenvectors are normalized
        eig_vec = eig_vec/la.norm(eig_vec, axis=0)
        
        # transformation from normal to original
        self.o_2_n = np.dot(eig_vec.T, dig_t)
        
        # trasnformation from original to normal
        self.n_2_o = la.inv(self.o_2_n)
        
    def original_2_normal(self, vec):
        """
        returns given coordinates into normal 
        
        :param vec: given cooordinates
        :return : normal coordinates
        :rtype : ndarray
        """
        
        return np.dot(self.o_2_n, vec)
        
    def normal_2_original(self, vec):
        """
        returns normal coordinates into given 
        
        :param vec: normal cooordinates
        :return : given coordinates
        :rtype : ndarray
        """
        
        return np.dot(self.n_2_o, vec)
     
