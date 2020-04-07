# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 21:16:35 2018

@author: Dinusha Herath Mudiyanselage
"""
import numpy as np

class EulerRotationMatrix(object):
    
    def rotate(self,angle):
        
        self.phi   = angle[0]
        self.theta = angle[1]
        self.psi   = angle[2]
        
        
        D = np.array([[ np.cos(self.phi),  np.sin(self.phi),    0],
                      [-np.sin(self.phi),  np.cos(self.phi),    0],
                      [                0,                 0,    1]
                     ])
          
        C = np.array([[1,                   0,                  0],
                      [0,  np.cos(self.theta), np.sin(self.theta)],
                      [0, -np.sin(self.theta), np.cos(self.theta)]
                     ])
    
        B = np.array([[ np.cos(self.psi),  np.sin(self.psi),    0],
                      [-np.sin(self.psi),  np.cos(self.psi),    0],
                      [                0,                 0,    1]
                     ])
    
        A = np.dot(B,np.dot(C,D))
        
        return A