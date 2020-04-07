# -*- coding: utf-8 -*-
"""
Created on Sat Feb 17 13:26:57 2018

@author: HerathMudiyanselage Dinusha


"""


import numpy as np

class euler_Rotation(object) :
      def _init_(self,angle) :
    
            self.phi   = angle[0]
            self.theta = angle[1]
            self.psi   = angle[2]
    
     
       
                     
            B =   np.array([[ np.cos(self.psi),  np.sin(self.psi)  ,    0],
                    [-np.sin(self.psi),  np.cos(self.theta),    0],
                    [ 0               ,  0                 ,    1]
                    ])
    
    
            C =    np.array([[1,   0                 , 0                 ],
                     [0,   np.cos(self.theta), np.sin(self.theta)],
                     [0,  -np.sin(self.theta), np.cos(self.theta)]
                     ])
     
                 
            D =    np.array([[ np.cos(self.phi),    np.sin(self.phi),   0],
                      [-np.sin(self.phi),    np.cos(self.phi),   0],
                     [0                ,    0               ,   1]
                     ])
                     
                     
            A = np.dot(B, np.dot( C, D ))
 
            return A
