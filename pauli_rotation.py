# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 23:17:01 2018

@author: Dinusha Herath Mudiyanselage
"""
import numpy as np

class PauliRotation(object):
    
    def rotate(self,n,theta):
        
        self.n_x = n[0]
        self.n_y = n[1]
        self.n_z = n[2]
        
        self.theta = theta
        
        if self.n_x != 0:
            
            self.n_x = np.array([[1,0],[0,1]])
        
        else :
            
            self.n_x = np.array([[0,0],[0,0]])
        
        
        if self.n_y != 0:
            
            self.n_y = np.array([[1,0],[0,1]])
        
        else :
            
            self.n_y = np.array([[0,0],[0,0]])
            
            
            
        if self.n_z != 0:
            
            self.n_x = np.array([[1,0],[0,1]])
        
        else :
            
            self.n_z = np.array([[0,0],[0,0]])    
            
            
            
            
        sigma_x = np.array([[0, 1],[ 1, 0]])    
        sigma_y = np.array([[0,-1j],[1j, 0]])  
        sigma_z = np.array([[1, 0],[0, -1]])    
        
        
        b = np.cos(self.theta/2)*np.array([[1,0],[0,1]])
        c = (np.dot(self.n_x,sigma_x)+np.dot(self.n_y,sigma_y)+np.dot(self.n_z,sigma_z))*np.sin(self.theta/2)
        a = b+c*1j
        
        return a
    