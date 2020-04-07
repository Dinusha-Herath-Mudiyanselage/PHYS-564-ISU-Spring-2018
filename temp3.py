# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
def eulerAnglesToRotationMatrix(self,angle) :
    
    self.phi= angle[0]
    self.theta=angle[1]
    self.psi=angle[2]
    
     
    R_x = np.array([[1,         0,                  0                   ],
                    [0,         np.cos(self.phi), -np.sin(self.phi) ],
                    [0,         np.sin(self.phi), np.cos(self.phi)  ]
                    ])
         
         
                     
    R_y = np.array([[np.cos(self.theta),    0,      np.sin(self.theta)  ],
                    [0,                     1,      0                   ],
                    [-np.sin(self.theta),   0,      np.cos(self.theta)  ]
                    ])
                 
    R_z = np.array([[np.cos(self.psi),    -np.sin(self.psi),    0],
                    [np.sin(self.psi),    np.cos(self.psi),     0],
                    [0,                     0,                      1]
                    ])
                     
                     
    R = np.dot(R_z, np.dot( R_y, R_x ))
 
    return R