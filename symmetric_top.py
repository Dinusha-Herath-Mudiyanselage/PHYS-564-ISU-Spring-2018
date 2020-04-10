# -*- coding: utf-8 -*-
"""
Created on Sun Jan 21 17:44:51 2018

@author: Dinusha Herath Mudiyanselage
"""

import numpy as np
from scipy.optimize import brenth
from scipy.integrate import quad, odeint

class SymTopIniCond(object):
    """Class defining the initial conditions for the symmetrical top
        
    """
    
    def __init__(self, alpha, beta, a_c , b_c, i_rat):
        """The constructor
        
        :param alpha: alpha constant defined in class (Energy)
        :param beta: beta constant defined in class (gravity)
        :param a_c: a constant, angular momentum around :math::\\Psi
        :param b_c: b constant, angular momentum around :math::\\phi
        :param i_rat: ratio of inertia moments: I/I_3
        """
        
        self.alpha = alpha
        self.beta = beta
        self.a_c = a_c
        self.b_c = b_c
        self.om = i_rat*a_c
       
        # first we will find if the specified conditions 
        # define valid equations of motion
        
        # The function f(u) is a polynomial of degree 3 with coefficients
        c_3 = beta
        c_2 = -(alpha+a_c**2)
        c_1 = 2*a_c*b_c-beta
        c_0 = alpha-b_c**2
        
        self.p_coeffs = np.array([c_3, c_2, c_1, c_0])
        # define the function f(u)
        self.f_top = np.poly1d(self.p_coeffs)
        
        # compute the only point u where the polynomial is a maximum 
        self.max = self.rmax()
        
        # If the maximum does not exist or it is not between -1,1 then stop
        if np.iscomplexobj(self.max) or (np.abs(self.max) > 1):
            print(self.max)
            raise ValueError('There are no solutions for these values')
        
        # If the maximum is NOT positive, then there is only 1 solution
        if self.f_top(self.max) < 0:
            print(self.max, self.f_top(self.max))
            raise ValueError('There is maxima but it is negative')
            
        # find the two points of return
        # they are between -1, max and max, 1
        self.u_l = brenth(self.f_top, -1.0, self.max)
        self.u_r = brenth(self.f_top, self.max, 1.0)
        
        # define u prime u_prime = b/a
        self.u_prime = self.b_c/self.a_c
        
        # compute the corresponding angle
        self.th_1 = np.arccos(self.u_l)
        self.th_2 = np.arccos(self.u_r)
        
        # compute the period of the motion
        def intv(x):
            """Integrand
            """
            return 1/np.sqrt(self.f_top(x))
    
        self.t_half = quad(intv, self.u_l, self.u_r)[0]
        
        # inform about the result
        print('Motion confined between ', self.th_2,' and ' , self.th_1)
        print('Half period is', self.t_half)
        
        # create the function to integrate the equations of motion
        def ode_int(y, t):
            """
            integrand
            """
            aux1 = 1.5*c_3*y[1]**2 + c_2*y[1] + 0.5*c_1
            aux2 = y[0]
            aux3 = (b_c-a_c*y[1])/(1-y[1]**2)
            aux4 = self.om-aux3*y[1]
            
            return np.array([aux1, aux2, aux3, aux4])

        self.ode_int = ode_int
        
        # The solution of the equation with full details will be stored here
        self.ode_sol = False
        
        # time over which the points are solved
        self.time = False
        
    def rmax(self):
        """Returns the root of a quadratic polynomial corresponding to maximum
        
        :return : roots
        :rtype: float
        """
            
        # derivative coefficients
        p = np.array([3, 2, 1])*self.p_coeffs[:-1]
        
        # comptue the radical
        rad = p[1]**2-4*p[0]*p[2]
        
        # if the radical is negative, then there is no solution we return i
        if rad < 0:
            return 1j
        # the free case p[0]=0.0 needs to be treated separately
        elif np.allclose(p[0], 0.0):
            return -p[2]/p[1]
        else:
            w = np.sqrt(rad)
            # roots of the derivative
            x_1 = (-p[1]+w)/(2*p[0])
            x_2 = (-p[1]-w)/(2*p[0])
                
            #This is the second derivative, sign tells wheter max or min
            val = 2*p[0]*x_1 + p[1]
            
            if val < 0:
                return x_1
                
            return x_2
    
    def solv_eqn(self, num_pnts, n_p):
        """
        :param num_pnts: number of points the solutions returns
        :param n_p: number of periods
        
        :return : solutions in the form of [\phi, \theta, \psi]
        :rtype : ndarray
        """
        
        t_end = 2*self.t_half*n_p
        self.time = np.linspace(0.0, t_end, num_pnts)
        # store solutions here
        sols = np.zeros([num_pnts, 3])
        
        # initial conditions
        y0 = np.array([0.0, self.u_r, 0.0, 0.0])
        
        # solve the equation
        solv =  odeint(self.ode_int, y0, self.time)
        
        self.solv_eqn = solv
        
        sols[:, 0] = solv[:, 2]
        sols[:, 1] = np.arccos(solv[:, 1])
        sols[:, 2] = solv[:, 3]

        return sols
