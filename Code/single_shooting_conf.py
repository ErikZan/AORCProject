# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 09:47:07 2019

@author: student
"""

import numpy as np
import os
from math import sqrt

np.set_printoptions(precision=3, linewidth=200, suppress=True)
LINE_WIDTH = 60
x0 = np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]) # Initial State
T = 3.0     # OCP horizon
dt = 0.15   # OCP time step

integration_scheme = 'RK-4'
use_finite_difference = False
weight_vel = 0.0001          # cost function weight for final velocity (for position it's implicitely 1)
    
weight_r= 1.0e-2             # weight of runnig cost
q_des = np.array([4.0,0.0,0.0,0.0,0.0,0.0])
v_des = np.array([0.0,0-.0,0.0,0.0,0.0,0.0])
weight_run_state = np.array([
        [0.0, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.], # weight of x cost
       [0., 0.0, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],  # weight of y cost
       [0., 0., 0.0, 0., 0., 0., 0., 0., 0., 0., 0., 0.],  # weight of z cost
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],   # weight of phi cost
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],   # weight of theta cost
       [0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.],   # weight of psi cost
       [0., 0., 0., 0., 0., 0., 0.0, 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0.0, 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.0, 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.0, 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,0.0 , 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.0]]) *1e-15

weight_const = 0.001 # cost for constraint : not used now
