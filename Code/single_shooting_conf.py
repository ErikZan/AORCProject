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
x0 = np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
T = 1.5     # 1,5                  # OCP horizon
dt = 0.05                     # OCP time step
integration_scheme = 'RK-4'
use_finite_difference = False
weight_vel = 0.5           #a=0.5              # cost function weight for final velocity (for position it's implicitely 1)
#weight_u = 1.0                                 # cost function weight for control
weight_run_state= 2300000.0 #a=2300000.0        # np.identity(12) # weight matrix for the state
weight_r= 0.1e-7           #a=0.1e-7           # weight of runnig cost
q_des = np.array([5.0,4.0,5.0,0.0,0.0,0.0])
v_des = np.array([0.0,0.0,0.0,0.0,0.0,0.0])

"""
'a' combination has good result when 

 """