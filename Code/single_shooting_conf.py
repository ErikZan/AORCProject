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
T = 1.5                       # OCP horizon
dt = 0.001                     # OCP time step
integration_scheme = 'RK-4'
use_finite_difference = False
weight_vel = 1e-1000 # cost function weight for final velocity (for position it's implicitely 1)
weight_u = 1e-3  # cost function weight for control
q_des = np.array([0.0,0.0,100.0,0.0,0.0,0.0])
v_des = np.array([0.0,0.0,0.0,0.0,0.0,0.0])
dp_des=np.array([0.0,0.0,0.0,0.0,0.0,0.0])



