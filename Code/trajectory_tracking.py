import numpy as np
from numpy.linalg import norm
from scipy.optimize import minimize
#from ode import ODERobot
from ode_quadcopter import Control_dynamics_1 as ode
from numerical_integration_quad import Integrator
import math
import single_shooting_conf as conf
from quadcopter_parameter import get_quad_data
import pandas as pd 
import matplotlib as matplot
import plot_utils as plut
import quadcopter_simu_model
import csv
from numpy import genfromtxt

if __name__=='__main__':
  
  T=conf.T
  dt=conf.dt
  N = int(T/dt) 
  
  # retrive control input
  my_data_u = genfromtxt('U_control.csv', delimiter=',')
  U_ctrl= my_data_u[1:N+1,1:5]
  #retrive state 
  my_data_x = genfromtxt('optimizationresult.csv', delimiter=',')
  X_state= my_data_x[1:N+1,1:13]
  
  

