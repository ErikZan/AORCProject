from __future__ import print_function
import numpy as np
import matplotlib as matplot 
import sys
import math
import numpy.linalg as la
from scipy.linalg import expm
from math import sqrt
from scipy.integrate import odeint



#import quadcopter_parameter 


class Empty:
    def __init__(self):
        pass

def  get_quad_data(name):
    param = Empty()
    if(name=='Quad1'):
        param.m = 0.468         # motor inertia
        param.l = 0.225
        param.k = 2.98e-06
        param.b = 1.140e-07 
        param.I__M = 3.357e-05
        param.Ixx = 4.85e-03
        param.Iyy = 4.85e-03
        param.Izz = 8.802e-03  
        param.A = 0.25
        
    return param
  
  