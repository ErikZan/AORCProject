from __future__ import print_function
import numpy as numpy
import matplotlib as matplot 
import sys
import math
import numpy.linalg as la
from scipy.linalg import expm
from math import sqrt



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
        
        
#complex parameter
#def get_quad_matrices(name):
#   if(name=='Quad1'):
'''         J = numpy.mat([[Ixx,0,-math.sin(theta) * Ixx],[0,(Iyy - Izz) * math.cos(phi) ** 2 + Izz,math.cos(phi) * math.sin(phi) * math.cos(theta) * (Iyy - Izz)],[-math.sin(theta) * Ixx,math.cos(phi) * math.sin(phi) * math.cos(theta) * (Iyy - Izz),((-Iyy + Izz) * math.cos(phi) ** 2 - Ixx + Iyy) * math.cos(theta) ** 2 + Ixx]])
        C_eta = numpy.mat([-psi_dot * theta_dot * math.cos(theta) * Ixx - (-theta_dot ** 2 * math.cos(phi) * math.sin(phi) + 2 * psi_dot * math.cos(theta) * (math.cos(phi) ** 2 - 0.1e1 / 0.2e1) * theta_dot + psi_dot ** 2 * math.cos(phi) * math.sin(phi) * math.cos(theta) ** 2) * (Iyy - Izz),-(((-2 * math.cos(phi) ** 2 * math.cos(theta) + math.cos(theta)) * psi_dot + 2 * theta_dot * math.cos(phi) * math.sin(phi)) * phi_dot + math.sin(theta) * psi_dot * theta_dot * math.cos(phi) * math.sin(phi)) * (Iyy - Izz) - psi_dot * (math.cos(theta) * math.sin(theta) * ((Iyy - Izz) * math.cos(phi) ** 2 + Ixx - Iyy) * psi_dot - math.sin(theta) * math.cos(phi) * math.sin(phi) * (Iyy - Izz) * theta_dot - Ixx * phi_dot * math.cos(theta)),-theta_dot ** 2 * math.sin(theta) * math.cos(phi) * math.sin(phi) * (Iyy - Izz) + 2 * math.cos(theta) * (((Iyy - Izz) * math.cos(phi) ** 2 - Ixx / 2 - Iyy / 2 + Izz / 2) * phi_dot + psi_dot * math.sin(theta) * ((Iyy - Izz) * math.cos(phi) ** 2 + Ixx - Iyy)) * theta_dot + 2 * math.cos(theta) ** 2 * phi_dot * psi_dot * math.cos(phi) * math.sin(phi) * (Iyy - Izz)])
        taub = numpy.mat([(-omega2 ** 2 + omega4 ** 2) * k * l,(-omega1 ** 2 + omega3 ** 2) * k * l,0.1e1 / (omega1 + 0.1e-7) * abs(omega1 + 0.1e-7) * omega1 ** 2 * b + sympy.diff(omega1, t) * I__M + 0.1e1 / (omega2 + 0.1e-7) * abs(omega2 + 0.1e-7) * omega2 ** 2 * b + sympy.diff(omega2, t) * I__M + 0.1e1 / (omega3 + 0.1e-7) * abs(omega3 + 0.1e-7) * omega3 ** 2 * b + sympy.diff(omega3, t) * I__M + 0.1e1 / (omega4 + 0.1e-7) * abs(omega4 + 0.1e-7) * omega4 ** 2 * b + sympy.diff(omega4, t) * I__M])
        eq2 = numpy.mat([(math.sin(psi) * math.sin(phi) + math.cos(psi) * math.sin(theta) * math.cos(phi)) * (0.1e1 / (omega1 + 0.1e-7) * abs(omega1 + 0.1e-7) * omega1 ** 2 * k - 0.1e1 / (omega2 + 0.1e-7) * abs(omega2 + 0.1e-7) * omega2 ** 2 * k + 0.1e1 / (omega3 + 0.1e-7) * abs(omega3 + 0.1e-7) * omega3 ** 2 * k - 0.1e1 / (omega4 + 0.1e-7) * abs(omega4 + 0.1e-7) * omega4 ** 2 * k) - A / m * x_dot,(-math.cos(psi) * math.sin(phi) + math.sin(psi) * math.sin(theta) * math.cos(phi)) * (0.1e1 / (omega1 + 0.1e-7) * abs(omega1 + 0.1e-7) * omega1 ** 2 * k - 0.1e1 / (omega2 + 0.1e-7) * abs(omega2 + 0.1e-7) * omega2 ** 2 * k + 0.1e1 / (omega3 + 0.1e-7) * abs(omega3 + 0.1e-7) * omega3 ** 2 * k - 0.1e1 / (omega4 + 0.1e-7) * abs(omega4 + 0.1e-7) * omega4 ** 2 * k) - A / m * y_dot,-g + math.cos(theta) * math.cos(phi) * (0.1e1 / (omega1 + 0.1e-7) * abs(omega1 + 0.1e-7) * omega1 ** 2 * k - 0.1e1 / (omega2 + 0.1e-7) * abs(omega2 + 0.1e-7) * omega2 ** 2 * k + 0.1e1 / (omega3 + 0.1e-7) * abs(omega3 + 0.1e-7) * omega3 ** 2 * k - 0.1e1 / (omega4 + 0.1e-7) * abs(omega4 + 0.1e-7) * omega4 ** 2 * k) - A / m * z_dot])
 '''
