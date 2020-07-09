from __future__ import print_function
import numpy as np
import matplotlib as matplot 
import plot_utils as plut
import sys
import math
import numpy.linalg as la
from scipy.linalg import expm
from math import sqrt
from scipy.integrate import odeint

from quadcopter_parameter import get_quad_data

class Quadcopter:
    def __init__(self, dt, params):
        # store motor parameters in member variables
        self.g = 9.81
        self.dt  = dt               # simulation time step
        self.m   = param.m                # mass of drone
        self.l   = param.l          # arm length
        self.k   = param.k          # air lift-force parameter
        self.b = param.b            # viscous friction parameter
        self.I__M = param.I__M      # rotor inertia
        self.Ixx   = param.Ixx      # X-axis inertia
        self.Iyy = param.Iyy        # Y-axis inertia
        self.Izz = param.Izz        # Z-axis inertia
        self.A = param.A            # areodynamic drag coefficient
        
        # set initial states
        self.x = np.zeros(6,dtype=float)        # position and orientation  of drone w.r.t absolute frame
        self.x_dot = np.zeros(6,dtype=float)    # linear and angular velocities of drone
        self.omega = np.zeros(4,dtype=float)     # rotors angular velocities
        self.torque = np.zeros(4,dtype=float)    # rotors input torques
        

    def set_state(self, x,x_dot,omega):
        self.x = np.copy(x)
        self.x_dot=np.copy(x_dot)
        self.omega=np.copy(omega)

    #Simulate rotors dynamics
    def simulate_rotors(self, torque):
        
        I__M = self.I__M
        b = self.b
        dt = self.dt
    
        self.torque = torque
        
        self.omega = self.omega + (torque -b*self.omega**2*np.sign(self.omega))/I__M*dt
        
    #Simulate drone dynamics
    def simulate_drone(self,omega):
        self.omega = omega
        
        g=self.g
        dt=self.dt
        m= self.m   
        l= self.l   
        k= self.k   
        Ixx= self.Ixx   
        Iyy= self.Iyy 
        Izz= self.Izz
        A = self.A         
        
        Ax = np.array([0.1e1 / m * (math.sin(self.x[5]) * math.sin(self.x[3]) + math.cos(self.x[5]) * math.sin(self.x[4]) * math.cos(self.x[3])) * (0.1e1 / (omega[0] + 0.1e-7) * abs(omega[0] + 0.1e-7) * omega[0] ** 2 * k - 0.1e1 / (omega[1] + 0.1e-7) * abs(omega[1] + 0.1e-7) * omega[1] ** 2 * k + 0.1e1 / (omega[2] + 0.1e-7) * abs(omega[2] + 0.1e-7) * omega[2] ** 2 * k - 0.1e1 / (omega[3] + 0.1e-7) * abs(omega[3] + 0.1e-7) * omega[3] ** 2 * k) - A / m * self.x_dot[0],0.1e1 / m * (-math.cos(self.x[5]) * math.sin(self.x[3]) + math.sin(self.x[5]) * math.sin(self.x[4]) * math.cos(self.x[3])) * (0.1e1 / (omega[0] + 0.1e-7) * abs(omega[0] + 0.1e-7) * omega[0] ** 2 * k - 0.1e1 / (omega[1] + 0.1e-7) * abs(omega[1] + 0.1e-7) * omega[1] ** 2 * k + 0.1e1 / (omega[2] + 0.1e-7) * abs(omega[2] + 0.1e-7) * omega[2] ** 2 * k - 0.1e1 / (omega[3] + 0.1e-7) * abs(omega[3] + 0.1e-7) * omega[3] ** 2 * k) - A / m * self.x_dot[1],-g + 0.1e1 / m * math.cos(self.x[4]) * math.cos(self.x[3]) * (0.1e1 / (omega[0] + 0.1e-7) * abs(omega[0] + 0.1e-7) * omega[0] ** 2 * k - 0.1e1 / (omega[1] + 0.1e-7) * abs(omega[1] + 0.1e-7) * omega[1] ** 2 * k + 0.1e1 / (omega[2] + 0.1e-7) * abs(omega[2] + 0.1e-7) * omega[2] ** 2 * k - 0.1e1 / (omega[3] + 0.1e-7) * abs(omega[3] + 0.1e-7) * omega[3] ** 2 * k) - A / m * self.x_dot[2]])
        
        self.x_dot[:3] += Ax*dt

        
        self.x[:3] = self.x[:3] + self.x_dot[:3]*dt   
        
        Jinv = np.mat([[(Iyy ** 2 * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 + Iyy ** 2 * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - 2 * Iyy * Izz * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 - 2 * Iyy * Izz * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + Izz ** 2 * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 + Izz ** 2 * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + Ixx * Iyy * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Ixx * Izz * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Iyy ** 2 * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + 2 * Iyy * Izz * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Izz ** 2 * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Ixx * Iyy * math.cos(self.x[3]) ** 2 + Ixx * Izz * math.cos(self.x[3]) ** 2 + Ixx * Izz * math.cos(self.x[4]) ** 2 - Iyy * Izz * math.cos(self.x[4]) ** 2 - Ixx * Izz) / Ixx / (Iyy ** 2 * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 + Iyy ** 2 * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - 2 * Iyy * Izz * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 - 2 * Iyy * Izz * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + Izz ** 2 * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 + Izz ** 2 * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + Ixx * Iyy * math.sin(self.x[4]) ** 2 * math.cos(self.x[3]) ** 2 + Ixx * Iyy * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Ixx * Izz * math.sin(self.x[4]) ** 2 * math.cos(self.x[3]) ** 2 - Ixx * Izz * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Iyy ** 2 * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + 2 * Iyy * Izz * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Izz ** 2 * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Ixx * Iyy * math.cos(self.x[3]) ** 2 + Ixx * Izz * math.sin(self.x[4]) ** 2 + Ixx * Izz * math.cos(self.x[3]) ** 2 + Ixx * Izz * math.cos(self.x[4]) ** 2 - Iyy * Izz * math.cos(self.x[4]) ** 2 - Ixx * Izz),math.sin(self.x[4]) * (Iyy - Izz) * math.cos(self.x[3]) * math.sin(self.x[3]) * math.cos(self.x[4]) / (Iyy ** 2 * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 + Iyy ** 2 * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - 2 * Iyy * Izz * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 - 2 * Iyy * Izz * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + Izz ** 2 * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 + Izz ** 2 * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + Ixx * Iyy * math.sin(self.x[4]) ** 2 * math.cos(self.x[3]) ** 2 + Ixx * Iyy * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Ixx * Izz * math.sin(self.x[4]) ** 2 * math.cos(self.x[3]) ** 2 - Ixx * Izz * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Iyy ** 2 * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + 2 * Iyy * Izz * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Izz ** 2 * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Ixx * Iyy * math.cos(self.x[3]) ** 2 + Ixx * Izz * math.sin(self.x[4]) ** 2 + Ixx * Izz * math.cos(self.x[3]) ** 2 + Ixx * Izz * math.cos(self.x[4]) ** 2 - Iyy * Izz * math.cos(self.x[4]) ** 2 - Ixx * Izz),-math.sin(self.x[4]) * (math.cos(self.x[3]) ** 2 * Iyy - math.cos(self.x[3]) ** 2 * Izz + Izz) / (Iyy ** 2 * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 + Iyy ** 2 * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - 2 * Iyy * Izz * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 - 2 * Iyy * Izz * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + Izz ** 2 * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 + Izz ** 2 * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + Ixx * Iyy * math.sin(self.x[4]) ** 2 * math.cos(self.x[3]) ** 2 + Ixx * Iyy * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Ixx * Izz * math.sin(self.x[4]) ** 2 * math.cos(self.x[3]) ** 2 - Ixx * Izz * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Iyy ** 2 * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + 2 * Iyy * Izz * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Izz ** 2 * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Ixx * Iyy * math.cos(self.x[3]) ** 2 + Ixx * Izz * math.sin(self.x[4]) ** 2 + Ixx * Izz * math.cos(self.x[3]) ** 2 + Ixx * Izz * math.cos(self.x[4]) ** 2 - Iyy * Izz * math.cos(self.x[4]) ** 2 - Ixx * Izz)],[math.sin(self.x[4]) * (Iyy - Izz) * math.cos(self.x[3]) * math.sin(self.x[3]) * math.cos(self.x[4]) / (Iyy ** 2 * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 + Iyy ** 2 * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - 2 * Iyy * Izz * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 - 2 * Iyy * Izz * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + Izz ** 2 * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 + Izz ** 2 * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + Ixx * Iyy * math.sin(self.x[4]) ** 2 * math.cos(self.x[3]) ** 2 + Ixx * Iyy * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Ixx * Izz * math.sin(self.x[4]) ** 2 * math.cos(self.x[3]) ** 2 - Ixx * Izz * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Iyy ** 2 * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + 2 * Iyy * Izz * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Izz ** 2 * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Ixx * Iyy * math.cos(self.x[3]) ** 2 + Ixx * Izz * math.sin(self.x[4]) ** 2 + Ixx * Izz * math.cos(self.x[3]) ** 2 + Ixx * Izz * math.cos(self.x[4]) ** 2 - Iyy * Izz * math.cos(self.x[4]) ** 2 - Ixx * Izz),(Iyy * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Izz * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + math.sin(self.x[4]) ** 2 * Ixx + Ixx * math.cos(self.x[4]) ** 2 - Iyy * math.cos(self.x[4]) ** 2 - Ixx) / (Iyy ** 2 * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 + Iyy ** 2 * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - 2 * Iyy * Izz * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 - 2 * Iyy * Izz * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + Izz ** 2 * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 + Izz ** 2 * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + Ixx * Iyy * math.sin(self.x[4]) ** 2 * math.cos(self.x[3]) ** 2 + Ixx * Iyy * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Ixx * Izz * math.sin(self.x[4]) ** 2 * math.cos(self.x[3]) ** 2 - Ixx * Izz * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Iyy ** 2 * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + 2 * Iyy * Izz * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Izz ** 2 * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Ixx * Iyy * math.cos(self.x[3]) ** 2 + Ixx * Izz * math.sin(self.x[4]) ** 2 + Ixx * Izz * math.cos(self.x[3]) ** 2 + Ixx * Izz * math.cos(self.x[4]) ** 2 - Iyy * Izz * math.cos(self.x[4]) ** 2 - Ixx * Izz),(Iyy - Izz) * math.cos(self.x[3]) * math.sin(self.x[3]) * math.cos(self.x[4]) / (Iyy ** 2 * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 + Iyy ** 2 * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - 2 * Iyy * Izz * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 - 2 * Iyy * Izz * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + Izz ** 2 * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 + Izz ** 2 * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + Ixx * Iyy * math.sin(self.x[4]) ** 2 * math.cos(self.x[3]) ** 2 + Ixx * Iyy * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Ixx * Izz * math.sin(self.x[4]) ** 2 * math.cos(self.x[3]) ** 2 - Ixx * Izz * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Iyy ** 2 * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + 2 * Iyy * Izz * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Izz ** 2 * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Ixx * Iyy * math.cos(self.x[3]) ** 2 + Ixx * Izz * math.sin(self.x[4]) ** 2 + Ixx * Izz * math.cos(self.x[3]) ** 2 + Ixx * Izz * math.cos(self.x[4]) ** 2 - Iyy * Izz * math.cos(self.x[4]) ** 2 - Ixx * Izz)],[-math.sin(self.x[4]) * (math.cos(self.x[3]) ** 2 * Iyy - math.cos(self.x[3]) ** 2 * Izz + Izz) / (Iyy ** 2 * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 + Iyy ** 2 * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - 2 * Iyy * Izz * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 - 2 * Iyy * Izz * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + Izz ** 2 * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 + Izz ** 2 * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + Ixx * Iyy * math.sin(self.x[4]) ** 2 * math.cos(self.x[3]) ** 2 + Ixx * Iyy * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Ixx * Izz * math.sin(self.x[4]) ** 2 * math.cos(self.x[3]) ** 2 - Ixx * Izz * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Iyy ** 2 * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + 2 * Iyy * Izz * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Izz ** 2 * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Ixx * Iyy * math.cos(self.x[3]) ** 2 + Ixx * Izz * math.sin(self.x[4]) ** 2 + Ixx * Izz * math.cos(self.x[3]) ** 2 + Ixx * Izz * math.cos(self.x[4]) ** 2 - Iyy * Izz * math.cos(self.x[4]) ** 2 - Ixx * Izz),(Iyy - Izz) * math.cos(self.x[3]) * math.sin(self.x[3]) * math.cos(self.x[4]) / (Iyy ** 2 * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 + Iyy ** 2 * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - 2 * Iyy * Izz * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 - 2 * Iyy * Izz * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + Izz ** 2 * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 + Izz ** 2 * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + Ixx * Iyy * math.sin(self.x[4]) ** 2 * math.cos(self.x[3]) ** 2 + Ixx * Iyy * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Ixx * Izz * math.sin(self.x[4]) ** 2 * math.cos(self.x[3]) ** 2 - Ixx * Izz * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Iyy ** 2 * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + 2 * Iyy * Izz * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Izz ** 2 * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Ixx * Iyy * math.cos(self.x[3]) ** 2 + Ixx * Izz * math.sin(self.x[4]) ** 2 + Ixx * Izz * math.cos(self.x[3]) ** 2 + Ixx * Izz * math.cos(self.x[4]) ** 2 - Iyy * Izz * math.cos(self.x[4]) ** 2 - Ixx * Izz),-(math.cos(self.x[3]) ** 2 * Iyy - math.cos(self.x[3]) ** 2 * Izz + Izz) / (Iyy ** 2 * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 + Iyy ** 2 * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - 2 * Iyy * Izz * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 - 2 * Iyy * Izz * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + Izz ** 2 * math.cos(self.x[3]) ** 4 * math.cos(self.x[4]) ** 2 + Izz ** 2 * math.cos(self.x[3]) ** 2 * math.sin(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + Ixx * Iyy * math.sin(self.x[4]) ** 2 * math.cos(self.x[3]) ** 2 + Ixx * Iyy * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Ixx * Izz * math.sin(self.x[4]) ** 2 * math.cos(self.x[3]) ** 2 - Ixx * Izz * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Iyy ** 2 * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 + 2 * Iyy * Izz * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Izz ** 2 * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) ** 2 - Ixx * Iyy * math.cos(self.x[3]) ** 2 + Ixx * Izz * math.sin(self.x[4]) ** 2 + Ixx * Izz * math.cos(self.x[3]) ** 2 + Ixx * Izz * math.cos(self.x[4]) ** 2 - Iyy * Izz * math.cos(self.x[4]) ** 2 - Ixx * Izz)]])

        tau_b = np.array([(0.1e1 / (self.omega[1] + 0.1e-7) * abs(self.omega[1] + 0.1e-7) * self.omega[1] ** 2 * k - 0.1e1 / (self.omega[3] + 0.1e-7) * abs(self.omega[3] + 0.1e-7) * self.omega[3] ** 2 * k) * l,(-0.1e1 / (self.omega[0] + 0.1e-7) * abs(self.omega[0] + 0.1e-7) * self.omega[0] ** 2 * k + 0.1e1 / (self.omega[2] + 0.1e-7) * abs(self.omega[2] + 0.1e-7) * self.omega[2] ** 2 * k) * l,self.torque[0] + self.torque[1] + self.torque[2] + self.torque[3]])
        
        C_eta = np.array([-self.x_dot[5] * self.x_dot[4] * math.cos(self.x[4]) * Ixx - (Iyy - Izz) * (-self.x_dot[4] ** 2 * math.cos(self.x[3]) * math.sin(self.x[3]) + 2 * (math.cos(self.x[3]) ** 2 - 0.1e1 / 0.2e1) * math.cos(self.x[4]) * self.x_dot[5] * self.x_dot[4] + self.x_dot[5] ** 2 * math.cos(self.x[3]) * math.sin(self.x[3]) * math.cos(self.x[4]) ** 2),-(Iyy - Izz) * (((-2 * math.cos(self.x[3]) ** 2 * math.cos(self.x[4]) + math.cos(self.x[4])) * self.x_dot[5] + 2 * self.x_dot[4] * math.cos(self.x[3]) * math.sin(self.x[3])) * self.x_dot[3] + math.sin(self.x[4]) * self.x_dot[5] * self.x_dot[4] * math.cos(self.x[3]) * math.sin(self.x[3])) - self.x_dot[5] * (math.sin(self.x[4]) * ((Iyy - Izz) * math.cos(self.x[3]) ** 2 + Ixx - Iyy) * math.cos(self.x[4]) * self.x_dot[5] - math.sin(self.x[4]) * math.cos(self.x[3]) * math.sin(self.x[3]) * (Iyy - Izz) * self.x_dot[4] - Ixx * self.x_dot[3] * math.cos(self.x[4])),-self.x_dot[4] ** 2 * math.sin(self.x[4]) * math.cos(self.x[3]) * math.sin(self.x[3]) * (Iyy - Izz) + 2 * (((Iyy - Izz) * math.cos(self.x[3]) ** 2 - Ixx / 2 - Iyy / 2 + Izz / 2) * self.x_dot[3] + self.x_dot[5] * math.sin(self.x[4]) * ((Iyy - Izz) * math.cos(self.x[3]) ** 2 + Ixx - Iyy)) * math.cos(self.x[4]) * self.x_dot[4] + 2 * math.cos(self.x[4]) ** 2 * self.x_dot[3] * self.x_dot[5] * math.cos(self.x[3]) * math.sin(self.x[3]) * (Iyy - Izz)])

        self.x_dot[3:6] = self.x_dot[3:6] + Jinv.dot(np.transpose(tau_b-C_eta))*dt
        self.x[3:6] = self.x[3:6] + self.x_dot[3:6]*dt
        
    def q(self):
        return self.x
    def p(self):
        return self.omega
    
if __name__=='__main__':  
    
    param = get_quad_data('Quad1')
    
    dt = 1e-4
    T = 1
    N = int(T/dt)  #number of time steps
    q = np.zeros((6,N+1))
    o = np.zeros((4,N+1))
    Torque = np.zeros((4,N+1))
    drone = Quadcopter(dt,param)
    
    #set initial conditions
    
    drone.set_state(np.array([0.0,0.0,1.0,0.0,0.0,0.0]),np.array([0.0,0.0,0.0,0.0,0.0,0.0]),np.array([620.6107625,-620.6107625,620.6107625,-620.6107625])) 
    q[:6,0] = np.array([0,0,1,0,0,0])
    
    # simulate rotor and drone dynamics
    
    for i in range(N+1):
        Torque[:4,i]= np.array([0.04390797988,-0.043907979880,0.04390797988,-0.04390797988])
        drone.simulate_rotors(Torque[:4,i])
        drone.simulate_drone(drone.p())
        q[:6,i] = drone.q()
     
        
    #Plot stuffs
    
    f, ax = plut.create_empty_figure(1)
    time = np.arange(0.0, T+dt, dt)
    time = time[:N+1]
    ax.plot(time, q[0,:N+1], label ='x')
    ax.legend()
    matplot.pyplot.xlabel('Time [s]')
    

    print("Final height", q[0,-1])
    
    f, ax = plut.create_empty_figure(1)
    time = np.arange(0.0, T+dt, dt)
    time = time[:N+1]
    ax.plot(time, q[1,:N+1], label ='y')
    ax.legend()
    matplot.pyplot.xlabel('Time [s]')
    

    print("Final height", q[1,-1])
    
    f, ax = plut.create_empty_figure(1)
    time = np.arange(0.0, T+dt, dt)
    time = time[:N+1]
    ax.plot(time, q[2,:N+1], label ='z')
    ax.legend()
    matplot.pyplot.xlabel('Time [s]')
    

    print("Final height", q[2,-1])
    
    f, ax = plut.create_empty_figure(1)
    time = np.arange(0.0, T+dt, dt)
    time = time[:N+1]
    ax.plot(time, q[3,:N+1], label ='phi')
    ax.legend()
    matplot.pyplot.xlabel('Time [s]')
    

    print("Final height", q[3,-1])
    
    f, ax = plut.create_empty_figure(1)
    time = np.arange(0.0, T+dt, dt)
    time = time[:N+1]
    ax.plot(time, q[4,:N+1], label ='theta')
    ax.legend()
    matplot.pyplot.xlabel('Time [s]')


    print("Final height", q[4,-1])
    
    f, ax = plut.create_empty_figure(1)
    time = np.arange(0.0, T+dt, dt)
    time = time[:N+1]
    ax.plot(time, q[5,:N+1], label ='psi')
    ax.legend()
    matplot.pyplot.xlabel('Time [s]')
    matplot.pyplot.show()

    print("Final height", q[5,-1])
    
    # f, ax = matplot.pyplot.subplots(2,1,sharex=True)
    # time = np.arange(0.0, T+dt, dt)
    # time = time[:N+1]
    # ax[0].plot(time, q[0], label ='x')
    # ax[0].plot(time, q[1], label ='y')
    # ax[0].plot(time, q[2], label ='z')
    # ax[1].plot(time, q[3], label ='phi')
    # ax[1].plot(time, q[4], label ='theta')
    # ax[1].plot(time, q[5], label ='psi')
    # for i in range(2): ax[i].legend()
    # matplot.pyplot.xlabel('Time [s]')
    # matplot.pyplot.show()
