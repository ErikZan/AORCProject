from __future__ import print_function
import numpy as np
import matplotlib as matplot 
import sys
import math
import numpy.linalg as la
from scipy.linalg import expm
from math import sqrt
from scipy.integrate import odeint



from quadcopter_parameter import get_quad_data




 
class Quadcopter
    def __init__(self, dt, params):
        # store motor parameters in member variables
        self.dt  = dt               # simulation time step
        self.l   = param.l          # arm length
        self.k   = param.k          # air lift-force parameter
        self.b = param.b            # viscous friction parameter
        self.I__M = param.I__M      # rotor inertia
        self.Ixx   = param.Ixx      # X-axis inertia
        self.Iyy = param.Iyy        # Y-axis inertia
        self.Izz = param.Izz        # Z-axis inertia
        self.A = param.A            # areodynamic drag coefficient
        
        # set initial states
        self.x = np.zeros(6)        # position and orientation  of drone w.r.t absolute frame
        self.x_dot = np.zeros(6)    # linear and angular velocities of drone
        self.omega = n.zeros(4)     # rotors angular velocities
        self.torque = n.zeros(4)    # rotors input torques
        

    def set_state(self, x):
        self.x = np.copy(x)


#Simulate rotors dynamics
    def simulate_rotors(self, torque):
        self.torque = torque
        
        omega = omega + (torque -b*omega**2)/I__M*dt
#Simulate drone dynamics
    def simulate_drone(self,)
        
    
        

if __name__=='__main__':  
    param = get_quad_data('Quad1')