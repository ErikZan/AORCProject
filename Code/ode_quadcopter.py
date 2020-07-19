import numpy as np
from quadcopter_parameter import get_quad_data
import math

class Control_dynamics_1:
    ''' An ordinary differential equation representing a drone system
    '''
    
    def __init__(self, name, param,dt):
        self.g = 9.81
        self.dt  = dt               # simulation time step
        self.m   = param.m                # mass of drone
        self.l   = param.l          # arm length
        self.k   = param.k          # air lift-force parameter
        self.b = param.b            # viscous friction parameter
        self.I__M = param.I__M      # rotor inertia
        self.Ixx   = param.Ixx      # X-axis inertia
        self.Iyy = param.Iyy        # Y-axis inertia
        self.Izz = param.Izz  
        self.nx =12   
        self.nu =4    # Z-axis inertia
        #self.A = param.A            # areodynamic drag coefficient
        
        
    ''' System dynamics
    generated and exported with Maple '''
    def f(self, x, u, t, jacobian=False):
        #nq = 6
        #nv = 6
        self.dx=np.copy(x)
        q = x[:6]
        v = x[6:]
        m=self.m
        Ixx=self.Ixx
        Izz=self.Izz
        Iyy=self.Iyy
        g=self.g
        
        f_dyn = np.array([v[2] * math.sin(q[4]) * math.cos(q[3]) * math.cos(q[5]) + v[1] * math.sin(q[4]) * math.cos(q[5]) * math.sin(q[3]) + v[2] * math.sin(q[5]) * math.sin(q[3]) - v[1] * math.sin(q[5]) * math.cos(q[3]) + math.cos(q[5]) * math.cos(q[4]) * v[0],v[2] * math.sin(q[4]) * math.sin(q[5]) * math.cos(q[3]) + v[1] * math.sin(q[4]) * math.sin(q[5]) * math.sin(q[3]) - v[2] * math.cos(q[5]) * math.sin(q[3]) + v[1] * math.cos(q[3]) * math.cos(q[5]) + math.sin(q[5]) * math.cos(q[4]) * v[0],-math.sin(q[4]) * v[0] + math.cos(q[4]) * math.sin(q[3]) * v[1] + math.cos(q[4]) * math.cos(q[3]) * v[2],(math.cos(q[3]) ** 2 * math.cos(q[4]) * v[3] + math.sin(q[3]) ** 2 * math.cos(q[4]) * v[3] + math.sin(q[4]) * math.cos(q[3]) * v[5] + math.sin(q[4]) * math.sin(q[3]) * v[4]) / math.cos(q[4]) / (math.cos(q[3]) ** 2 + math.sin(q[3]) ** 2),(math.cos(q[3]) * v[4] - math.sin(q[3]) * v[5]) / (math.cos(q[3]) ** 2 + math.sin(q[3]) ** 2),(math.cos(q[3]) * v[5] + math.sin(q[3]) * v[4]) / math.cos(q[4]) / (math.cos(q[3]) ** 2 + math.sin(q[3]) ** 2),-0.1e1 / m * (-g * math.sin(q[4]) * m + v[4] * v[2] - v[5] * v[1]),(-g * math.cos(q[4]) * math.sin(q[3]) * m + v[3] * v[2] - v[5] * v[0]) / m,-0.1e1 / m * (g * math.cos(q[4]) * math.cos(q[3]) * m + v[3] * v[1] - v[4] * v[0] - u[0]),0.1e1 / Ixx * (v[5] * Iyy * v[4] - v[4] * Izz * v[5] + u[1]),-0.1e1 / Iyy * (v[5] * Ixx * v[3] - v[3] * Izz * v[5] - u[2]),0.1e1 / Izz * (v[4] * Ixx * v[3] - v[3] * Iyy * v[4] + u[3])])

        self.dx[:6]  = f_dyn[:6] 
        self.dx[6:] = f_dyn[6:] 
        
        if(jacobian):
            
            self.Fx = np.array([[0,0,0,-v[2] * math.sin(q[4]) * math.sin(q[3]) * math.cos(q[5]) + v[1] * math.sin(q[4]) * math.cos(q[5]) * math.cos(q[3]) + v[2] * math.sin(q[5]) * math.cos(q[3]) + v[1] * math.sin(q[5]) * math.sin(q[3]),v[2] * math.cos(q[4]) * math.cos(q[3]) * math.cos(q[5]) + v[1] * math.cos(q[4]) * math.cos(q[5]) * math.sin(q[3]) - math.cos(q[5]) * math.sin(q[4]) * v[0],-v[2] * math.sin(q[4]) * math.sin(q[5]) * math.cos(q[3]) - v[1] * math.sin(q[4]) * math.sin(q[5]) * math.sin(q[3]) + v[2] * math.cos(q[5]) * math.sin(q[3]) - v[1] * math.cos(q[3]) * math.cos(q[5]) - math.sin(q[5]) * math.cos(q[4]) * v[0],math.cos(q[5]) * math.cos(q[4]),math.sin(q[4]) * math.cos(q[5]) * math.sin(q[3]) - math.sin(q[5]) * math.cos(q[3]),math.sin(q[4]) * math.cos(q[3]) * math.cos(q[5]) + math.sin(q[5]) * math.sin(q[3]),0,0,0],[0,0,0,-v[2] * math.sin(q[4]) * math.sin(q[5]) * math.sin(q[3]) + v[1] * math.sin(q[4]) * math.sin(q[5]) * math.cos(q[3]) - v[2] * math.cos(q[5]) * math.cos(q[3]) - v[1] * math.sin(q[3]) * math.cos(q[5]),v[2] * math.cos(q[4]) * math.sin(q[5]) * math.cos(q[3]) + v[1] * math.cos(q[4]) * math.sin(q[5]) * math.sin(q[3]) - math.sin(q[5]) * math.sin(q[4]) * v[0],v[2] * math.sin(q[4]) * math.cos(q[3]) * math.cos(q[5]) + v[1] * math.sin(q[4]) * math.cos(q[5]) * math.sin(q[3]) + v[2] * math.sin(q[5]) * math.sin(q[3]) - v[1] * math.sin(q[5]) * math.cos(q[3]) + math.cos(q[5]) * math.cos(q[4]) * v[0],math.sin(q[5]) * math.cos(q[4]),math.sin(q[4]) * math.sin(q[5]) * math.sin(q[3]) + math.cos(q[3]) * math.cos(q[5]),math.sin(q[4]) * math.sin(q[5]) * math.cos(q[3]) - math.cos(q[5]) * math.sin(q[3]),0,0,0],[0,0,0,math.cos(q[4]) * math.cos(q[3]) * v[1] - math.cos(q[4]) * math.sin(q[3]) * v[2],-math.cos(q[4]) * v[0] - math.sin(q[4]) * math.sin(q[3]) * v[1] - math.sin(q[4]) * math.cos(q[3]) * v[2],0,-math.sin(q[4]),math.cos(q[4]) * math.sin(q[3]),math.cos(q[4]) * math.cos(q[3]),0,0,0],[0,0,0,(-math.sin(q[4]) * math.sin(q[3]) * v[5] + math.sin(q[4]) * math.cos(q[3]) * v[4]) / math.cos(q[4]) / (math.cos(q[3]) ** 2 + math.sin(q[3]) ** 2),(-math.cos(q[3]) ** 2 * math.sin(q[4]) * v[3] - math.sin(q[3]) ** 2 * math.sin(q[4]) * v[3] + math.cos(q[4]) * math.cos(q[3]) * v[5] + math.cos(q[4]) * math.sin(q[3]) * v[4]) / math.cos(q[4]) / (math.cos(q[3]) ** 2 + math.sin(q[3]) ** 2) + (math.cos(q[3]) ** 2 * math.cos(q[4]) * v[3] + math.sin(q[3]) ** 2 * math.cos(q[4]) * v[3] + math.sin(q[4]) * math.cos(q[3]) * v[5] + math.sin(q[4]) * math.sin(q[3]) * v[4]) / math.cos(q[4]) ** 2 / (math.cos(q[3]) ** 2 + math.sin(q[3]) ** 2) * math.sin(q[4]),0,0,0,0,(math.cos(q[3]) ** 2 * math.cos(q[4]) + math.sin(q[3]) ** 2 * math.cos(q[4])) / math.cos(q[4]) / (math.cos(q[3]) ** 2 + math.sin(q[3]) ** 2),math.sin(q[4]) * math.sin(q[3]) / math.cos(q[4]) / (math.cos(q[3]) ** 2 + math.sin(q[3]) ** 2),math.sin(q[4]) * math.cos(q[3]) / math.cos(q[4]) / (math.cos(q[3]) ** 2 + math.sin(q[3]) ** 2)],[0,0,0,(-math.sin(q[3]) * v[4] - math.cos(q[3]) * v[5]) / (math.cos(q[3]) ** 2 + math.sin(q[3]) ** 2),0,0,0,0,0,0,math.cos(q[3]) / (math.cos(q[3]) ** 2 + math.sin(q[3]) ** 2),-math.sin(q[3]) / (math.cos(q[3]) ** 2 + math.sin(q[3]) ** 2)],[0,0,0,(math.cos(q[3]) * v[4] - math.sin(q[3]) * v[5]) / math.cos(q[4]) / (math.cos(q[3]) ** 2 + math.sin(q[3]) ** 2),(math.cos(q[3]) * v[5] + math.sin(q[3]) * v[4]) / math.cos(q[4]) ** 2 / (math.cos(q[3]) ** 2 + math.sin(q[3]) ** 2) * math.sin(q[4]),0,0,0,0,0,math.sin(q[3]) / math.cos(q[4]) / (math.cos(q[3]) ** 2 + math.sin(q[3]) ** 2),math.cos(q[3]) / math.cos(q[4]) / (math.cos(q[3]) ** 2 + math.sin(q[3]) ** 2)],[0,0,0,0,g * math.cos(q[4]),0,0,0.1e1 / m * v[5],-0.1e1 / m * v[4],0,-0.1e1 / m * v[2],0.1e1 / m * v[1]],[0,0,0,-g * math.cos(q[4]) * math.cos(q[3]),g * math.sin(q[4]) * math.sin(q[3]),0,-0.1e1 / m * v[5],0,v[3] / m,0.1e1 / m * v[2],0,-v[0] / m],[0,0,0,g * math.cos(q[4]) * math.sin(q[3]),g * math.sin(q[4]) * math.cos(q[3]),0,0.1e1 / m * v[4],-v[3] / m,0,-0.1e1 / m * v[1],v[0] / m,0],[0,0,0,0,0,0,0,0,0,0,0.1e1 / Ixx * (Iyy * v[5] - Izz * v[5]),0.1e1 / Ixx * (Iyy * v[4] - Izz * v[4])],[0,0,0,0,0,0,0,0,0,-0.1e1 / Iyy * (Ixx * v[5] - Izz * v[5]),0,-0.1e1 / Iyy * (Ixx * v[3] - Izz * v[3])],[0,0,0,0,0,0,0,0,0,0.1e1 / Izz * (Ixx * v[4] - Iyy * v[4]),0.1e1 / Izz * (Ixx * v[3] - Iyy * v[3]),0]])
            self.Fu = np.array([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0.1e1 / m,0,0,0],[0,0.1e1 / Ixx,0,0],[0,0,0.1e1 / Iyy,0],[0,0,0,0.1e1 / Izz]])
            return (np.copy(self.dx), np.copy(self.Fx), np.copy(self.Fu))
        
        return np.copy(self.dx)
        
                
    def f_x_fin_diff(self, x, u, t, delta=1e-8):
        ''' Partial derivatives of system dynamics w.r.t. x computed via finite differences '''
        f0 = self.f(x, u, t)
        Fx = np.zeros((self.nx, self.nx))
        for i in range(self.nx):
            xp = np.copy(x)
            xp[i] += delta
            fp = self.f(xp, u, t)
            Fx[:,i] = (fp-f0)/delta
        return Fx
        
        
    def f_u_fin_diff(self, x, u, t, delta=1e-8):
        ''' Partial derivatives of system dynamics w.r.t. u computed via finite differences '''
        f0 = self.f(x, u, t)
        Fu = np.zeros((self.nx, self.nu))
        for i in range(self.nu):
            up = np.copy(u)
            up[i] += delta
            fp = self.f(x, up, t)
            Fu[:,i] = (fp-f0)/delta
        return Fu