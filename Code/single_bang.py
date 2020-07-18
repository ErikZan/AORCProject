# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 18:12:04 2020

@author: student
"""

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
import os
class Empty:
    def __init__(self):
        pass
    

class SingleShootingProblem:
    ''' A simple solver for a single shooting OCP.
        In the current version, the solver considers only a cost function (no path constraints).
    '''
    
    def __init__(self, name, ode, x0, dt, N, integration_scheme):
        self.name = name
        self.ode = ode
        self.integrator = Integrator('integrator')
        self.x0 = x0
        self.dt = dt
        self.N = N
        self.integration_scheme = integration_scheme
               
        self.nq = int(x0.shape[0]/2)
        self.nx = x0.shape[0]
        self.nu = 4
        self.X = np.zeros((N, self.x0.shape[0]))
        self.U = np.zeros((N, 6))
        self.last_cost = 0.0
        self.running_costs = []
        self.final_costs = []
       #Constraints -------------------------------------------------------------------- 
        self.U_cons = np.zeros((N, 6))
        self.X_cons = np.zeros((N, self.x0.shape[0]))
        self.dXdU_cons = np.array([])

        self.bound_phi = 1.05 # 0.78 -> 45 deg
        self.grad_phi = np.array([])

        self.bound_theta = 1.05 # 0.78 -> 45 deg
        self.grad_theta = np.array([])
        
        self.dist_wind =2.0
        self.offset = 0.1
        
        self.position_w_y = 0.0
        self.size_y =1.0
        
        self.position_w_z = 3.0
        self.size_z =1.0
        self.grad_line_dwn = np.array([])
        
        # Constraint cost 
        


        
    def add_running_cost(self, c, weight=1):
        self.running_costs += [(weight,c)]
    
    def add_final_cost(self, c, weight=1):
        self.final_costs += [(weight,c)]
        
        
    def running_cost(self, X, U):
        ''' Compute the running cost integral '''
        cost = 0.0
        t = 0.0
        for i in range(U.shape[0]):
            for (w,c) in self.running_costs:
                cost += w * dt * c.compute(X[i,:], U[i,:], t, recompute=True)
                t += self.dt
        return cost
        
    def running_cost_w_gradient(self, X, U, dXdU):
        ''' Compute the running cost integral and its gradient w.r.t. U'''
        cost = 0.0
        grad = np.zeros(self.N*self.nu)
        t = 0.0
        nx, nu = self.nx, self.nu
        for i in range(U.shape[0]):
            for (w,c) in self.running_costs:
                ci, ci_x, ci_u = c.compute_w_gradient(X[i,:], U[i,:], t, recompute=True)
                dci = ci_x.dot(dXdU[i*nx:(i+1)*nx,:]) 
                dci[i*nu:(i+1)*nu] += ci_u
                
                cost += w * self.dt * ci
                grad += w * self.dt * dci
                t += self.dt
        return (cost, grad)
        
    def final_cost(self, x_N):
        ''' Compute the final cost '''
        cost = 0.0
        for (w,c) in self.final_costs:
            cost += w * c.compute(x_N, recompute=True)
        return cost
        
    def final_cost_w_gradient(self, x_N, dxN_dU):
        ''' Compute the final cost and its gradient w.r.t. U'''
        cost = 0.0
        grad = np.zeros(self.N*self.nu)
        for (w,c) in self.final_costs:
            ci, ci_x = c.compute_w_gradient(x_N, recompute=True)
            dci = ci_x.dot(dxN_dU)
            cost += w * ci
            grad += w * dci
        return (cost, grad)
        
    def compute_cost(self, y):
        ''' Compute cost function '''
        # compute state trajectory X from control y
        U = y.reshape((self.N, self.nu))
        t0, ndt = 0.0, 1
        X = self.integrator.integrate(self.ode, self.x0, U, t0, self.dt, ndt, 
                                      self.N, self.integration_scheme)
        
        # compute cost
        run_cost = self.running_cost(X, U)
        fin_cost = self.final_cost(X[-1,:])
       
        cost = run_cost + fin_cost
        # store X, U and cost
        self.X, self.U = X, U
        self.last_cost = cost        
        return cost
        
    def compute_cost_w_gradient_fd(self, y):
        ''' Compute both the cost function and its gradient using finite differences '''
        eps = 1e-8
        y_eps = np.copy(y)
        grad = np.zeros_like(y)
        cost = self.compute_cost(y)
        for i in range(y.shape[0]):
            y_eps[i] += eps
            cost_eps = self.compute_cost(y_eps)
            y_eps[i] = y[i]
            grad[i] = (cost_eps - cost) / eps
        return (cost, grad)
        
    def compute_cost_w_gradient(self, y):
        ''' Compute cost function and its gradient '''
        # compute state trajectory X from control y
        U = y.reshape((self.N, self.nu))
        t0 = 0.0
        X, dXdU = self.integrator.integrate_w_sensitivities_u(self.ode, self.x0, U, t0, 
                                                        self.dt, self.N, 
                                                        self.integration_scheme)
        
        # compute cost
        (run_cost, grad_run) = self.running_cost_w_gradient(X, U, dXdU)
        (fin_cost, grad_fin) = self.final_cost_w_gradient(X[-1,:], dXdU[-self.nx:,:])
       
        cost = run_cost + fin_cost 
        grad = grad_run + grad_fin 
        
        # store X, U and cost
        self.X, self.U = X, U
        self.last_cost = cost        
        return (cost, grad)
        
    
    def solve_bounds(self, y0=None, method='BFGS', use_finite_difference=False,bnds=None):
        ''' Solve the optimal control problem '''
        # if no initial guess is given => initialize with zeros
        if(y0 is None):
            y0 = np.zeros(self.N*self.nu)
        
        self.iter = 0
        print('Start optimizing')
        if(use_finite_difference):
            r = minimize(self.compute_cost_w_gradient_fd, y0, jac=True, method=method, # 
                     callback=self.clbk, options={'maxiter': 20, 'disp': True},bounds=bnds) # cons not implemented 
        else:
            r = minimize(self.compute_cost_w_gradient, y0, jac=True, method=method, 
                     callback=self.clbk, options={'maxiter': 50, 'disp': True },bounds=bnds,
                     constraints=(
                         {'type':'ineq','fun': self.fun_cons_phi}, # ,'jac':self.jac_cons_phi
                                  {'type':'ineq','fun': self.fun_cons_theta},
                                   {'type':'ineq','fun': self.fun_cons_phi}, # ,'jac':self.jac_cons_theta
                                   {'type':'ineq','fun': self.cons_line_up},
                                   {'type':'ineq','fun': self.cons_line_dwn}  #  ,'jac':self.jac_cons_line_dwn                  np.tile( ,(self.N)) ,,'jac':lambda x : np.array([1.0,0.0,0.0,0.0])
                                  ))
        return r
    
    def solve(self, y0=None, method='BFGS', use_finite_difference=False):
        ''' Solve the optimal control problem '''
        # if no initial guess is given => initialize with zeros
        if(y0 is None):
            y0 = np.zeros(self.N*self.nu)
            
        self.iter = 0
        print('Start optimizing')
        if(use_finite_difference):
            r = minimize(self.compute_cost_w_gradient_fd, y0, jac=True, method=method, 
                     callback=self.clbk, options={'maxiter': 200, 'disp': True}) 
        else:
            r = minimize(self.compute_cost_w_gradient, y0, jac=True, method=method, 
                     callback=self.clbk, options={'maxiter': 200, 'disp': True})
        return r
        
    def sanity_check_cost_gradient(self, N_TESTS=10):
        ''' Compare the gradient computed with finite differences with the one
            computed by deriving the integrator
        '''
        for i in range(N_TESTS):
            y = np.random.rand(self.N*self.nu)
            (cost, grad_fd) = self.compute_cost_w_gradient_fd(y)
            (cost, grad) = self.compute_cost_w_gradient_fd(y)
            grad_err = grad-grad_fd
            if(np.max(np.abs(grad_err))>1e-4):
                print('Grad:   ', grad)
                print('Grad FD:', grad_fd)
            else:
                print('Everything is fine', np.max(np.abs(grad_err)))
        
    def clbk(self, xk):
        print('Iter %3d, cost %5f'%(self.iter, self.last_cost))
        self.iter += 1
       
        file = pd.DataFrame(self.X).to_csv(f"/home/test/Desktop/Desktop/GitAORC/AORCProject/Code/stored_trajectory/file{self.iter}.csv")
        np.savetxt(f"/home/test/Desktop/Desktop/GitAORC/AORCProject/Code/stored_trajectory/file{self.iter}.txt", pd.DataFrame(self.X).values, fmt='%f')
        #os.system("./home/test/Desktop/Desktop/GitAORC/AORCProject/Code/stored_trajectory/plot_stuff.sh")
        #file = True
        return file

   #Bounds on phi ----------------
    def fun_cons_phi(self, y ):

        self.compute_dyn_cons(y) 
        X = self.X_cons
        U = self.U_cons
        dXdU =self.dXdU_cons                                              

        nx = self.nx
        nu = self.nu
        N = self.N

        cons = np.array([])
        grad = np.zeros(self.N*self.nu)
        
        for i in range(U.shape[0]):

            ci =  self.bound_phi - np.absolute(X[i,3])
            grad += np.array([0.0,0.0,0.0,-1.0*np.sign(X[i,3]),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]).dot(dXdU[i*nx:(i+1)*nx,:])
            cons = np.append(cons,ci)
        self.grad_phi = grad
        return cons

    def jac_cons_phi(self,y):

        return self.grad_phi


    #Bounds on theta ----------------
    def fun_cons_theta(self, y ):
        
        self.compute_dyn_cons(y) 
        X = self.X_cons
        U = self.U_cons
        dXdU =self.dXdU_cons  

        nx = self.nx
        nu = self.nu
        N = self.N                                            
    
        cons = np.array([])
        grad = np.zeros(self.N*self.nu)

        for i in range(U.shape[0]):
            ci =  self.bound_theta - np.absolute(X[i,4])
            grad += np.array([0.0,0.0,0.0,0.0,-1.0*np.sign(X[i,4]),0.0,0.0,0.0,0.0,0.0,0.0,0.0]).dot(dXdU[i*nx:(i+1)*nx,:])
            cons = np.append(cons,ci)
        self.grad_theta = grad
        return cons

    def fun_cons_psi(self, y ):
        
        self.compute_dyn_cons(y) 
        X = self.X_cons
        U = self.U_cons
        dXdU =self.dXdU_cons  

        nx = self.nx
        nu = self.nu
        N = self.N                                            
    
        cons = np.array([])
        grad = np.zeros(self.N*self.nu)

        for i in range(U.shape[0]):
            ci =  self.bound_theta - np.absolute(X[i,5])
            grad += np.array([0.0,0.0,0.0,0.0,0.0,-1.0*np.sign(X[i,5]),0.0,0.0,0.0,0.0,0.0,0.0]).dot(dXdU[i*nx:(i+1)*nx,:])
            cons = np.append(cons,ci)
        self.grad_theta = grad
        return cons
    
    def jac_cons_theta(self,y):

        return self.grad_theta
    
    def ywindows(self,y):
        U = y.reshape((self.N, self.nu))
        t0 = 0.0
        X, dXdU = self.integrator.integrate_w_sensitivities_u(self.ode, self.x0, U, t0, 
                                                        self.dt, self.N, 
                                                        self.integration_scheme)
        cons = np.array([])
        grad = np.array([])
        
        for i in range(U.shape[0]):
            if X[i,0] > self.dist_wind-self.offset and X[i,0] < self.dist_wind+self.offset:
            #if X[i,0] > self.dist_wind-self.offset:
                cry = -np.absolute(X[i,1] - self.position_w_y) + self.size_y/2
                cons = np.append(cons,cry)
        return cons
    
    def zwindows(self,y):
        U = y.reshape((self.N, self.nu))
        t0 = 0.0
        X, dXdU = self.integrator.integrate_w_sensitivities_u(self.ode, self.x0, U, t0, 
                                                        self.dt, self.N, 
                                                        self.integration_scheme)
        cons = np.array([])
        grad = np.array([])
        grad1 = np.array([1.0,0.0,0.0,0.0])
        for i in range(U.shape[0]):
            if X[i,0] > self.dist_wind-self.offset and X[i,0] < self.dist_wind+self.offset:
            #if X[i,0] > self.dist_wind-self.offset:
                cuz = -np.absolute(X[i,2] - self.position_w_z) + self.size_z/2
                cons = np.append(cons,cuz)
                grad = np.append(grad,grad1)
                
        return cons
    
    def compute_dyn_cons(self,y):
        if np.array_equal(y,self.U_cons):

            return False
        else:
            self.U_cons = y.reshape((self.N, self.nu))
            t0 = 0.0
            self.X_cons, self.dXdU_cons = self.integrator.integrate_w_sensitivities_u(self.ode, self.x0, self.U_cons, t0, 
                                                        self.dt, self.N, 
                                                        self.integration_scheme)
        return 0
     
    def cons_line_dwn(self,y):
        
        self.compute_dyn_cons(y)
        X = self.X_cons
        U = self.U_cons
        dXdU = self.dXdU_cons 

        nx = self.nx
        nu = self.nu
        N = self.N

        mz = ((self.position_w_z-self.size_z/2.0)-(self.x0[2]-5.0))/(self.dist_wind-self.x0[0])
        
        cons = np.array([])
        grad = np.zeros(self.N*self.nu)
        
        for i in range(U.shape[0]):
            if X[i,0] < self.dist_wind:
                cuz = X[i,2] - (mz*X[i,0] + self.x0[2]-5.0) # prova ad aggiungere un uno
                cons = np.append(cons,cuz)
                grad += np.array([-mz,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]).dot(dXdU[i*nx:(i+1)*nx,:])
            else:
                #cuz = X[i,2] - (mz*X[i,0] + x0[2])
                #cuz = X[i,2]-((self.position_w_z-self.size_z/2) -100*(X[i,0])-self.dist_wind)
                cuz = 0
                cons = np.append(cons,cuz)
                grad += np.array([-100.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]).dot(dXdU[i*nx:(i+1)*nx,:])
        self.grad_line_dwn = grad  
        return cons
    def jac_cons_line_dwn(self,y):

        return self.grad_line_dwn

    def cons_line_up(self,y):
         
        self.compute_dyn_cons(y)
        X = self.X_cons
        U = self.U_cons
        
        mz = ((self.position_w_z+self.size_z/2)-(self.x0[2]+5.0))/(self.dist_wind-self.x0[0])
        
        cons = np.array([])
        grad = np.array([])
        
        for i in range(U.shape[0]):
            if X[i,0] < self.dist_wind:
                cuz = -X[i,2] + mz*X[i,0] + (self.x0[2]+5.0)
                cons = np.append(cons,cuz)
            else:
                cuz = 0
                cons = np.append(cons,cuz)
                
        return cons

    '''     def cost_on_cons(self,X):
        cost = 0.0
        grad = np.zeros(self.N*self.nu)
        t = 0.0
        line = self.create_line(X)
        w = 1.0
        nu =4
        nx=12
        for i in range(X.shape[0]):
            ci, ci_x, ci_u = c.compute_cost_constraint(line[i])
            dci = ci_x.dot(self.dXdU_cons[i*nx:(i+1)*nx,:]) 
            dci[i*nu:(i+1)*nu] += ci_u
                
            cost += w * self.dt * ci
            grad += w * self.dt * dci
            t += self.dt
        return (cost, grad)
     '''
   

        
if __name__=='__main__':
    import plot_utils as plut
    import matplotlib.pyplot as plt
    import time
    from cost_functions_quad import OCPFinalCostState, OCPRunningCostQuadraticControl , OCPRunningConstraint , OCPFinalCostLength
    import single_shooting_conf as conf
    import os
    import numpy as np
    # delete previous csv 
    
    os.system('./home/test/Desktop/Desktop/GitAORC/AORCProject/Code/stored_trajectory/delete_csv.sh')
   # os.system('./stored_trajectory/delete_csv.sh')
   
    np.set_printoptions(precision=3, linewidth=200, suppress=True)
    
    dt = conf.dt                 # time step
    T = conf.T
    N = int(T/dt)        # horizon size
    PLOT_STUFF = 1
    thrust_guess = np.tile(np.array([4.6,0.0,0.0,0.0]),(N))
    linestyles = ['-*', '--*', ':*', '-.*']
    ode = ode('ode',get_quad_data('Quad1'),conf.dt)
    # create OCP
    problem = SingleShootingProblem('ssp', ode, conf.x0, dt, N, conf.integration_scheme)
    
    # simulate motion with initial guess    
        #Noooope meme
    # create cost function terms
    final_cost_state = OCPFinalCostState( conf.q_des, conf.v_des, conf.weight_vel)
    problem.add_final_cost(final_cost_state)
    effort_cost = OCPRunningCostQuadraticControl(dt,conf.weight_run_state) # effort cost modificato
    problem.add_running_cost(effort_cost, conf.weight_r) 
    
    #traj_cost = OCPFinalCostLength(conf.q_des, conf.v_des, conf.weight_vel,dt,N)
    #problem.add_final_cost(traj_cost)
    #constraint = OCPRunningConstraint(dt,conf.weight_const,conf.x0,problem.position_w_z,problem.size_z,problem.dist_wind,N) # da cambiare
    #problem.add_running_cost(constraint,conf.weight_r)   
#    problem.sanity_check_cost_gradient()
   
    """
    # solve OCP
#    problem.solve(method='Nelder-Mead') # solve with non-gradient based solver
    problem.solve(use_finite_difference=conf.use_finite_difference)
    print('U norm:', norm(problem.U))
    print('X_N\n', problem.X[-1,:].T)
    
     
    """    
    # Buonds on u value # sono bound su gli input, il problema è che non riesco a fare bound sugli stati, non so come passarli alla funzione
    a = (-0.05,0.05) # quindi questi bound sono in teoria su manovra phi,theta,psi quelli (None,NOne ) e quello solo positivi il trust
    b = (0.0,10.0) # de ve essere superiore a zero
    bnds = (b,a,a,a)*N
    #print('bounds:', bnds)
    
    # bounds on X state value, implemented as constraints ?
    """cons1 = {'type': 'ineq', 'fun': compute_cost_w_gradient_fd.X[3] +0.1 } 
    cons2 = {'type': 'ineq', 'fun': -compute_cost_w_gradient_fd.X[3] -0.1 }
    
    all_cons = [cons1,cons2] # 
    """
    # function that use bounds
    
    problem.solve_bounds(method='slsqp',y0=thrust_guess,use_finite_difference=conf.use_finite_difference,bnds=bnds) # l-bfgs-b -> sucks // slsqp -> several runtime error or NaN
    print('U norm:', norm(problem.U))
    print('X_N\n', problem.X[-1,:].T)                                       
    
    import datetime

    datetime_object = datetime.datetime.minute
    pd.DataFrame(problem.X).to_csv(f"/home/test/Desktop/Desktop/GitAORC/AORCProject/Code/optimizationresult.csv")
    pd.DataFrame(problem.U).to_csv(f"/home/test/Desktop/Desktop/GitAORC/AORCProject/Code/U_control.csv")
    
    mz = ((problem.position_w_z-problem.size_z/2)-(problem.X[2]-5.0))/(problem.dist_wind-problem.x0[0])
    
    Y=np.zeros(1)
    for i in range(N):
        Y=np.append(Y,[mz*problem.X[i,0] + problem.x0[2]-5.0])

    print(Y)
    # all in one plot
    
    f, ax = plut.create_empty_figure(1)
    time = np.arange(0.0, T+dt, dt)
    time = time[:N]
    ax.plot(time, problem.U[:N,1], label ='U manouvre phi')
    ax.plot(time, problem.U[:N,2], label ='U manouvre theta')
    ax.plot(time, problem.U[:N,3], label ='U manouvre psi')
    ax.legend()
    matplot.pyplot.xlabel('Time [s]')
    
    f, ax = plut.create_empty_figure(1)
    time = np.arange(0.0, T+dt, dt)
    time = time[:N]
    ax.plot(time, problem.U[:N,0], label ='Thrust')
    ax.legend()
    matplot.pyplot.xlabel('Time [s]')
    
    f, ax = plut.create_empty_figure(1)
    time = np.arange(0.0, T+dt, dt)
    time = time[:N]
    ax.plot(time, problem.X[:N,3]*180/3.14, label ='angle phi')
    ax.plot(time, problem.X[:N,4]*180/3.14, label ='angle theta')
    ax.plot(time, problem.X[:N,5]*180/3.14, label ='angle psi')
    ax.legend()
    matplot.pyplot.xlabel('Time [s]')
    
    
    f, ax = plut.create_empty_figure(1)
    time = np.arange(0.0, T+dt, dt)
    time = time[:N]
    ax.plot(problem.X[:N,0], problem.X[:N,1], label ='path y')
    ax.plot(problem.X[:N,0], problem.X[:N,2], label ='path z')
    #ax.plot(problem.X[:N,0], Y[:N], label ='cons')
    ax.legend()
    matplot.pyplot.xlabel('X-coord [m]')
    plt.show() 