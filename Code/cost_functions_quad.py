# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 18:12:04 2020

@author: student
"""

import numpy as np
from numpy.linalg import norm


class Empty:
    def __init__(self):
        pass
        
        
class OCPFinalCostState:
    ''' Cost function for reaching a desired state of the robot
    '''
    def __init__(self, q_des, v_des, weight_vel):
        
        self.nq = 6
        self.q_des = q_des   # desired joint angles
        self.v_des = v_des  # desired joint velocities
        self.weight_vel = weight_vel
        
    def compute(self, x, recompute=True):
        ''' Compute the cost given the final state x '''
        q = x[:self.nq]
        v = x[self.nq:]
        e = q-self.q_des
        de = v - self.v_des
        cost = 0.5*e.dot(e) + 0.5*self.weight_vel*de.dot(de)
        return cost
        
    def compute_w_gradient(self, x, recompute=True):
        ''' Compute the cost and its gradient given the final state x '''
        q = x[:self.nq]
        v = x[self.nq:]
        e = q-self.q_des
        de = v - self.v_des
        cost = 0.5*e.dot(e) + 0.5*self.weight_vel*de.dot(de)
        grad =  np.concatenate((e, self.weight_vel*de))
        return (cost, grad)
        
        
class OCPRunningCostQuadraticControl:
    ''' Quadratic cost function for penalizing control inputs '''
    def __init__(self, dt, weight_run_state):
        self.dt = dt
        self.weight_run_state = weight_run_state
    def compute(self, x, u, t, recompute=True):
        ''' Compute the cost for a single time instant'''
        cost = 0.5*u.dot(u) + 0.5*self.weight_run_state*x.dot(x)
        return cost
        
    def compute_w_gradient(self, x, u, t, recompute=True):
        ''' Compute the cost for a single time instant and its gradient w.r.t. x and u '''
        cost = 0.5*u.dot(u)+ 0.5*self.weight_run_state*x.dot(x)
        grad_x = self.weight_run_state*x
        grad_u = u
        return (cost, grad_x, grad_u)


""" class OCPFinalCostFrame:
    ''' Cost function for reaching a desired position-velocity with a frame of the robot
        (typically the end-effector).
    '''
    def __init__(self, p_des, dp_des, weight_vel):
        self.robot = robot
        self.nq = robot.model.nq
        self.frame_id = robot.model.getFrameId(frame_name)
        assert(robot.model.existFrame(frame_name))
        self.p_des  = p_des   # desired 3d position of the frame
        self.dp_des = dp_des  # desired 3d velocity of the frame
        self.weight_vel = weight_vel
        
    def compute(self, x, recompute=True):
        ''' Compute the cost given the final state x '''
        q = x[:self.nq]
        v = x[self.nq:]

        H = self.robot.framePlacement(q, self.frame_id, recompute)
        p = H.translation # take the 3d position of the end-effector
        v_frame = self.robot.frameVelocity(q, v, self.frame_id, recompute)
        dp = v_frame.linear # take linear part of 6d velocity
        
        cost = norm(p-self.p_des) + self.weight_vel*norm(dp - self.dp_des)
        
        return cost
         """
    # gradient not implemented yet