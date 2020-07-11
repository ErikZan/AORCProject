# Functions 
Define instance of class OCPfinalCostState
#### OCPfinalCostState
Class that define the final cost state with argument (desired position,desired velocity,weight on velocity). 
This class has two methods :
* compute()
takes as input state and define the cost as a quadratic function of position error and weighted velocity error 
* compute_w_gradient()
takes as input state and define the cost as a quadratic function of position error and weighted velocity error and compute the gradient : return an array with position error and weighted velocity error
# 
calls the method add_final_cost() on the instance of class SingleShootigProblem with argument the instance of the class OCPfinalCostState
##### SingleShootigProblem.add_final_cost()
this  method takes as input a class (final cost instance) and a optional weight and increases the class variable final_costs(list) with an array containig weight and class
# 
define instance of class OCPRunnigCostQuadraticControl. 
#### OCPRunnigCostQuadraticControl
class that define the running cost and has a class variable dt. It has teo methods :
* compute() : which takes as input state,control,t,recompute = True(boolean) and compute the quadratic cost
* compute_w_gradient() : same argument as compute() and compute the quadratic cost return also its gradient wrt x,u 
If we want to add a running cost on state (trajectory) we should change the cost variable adding quadratic cost of states(x)

