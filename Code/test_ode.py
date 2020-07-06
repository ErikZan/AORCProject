from __future__ import print_function
import numpy as np
import matplotlib as matplot 
import sys
import math
import numpy.linalg as la
import matplotlib.pyplot as plt
from math import sqrt
from scipy.integrate import odeint

def model(y,t):
    k=0.3
    dydt = -k*y
    return dydt
  
ys = 5

t=np.linspace(0,20)

y=odeint(model,ys,t)

plt.plot(t,y)
plt.show()