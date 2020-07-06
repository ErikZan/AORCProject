from __future__ import print_function
import numpy as np
import matplotlib as matplot 
import sys
import math
import numpy.linalg as la

from math import sqrt
from scipy.integrate import odeint

def model(y,t):
    k=0.3
    dydt = -k*y
    return dydt
  
y0 = 5

t=np.linspace(0,20)

y=odeint(model,y0,t)

matplot.pyplot.plot(t.y)