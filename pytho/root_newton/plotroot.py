from numpy import *
from os import *

import matplotlib.pyplot as plt
import numpy as np
import math


# This is an example use of the plot function 
# used in association with the finding of a root
# problem solved earlier.
# % a is the lower value of the range
# b is the upper value of the range
eps=0.0000001
inf=9.9e30
a = 0.0
fa = -inf
x = 3.0
fb = inf
b=x
delta = (b-a)/100.0

xx = np.arange(a,b,delta)
asize=xx.size
yy = np.zeros(asize)
for i in range(1,asize):
    yy[i]=0.5*sin(2*(xx[i]-(math.pi/4)))+0.5*sin(xx[i]);



while (sqrt(math.pow(x-a,2)))>eps:
  a=x;
  fx = 0.5*sin(2*(x-(math.pi/4)))+0.5*sin(x);
  dfx= cos(2*(x-(math.pi/4)))+0.5*cos(x);
  x=a-(fx/dfx);
  print x-a
 


print ' The root is :',x
plt.plot(xx,yy)