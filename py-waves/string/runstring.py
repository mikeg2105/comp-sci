#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 16:26:59 2020

@author: mikegriffiths


"""


from numpy import *
from os import *

import matplotlib.pyplot as plt
import numpy as np
import math



# Constants
g = 9.81;
u0 = 0;
v0 = 0;
b = 2;
h0 = 5030;
damp=0.0;
force=1;
forcefreq=1.58114;
k=2.5;

# Define the x domain
ni = 51;
xmax = 1.0;
dx = xmax/(ni-1);

x = np.arange(0,xmax,dx);


# Define the y domain
nj = 51;
ymax = 1.0;
dy = ymax/(nj-1);

y = np.arange(0,ymax,dy);



# Define the wavespeed
wavespeed = u0+math.sqrt(k);

# Define time-domain
dt = 10*(0.68*dx)/wavespeed;
dt=0.025;
tmax = 1;
#t = [0:dt:tdomain];
t = np.arange(0,tmax,dt)
nt=(np.size(t));
courant = (wavespeed*dt)/dx;
nt=32000;
# Build empty u, v, b matrices
u = np.zeros(  ((np.size(x)),(np.size(y)))  );
v = np.zeros(  ((np.size(x)),(np.size(y)))  );
tv = np.zeros( ((np.size(x)),(np.size(y)))  );

