#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 16:42:18 2020

@author: mikegriffiths
"""


"""
% Example of solving hydrodynamics equation to advect fluid
% Illustrates the problem of numerical instability
%
% Higher order terms for the finite element derivatives are n0t sufficient
% to remove discontinuities
%
%
% Use two step Lax-Wendroff method to stabilise solution
% http://en.wikipedia.org/wiki/Lax-Wendroff_method


%    The  time-dependent advection equation is:
%
%      du/dt +  du/dx = 0
%

%
%      du/dt + dF/dx = 0
%
%    For the advection equation,  define 
%
%      F(x,t) =  u,
%      A(x,t) = dF/dx = u
%
%    and then the Lax-Wendroff method approximates the solution 
%    using the iteration:
%
%      u(x,t+dt) = u(t) - dt dF/dx + 1/2 dt^2 d/dx A dF/dx
%
%    which can be written:
%
%      u(x,t+dt) = u(x,t) - dt ( F(x+dx,t) - F(x-dx,t) ) / ( 2 * dx )
%        + 1/2 dt^2/dx^2 ( A(x+dx/2,t) * ( F(x+dx,t) - F(x,t) )
%                        - A(x-dx/2,t) * ( F(x,t) - F(x-dx,t) )
%
%    Use the approximation:
%
%      A(x+dx/2,t) = 1/2 ( u(x+dx,t) + u(x,t) )
%      A(x-dx/2,t) = 1/2 ( u(x,t) + u(x-dx,t) )
%
%    There is a stability condition that applies here, which requires that
%
%      dt * max ( abs ( u ) ) / dx <= 1
"""

from numpy import *
from os import *


import numpy as np
import math


import matplotlib.pyplot as plt
import matplotlib.animation as animation



# Constants
g = 9.81;
u0 = 0;
v0 = 0;
b = 2;
h0 = 5030;
damp=1.0;
force=0.0;
forcefreq=0.01;
k=2.5;


damp=0.5;
force=-0.0;


# Define the x domain
ni = int(151);
xmax = 1.0;
dx = xmax/(ni-1);

x = np.arange(0,xmax,dx);


# Define the y domain
nj = int(1);
ymax = 1.0;
dy = ymax/(nj);

y = np.arange(0,ymax,dy);

ni = int(np.size(x));
nj = int(np.size(y));

# Define the wavespeed
wavespeed = u0+math.sqrt(k);
wavespeed=0.01;

# Define time-domain
dt = 50*(0.68*dx)/wavespeed;
#dt=0.025;
#dt=0.25;
tmax = 1;
#t = [0:dt:tdomain];
t = np.arange(0,tmax,dt)
nt=(np.size(t));
courant = (wavespeed*dt)/dx;
nt=151;
# Build empty u, v, b matrices
u = np.zeros(  ((np.size(x)))  );
v = np.zeros(  ((np.size(x)))  );
uh = np.zeros(  ((np.size(x)))  );
vh = np.zeros(  ((np.size(x)))  );
tv = np.zeros( ((np.size(x))) );



for i in range( 0,ni):
    #if i>75:
    #    u[i]=1.0
    #else:
    #    u[i]=0.0
    u[i]  =  b*math.sin(2*i*math.pi/(ni-1)); 



#1st order periodic boundary condition
u[ni-1]=u[2];
u[1]=u[ni-2];

        
"""
u[0]  = u[nj-4] ;
u[1]  = u[nj-3] ;
u[nj-1]  = u[2] ;
u[nj-2]  = u[3] ;
v[0]  = v[nj-4] ;
v[1]  = v[nj-3] ;
v[nj-1]  = v[3] ;
v[nj-2]  = v[2] ;
"""


fig, ax = plt.subplots()

xv=x
av=u

print(np.shape(xv))
print(np.shape(av))
print(np.shape(v))
#x = np.arange(0, 2*np.pi, 2*np.pi/50)
line, = ax.plot(xv, av)



def init():  # only required for blitting to give a clean slate.
    line.set_ydata([np.nan] * len(xv))
    return line,


def animate(n):
    #it=1
    #if i<120 :
    #    it=i
    #it=i
    #print(it)
    #outfile='out/outfile_'+str(it)+'.out.npz';
    #np.load(outfile)
    #uall=data['v']


    #first order 
    starti=2;
    finishi=ni-1;
    
    #second order %uncomment these to test second order differencing
    #starti=3;
    #finishi=ni-2;



    
    #dt=0.05; 
    dt=0.00287;
    t=n*dt;

    c=-wavespeed*dt/dx;
    
    for i in range(starti,finishi):        
    #comment the above     
    #for i in range(3,ni-2)
    #Lax-Wendroff   
        v[i] = u[i]+c*(u[i+1]-u[i-1])/2+math.pow(c,2)*(u[i+1]-2*u[i]+u[i-1])/2;
    #first order central differencing
        #v(i) = u(i)+c*(u(i+1)-u(i-1))/2;%+c*(uold(i+1)-uold(i-1))/2;
    #second order central differencing %includes Lax-Wedroff correction term
       #v(i) = u(i)+c*(u(i-2)+8*u(i+1)-8*u(i-1)-u(i+2))/12;%+c^2*(u(i+1)-2*u(i)+u(i-1))/2;
   
    uold=u;
    for i in range(starti,finishi+1):        
       u[i]=v[i];
    #1st order periodic boundary condition
    u[ni-1]=u[2];
    u[1]=u[ni-2];
  
  #2nd order periodic boundary condition
  #u(ni-1)=u(3);
  #u(2)=u(ni-2);
  #u(ni)=u(4);
  #u(1)=u(ni-3);
    
   
    
   
    
   
    av=u[:]
    #print(u)
    line.set_ydata(av)
    #line.set_ydata(np.sin(x + i / 100))  # update the data.
    #print(av[60])
    return line,


ani = animation.FuncAnimation(
    fig, animate, init_func=init, interval=2, blit=True, save_count=50)

# To save the animation, use e.g.
#
# ani.save("movie.mp4")
#
# or
#
# from matplotlib.animation import FFMpegWriter
# writer = FFMpegWriter(fps=15, metadata=dict(artist='Me'), bitrate=1800)
# ani.save("movie.mp4", writer=writer)

#plt.plot(x)
plt.show()




