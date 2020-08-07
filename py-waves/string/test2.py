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
nj = int(6);
ymax = 1.0;
dy = ymax/(nj-1);

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
u = np.zeros(  ((np.size(x)),(np.size(y)))  );
v = np.zeros(  ((np.size(x)),(np.size(y)))  );
uh = np.zeros(  ((np.size(x)),(np.size(y)))  );
vh = np.zeros(  ((np.size(x)),(np.size(y)))  );
tv = np.zeros( ((np.size(x)),(np.size(y)))  );
u[int((ni+1)/2)][int((nj+1)/2)]=b;


for i in range( 0,ni):
    for j in  range(0,nj):
        if i>75:
            u[i][j]=1.0
        else:
            u[i][j]=0.0
        #u[i][j]  =  b*math.sin(2*i*math.pi/(ni-1)); 




for j in  range(0,nj):
        u[0][j]  = u[ni-4][j] ;
        u[1][j]  = u[ni-3][j] ;
        u[ni-1][j]  = u[2][j] ;
        u[ni-2][j]  = u[3][j] ;
        v[0][j]  = v[ni-4][j] ;
        v[1][j]  = v[ni-3][j] ;
        v[ni-1][j]  = v[3][j] ;
        v[ni-2][j]  = v[2][j] ;
        
for i in  range(0,ni):
        u[i][0]  = u[i][nj-4] ;
        u[i][1]  = u[i][nj-3] ;
        u[i][nj-1]  = u[i][2] ;
        u[i][nj-2]  = u[i][3] ;
        v[i][0]  = v[i][nj-4] ;
        v[i][1]  = v[i][nj-3] ;
        v[i][nj-1]  = v[i][3] ;
        v[i][nj-2]  = v[i][2] ;



fig, ax = plt.subplots()

xv=x
av=u[:,1]

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
    
    
    #dt=0.05; 
    dt=0.00287;
    t=n*dt;
    for i in range( 1,ni-1):
        for j in  range(1,nj-1):  
            tv[i][j] = (1.0-damp)*v[i][j]-(dt/2)*dx*k*(4*u[i][j]-u[i+1][j]-u[i-1][j]-u[i][j+1]-u[i][j-1])+force*dt*math.sin(forcefreq*t*math.pi);
            
    #print(np.max(tv))           
    #tv[int((ni+1)/2)][int((nj+1)/2)] = tv[int((ni+1)/2)][int((nj+1)/2)]+(force*dt*math.sin(forcefreq*t*math.pi));            
    for i in range( 0,ni-1):
        for j in  range(0,nj-1):
            u[i][j]  = u[i][j] +dt*tv[i][j] ;
    for i in range( 0,ni-1):
        for j in  range(0,nj-1):
            v[i][j]  = v[i][j] +dt*tv[i][j] ;  
    for i in range( 1,ni-1):
        for j in  range(1,nj-1):
            v[i][j] = v[i][j]-(dt)*dx*k*(4*u[i][j]-u[i+1][j]-u[i-1][j]-u[i][j+1]-u[i][j-1])+force*dt*math.sin(forcefreq*t*math.pi);
   # v[int((ni+1)/2)][int((nj+1)/2)] = v[int((ni+1)/2)][int((nj+1)/2)]+force*dt*math.sin(forcefreq*t*math.pi);  


    for j in  range(0,nj):
            u[0][j]  = u[ni-4][j] ;
            u[1][j]  = u[ni-3][j] ;
            u[ni-1][j]  = u[2][j] ;
            u[ni-2][j]  = u[3][j] ;
            v[0][j]  = v[ni-4][j] ;
            v[1][j]  = v[ni-3][j] ;
            v[ni-1][j]  = v[3][j] ;
            v[ni-2][j]  = v[2][j] ;
            
    for i in  range(0,ni):
            u[i][0]  = u[i][nj-4] ;
            u[i][1]  = u[i][nj-3] ;
            u[i][nj-1]  = u[i][2] ;
            u[i][nj-2]  = u[i][3] ;
            v[i][0]  = v[i][nj-4] ;
            v[i][1]  = v[i][nj-3] ;
            v[i][nj-1]  = v[i][3] ;
            v[i][nj-2]  = v[i][2] ;

    #set boundaries
# Define Boundary Conditions
#  u[:,ni] = 2.5*u[]:,n+1)-2*u(3,:,n+1)+0.5*u(4,:,n+1);
    
#  //u(1,:,n+1) = 2.5*u(2,:,n+1)-2*u(3,:,n+1)+0.5*u(4,:,n+1);  
#  //u(max(size(x)),:,n+1) = 2.5*u(ni-1,:,n+1)-2*u(ni-2,:,n+1)+0.5*u(ni-3,:,n+1);
#  //u(:,1,n+1) = 2.5*u(:,2,n+1)-2*u(:,3,n+1)+0.5*u(:,4,n+1);
#  //u(:,max(size(y)),n+1) = 2.5*u(:,nj-1,n+1)-2*u(:,nj-2,n+1)+0.5*u(:,nj-3,n+1);           
    
    av=u[:,1]
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




