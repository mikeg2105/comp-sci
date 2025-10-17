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
ni = int(51);
xmax = 1.0;
dx = xmax/(ni-1);

x = np.arange(0,xmax,dx);


# Define the y domain
nj = int(51);
ymax = 1.0;
dy = ymax/(nj-1);

y = np.arange(0,ymax,dy);



# Define the wavespeed
wavespeed = u0+math.sqrt(k);

# Define time-domain
dt = 10*(0.68*dx)/wavespeed;
dt=0.025;
dt=0.25;
tmax = 1;
#t = [0:dt:tdomain];
t = np.arange(0,tmax,dt)
nt=(np.size(t));
courant = (wavespeed*dt)/dx;
nt=128;
# Build empty u, v, b matrices
u = np.zeros(  ((np.size(x)),(np.size(y)))  );
v = np.zeros(  ((np.size(x)),(np.size(y)))  );
tv = np.zeros( ((np.size(x)),(np.size(y)))  );
u[int((ni+1)/2)][int((nj+1)/2)]=b;
#for i in range(2,ni):
#    for j in range(2,nj):       
#      #u(i,j) = -b*sin(i*dx*%pi)*sin(j*dx*%pi);
#      #v(i,j,1) = b*cos(i*dx*%pi)*cos(j*dx*%pi);    
     
# Employ Lax
for n in range(1,nt+1):
    dt=0.05; 
    t=n*dt;
    for i in range( 2,ni-2):
        for j in  range(2,nj-2):  
            tv[i][j] = (1.0-damp)*v[i][j]-(dt/2)*dx*k*(4*u[i][j]-u[i+1][j]-u[i-1][j]-u[i][j+1]-u[i][j-1]);#+force*dt*sin(forcefreq*t*%pi);
    tv[int((ni+1)/2)][int((nj+1)/2)] = tv[int((ni+1)/2)][int((nj+1)/2)]+(force*dt*math.sin(forcefreq*t*math.pi));            
    for i in range( 2,ni-1):
        for j in  range(2,nj-1):
            u[i][j]  = u[i][j] +dt*tv[i][j] ;
    for i in range( 2,ni-1):
        for j in  range(2,nj-1):
            v[i][j]  = v[i][j] +dt*tv[i][j] ;  
    for i in range( 2,ni-2):
        for j in  range(2,nj-2):
            v[i][j] = v[i][j]-(dt)*dx*k*(4*u[i][j]-u[i+1][j]-u[i-1][j]-u[i][j+1]-u[i][j-1]);  #//+force*dt*sin(forcefreq*t*%pi);
    v[int((ni+1)/2)][int((nj+1)/2)] = v[int((ni+1)/2)][int((nj+1)/2)]+force*dt*math.sin(forcefreq*t*math.pi);                
   # Define Boundary Conditions
   #u(1,:,n+1) = 2.5*u(2,:,n+1)-2*u(3,:,n+1)+0.5*u(4,:,n+1);
   #u(max(size(x)),:,n+1) = 2.5*u(ni-1,:,n+1)-2*u(ni-2,:,n+1)+0.5*u(ni-3,:,n+1);
   #u(:,1,n+1) = 2.5*u(:,2,n+1)-2*u(:,3,n+1)+0.5*u(:,4,n+1);
   #u(:,max(size(y)),n+1) = 2.5*u(:,nj-1,n+1)-2*u(:,nj-2,n+1)+0.5*u(:,nj-3,n+1);
    outfile='out/outfile_'+str(n)+'.out';
    np.savez(outfile,u=u,v=v);
   #realtime(i); //wait till date 0.1*i seconds
   #s.data.z = (sin((I(i)/10)*x)'*cos((I(i)/10)*y))';
#end of lax timeloop
