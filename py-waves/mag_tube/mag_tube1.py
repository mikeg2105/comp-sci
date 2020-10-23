#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 21:13:54 2020

@author: mikegriffiths
"""


import numpy as np
from numpy import *
import scipy.io
from scipy import special
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import struct

import csv
import numpy

from scipy import interpolate




def deriv1(f,x):
    #deriv1(f,x, dim)
    nel=np.size(f);
    nel1=np.size(x);
    
    f1=np.reshape(f,(nel1,1));
    num=np.shape(f);
    
    # if nel(1)>nel(2)
    #     num=nel(1);
    # else
    #     num=nel(2);
    # end
    
    # if (nel ne nel1) then begin
    #  print,'Inconsistant input, stop.'
    #  stop
    # endif
    
    res=np.zeros(nel);
    for i in range(3,nel-1):
        res[i]=(1/12.0/(x[i+1]-x[i]))*(8.0*f1[i+1]-8.0*f1[i-1]-f1[i+2]+f1[i-2]);
    
    #for i=1,nel-2 do res(i)=(1.d0/2.d0/(x(i+1)-x(i)))*(f(i+1)-f(i-1))
    res[1]=res[3];
    res[2]=res[3];
    res[nel]=res[nel-2];
    res[nel-1]=res[nel-2];
    return res,






def par4(x,x0,A):
    res=A*np.exp(-x/(x0));
    return res,


#%matplotlib inline  



#xmin=133333.33;
xmin=199219.0+6478.25;
#ymin=1953.1;
ymin=39687.5+1292.93;
#zmin=1953.1;
zmin=39687.5+1292.93;
xmax=5955555.6e0;
#xmax=12.8496e6;
#xmax=5.9e6;
#ymax=4.0e6;
#zmax=4.0e6;
ymax=2.559984e6;
zmax=2.559984e6;

#xmin=5.0e3

# Define the x domain
ni = int(51);
#xmax = 1.0;
dx = xmax/(ni-1);
x = np.arange(0,xmax,dx);


# Define the y domain
nj = int(51);
#ymax = 1.0;
dy = ymax/(nj-1);
y = np.arange(0,ymax,dy);



# Define the z domain
nk = int(51);
#zmax = 1.0;
dz = zmax/(nk-1);
z = np.arange(0,zmax,dz);

mu=0.6e0; #magnetic permeability
R=8.31e3;
ggg= -274.0;  #-274.0e0 % acceleration due to gravity on the sun


gamma=2.0/3.0;  # adiabatic constant
mu=4*math.pi/10000000; #magnetic permeability
p0=0;  #pressure in tube
pe=0;  #presuure outside tube
b0=0;  #mag field in tube
be=0;  #mag field outside tube
rho0=0;
rhoe=0;

c0=0;
va=0;
vae=0;
ce=0;


kw=0; #wavenumber
om=0; #wave frequency


mu_thermal=0.6e0;
R=8.31e3;






tempg=5109.4*np.ones(ni);
presg=7160.2*np.ones(ni);
densg=0.0001*np.ones(ni);
energ=np.zeros(ni);
rheight=np.zeros(ni);

for i in range(0,ni):
    rheight[i]=xmin+dx*i;

p = np.zeros(  ((np.size(x)),(np.size(y)),(np.size(z)))  ); #pressure 
en = np.zeros(  ((np.size(x)),(np.size(y)),(np.size(z)))  ); #energy
b = np.zeros(  ((np.size(x)),(np.size(y)),(np.size(z)))  ); #b field
rho = np.zeros(  ((np.size(x)),(np.size(y)),(np.size(z)))  ); #density



#load VAL IIIc data from text file
asz=0
with open('atmos.csv','rt')as f:
  data = csv.reader(f)
  for row in data:
        print(row)
        asz=asz+1
        
hval3c = np.zeros(  ((asz))  ); #density
rhoval3c = np.zeros(  ((asz))  ); #density
eval3c = np.zeros(  ((asz))  ); #energy
#t1a=np.zeros(asz,4)

#t1a=np.zeros(asz,4);
i=0
t1a=[]

#height, temp, density, energy
with open('atmos.csv','rt')as f:
  data = csv.reader(f)
  for row in data:
        t1a.append(row)



for i in range(0,asz):
    hval3c[i]=float(t1a[2047-i][0])
    rhoval3c[i]=float(t1a[2047-i][2])
    
f = interpolate.interp1d(hval3c, rhoval3c,fill_value="extrapolate")
densg = f(rheight)

for i in range(0,ni-1):
    if densg[i]<0:
        densg[i]=min(rhoval3c)

#now read through each row of t1a and set the height, density,  values

#reader = csv.reader(open("atmos.csv", "rb"), delimiter=",")
#tx = list(reader)
#t1a = numpy.array(tx).astype("float")

#t1a=np.loadtxt('atmos.csv')
#asz=data.size
#t1=t1a[1:asz,0]
#t2=t1a[1:asz,1]


#set up the atmosphere
iniene=731191.34e0*R*(0.0001)/mu_thermal/(gamma-1.0);

for i in range(0,ni-1):
    for j in range(0,nj-1):
        for k in range(0,nk-1):
            rho[i][j][k]=densg[i];  #density
            en[i][j][k]=iniene;  #energy
                                    







for i in range(0,ni-1):
    presg[i]=(gamma-1)*iniene;


presg1=presg;


#for i=nx1-1:-1:1
for i in range(1,ni-2):
    comi=-abs(rheight[i+1]-rheight[i]);
    presg[i]=presg[i+1]-densg[i]*comi*ggg;




#for i=3:nx1-2
for i in range(3,ni-3):
     comi=-abs(rheight[i+1]-rheight[i]);
     #densg(i)=densg(i)-(1.0/consts.ggg)*(  (1.0/(12*(rheight(i+1)-rheight(i)))) *(presg(i+2)-8*presg(i+1)+8*presg(i-1)-presg(i-2))     );
     densg[i]=(1.0/ggg)*(  (1.0/(12*(rheight[i+1]-rheight[i]))) *(presg[i+2]-8*presg[i+1]+8*presg[i-1]-presg[i-2])     );


for i in range(0,ni-1):
    if densg[i]<0:
        densg[i]=min(rhoval3c)


#for i=nx1-4:nx1-2
for i in range(ni-2,ni-5):
  p_1=presg[i-2]-8*presg[i-1]+8*presg[i+1];
  p_2=-densg[i]*ggg;
  presg[i+2]= p_1-12.0*(rheight[i]-rheight[i-1])*p_2;


presg2=presg;

#for i=1:nx1
for i in range(0,ni-1):
    energ[i]=presg[i]/(gamma -1);
    

#rebuild array    
for i in range(0,ni-1):
    for j in range(0,nj-1):
        for k in range(0,nk-1):
            rho[i][j][k]=densg[i];  #density
            en[i][j][k]=energ[i];  #energy    



#generate the field

bx=np.zeros(  ((np.size(x)),(np.size(y)),(np.size(z)))  );
by=np.zeros(  ((np.size(x)),(np.size(y)),(np.size(z)))  );
bz=np.zeros(  ((np.size(x)),(np.size(y)),(np.size(z)))  );



#intermediate fields to generate b field
b0z=np.zeros(  ((np.size(x)) )  );
    
xf=np.zeros(  ((np.size(x)),(np.size(y)),(np.size(z)))  );


    
    
b0z=np.zeros(ni);   #zeros(nx1,1);
    
xf=np.zeros(  ((np.size(x)),(np.size(y)),(np.size(z)))  );
    
Bmax=0.15  ; #mag field Tesla
#Bmin=0.0006d0  ; #mag field Tesla
Bmin=0.0002  ; #mag field Tesla

    
d_z=0.5; # width of Gaussian in Mm
z_shift= 0.0; # shift in Mm
A=0.45; # amplitude
scale=1.0e6;
b0z_top=0.08;
   
f0=2.0e6; #tube opening factor

Ab0z=20.e0; # bz - amplitude


xr=0.1e6;
yr=0.1e6;

R2=(xr**2+yr**2);

A=R2/2;


#rebuild array    
for i in range(0,ni-1):
    for j in range(0,nj-1):
        for k in range(0,nk-1):
            f=(y[j]**2+z[k]**2)/R2;
            xf[i][j][k]=np.exp(-f);

maxxf=np.max(xf);
b0zz=0.0001;

for i in range(0,ni-1):
    for j in range(0,nj-1):
        for k in range(0,nk-1):

            xf[i][j][k]=xf[i][j][k]/maxxf;





for i in range(0,ni-1):
        #b0z(i)=(par4((x(i)/scale-z_shift),d_z,A)).^2;
        #b0z(i)=(par3((x(i)/scale-z_shift),d_z,A));
        #b0z[i]=par4((x[i]/scale-z_shift),d_z,A);
        b0z[i]=A*np.exp(-((x[i]/scale-z_shift))/(d_z));


bnmin=np.min(b0z);
bnmax=np.max(b0z);




for i in range(0,ni-1):
	b0z[i]=((Bmax-Bmin)/(bnmax-bnmin))*(b0z[i]-bnmin)+Bmin;




for i in range(0,ni-1):
    for j in range(0,nj-1):
        for k in range(0,nk-1):
            bz[i][j][k]=bz[i][j][k]+b0zz*xf[i][j][k];
            bx[i][j][k]=bx[i][j][k];
            by[i][j][k]=by[i][j][k];

# bz(i,j,k)=bz(i,j,k)+(b0z(i)/sqrt((x(j)-ybp).^2+(y(k)-zbp).^2)*xf(i,j,k));
# bx(i,j,k)=bx(i,j,k)-(dbz(i)*(x(j)-ybp)/sqrt((x(j)-ybp).^2+(y(k)-zbp).^2)*xf(i,j,k));
# by(i,j,k)=by(i,j,k)-dbz(i)*(y(k)-zbp)/sqrt((x(j)-ybp).^2+(y(k)-zbp).^2)*xf(i,j,k);
#bz(i,j,k)=bz(i,j,k)+(b0z(i)/sqrt((x(j)).^2+(y(k)).^2)*xf(i,j,k));
#bx(i,j,k)=bx(i,j,k)-(dbz(i)*(x(j))/sqrt((x(j)).^2+(y(k)).^2)*xf(i,j,k));
#by(i,j,k)=by(i,j,k)-dbz(i)*(y(k))/sqrt((x(j)).^2+(y(k)).^2)*xf(i,j,k);

#step through the cycles set up the oscillations




#save oscillation configuration at each step










"""
fig = plt.figure()
ax = fig.gca(projection='3d')
X = alldat[0,:,:]
Y = alldat[1,:,:]
Z = alldat[2,:,:]

dens = alldat[2,:,:]
#bsq=alldat[6,:,:]*alldat[6,:,:]+alldat[7,:,:]*alldat[7,:,:]
#bmag=sqrt(bsq)

surf = ax.plot_surface(X, Y, dens, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
#ax.set_zlim(-1.01, 1.01)

#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
ax.view_init(elev=90,azim=0) 
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
"""