import os,sys

from sys import *
from numpy import *
#from os import *
from math import *


import numpy as np


#shooting method example using python

#      PROGRAM STARTSHOOTING
#  THIS PROVIDES A SOLUTION FOR A BOUNDARY-VALUE PROBLEM FOR A
#  FIRST-ORDER ODE WITH ONE UNKNOWN PARAMETER TO BE DETERMINED.
#  THE ODE IS OF THE FORM     Dy/Dx = f(x,Q)  WHERE Q IS THE 
#  UNKNOWN PARAMETER.   THE BOUNDARY CONDITIONS ARE y = Ya FOR
#  x = a AND y = Yb FOR X = b.   
#  THE FUNCTION f(x,Q) IS PROVIDED AS A FUNCTION STATEMENT.
def f(x,q):
    f=-15.915494*q/(math.pow(2-x,2))
    return f



#main routine starts here
y=zeros(1000)
w=zeros(4)
dy=zeros(4)

w=[0.0,0.5,0.5,1.0]

print 'INPUT THE FIRST BOUNDARY CONDITION a,Ya'
print 'WHERE y = Ya WHEN x = a.'

xi=float(input('a'))
y[0]=float(input('Ya'))

print 'INPUT THE FINAL BOUNDARY CONDITION  b,Yb'
print 'WHERE y = YbWHEN x = b.'
xf=float(input('b'))
yf=float(input('Yb'))

print' INPUT THE INTEGRATION STEP LENGTH h IN THE '
print' FORM OF AN INTEGER N WHERE h = [b - a]/N. '
print' IF YOU WANT OUTPUT OF Y AT INTERVALS OF [B-A]/M'
print' THEN MAKE N A MULTIPLE OF M.'


n=int(input('n'))
h=(xf-xi)/n

q=float(input('INPUT ESTIMATE OF THE UNKNOWN PARAMETER, Q.'))


#begin runge-kutta integration
for i in range(1,n):
    x=(i-1)*h
    for j in range(1,4):
        xx=x+w[j]*h
        dy[j]=h*f(xx,q)
    y[i]=y[i-1]+(dy[0]+dy[3]+2*(dy[1]+dy[2]))/6.0
print'The value of Yb is ',y[n]

m=input('Output interval?')
#printed output?
fd = open("mydoc.txt", "w")
k=n/m
for i in range(0,n,k):
    xx=i*h
    fouts=repr(xx)+' '+repr(y[i])+'\n'
    print(fouts)
    fd.write(fouts)
fd.close()

