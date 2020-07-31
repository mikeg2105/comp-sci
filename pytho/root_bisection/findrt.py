
import math

# here is an example use of the while statement
# which is used for finding the root of a polynomial 
# which is known to lie within a certain interval.
# a is the lower value of the range
# b is the upper value of the range
eps=0.0000001
inf=9.9e30
a = 0.0
fa = -inf
b = 4.0
fb = inf
test=1.0
while abs(b-a)>(eps*b):
  x = (a+b)/2.0
  fx = math.pow(x,3)-2*x-5
  print a,b,(a+b)/2,x,fx,fa
  if math.copysign(test,fx) == math.copysign(test,fa):
    a = x
    fa = fx
  else:
    b = x
    fb = fx


print ' The root is :',x

