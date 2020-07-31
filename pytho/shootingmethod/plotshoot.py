from numpy import *
from os import *

import matplotlib.pyplot as plt
import numpy as np
#application which reads columnar data from file
#using matplotlib,scipy and numpy


#load data from text file
t1a=np.loadtxt('mydoc.txt')
asz=t1a.size
t1=t1a[1:asz,0]
t2=t1a[1:asz,1]

#data = zeros((3,3))
#t1 = np.arange(0.0, 5.0, 0.1)
#t2 = np.arange(0.0, 5.0, 0.1)

plt.plot(t1, t2)

plt.ylabel('some numbers')
plt.show()
