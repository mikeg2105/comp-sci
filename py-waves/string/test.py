#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 12:37:28 2020

@author: mikegriffiths
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig, ax = plt.subplots()

x = np.arange(0, 2*np.pi, 2*np.pi/50)
line, = ax.plot(x, np.sin(x))
n=120

outfile='out/outfile_'+str(n)+'.out.npz';
data=np.load(outfile)  
print(data['v'])
uall=data['v']
yt=uall[25]
print(yt)
print(np.shape(yt))
it=120

#print(it)
#outfile='out/outfile_'+str(it)+'.out.npz';
#np.load(outfile)
#uall=data['v']
#yt=uall[25]
