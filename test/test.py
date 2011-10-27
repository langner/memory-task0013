import random

import numpy as np
from scipy.spatial import distance
import pylab

population = 1000
size = 64

C = []
Nx = int(np.sqrt(population))
Ny = int(population // Nx)
n = population - Nx*Ny
dX = size * 1.0 / Nx
dY = size * 1.0 / (Ny+1)
for i in range(population):
    dx = dX*random.random()/3.25
    dy = dY*random.random()/3.25
    if i < Nx*Ny:
        X = dX * (i%Nx + dx)
        Y = dY * (i//Nx + dy)
    else:
        X = (i - Nx*Ny + dx) * size * 1.0 / n
        Y = dY * (Ny + dy)
    C.append([X,Y])
    print X,Y
pds = distance.pdist(C)
hist,bins = np.histogram(pds, bins=256, range=(0.0,32.0))
print hist
pylab.plot(bins[:-1],hist)
pylab.show()