'''
Spectral methods in MATLAB.
compare to Trefethen p19.m
'''

# 2nd order wave eq. on Chebyshev grid

from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from chebfft import *
from mayavi import mlab

N = 80
x = cos(pi*arange(0,N+1)/N)
dt = 8.0/N**2
v = exp(-200*x**2)
vold = exp(-200*(x-dt)**2)
tmax = 4
tplot = 0.075
plotgap = int(round(tplot/dt))
dt = tplot/plotgap
nplots = int(round(tmax/tplot))
plotdata = vstack((v, zeros((nplots,N+1))))
tdata = 0
for i in range(0,nplots):
    for n in range(0,plotgap):
        w = chebfft(chebfft(v)).T
        w[0] = 0
        w[N] = 0
        vnew = 2*v - vold +dt**2*w
        vold = v
        v = vnew
    plotdata[i+1,:] = v
    tdata = vstack((tdata, dt*i*plotgap))

# Plot results

X, Y = meshgrid(x, tdata)

mlab.surf(plotdata,warp_scale="auto")
mlab.show()
