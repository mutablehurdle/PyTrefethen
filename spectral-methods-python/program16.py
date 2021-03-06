'''
Spectral methods in MATLAB.
compare to Trefethen p16.m
'''

# Poisson equation on [-1,1]x[-1,1] with u=o on boundary

from mayavi import mlab
from numpy import *
import time
from cheb import cheb
from numpy.linalg import matrix_power,solve
from matplotlib import pyplot as plt
from scipy.interpolate import interp2d,griddata,bisplev,bisplrep
# from mpl_toolkits.mplot3d import axes3d

# Set up grids and tensor product Laplacian and solve for u:

N = 24
D,x = cheb(N)
y = x
xx,yy = meshgrid(x[1:N],y[1:N])
xx = hstack(xx[:]); yy = hstack(yy[:]);                   # strech 2D grids to 1D vectors
f = 10*sin(8*xx*(yy-1))
D2 = matrix_power(D, 2)
D2 = D2[1:N,1:N]
I = identity(N-1)
L = kron(I,D2) + kron(D2,I)          # Laplacian
plt.figure(1)
plt.spy(L)
tic = time.time(); u = solve(L,f); toc = time.time() - tic;      # Solve problem and watch the clock

# Reshape long 1D results onto 2D grid:
uu = zeros((N+1,N+1))
uu[1:N,1:N]=u.reshape(N-1,N-1)
xx,yy = meshgrid(x,y)
value = uu[int(N/4.+1),int(N/4.+1)]

# Interpolate to finer grid and plot:
uuu = interp2d(x, y, uu, kind='cubic')

mlab.surf(uuu(arange(-1,1,0.04),arange(-1,1,0.04)),warp_scale="auto")
mlab.show()
