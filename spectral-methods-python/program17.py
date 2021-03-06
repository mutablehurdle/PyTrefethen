'''
Spectral methods in MATLAB.
compare to Trefethen p17.m
'''

# Helmholtz equation u_xx + u_yy + (k^2)u = f
# on [-1,1]x[-1,1]

from mayavi import mlab
from numpy import *
from cheb import *
import time
from numpy.linalg import matrix_power,solve
from matplotlib import pyplot as plt
from scipy.interpolate import interp2d,griddata,bisplev,bisplrep
from mpl_toolkits.mplot3d import axes3d

# Set up spectral grid and tensor product Helmholtz operator:

N = 24
D,x = cheb(N)
y = x
xx,yy = meshgrid(x[1:N],y[1:N])
xx = hstack(xx[:]); yy = hstack(yy[:]);                   # strech 2D grids to 1D vectors
f = exp(-10*((yy-1)**2+(xx-0.5)**2))
D2 = matrix_power(D, 2)
D2 = D2[1:N,1:N]
I = identity(N-1)
k = 9
L = kron(I,D2) + kron(D2,I) + k**2*identity((N-1)**2)

# Solve for u, reshape to 2D grid, and plot:
u = solve(L,f);
uu = zeros((N+1,N+1))
uu[1:N,1:N]=u.reshape(N-1,N-1)
xx,yy = meshgrid(x,y)
ar = arange(-1,1,0.0333)
xxx,yyy = meshgrid(ar,ar)
uuu = interp2d(x, y, uu, kind='cubic')

mlab.surf(uuu(ar,ar),warp_scale="auto")
mlab.show()

#plot contours

plt.figure()
plt.contour(xxx, yyy, uuu(ar,ar))
plt.axis('equal')

plt.show()
