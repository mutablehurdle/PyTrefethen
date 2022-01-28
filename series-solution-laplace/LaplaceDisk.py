


#import necessary modules

import numpy as np
import matplotlib.pyplot as plt

#define center of disk using complex coordinates (real numbers are x, imaginary are y)

c = 3 +1j #center of disk
r = 1    #radius of disk

#use power series expansion to solve Laplace's equation

N =10      #number of expansion terms
npts = 3*N #number of sample points

npts_range=np.linspace(0,npts,npts)
z = c + r*np.exp(2j*np.pi*npts_range/npts) #sample points along disk boundary

rhs = -np.log(np.abs(z))+np.log(np.abs(z-c)) #right-hand side
A = np.ones((npts,2*N+1))

for k in np.arange(1,N+1):
    A[:,2*k-1] = np.real((z-c)**(-k))
    A[:,2*k] = np.imag((z-c)**(-k))
    
a = np.linalg.lstsq(A,rhs)[0]


#plot disk to check that it is centered where you want it, with the radius given and sample points
plt.plot(np.real(z-c),np.imag(z-c),'.')


#Define power series solution to Laplace's equation in complex numbers
def disk1fun(z):
    u = np.log(np.abs(z)) - np.log(np.abs(z-c)) + a[0]
    
    for k in np.arange(1,N+1):
        u+=a[2*k-1]*np.real((z-c)**(-k))+a[2*k]*np.imag((z-c)**(-k))
    u[abs(z-c)<=r]=np.NaN
    
    return u


#contour plot the solution on a grid
xx,yy = np.linspace(-5,5,145),np.linspace(-4,4,115)
[xx,yy] = np.meshgrid(xx,yy)
zz = xx + 1j*yy
uu = disk1fun(zz)
levels = np.linspace(-3,-.25,10)
plt.contourf(xx,yy,uu)