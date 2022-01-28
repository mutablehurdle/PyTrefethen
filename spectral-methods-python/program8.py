'''
Spectral methods in MATLAB.
compare to Trefethen p8.m
'''

#eigenvalues of harmonic oscillator -u"+x^2 u on R

from numpy import *
from scipy.linalg import toeplitz
set_printoptions(linewidth=200,precision=14)

L = 8.                             # domain is [-L L], periodic
for N in range(6,42,6):
	h = 2*pi/N; x = h*arange(1,N+1); x = L*(x-pi)/pi
	i = array( range(1,N) )
	column = hstack(( [-pi**2/(3*h**2)-1./6.] ,\
                     -.5*(-1)**i/sin(h*i/2.)**2 ))
	D2 = (pi/L)**2*toeplitz(column)  # 2nd-order differentiation
	eigenvalues = sort(linalg.eigvals(-D2 + diag(x**2)))
	print(N, eigenvalues[0:4])
