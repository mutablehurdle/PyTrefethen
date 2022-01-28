
'''
Spectral methods in MATLAB.
compare to Trefethen p9.m
'''

# polynomial interpolation in equispaced and Chebyshev pts

from numpy import *
from pylab import *
N = 16
xx = arange(-1.01,1.015,.005) #-1.01:.005:1.01
for i in [1,2]:
	if i==1: s = 'equispaced points'; x = -1. + 2.*arange(0,N+1)/N
	if i==2: s =  'Chebyshev points'; x = cos(pi*arange(0,N+1)/N)
	subplot(1,2,i)
	u  = 1/(1+16*x**2 )
	uu = 1/(1+16*xx**2)
	p = polyfit(x,u,N)              # interpolation
	pp = polyval(p,xx)              # evaluation of interpolant
	plot(x,u,'.',markersize=8)
	plot(xx,pp)
	xlim(-1.1,1.1); ylim(-1,1.5); title(s)
	error = linalg.norm(uu-pp,inf)
	text(-.5,-.5,'max error = '+str('%.4f' % error) )
show()
