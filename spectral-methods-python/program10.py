'''
Spectral methods in MATLAB.
compare to Trefethen p10.m
'''

#polynomials and corresponding equipotential curves

from numpy import *
from pylab import *
N = 16
for i in [0,1]:
	if i==0: s = 'equispaced points'; x = linspace(-1.,1.,N+1)
	if i==1: s = 'Chebyshev points' ; x = cos(linspace(0,pi,N+1))
	p = poly(x)
	# Plot p(x) over [-1,1]:
	xx = arange(-1.,1.005,.005); pp = polyval(p,xx)
	subplot(2,2,2*i+1)
	plot(x,0*x,'.',markersize=12)
	plot(xx,pp,linewidth=2); grid(); xlim(-1.,1.)
	xticks([-1.,-.5,0,.5,1]); title(s,fontsize=8)
	# Plot equipotential curves:
	subplot(2,2,2*i+2)
	plot(real(x),imag(x),'.',markersize=12)
	xgrid = arange(-1.4,1.42,.02); ygrid = arange(-1.12,1.14,.02)
	xx,yy = meshgrid(xgrid,ygrid); zz = xx+1j*yy
	pp = polyval(p,zz); levels = 10.**arange(-4,1)
	contour(xx,yy,abs(pp),levels,colors=['m']); title(s,fontsize=8);
	
show()
