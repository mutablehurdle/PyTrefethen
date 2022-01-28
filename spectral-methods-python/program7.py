'''
Spectral methods in MATLAB.
compare to Trefethen p7.m
'''

from numpy import *
from scipy.linalg import toeplitz
from pylab import *
set_printoptions(linewidth=200,precision=5)
# Compute derivatives for various values of N:
Nmax = 50; E = zeros((4,Nmax//2-2))
Ns = range(6,Nmax+2,2)
for N in Ns:
	h = 2*pi/N; x = h*arange(1,N+1)
	i = array( range(1,N) )
	column   = hstack(( [0.], .5*(-1)**i/tan(i*h/2.) ))
	row      = hstack(( [0.], column[N-1:0:-1]       ))
	D = toeplitz(column,row)
	#print D
	v = abs(sin(x))**3								# 3rd deriv in BV
	vprime = 3*sin(x)*cos(x)*abs(sin(x))
	E[0,N//2-3] = linalg.norm(dot(D,v)-vprime,inf)
	v = exp(-sin(x/2)**(-2))						# C-infinity
	vprime = .5*v*sin(x)/sin(x/2)**4
	E[1,N//2-3] = linalg.norm(dot(D,v)-vprime,inf)
	v = 1/(1+sin(x/2)**2)							# analytic in a strip
	vprime = -sin(x/2)*cos(x/2)*v**2
	E[2,N//2-3] = linalg.norm(dot(D,v)-vprime,inf)
	v = sin(10*x); vprime = 10*cos(10*x)		# band-limited
	E[3,N//2-3] = linalg.norm(dot(D,v)-vprime,inf)

# Plot results:
titles = ['|sin(x)|^3','exp(-sin^{-2}(x/2))',\
          '1/(1+sin^2(x/2))','sin(10x)'      ]
for iplot in [0,1,2,3]:
	subplot(2,2,iplot+1)
	semilogy(Ns,E[iplot,:],'.',markersize=12)
	plot(Ns,E[iplot,:])
	ylim(1e-16, 1e3)
	xticks(arange(0,Nmax+10,10));
	if iplot>1: xlabel('N')
	yticks(10.**(arange(-15,5,5)));
	if iplot%2==0: ylabel('error')
	grid()
	title(titles[iplot])
show()
