'''
Spectral methods in MATLAB.
compare to Trefethen p30.m
'''


# Spectral integration, ODE style

from numpy import *
from numpy.linalg import inv
from numpy.polynomial import legendre
from scipy.special import erf
from matplotlib import pyplot as plt
from cheb import *
from clencurt import *
from gauss import *

#Computation: various values of N, four functions:
Nmax = 50
E = zeros((4,Nmax))
for N in range(1,Nmax+1):
    '''
    CHEB
    '''

    # i = slice(0,N)
    # D,x = cheb(N)
    # x = x[i]
    # Di = inv(D[i,i])
    # w = Di[0,:]


    '''
    CLENCURT
    '''
    # x, w = clencurt(N)
    

    '''GAUSS'''
    x,w = gauss(N)
    x2, w2 = legendre.leggauss(N)

    assert allclose(x, x2)
    assert allclose(w2, w)
    f = abs(x)**3
    E[0,N-1] = abs(w.dot(f) - 0.5)
    f = exp(-x**(-2))
    E[1,N-1] = abs(w.dot(f) - 2*(exp(-1)+sqrt(pi)*(erf(1)-1)))
    f = 1.0/(1+x**2)
    E[2,N-1] = abs(w.dot(f) - pi/2)
    f = x**2
    E[3,N-1] = abs(w.dot(f) - 2.0/11)

titles=['|x^3|', 'exp(-x^{-2})', '1/(1+x^2)', 'x^(10)']
for iplot in range(0,4):
    plt.subplot(220+iplot+1)
    plt.semilogy(E[iplot,:]+1e-100,'o')
    plt.plot(E[iplot,:]+1e-100, linestyle='-')
    plt.title(titles[iplot])
    plt.ylabel('error')


plt.show()
