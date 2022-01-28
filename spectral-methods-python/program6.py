'''
Spectral methods in MATLAB.
compare to Trefethen p6.m
'''

# Variable coefficient wave equation


from numpy import *
from mayavi import mlab
from matplotlib import pyplot as plt


#Grid, variable coefficient, and initial data:
N = 128
h = 2*pi/N
x = arange(0,N)*h
t = 0
dt = h/4.0
c = 0.2 + sin(x-1)**2
v = exp(-100*(x-1)**2)
vold = exp(-100*(x-0.2*dt-1)**2)

# Time stepping by leap frog formula:
tmax = 8
tplot = 0.15
plotgap = int(round(tplot/dt))
dt = tplot/plotgap
nplots = int(round(tmax/tplot))
data = vstack((v, zeros((nplots,N))))
tdata = t
for i in range(0,nplots):
	for n in range(0,plotgap):
		t += dt
		v_hat = fft.fft(v)
		w_hat = 1j*fft.fftfreq(N, 1./N)*v_hat#hstack((arange(0,N/2.0),0,arange(-N/2.0,-1)))*(v_hat)
		w = real(fft.ifft(w_hat))
		vnew = vold - 2*dt*c*(w)
		vold = v; v = vnew;
	data[i+1,:] = v
	tdata = vstack((tdata, t))

X, Y = meshgrid(x, tdata)
mlab.mesh(X,Y,data)
mlab.show()
