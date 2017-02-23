import numpy as np 
import scipy as sp
from matplotlib import pyplot as plt 
from scipy.integrate import quad 

def integrandx(xi, n, l, x0, lsx):
	return  np.exp(-((xi-x0)**2 ) / (2*lsx**2) ) * np.sin(np.pi*n*xi/l)

def integrandt(tau, n, l, t0, lst, kappa, t, a):
	lambdan = np.sqrt((a*np.pi*n/l)**2 - ((kappa**2) / 4.))
	return  np.exp(-((tau-t0)**2 ) / (2*lst**2) ) * \
	        np.exp(0.5*kappa*tau)                 * \
	        np.sin(lambdan*(t-tau))

n = 1.
l = 0.75
tmax = 10.

x0 = 0.1*l
t0 = 0.33*tmax

lsx = 0.05
lst = 0.05

a = 2.
kappa = 1.	

Nx = 10
Nt = 100
x = np.linspace(0.0*l,1.0*l,Nx).reshape(-1,1)
t = np.linspace(0,tmax,Nt).reshape(-1,1)

K = 3

W = np.zeros((t.size ,x.size))	

plt.figure(1)
for n in range(1, K+1):
	Wn = np.zeros((t.size ,x.size))	
	Ix = quad(integrandx, 0., l, args=(n, l, x0, lsx))
	for c in range (0, x.size):
		for r in range(0, t.size):
			It = quad(integrandt, 0., t[r], args=(n, l, t0, lst, kappa, t[r], a))
			Wn[r,c] = np.exp(-0.5*kappa*t[r]) * np.sin(np.pi*n*x[c]/l) * Ix[0] * It[0]
	lambdan = np.sqrt((a*np.pi*n/l)**2 - ((kappa**2) / 4.))
	Wn =  (1. / lambdan)* Wn 
	W  = W + Wn
	plt.plot(t,Wn[:,0])
W = (2./l) * W 

plt.figure(2)
plt.plot(t,np.sum(W,1))


plt.figure(3)
for j in range(0,t.size):
	plt.plot(x, W[j,:],'b')
	plt.xlim([0, l])
	plt.ylim([-0.002, 0.002])
	plt.pause(0.05)
	plt.hold(False)



plt.show()
