# Author : Koshvendra Singh
# Date   : 12/06/2020
# Email  : koshvendra.singh@tifr.res.in
# Description : Numerical solution of schrodinger equation for V = (x**4-3)*exp(-x*x/2) in a box  [-5,5]
# So  b.c are    w.f.(-5) = w.f.(5) = 0   , so only N-2 grid points would be analysed
# Method : finite differnece method , writing the matrix of hamiltonian and then finding it's eigenvalue

import numpy as np
import matplotlib.pyplot as plt

def kin(xmin,xmax,N=500):
    K = np.zeros((N-2)**2).reshape(N-2,N-2)
    for i in range(N-2):
        for j in range(N-2):
            if i==j:
                K[i,j]=-2
            elif np.absolute(i-j) == 1:
                K[i,j] = 1
            else :
                K[i,j] =0
    return K

def potential(x):
    return (x**4 - 3)*np.exp(-x*x/2)

def pot(xmin,xmax,N=500):
    x=np.linspace(xmin,xmax,N)
    dx=x[1]-x[0]
    V = np.zeros((N-2)**2).reshape(N-2,N-2)
    for i in range(N-2):
        for j in range(N-2):
            if i==j:
                V[i,j] = potential(x[i+1])
            else:
                V[i,j] = 0
    return V

#hamiltonian and eigenvalue
def H(xmin,xmax,N=500):
    dx=(xmax-xmin)/(N-1)
    H=-0.5*kin(xmin,xmax)/(dx**2) + pot(xmin,xmax)
    eigenval =np.sort( np.linalg.eigvals(H) )
    return (eigenval)

xmin=-5
xmax=5
k=H(xmin,xmax)
print('the first five energy levels are',k[0:5])

energy_mat=np.zeros([20,40])
for i in range(40):
    xmax = xmax + i*0.50
    xmin = xmin - i*0.50
    energy=H(xmin,xmax)
    energy_mat[:,i]=energy[0:20]
   

for i in range(20):
    p=np.linspace(10,50,40)
    plt.plot(p,energy_mat[i,:],label=i)

plt.legend()
plt.xlabel('size of box')
plt.ylabel('energy')
plt.title('20 energy levels for various size of box')
plt.show()

for i in range(10):
    p=np.linspace(10,50,40)
    plt.plot(p,energy_mat[2*i,:],label=i*2)

plt.legend()
plt.xlabel('size of box')
plt.ylabel('energy')
plt.title('even parity levels')
plt.show()

for i in range(10):
    p=np.linspace(10,50,40)
    plt.plot(p,energy_mat[2*i+1,:],label=2*i+1)

plt.legend()
plt.xlabel('size of box')
plt.ylabel('energy')
plt.title('Odd parity levels')
plt.show()

