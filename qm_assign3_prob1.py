# Author : Koshvendra Singh
# Date   : 10/06/2020
# Description : Solving schrodinger equation for potential V=mga(1-cos(phi)) for quantum pendulum

import numpy as np
import matplotlib.pyplot as plt

# hamiltonian  H = (h^2* L^2)/(8*I*pi^2)   L=-id/dphi
#h_bar = 6.625 * ( 10 **(-34))/(4*(np.pi**2))
d = (1/600)                   # 1/(2I)

num = 1500                         # no.of grid points
h= (2*np.pi -0)/(num-1)              # step size
phi = np.arange(0,2*np.pi,h)

# boundary condition  ,  w.f(0)= w.f(2*np.pi)

# defining kinetic energy matrix
K = np.zeros((num-1)**2).reshape(num-1,num-1)
for i in range(num-1):
    for j in range(num-1):
        
        if i==j :
            K[i,j]=-2
        elif np.absolute(i-j)==1:
            K[i,j]=1
        else:
            K[i,j]=0
K[0,num-2]=1
K[num-2,0]=1

# potential matrix
def potential(phi):
    return(1-np.cos(phi))               

V=np.zeros((num-1)**2).reshape(num-1,num-1)
for i in range(num-1):
    for j in range(num-1):
        if i==j:
            V[i,j]=potential(phi[i])
        else:
            V[i,j]=0

# Hamiltonian matrix
H = -K*d/(h**2)  + V
            
# solving the hamiltonian matrix to get it's eigenvalues
val =np.linalg.eigvals(H)
val_1=np.sort(val)
val_2=val_1/(val_1[0])           # energies normalised by grount state energy

eig_val_300=val_2[0:300]        # first 300 eigen energies
x=np.arange(0,300,1)

# function for finding starting point of degenracy
def degenracy(arr)
    i=0
    a=[]
    while arr[i] != arr[i+1]:
        a.append(i)
        i=i+1
    return (len(a))

# plotting
plt.plot(x,eig_val_300)
plt.xlabel('excited states')
plt.ylabel('energy of states(normalized by ground state)')
plt.title('energy eigenvalues of quantum pendulum for fist 300 states')
plt.show()

print(val_1[0],val_1[49],val_1[200])
print(degenracy(eig_val_300))
