# Author : Koshvendra Singh
# Date   : 12/06/2020
# Description : Solving Schrodinger equation for periodic boundary condition
# Email  : koshvendra.singh@tifr.res.in

import numpy as np
import matplotlib.pyplot as plt

#defining potential fourier transform matrix
def pot_ft(xmin,xmax,n=256):
   # def pot(x):
        #return (x**4 - 3)*np.exp(-x*x/2)
    x=np.linspace(xmin,xmax,2*n)
    dx=x[1]-x[0]
    v=(x**4 - 3)*np.exp(-x*x/2)
    #for i in x:
        #v.append(pot(x[i]))
    v=np.fft.fftshift(np.fft.fft(v))
    freq=2*np.pi*np.fft.fftshift(np.fft.fftfreq(2*n,dx))
    factor = ( dx/np.sqrt(2*np.pi))*np.exp(-1j*freq*xmin)
    v=factor*v
    #print(len(v))
    v_mat=np.zeros([n-1,n-1])
    for i in range(n-1):
        for j in range(n-1):
            v_mat[i,j]=np.real(v[i-j+n])
    return (1/np.sqrt(xmax-xmin))*v_mat

# kinetic energy matrix h_bar=1  ,m=0.1
def kin(xmin,xmax,k,nk=20,n=256):
    x=np.linspace(xmin,xmax,2*n)
    dx=x[1]-x[0]
    freq=2*np.pi*np.fft.fftshift(np.fft.fftfreq(2*n,dx))
    g=freq[int(n/2)+1:int(3*n/2)]
    k_mat=5*np.diag((k+g)**2)
    return k_mat

# dispersion relation
def dispersion(xmin,xmax,n=30,m=256):
    k=np.linspace(-np.pi/(xmax-xmin),np.pi/(xmax-xmin),n)
    energy=np.zeros([n,m-1])
    for i in range(n):
        hamiltonian=kin(xmin,xmax,k[i])+pot_ft(xmin,xmax)
        e=np.linalg.eigvalsh(hamiltonian)
        energy[i,:]=np.sort(e)
        
    for i in range(7):
        plot= plt.plot(k,energy[:,i],label=i)
    return (energy,plot)
        
E,plot=dispersion(-5,5)
plt.xlabel('k')
plt.ylabel('energy')
plt.title('dispersion relation plot')
plt.legend()
plt.show()

#first five band gap
for i in range(5):
    band_gap=[]
    band=np.min(E[:,i+1])-np.max(E[:,i])
    #band_gap.append(band)
    print(band)
#print('first five band gap',band_gap)
        
    


    
    
    
    
