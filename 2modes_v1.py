from qutip import * 
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from math import factorial as fac
from IPython.display import HTML
from matplotlib import cm,colors
from matplotlib import animation
import time as time
starttime = time.time()

N=30
xvec=np.linspace(-5,5,500)
steps=10
alpha=2.
beta=2.
signal_dm=coherent_dm(N,beta)
p_array=np.zeros(steps+1)
phi_array=np.zeros(steps+1)
n_array=np.zeros(steps+1)

for i in range(steps+1):
	phi0=(1/alpha)*(1.*i/steps)
	phi_array[i]=phi0
	psi_f=np.exp(-1*beta**2/2.)*sum(beta**n/np.sqrt(fac(n))*tensor(basis(N,n),coherent(N,np.exp(-1j*phi0*n)*alpha)) for n in range(20))
	rho_f=psi_f*psi_f.dag()
	n=negativity(rho_f,0,method='eigenvalues')
	n_array[i]=n	
	rho_p=rho_f.ptrace(1)
	p=(rho_p**2).tr()
	p_array[i]=p

np.savetxt('2modespurity.txt',np.array([phi_array,p_array]).T)
np.savetxt('2modesnegativity.txt',np.array([phi_array,n_array]).T)
fig,axes = plt.subplots(1,2,figsize=(10,5))
axes[0].plot(phi_array,p_array,'bo')
axes[1].plot(phi_array,n_array,'bo')
plt.savefig('2mode.png', dpi=fig.dpi)
stoptime = time.time()
print("Program took", int(stoptime-starttime),"seconds")
plt.show()

# psi_f=np.exp(-1*beta**2/2.)*sum(beta**n/np.sqrt(fac(n))*tensor(basis(N,n),coherent(N,np.exp(-1j*0.0*n)*alpha)) for n in range(20))
# rho_f=psi_f*psi_f.dag()
# rho_p=rho_f.ptrace(1)
# p=(rho_f).tr()
# print p
#plt.bar(np.arange(0,N)-0.5,signal_dm.diag())
#plt.figure()
# W=wigner(rho_p,xvec,xvec,'iterative',2)
# wmap=wigner_cmap(W)
# nrm=colors.Normalize(-W.max(),W.max())
# color=plt.contourf(xvec,xvec,W,100,cmap=cm.RdBu,norm=nrm)
# clb=plt.colorbar(color)
# plt.show()