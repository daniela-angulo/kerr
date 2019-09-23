from qutip import * 
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from math import factorial as fac
from IPython.display import HTML
from matplotlib import cm,colors
from matplotlib import animation

N=30
xvec=np.linspace(-5,5,500)
alpha=2.
steps=10

r=0.5
phi2=(1/alpha)
beta=4.*r
gamma=np.sqrt(1-r**2)*2
signal_dm=coherent_dm(N,beta)
p_array=np.zeros(steps+1)
phi_array=np.zeros(steps+1)
n_array=np.zeros(steps+1)

for i in range(steps+1):
	phi1=(1/alpha)*(1.*i/steps)
	phi_array[i]=phi1
	psi_f=np.exp(-1*beta**2/2.)*np.exp(-1*gamma**2/2.)*sum(beta**n/np.sqrt(fac(n))*gamma**m/np.sqrt(fac(m))*tensor(basis(N,n),basis(N,m),coherent(N,np.exp(-1j*phi1*n-1j*phi2*m)*alpha)) for n in range(15) for m in range(15))
	#psi_f=np.exp(-1*beta**2/2.)*sum(beta**n/np.sqrt(fac(n))*tensor(basis(N,n),basis(N,2),coherent(N,np.exp(-1j*phi1*n-1j*phi2*2)*alpha)) for n in range(20))
	rho_f=psi_f*psi_f.dag()
	rho_p=rho_f.ptrace([1,2])
	n=negativity(rho_p,0,method='eigenvalues')
	n_array[i]=n
	p=(rho_p**2).tr()
	p_array[i]=p

np.savetxt('3modespurity.txt',np.array([phi_array,p_array]).T)
np.savetxt('3modesnegativity.txt',np.array([phi_array,n_array]).T)
fig,axes = plt.subplots(1,2,figsize=(10,5))
axes[0].plot(phi_array,p_array,'ro')
axes[1].plot(phi_array,n_array,'ro')
plt.savefig('3mode.png', dpi=fig.dpi)
plt.show()

#print("Trace squared of probe state is %1.3f" %p)
# plt.bar(np.arange(0,N)-0.5,signal_dm.diag())
# plt.figure()
# W=wigner(rho_p,xvec,xvec,'iterative',2)
# wmap=wigner_cmap(W)
# nrm=colors.Normalize(-W.max(),W.max())
# color=plt.contourf(xvec,xvec,W,100,cmap=cm.RdBu,norm=nrm)
# clb=plt.colorbar(color)
# plt.show()