from qutip import * 
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from IPython.display import HTML
from matplotlib import cm,colors
from matplotlib import animation
import time as time
starttime = time.time()

N=30
xvec=np.linspace(-5,5,500)
t_steps=240
tlist=np.linspace(0.0,240.0, t_steps)
alpha=2.0
probe0=coherent_dm(N, alpha)
beta=2.
signal0=coherent_dm(N,beta)
#Field operators
a2=tensor(identity(N),destroy(N))
a1=tensor(destroy(N),identity(N))
#ap=tensor(identity(N),identity(N),destroy(N))
n1=a1.dag()*a1
n2=a2.dag()*a2

#Hamiltonian
phi=(1/alpha)*np.arange(10)/10.
H=[phi[i]*n1*n2 for i in range(10)]
U=[(-1j*H[i]).expm() for i in range(10)]
rho0=tensor(signal0,probe0)
rho_f=[U[i]*rho0*U[i].dag() for i in range(10)]
n=[negativity(rho_f[i],0,method='eigenvalues') for i in range(10)]
p=[((rho_f[i].ptrace(1))**2).tr() for i in range(10)]


fig,axes = plt.subplots(1,2,figsize=(10,5))
axes[0].plot(phi,p,'bo')
axes[1].plot(phi,n,'bo')
plt.savefig('2modeother.png', dpi=fig.dpi)
stoptime = time.time()
print("Program took", int(stoptime-starttime),"seconds")
plt.show()
# W=wigner(rho_p,xvec,xvec)
# wmap=wigner_cmap(W)
# nrm=colors.Normalize(-W.max(),W.max())
# color=plt.contourf(xvec,xvec,W,100,cmap=cm.RdBu,norm=nrm)
# clb=plt.colorbar(color)
# plt.show()