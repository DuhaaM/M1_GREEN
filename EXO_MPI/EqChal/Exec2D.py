import subprocess
import matplotlib.pyplot as plt
import numpy as np






subprocess.run(["g++","EulExp2D.cpp"])


Nx=10
Ny=10
Tfinal=0.1
K=1
kprint=1
Lx=1
Ly=1

TabU=[]

Commande= ["./a.out",str(Nx),str(Ny),str(Tfinal),str(K),str(kprint),str(Lx),str(Ly)]
result = subprocess.run(Commande,capture_output=True,text=True)

f=open('Data2D/ChalFinal.dat','r')
lines=f.readlines()
for line in lines :
    TabU.append(float(line))
f.close



uu = np.array(TabU).reshape((Ny+1,Nx+1))[::-1,:]
u = np.zeros((2*(Ny), 2*(Nx)))
u[0::2, 0::2] = uu[:-1,:-1]
u[1::2, 1::2] = uu[1:,1:]
u[1::2, 0::2] = uu[1:,:-1]
u[0::2, 1::2] = uu[:-1,1:]

plt.clf()
plt.imshow(u, extent=(0,Lx,0,Ly))
#plt.legend()
plt.show()
