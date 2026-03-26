import subprocess
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd






subprocess.run(["mpic++","EulerExp2D.cpp"])

core=6
Nx=100
Ny=100
Tfinal=0.1
K=1
kprint=1
Lx=1
Ly=1
methode= 2
function = 1



Commande= ["mpirun","-np",str(core),"./a.out",str(Nx),str(Ny),str(Tfinal),str(K),str(kprint),str(Lx),str(Ly),str(methode),str(function)]
result = subprocess.run(Commande,capture_output=True,text=True)

df=pd.read_csv("Data2D/Chal20.csv")
TabU=df["#Solution_obtenue"].to_numpy()



uu =TabU.reshape((Ny+1,Nx+1))[::-1,:]
u = np.zeros((2*(Ny), 2*(Nx)))
u[0::2, 0::2] = uu[:-1,:-1]
u[1::2, 1::2] = uu[1:,1:]
u[1::2, 0::2] = uu[1:,:-1]
u[0::2, 1::2] = uu[:-1,1:]

plt.clf()
plt.imshow(u, extent=(0,Lx,0,Ly))
#plt.legend()
plt.imsave("schema.png",u)
plt.show()
