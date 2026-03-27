import subprocess
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd






subprocess.run(["mpic++","EulerExp2D.cpp"])

core=6
#Nx=(10,20,40,80)
#Ny=(10,20,40,80)
N=(10,20,40,80)
Tfinal=0.1
K=1
kprint=1
Lx=1
Ly=1
methode= 2
function = 1

TabErrInf=[]
TabCell=[]

for NX in N :
    for NY in Ny:
        Commande= ["mpirun","-np",str(core),"./a.out",str(NX),str(NY),str(Tfinal),str(K),str(kprint),str(Lx),str(Ly),str(methode),str(function)]
        result = subprocess.run(Commande,capture_output=True,text=True)
        Erreur = [line for line in result.stdout.splitlines() if "Erreur" in line]
        TabErrInf.append(float(Erreur[0][10:]))
        TabCell.append(Lx*Ly/(NX*NY))

        print(NX,"\t",NY,"\n",result.stdout)

df=pd.read_csv("Data2D/Chal450.csv")
TabU=df["#Solution_obtenue"].to_numpy()



#uu =TabU.reshape((Ny+1,Nx+1))[::-1,:]
#u = np.zeros((2*(Ny), 2*(Nx)))
#u[0::2, 0::2] = uu[:-1,:-1]
#u[1::2, 1::2] = uu[1:,1:]
#u[1::2, 0::2] = uu[1:,:-1]
#u[0::2, 1::2] = uu[:-1,1:]

TabCell2=np.array(TabCell)**2*1e-3


#plt.clf()
#plt.imshow(u, extent=(0,Lx,0,Ly))
#plt.legend()
#plt.imsave("Picture2D/Solution.png",u)
#plt.show()


plt.clf()
plt.loglog(TabCell,TabErrInf,label = "Erreur infinie")
plt.loglog(TabCell,TabCell2,label="Courbe d'ordre 2")
plt.title("Étude d'erreur du schéma de différences finies")
plt.ylabel("||uh-uexact||Inf")
plt.xlabel("dx")
plt.legend()
plt.savefig("Picture2D/Erreur Infinie")
plt.show()
