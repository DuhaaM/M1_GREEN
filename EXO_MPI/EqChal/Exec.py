import subprocess
import matplotlib.pyplot as plt
import numpy as np






subprocess.run(["mpic++","EulerExp1D.cpp"])




core=6                      ## Nombre de coeur utilisés
cells=(100,200,400,800)      ## Nombre de cellules utilisées
Tfinal = 1                  ## Temps maximal
K=1                         ## Un entier "sin(K M_PI x)
kprint=100                  ## Pas d'écriture sur fichier
L=1                         ## Borne supèrieure du domaine 
methode=2                   ## Type de communication 1=Asynchrone 2=SendRecv


TabErr=[]
TabErrInf=[]
TabCell=[]

for cell in cells:
      
    

    Commande=["mpirun", "-np",str(core),"a.out",str(cell),str(Tfinal),str(K),str(kprint),str(L),str(methode)]

    result = subprocess.run(Commande,capture_output=True,text=True)
    Erreur = [line for line in result.stdout.splitlines() if "Erreur" in line]
        

            
    TabErr.append(float(Erreur[2][7:]))
    TabErrInf.append(float(Erreur[1][10:]))
    TabCell.append(L/cell)
   
    print(cell,"\n",result.stdout)


    


TabCell2=np.array(TabCell)**2*1e-3

plt.clf()
plt.loglog(TabCell,TabErr,label="Erreur Legendre-Gauss")
plt.loglog(TabCell,TabCell2,label="Courbe d'ordre 2")
plt.title("Étude d'erreur du schéma de différences finies")
plt.ylabel("||uh-uexact||L2")
plt.xlabel("dx")
plt.legend()
plt.savefig("Picture/Erreur Legendre-Gauss")
plt.show()
 
plt.clf()
plt.loglog(TabCell,TabErrInf,label = "Erreur infinie")
plt.loglog(TabCell,TabCell2,label="Courbe d'ordre 2")
plt.title("Étude d'erreur du schéma de différences finies")
plt.ylabel("||uh-uexact||Inf")
plt.xlabel("dx")
plt.legend()
plt.savefig("Picture/Erreur Infinie")
plt.show()
