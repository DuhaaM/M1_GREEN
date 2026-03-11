import subprocess
import matplotlib.pyplot as plt






subprocess.run(["mpic++","EulExp1D.cpp"])




core=8                      ## Nombre de coeur utilisés
cells=(100,200,400,800)      ## Nombre de cellules utilisées
Tfinal = 1                  ## Temps maximal
K=1                         ## Un entier "sin(K M_PI x)
kprint=100                  ## Pas d'écriture sur fichier
L=1                         ## Borne supèrieure du domaine 
methode=2                   ## Type de communication 1=Asynchrone 2=SendRecv


TabErr=[]
TabErrInf=[]
TabCell=[]


TabU=[]
Flag = True 
for cell in cells:
      
    

    Commande=["mpirun", "-np",str(core),"a.out",str(cell),str(Tfinal),str(K),str(kprint),str(L),str(methode)]

    result = subprocess.run(Commande,capture_output=True,text=True)
    Erreur = [line for line in result.stdout.splitlines() if "Erreur" in line]
        

    #Assertion de l'erreur si la fonction initiale est affine ou constante 

    if (cell == cells[0] ):

        f = open('Data/Chal0.dat','r')
        lines=f.readlines()
        for line in lines :
            TabU.append(float(line.split()[1]))
        f.close

        for i in range (1,len(TabU)-3):
            if round(TabU[i+1] - TabU[i],15) != round(TabU[i+2] - TabU[i+1],15):
                Flag=False
                break
        if (Flag == True):
            assert round(float(Erreur[0][7:]),15) == 0



            
    TabErr.append(float(Erreur[0][7:]))
    TabErrInf.append(float(Erreur[1][10:]))
    TabCell.append(L/cell)
    
    print(cell,"\n","Resulat:",result.stdout)


    


plt.clf()
plt.plot(TabCell,TabErr,label="Erreur Legendre-Gauss")
plt.title("Étude d'erreur du schéma de différences finies")
plt.ylabel("||uh-uexact||L2")
plt.xlabel("dx")
plt.legend()
plt.savefig("Picture/Erreur Legendre-Gauss")
#plt.show()
 
plt.clf()
plt.plot(TabCell,TabErrInf,label = "Erreur infinie")
plt.title("Étude d'erreur du schéma de différences finies")
plt.ylabel("||uh-uexact||Inf")
plt.xlabel("dx")
plt.legend()
plt.savefig("Picture/Erreur Infinie")
#plt.show()
