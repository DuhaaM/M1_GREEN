import subprocess
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd



NbExec=100


#subprocess.run(["mpic++","Eulerexp2D.cpp"])
subprocess.run(["mpic++","-O3","HeatEuler3DScale.cpp"])

cores=(1,2,4,8,16,32)
Nx=(32,64,128,256)
Ny=(64,128,256,512)
Tfinal=0.1
K=1
kprint=100
Lx=1
Ly=1
Lz=1
methode= 1
function = 1

for i in (1,2):
    plt.clf()
    for NX in Nx:
    
        TabTime=np.zeros(len(cores))
        TabRef=[]
        if i==1:
            TabCore=[]
            TabMax1=np.zeros(len(cores))
            TabMin1=[float('inf') for _ in range (len(cores))]
            for ie,core in enumerate(cores) :

                cell = NX*np.sqrt(core)

                for k in range(NbExec):

                    Commande= ["mpirun","-np",str(core),"--map-by","socket:PE=1","--bind-to","core","./a.out",str(cell),str(cell),str(cell),str(Tfinal),str(K),str(kprint),str(Lx),str(Ly),str(Lz),str(methode),str(function),str(i)]
                   
                    result = subprocess.run(Commande,capture_output=True,text=True)
                    Erreur = [line for line in result.stdout.splitlines() if "Erreur" in line]
                    Time = [line for line in result.stdout.splitlines() if "time" in line]
                    
                    t= float(Time[0][17:])
                    TabTime[ie]+=t
                    if t > TabMax1[ie]:
                        TabMax1[ie] = t
                    if t < TabMin1[ie] :
                        TabMin1[ie] = t
                    
                    print(NX,"\n",result.stdout)
                
                TabCore.append(core)
            Tabtime= [T / NbExec for T in TabTime]
            time = [ T / Tabtime[0] for T in Tabtime]
            TabMarge1=[abs(TabMax1[j] - TabMin1[j]) for j in range(len(cores))] 
            
            
            plt.errorbar(TabCore,time,yerr=TabMarge1, label = str(NX),capsize = 4)
            
        if i==2 :
            TabCore=[]
            TabMax2=np.zeros(len(cores))
            TabMin2=[float('inf') for _ in range (len(cores))]
            
            for ie,core in enumerate(cores) :
                cell = NX
                for k in range(NbExec):
                            
                    Commande= ["mpirun","-np",str(core),"--map-by","socket:PE=1","--bind-to","core","./a.out",str(cell),str(cell),str(cell),str(Tfinal),str(K),str(kprint),str(Lx),str(Ly),str(Lz),str(methode),str(function),str(i)]
                    
                    result = subprocess.run(Commande,capture_output=True,text=True)
        
                    Erreur = [line for line in result.stdout.splitlines() if "Erreur" in line]
                    Time = [line for line in result.stdout.splitlines() if "time" in line]
        
                    t= float(Time[0][17:])
                    TabTime[ie]+=t
                    if t > TabMax2[ie]:
                        TabMax2[ie] = t
                    if t < TabMin2[ie] :
                        TabMin2[ie] = t
        
                    print(NX,"\n",result.stdout)
                TabCore.append(core)
            
            Tabtime= [T / NbExec for T in TabTime]
            time = [ T / Tabtime[0] for T in Tabtime]
            TabMarge2=[abs(TabMax2[j] - TabMin2[j]) for j in range(len(cores))] 
        
            plt.errorbar(TabCore,time,yerr=TabMarge2, label = str(NX),capsize = 4)
            
    if i==1 :
        plt.axhline(y=1.0, linestyle='--',color='black',label="Scalabilité idéale")
        plt.title("Étude de scalabilité faible")
        plt.ylabel("Temps normalisé")
        plt.xlabel("Cores")
        plt.legend()
        plt.savefig("CourbeFaible3D")
                
        
    if i == 2 :
        TabRef=[cores[0]/c for c in cores ]
        plt.loglog(TabCore,TabRef,label="Scalabilité idéale")
        plt.title("Étude de scalabilité forte")
        plt.ylabel("Temps normalisé")
        plt.xlabel("Cores")
        plt.legend()
        plt.savefig("CourbeForte3D")
            

