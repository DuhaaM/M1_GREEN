import subprocess
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


subprocess.run(["nvcc","HeatEuler2D.cu","-o","HEcuda"])
subprocess.run(["mpic++","-O3","../EXO_MPI/EqChal/2D/HeatEuler2D.cpp","-o","HEcpp"])
        
N=500
Tfinal=0.1
k=1
kprint=100
Lx=1
Ly=1
methode=1
function=1
Bdim=(1,2,4,8)

TabTimeCuda=[]
TabTimeCpp=[]
TabB=[]


for B in Bdim:

    Commande=["./HEcuda",str(N),str(N),str(Tfinal),str(k),str(kprint),str(Lx),str(Ly),str(function),str(B)]
    result=subprocess.run(Commande,capture_output=True,text=True)
    Time=[line for line in result.stdout.splitlines() if "time" in line]
    TabTimeCuda.append(float(Time[0][16:]))        
    
    
    print("CUDA",result.stdout)

    Commande=["mpirun","-np",str(B),"./HEcpp",str(N),str(N),str(Tfinal),str(k),str(kprint),str(Lx),str(Ly),str(methode),str(function)]
    result=subprocess.run(Commande,capture_output=True,text=True)
    Time=[line for line in result.stdout.splitlines() if "time" in line]
    TabTimeCpp.append(float(Time[0][17:]))

    print("Cpp",result.stdout)

    TabB.append(B)


plt.clf()
plt.plot(TabB,TabTimeCpp,label="CPU")
plt.plot(TabB,TabTimeCuda,label="GPU")
plt.legend()
plt.title("Temps de calcul")
plt.xlabel("sqrt( Nombre de thread par block) / core")
plt.ylabel("Temps en secondes")
plt.savefig("TimeHeat")
plt.show()
