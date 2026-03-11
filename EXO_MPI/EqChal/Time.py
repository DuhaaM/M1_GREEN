import matplotlib.pyplot as plt
import numpy as np



X=[100,200,400,800]
z=[3.5e-10,3e-10,2e-10,0]
print(X)
Z=[4.5e-2,3.1e-1,8.4e-1,1.5,2.6,4,7.3,9.5,13.6,17]
Y=[4.5e-2,3.5e-1,8.9e-1,1.8,3.2,5,7.8,11.2,15.7,19.8]
x=[3.3687e-10,1.19081e-10,4.2126e-11,1.49006e-11]

plt.clf()
plt.plot(X,x,label="Diff_finies")
plt.plot(X,z,label="Courbe de pente -1/2")
plt.ylabel("||uh-uexact||")
plt.xlabel("Nb de cellulles")
plt.title("Erreur du schéma")
plt.legend()
#plt.show()
plt.savefig("Error.png")
