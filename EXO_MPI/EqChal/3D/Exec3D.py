import subprocess
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import meshio

subprocess.run(["mpic++","-O3", "HeatEuler3D.cpp"])

core   = 8
Nx     = (10, 20, 40)
Nz     = (10, 20, 40)
Tfinal = 0.5
K      = 1
kprint = 10
Lx     = 1
Ly     = 1
Lz     = 1
methode  = 1
function = 2

TabErrInf = []
TabCell   = []
TabErrL2  = []

for NX in Nx:
    Commande = ["mpirun", "-np", str(core), "./a.out",
                str(NX), str(NX), str(NX),
                str(Tfinal), str(K), str(kprint),
                str(Lx), str(Ly), str(Lz),
                str(methode), str(function)]
    result = subprocess.run(Commande, capture_output=True, text=True)
    Erreur = [line for line in result.stdout.splitlines() if "Erreur" in line]
    TabErrInf.append(float(Erreur[0].split('= ')[-1]))
    TabErrL2.append(float(Erreur[1].split('= ')[-1]))
    TabCell.append(Lx / NX)

    print(NX, "\n", result.stdout)

df   = pd.read_csv("Data3D/Chal300.csv")
TabU = df["#Solution_obtenue"].to_numpy()
coords_x=  df["#Point_x"].to_numpy()
coords_y=  df["#Point_y"].to_numpy()
coords_z=  df["#Point_z"].to_numpy()

coords = np.stack([coords_x, coords_y, coords_z], axis=1)

N     = Nx[-1]
k_mid = N // 2
slice_2d = TabU[k_mid * (N+1)**2 : (k_mid+1) * (N+1)**2].reshape((N+1, N+1))[::-1, :]

TabCell2 = np.array(TabCell) ** 2

plt.clf()
plt.imshow(slice_2d, extent=(0, Lx, 0, Ly))
plt.colorbar(label="u")
plt.title(f"Solution 3D — coupe z = {k_mid * Lz / N:.2f}")
plt.imsave("Picture3D/Solution.png", slice_2d)
plt.show()

plt.clf()
plt.loglog(TabCell, TabErrInf, label="Erreur infinie")
plt.loglog(TabCell, TabErrL2,  label="Erreur L2")
plt.loglog(TabCell, TabCell2,  label="Courbe d'ordre 2")
plt.title("Étude d'erreur du schéma de différences finies 3D")
plt.ylabel("||uh-uexact||")
plt.xlabel("dx")
plt.legend()
plt.savefig("Picture3D/Erreur_Infinie_L2")
plt.show()

n_nodes = len(TabU)
N1      = round(n_nodes ** (1/3))
N_elem  = N1 - 1
IZ, IY, IX = np.meshgrid(np.arange(N_elem), np.arange(N_elem), np.arange(N_elem), indexing='ij')
IZ, IY, IX = IZ.ravel(), IY.ravel(), IX.ravel()
cells_hex = np.stack([
    IX     + IY    *N1 + IZ    *N1*N1,
    (IX+1) + IY    *N1 + IZ    *N1*N1,
    (IX+1) + (IY+1)*N1 + IZ    *N1*N1,
    IX     + (IY+1)*N1 + IZ    *N1*N1,
    IX     + IY    *N1 + (IZ+1)*N1*N1,
    (IX+1) + IY    *N1 + (IZ+1)*N1*N1,
    (IX+1) + (IY+1)*N1 + (IZ+1)*N1*N1,
    IX     + (IY+1)*N1 + (IZ+1)*N1*N1,
], axis=1)

meshio.write("EulExp3D.xdmf",
             meshio.Mesh(
                 points=coords,
                 cells=[("hexahedron", cells_hex)],
                 point_data={'Solution3D': TabU}))
