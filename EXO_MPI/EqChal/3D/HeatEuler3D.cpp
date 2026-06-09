#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string>
#include <cmath>
#include <mpi.h>
#include "Solution3D.hpp"
#include "Error3D.hpp"
#include "Mesh3D.hpp"


int main(int argc, char *argv[]) {
  int Nx, Ny, Nz, K, kprint, method, Function, rang;
  double Tfinal, Lx, Ly, Lz;
  double temps_debut, temps_fin;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rang);

  if (argc != 12) {
    std::cerr << "usage is : " << argv[0]
              << " <Nx> <Ny> <Nz> <t=final time> <K=integer> <kprint=print step> <Lx> <Ly> <Lz> <Comm_Method : Non-Blocking=1 or Sendrecv=2> <Function : 1=Sinus , 2=Linear , 3=Constant>" << std::endl;
    return 0;
  }

  Nx = std::atoi(argv[1]);
  Ny = std::atoi(argv[2]);
  Nz = std::atoi(argv[3]);
  Tfinal = std::atof(argv[4]);
  K = std::atoi(argv[5]);
  kprint = std::atoi(argv[6]);
  Lx = std::atof(argv[7]);
  Ly = std::atof(argv[8]);
  Lz = std::atof(argv[9]);
  method = std::atoi(argv[10]);
  Function = std::atoi(argv[11]);

  Solution* Sol = nullptr;
  if (Function == 1) Sol = new Sin(K, Lx, Ly, Lz);
  if (Function == 2) Sol = new Lin();
  if (Function == 3) Sol = new Const();

  Mesh Mesh(Nx, Ny, Nz, Lx, Ly, Lz);

  MPI_Comm comm_3D = Mesh.DefComm();
  Mesh.Defneighbor();
  Mesh.computeNbElemtsLoc();
  Mesh.NbLocGhost();
  Mesh.computeDeplacement();
  Mesh.setDebut();

  int* NxLoc = Mesh.NxLoc();
  int* NyLoc = Mesh.NyLoc();
  int* NzLoc = Mesh.NzLoc();

  double* x = Mesh.DefMeshx();
  double* y = Mesh.DefMeshy();
  double* z = Mesh.DefMeshz();

  double dx    = Lx / double(Nx);
  double dy    = Ly / double(Ny);
  double dz    = Lz / double(Nz);
  double hmin  = std::min(dx, std::min(dy, dz));
  double dtloc = 0.1 * hmin * hmin;   // condition CFL 3D : dt*(1/dx²+1/dy²+1/dz²) ≤ 1/2
  double t = 0, time = 0;
  int kT = 0, Ind = 0;

  double* Utemp = Mesh.AllocateMemory();

  MPI_Request request[12];

  double* U = Mesh.Initialize(Sol);

  Comm Comm(NxLoc, NyLoc, NzLoc, comm_3D, Mesh.Neighbor());
  Comm.DefMess();

  while (t < Tfinal) {

    if (t + dtloc > Tfinal) dtloc = Tfinal - t;

    if (rang == 0) temps_debut = MPI_Wtime();

    if (method == 1) {   // Communications asynchrones
      int nb_req = Comm.AsynchroneComm(U, request);
      Mesh.Next(Utemp, U, dtloc);
      Comm.Wait(request, nb_req);
      Mesh.BoundNext(Utemp, U, dtloc);
      Mesh.Copy(Utemp, U);
    }

    if (method == 2) {   // Communications synchrones
      Comm.SynchroneComm(U);
      Mesh.Next(Utemp, U, dtloc);
      Mesh.BoundNext(Utemp, U, dtloc);
      Mesh.Copy(Utemp, U);
    }

    if (rang == 0) {
      temps_fin = MPI_Wtime();
      time += temps_fin - temps_debut;
    }

    if (kT % kprint == 0) {
      Mesh.Write(Ind, U, Sol, t);
      Ind++;
    }

    t += dtloc;
    kT++;
  }

  if (rang == 0) std::cout << "Computing time = " << time << "s" << std::endl;

  Error Error(Lx, Ly, Lz, Nx, Ny, Nz, K, Tfinal, x, y, z,Mesh.Debutx(), Mesh.Debuty(), Mesh.Debutz(),NxLoc, NyLoc, NzLoc);
  Error.ErrorInf(U, Sol);
  Error.ErrorL2(U, Sol);

  MPI_Finalize();
  return 0;
}
