#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string>
#include <cmath>
#include <mpi.h>

class Error {
private:
  double _Lx, _Ly, _Lz;
  int _NbxElemtsGlob, _NbyElemtsGlob, _NbzElemtsGlob;
  int _k;
  double _Tfinal;
  double* _Mesh0x;
  double* _Mesh0y;
  double* _Mesh0z;
  int _rank;
  int* _NxLoc;
  int* _NyLoc;
  int* _NzLoc;
  int* _debutx;
  int* _debuty;
  int* _debutz;

public:
  Error(const double& Lx, const double& Ly, const double& Lz, const int& Nx, const int& Ny, const int& Nz, const int& k, const double& Tfinal, double* Initx, double* Inity, double* Initz, int* debutx, int* debuty, int* debutz, int* NxLoc, int* NyLoc, int* NzLoc)
    : _Lx(Lx), _Ly(Ly), _Lz(Lz), _NbxElemtsGlob(Nx), _NbyElemtsGlob(Ny), _NbzElemtsGlob(Nz), _k(k), _Tfinal(Tfinal), _Mesh0x(Initx), _Mesh0y(Inity), _Mesh0z(Initz), _debutx(debutx), _debuty(debuty), _debutz(debutz), _NxLoc(NxLoc), _NyLoc(NyLoc), _NzLoc(NzLoc) {
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
  }

  void ErrorInf(const double* Vect, Solution* Sol) {
    double Errloc = 0, Errglob = 0;
    const int nx = _NxLoc[_rank], ny = _NyLoc[_rank], nz = _NzLoc[_rank];
    const int sy = nx + 2;
    const int sz = (nx + 2) * (ny + 2);

    for (int k = 1; k <= nz; k++) {
      for (int j = 1; j <= ny; j++) {
        for (int i = 1; i <= nx; i++) {
          int I = i + j*sy + k*sz;
          double err = std::abs(Vect[I] - Sol->Exact( _Mesh0x[_debutx[_rank]+i-1], _Mesh0y[_debuty[_rank]+j-1], _Mesh0z[_debutz[_rank]+k-1], _Tfinal));
          Errloc = std::max(Errloc, err);
        }
      }
    }
    MPI_Reduce(&Errloc, &Errglob, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (_rank == 0) std::cout << "ErreurInf = " << Errglob << std::endl;
  }

  void ErrorL2(const double* Vect, Solution* Sol) {
    const double xi[2] = { -1.0/sqrt(3.0), 1.0/sqrt(3.0) };
    const double w[2]  = { 1.0, 1.0 };

    const double dx = _Mesh0x[1] - _Mesh0x[0];
    const double dy = _Mesh0y[1] - _Mesh0y[0];
    const double dz = _Mesh0z[1] - _Mesh0z[0];

    const int nx = _NxLoc[_rank], ny = _NyLoc[_rank], nz = _NzLoc[_rank];
    const int sy = nx + 2;
    const int sz = (nx + 2) * (ny + 2);

    double ErrlocL2 = 0.0, ErrglobL2 = 0.0;

    for (int k = 1; k < nz; k++) {
      for (int j = 1; j < ny; j++) {
        for (int i = 1; i < nx; i++) {
          double x0 = _Mesh0x[_debutx[_rank]+i-1], x1 = _Mesh0x[_debutx[_rank]+i];
          double y0 = _Mesh0y[_debuty[_rank]+j-1], y1 = _Mesh0y[_debuty[_rank]+j];
          double z0 = _Mesh0z[_debutz[_rank]+k-1], z1 = _Mesh0z[_debutz[_rank]+k];

          double cell_sum = 0.0;
          for (int p = 0; p < 2; p++) {
            for (int q = 0; q < 2; q++) {
              for (int s = 0; s < 2; s++) {
                double xg = 0.5*(x0+x1) + 0.5*dx*xi[p];
                double yg = 0.5*(y0+y1) + 0.5*dy*xi[q];
                double zg = 0.5*(z0+z1) + 0.5*dz*xi[s];

                double r  = (xg - x0) / dx;
                double rv = (yg - y0) / dy;
                double t  = (zg - z0) / dz;

                double uh = (1-r)*(1-rv)*(1-t) * Vect[ i + j *sy +  k *sz] + r *(1-rv)*(1-t) * Vect[(i+1) + j * sy + k *sz] + (1-r)* rv *(1-t) * Vect[ i + (j+1) * sy + k *sz] + r * rv *(1-t) * Vect[(i+1)  + (j+1)*sy +  k *sz] + (1-r)*(1-rv)* t * Vect[ i + j *sy + (k+1) *sz] + r *(1-rv)* t * Vect[(i+1) + j *sy + (k+1) *sz] + (1-r)* rv * t * Vect[ i + (j+1) *sy + (k+1) *sz] + r * rv * t * Vect[(i+1) + (j+1)*sy + (k+1) *sz];

                double uex  = Sol->Exact(xg, yg, zg, _Tfinal);
                double diff = uh - uex;
                cell_sum += w[p] * w[q] * w[s] * diff * diff;
              }
            }
          }
          ErrlocL2 += cell_sum * dx * dy * dz / 8.0;
        }
      }
    }

    MPI_Reduce(&ErrlocL2, &ErrglobL2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (_rank == 0) std::cout << "ErreurL2 = " << sqrt(ErrglobL2) << std::endl;
  }
};
