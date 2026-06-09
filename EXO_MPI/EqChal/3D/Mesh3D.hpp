#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string>
#include <cmath>
#include <mpi.h>


class Comm {
private:
  int* _NyLoc;
  int* _NxLoc;
  int* _NzLoc;                                   
  int _rank;
  MPI_Comm _comm_3D;
  int* _neighbor;
  MPI_Datatype _MPI_XFACE, _MPI_YFACE, _MPI_ZFACE;

public:
  Comm(int* NxLoc, int* NyLoc, int* NzLoc, MPI_Comm Comm_3D, int* neighbor)
    : _NxLoc(NxLoc), _NyLoc(NyLoc), _NzLoc(NzLoc), _comm_3D(Comm_3D), _neighbor(neighbor) {
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
  }

  void DefMess() {
    int nx = _NxLoc[_rank], ny = _NyLoc[_rank], nz = _NzLoc[_rank];
    int sy = nx + 2;
    int sz = (nx + 2) * (ny + 2);

    MPI_Datatype col_type;
    MPI_Type_vector(ny, 1, sy, MPI_DOUBLE, &col_type);
    MPI_Type_commit(&col_type);
    MPI_Aint lb, extent;
    MPI_Type_get_extent(MPI_DOUBLE, &lb, &extent);
    MPI_Type_create_hvector(nz, 1, (MPI_Aint)sz * extent, col_type, &_MPI_XFACE);
    MPI_Type_commit(&_MPI_XFACE);
    MPI_Type_free(&col_type);

    MPI_Type_vector(nz, nx, sz, MPI_DOUBLE, &_MPI_YFACE);
    MPI_Type_commit(&_MPI_YFACE);

    MPI_Type_vector(ny, nx, sy, MPI_DOUBLE, &_MPI_ZFACE);
    MPI_Type_commit(&_MPI_ZFACE);
  }

  int AsynchroneComm(double* x, MPI_Request* request) {
    int nx = _NxLoc[_rank], ny = _NyLoc[_rank], nz = _NzLoc[_rank];
    int sy = nx + 2;
    int sz = (nx + 2) * (ny + 2);
    int etiquette = 100;
    int nb_req = 0;

    for (int i = 0; i < 12; i++) request[i] = MPI_REQUEST_NULL;

    // Direction x (Left=0, Right=1)
    MPI_Irecv(&x[0 + 1*sy + 1*sz], 1, _MPI_XFACE, _neighbor[0], etiquette, _comm_3D, &request[nb_req++]);
    MPI_Irecv(&x[(nx+1) + 1*sy + 1*sz], 1, _MPI_XFACE, _neighbor[1], etiquette, _comm_3D, &request[nb_req++]);
    MPI_Isend(&x[1 + 1*sy + 1*sz], 1, _MPI_XFACE, _neighbor[0], etiquette, _comm_3D, &request[nb_req++]);
    MPI_Isend(&x[nx + 1*sy + 1*sz], 1, _MPI_XFACE, _neighbor[1], etiquette, _comm_3D, &request[nb_req++]);

    // Direction y (Bottom=2, Top=3)
    MPI_Irecv(&x[1 + 0*sy + 1*sz], 1, _MPI_YFACE, _neighbor[2], etiquette, _comm_3D, &request[nb_req++]);
    MPI_Irecv(&x[1 + (ny+1)*sy + 1*sz], 1, _MPI_YFACE, _neighbor[3], etiquette, _comm_3D, &request[nb_req++]);
    MPI_Isend(&x[1 + 1*sy + 1*sz], 1, _MPI_YFACE, _neighbor[2], etiquette, _comm_3D, &request[nb_req++]);
    MPI_Isend(&x[1 + ny*sy + 1*sz], 1, _MPI_YFACE, _neighbor[3], etiquette, _comm_3D, &request[nb_req++]);

    // Direction z (Front=4, Back=5)
    MPI_Irecv(&x[1 + 1*sy + 0*sz ], 1, _MPI_ZFACE, _neighbor[4], etiquette, _comm_3D, &request[nb_req++]);
    MPI_Irecv(&x[1 + 1*sy + (nz+1)*sz ], 1, _MPI_ZFACE, _neighbor[5], etiquette, _comm_3D, &request[nb_req++]);
    MPI_Isend(&x[1 + 1*sy + 1*sz ], 1, _MPI_ZFACE, _neighbor[4], etiquette, _comm_3D, &request[nb_req++]);
    MPI_Isend(&x[1 + 1*sy + nz*sz ], 1, _MPI_ZFACE, _neighbor[5], etiquette, _comm_3D, &request[nb_req++]);

    return nb_req;
  }

  void Wait(MPI_Request* request, int nb_req) {
    MPI_Waitall(nb_req, request, MPI_STATUS_IGNORE);
  }

  void SynchroneComm(double* x) {
    int nx = _NxLoc[_rank], ny = _NyLoc[_rank], nz = _NzLoc[_rank];
    int sy = nx + 2;
    int sz = (nx + 2) * (ny + 2);
    int etiquette = 100;

    // Direction x
    MPI_Sendrecv(&x[1 + 1*sy + 1*sz], 1, _MPI_XFACE, _neighbor[0], etiquette,&x[(nx+1) + 1*sy + 1*sz], 1, _MPI_XFACE, _neighbor[1], etiquette, _comm_3D, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&x[nx + 1*sy + 1*sz], 1, _MPI_XFACE, _neighbor[1], etiquette,&x[0 + 1*sy + 1*sz], 1, _MPI_XFACE, _neighbor[0], etiquette, _comm_3D, MPI_STATUS_IGNORE);

    // Direction y
    MPI_Sendrecv(&x[1 + 1*sy + 1*sz], 1, _MPI_YFACE, _neighbor[2], etiquette, &x[1 + (ny+1)*sy + 1*sz], 1, _MPI_YFACE, _neighbor[3], etiquette, _comm_3D, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&x[1 + ny*sy + 1*sz], 1, _MPI_YFACE, _neighbor[3], etiquette, &x[1 + 0*sy + 1*sz], 1, _MPI_YFACE, _neighbor[2], etiquette, _comm_3D, MPI_STATUS_IGNORE);

    // Direction z
    MPI_Sendrecv(&x[1 + 1*sy + 1*sz ], 1, _MPI_ZFACE, _neighbor[4], etiquette, &x[1 + 1*sy + (nz+1)*sz], 1, _MPI_ZFACE, _neighbor[5], etiquette, _comm_3D, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&x[1 + 1*sy + nz*sz ], 1, _MPI_ZFACE, _neighbor[5], etiquette, &x[1 + 1*sy + 0*sz ], 1, _MPI_ZFACE, _neighbor[4], etiquette, _comm_3D, MPI_STATUS_IGNORE);
  }
  
};


class Mesh {
private:
  double _Lx;
  double _Ly;
  double _Lz;
  int _NbxElemtsGlob;
  int _nx;
  int _NbyElemtsGlob;
  int _ny;
  int _NbzElemtsGlob;
  int _nz;
  int* _debutx;
  int* _debuty;
  int* _debutz;
  int* _NbElemtsLoc;
  int* _NyLoc;
  int* _NxLoc;
  int* _NzLoc;
  int* _deplacement;
  int _nbp;
  int _rank;
  double* _Init;
  double* _Mesh0x;
  double* _Mesh0y;
  double* _Mesh0z;
  int* _dims;
  int _neighbor[6];
  MPI_Comm comm_3D;
  int* _coords;
  int* _NbLocGhost;

  // Retourne vrai si la cellule (i,j,k) touche une frontière physique du domaine global
  bool isPhysicalBoundary(int i, int j, int k) const {
    if (i == 1 && _neighbor[0] == MPI_PROC_NULL) return true;
    if (i == _NxLoc[_rank] && _neighbor[1] == MPI_PROC_NULL) return true;
    if (j == 1 && _neighbor[2] == MPI_PROC_NULL) return true;
    if (j == _NyLoc[_rank] && _neighbor[3] == MPI_PROC_NULL) return true;
    if (k == 1 && _neighbor[4] == MPI_PROC_NULL) return true;
    if (k == _NzLoc[_rank] && _neighbor[5] == MPI_PROC_NULL) return true;
    return false;
  }

public:
  Mesh(const int& Nx, const int& Ny, const int& Nz, const double& Lx, const double& Ly, const double& Lz)
    : _NbxElemtsGlob(Nx), _NbyElemtsGlob(Ny), _NbzElemtsGlob(Nz), _Lx(Lx), _Ly(Ly), _Lz(Lz) {
    MPI_Comm_size(MPI_COMM_WORLD, &_nbp);
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    _ny = _NbyElemtsGlob + 1;
    _nx = _NbxElemtsGlob + 1;
    _nz = _NbzElemtsGlob + 1;
  }

  MPI_Comm DefComm() {
    int ndims = 3;
    int periode[3]  = {false, false, false};
    int reorganisation = false;
    _coords = new int[ndims];
    _dims   = new int[ndims];
    _dims[0] = 0; _dims[1] = 0; _dims[2] = 0;

    MPI_Dims_create(_nbp, 3, _dims);
    if (_rank == 0) std::cout << "dims[x]:" << _dims[0] << "\t dims[y]:" << _dims[1] << "\t dims[z]:" << _dims[2] << std::endl;

    MPI_Cart_create(MPI_COMM_WORLD, ndims, _dims, periode, reorganisation, &comm_3D);
    MPI_Cart_coords(comm_3D, _rank, ndims, _coords);
    return comm_3D;
  }

  void Defneighbor() {
    int Left=0, Right=1, Bottom=2, Top=3, Front=4, Back=5;
    MPI_Cart_shift(comm_3D, 0, 1, &_neighbor[Left],   &_neighbor[Right]);
    MPI_Cart_shift(comm_3D, 1, 1, &_neighbor[Bottom], &_neighbor[Top]);
    MPI_Cart_shift(comm_3D, 2, 1, &_neighbor[Front],  &_neighbor[Back]);  
  }

  int* Neighbor() { return _neighbor; }

  void computeNbElemtsLoc() {
    _NbElemtsLoc = new int[_nbp];
    _NyLoc = new int[_nbp];
    _NxLoc = new int[_nbp];
    _NzLoc = new int[_nbp];              

    _NxLoc[_rank] = _nx / _dims[0];
    if (_coords[0] < _nx % _dims[0] && _nx % _dims[0] > 0) _NxLoc[_rank]++;

    _NyLoc[_rank] = _ny / _dims[1];
    if (_coords[1] < _ny % _dims[1] && _ny % _dims[1] > 0) _NyLoc[_rank]++;

    _NzLoc[_rank] = _nz / _dims[2];           
    if (_coords[2] < _nz % _dims[2] && _nz % _dims[2] > 0) _NzLoc[_rank]++;

    _NbElemtsLoc[_rank] = _NxLoc[_rank] * _NyLoc[_rank] * _NzLoc[_rank];

    MPI_Allgather(&_NbElemtsLoc[_rank], 1, MPI_INT, _NbElemtsLoc, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&_NxLoc[_rank], 1, MPI_INT, _NxLoc, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&_NyLoc[_rank], 1, MPI_INT, _NyLoc, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&_NzLoc[_rank], 1, MPI_INT, _NzLoc, 1, MPI_INT, MPI_COMM_WORLD);

  }
  
  int* NyLoc() { return _NyLoc; }
  int* NxLoc() { return _NxLoc; }
  int* NzLoc() { return _NzLoc; }

  void NbLocGhost() {
    _NbLocGhost = new int[_nbp];
    _NbLocGhost[_rank] = (_NxLoc[_rank]+2) * (_NyLoc[_rank]+2) * (_NzLoc[_rank]+2);
  }

  void computeDeplacement() {
    _deplacement    = new int[_nbp];
    _deplacement[0] = 0;
    for (int j = 1; j < _nbp; j++)   _deplacement[j] = _deplacement[j-1] + _NbElemtsLoc[j-1];
  }

  void setDebut() {
    _debuty = new int[_nbp];
    _debutx = new int[_nbp];
    _debutz = new int[_nbp];

    _debuty[_rank] = _coords[1] * (_ny / _dims[1]) + std::min(_coords[1], _ny % _dims[1]);
    _debutx[_rank] = _coords[0] * (_nx / _dims[0]) + std::min(_coords[0], _nx % _dims[0]);
    _debutz[_rank] = _coords[2] * (_nz / _dims[2]) + std::min(_coords[2], _nz % _dims[2]);  

    MPI_Allgather(&_debuty[_rank], 1, MPI_INT, _debuty, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&_debutx[_rank], 1, MPI_INT, _debutx, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&_debutz[_rank], 1, MPI_INT, _debutz, 1, MPI_INT, MPI_COMM_WORLD);
  }

  int* Debuty() { return _debuty; }
  int* Debutx() { return _debutx; }
  int* Debutz() { return _debutz; }

  double* AllocateMemory() {
    return new double[_NbLocGhost[_rank]];
  }

  double* DefMeshx() {
    _Mesh0x = new double[_nx];
    for (int i = 0; i <= _NbxElemtsGlob; i++) _Mesh0x[i] = i * _Lx / double(_NbxElemtsGlob);
    return _Mesh0x;
  }

  double* DefMeshy() {
    _Mesh0y = new double[_ny];
    for (int i = 0; i <= _NbyElemtsGlob; i++) _Mesh0y[i] = i * _Ly / double(_NbyElemtsGlob);
    return _Mesh0y;
  }

  double* DefMeshz() {
    _Mesh0z = new double[_nz];
    for (int i = 0; i <= _NbzElemtsGlob; i++) _Mesh0z[i] = i * _Lz / double(_NbzElemtsGlob);
    return _Mesh0z;
  }

  double* Initialize(Solution* Sol) {
    _Init = new double[_NbLocGhost[_rank]];
    const int nx = _NxLoc[_rank], ny = _NyLoc[_rank], nz = _NzLoc[_rank];
    const int sy = nx + 2;
    const int sz = (nx + 2) * (ny + 2);
    for (int k = 1; k <= nz; k++) {
      for (int j = 1; j <= ny; j++) {
        for (int i = 1; i <= nx; i++) {
          int I = i + j*sy + k*sz;
          _Init[I] = Sol->Initial(_Mesh0x[_debutx[_rank]+i-1], _Mesh0y[_debuty[_rank]+j-1],  _Mesh0z[_debutz[_rank]+k-1]);
        }
      }
    }
    return _Init;
  }

  double* AllocateCommMemory() {
    return new double[_NyLoc[_rank] * _NzLoc[_rank]];
  }

  double EulerExp(double* Utemp, double* U, const double& dtloc,const double& dx, const double& dy, const double& dz, const int& I) {
    const int sy = _NxLoc[_rank] + 2;
    const int sz = (_NxLoc[_rank] + 2) * (_NyLoc[_rank] + 2);
    Utemp[I] = U[I] + dtloc/dx/dx * (U[I+1]  - 2*U[I] + U[I-1]) + dtloc/dy/dy * (U[I+sy] - 2*U[I] + U[I-sy]) + dtloc/dz/dz * (U[I+sz] - 2*U[I] + U[I-sz]);
    return Utemp[I];
  }

  void Next(double* Utemp, double* U, const double& dtloc) {
    const double dx = _Lx / double(_NbxElemtsGlob);
    const double dy = _Ly / double(_NbyElemtsGlob);
    const double dz = _Lz / double(_NbzElemtsGlob);
    const int nx = _NxLoc[_rank], ny = _NyLoc[_rank], nz = _NzLoc[_rank];
    const int sy = nx + 2;
    const int sz = (nx + 2) * (ny + 2);
    for (int k = 2; k < nz; k++) {
      for (int j = 2; j < ny; j++) {
        for (int i = 2; i < nx; i++) {
          int I = i + j*sy + k*sz;
          EulerExp(Utemp, U, dtloc, dx, dy, dz, I);
        }
      }
    }
  }

  void BoundNext(double* Utemp, double* U, const double& dtloc) {
    const double dx = _Lx / double(_NbxElemtsGlob);
    const double dy = _Ly / double(_NbyElemtsGlob);
    const double dz = _Lz / double(_NbzElemtsGlob);
    const int nx = _NxLoc[_rank], ny = _NyLoc[_rank], nz = _NzLoc[_rank];
    const int sy = nx + 2;
    const int sz = (nx + 2) * (ny + 2);

    // Face j=1 (Bottom) : tous les k, tous les i
    for (int k = 1; k <= nz; k++) {
      for (int i = 1; i <= nx; i++) {
        int j = 1;
        int I = i + j*sy + k*sz;
        if (isPhysicalBoundary(i, j, k)) Utemp[I] = U[I];
        else EulerExp(Utemp, U, dtloc, dx, dy, dz, I);
      }
    }

    // Face j=NyLoc (Top) : tous les k, tous les i
    for (int k = 1; k <= nz; k++) {
      for (int i = 1; i <= nx; i++) {
        int j = ny;
        int I = i + j*sy + k*sz;
        if (isPhysicalBoundary(i, j, k)) Utemp[I] = U[I];
        else EulerExp(Utemp, U, dtloc, dx, dy, dz, I);
      }
    }

    // Face i=1 (Left) : tous les k, j intérieur
    for (int k = 1; k <= nz; k++) {
      for (int j = 2; j <= ny-1; j++) {
        int i = 1;
        int I = i + j*sy + k*sz;
        if (isPhysicalBoundary(i, j, k)) Utemp[I] = U[I];
        else EulerExp(Utemp, U, dtloc, dx, dy, dz, I);
      }
    }

    // Face i=NxLoc (Right) : tous les k, j intérieur
    for (int k = 1; k <= nz; k++) {
      for (int j = 2; j <= ny-1; j++) {
        int i = nx;
        int I = i + j*sy + k*sz;
        if (isPhysicalBoundary(i, j, k)) Utemp[I] = U[I];
        else EulerExp(Utemp, U, dtloc, dx, dy, dz, I);
      }
    }

    // Face k=1 (Front) : j intérieur, i intérieur (arêtes et coins déjà traités)
    for (int j = 2; j <= ny-1; j++) {
      for (int i = 2; i <= nx-1; i++) {
        int k = 1;
        int I = i + j*sy + k*sz;
        if (isPhysicalBoundary(i, j, k)) Utemp[I] = U[I];
        else EulerExp(Utemp, U, dtloc, dx, dy, dz, I);
      }
    }

    // Face k=NzLoc (Back) : j intérieur, i intérieur
    for (int j = 2; j <= ny-1; j++) {
      for (int i = 2; i <= nx-1; i++) {
        int k = nz;
        int I = i + j*sy + k*sz;
        if (isPhysicalBoundary(i, j, k)) Utemp[I] = U[I];
        else EulerExp(Utemp, U, dtloc, dx, dy, dz, I);
      }
    }
  }

  double* Copy(double* Utemp, double* U) {  
    const int nx = _NxLoc[_rank], ny = _NyLoc[_rank], nz = _NzLoc[_rank];
    const int sy = nx + 2;
    const int sz = (nx + 2) * (ny + 2);
    for (int k = 1; k <= nz; k++) {
      for (int j = 1; j <= ny; j++) {
        for (int i = 1; i <= nx; i++) {
          int I = i + j*sy + k*sz;
          U[I] = Utemp[I];
        }
      }
    }
    return U;
  }

  void Write(const int& Ind, const double* const x, Solution* Sol, double& T) {
    const int nx = _NxLoc[_rank], ny = _NyLoc[_rank], nz = _NzLoc[_rank];
    const int sy = nx + 2;
    const int sz = (nx + 2) * (ny + 2);

    // Extraction locale sans cellules fantômes
    double* TempLoc = new double[nx * ny * nz];
    for (int k = 1; k <= nz; k++) {
      for (int j = 1; j <= ny; j++) {
        for (int i = 1; i <= nx; i++) {
          int I    = i + j*sy + k*sz;
          int Iloc = (i-1) + (j-1)*nx + (k-1)*nx*ny;
          TempLoc[Iloc] = x[I];
        }
      }
    }

    double* Xtmp = nullptr;
    double* X    = nullptr;
    if (_rank == 0) {
      X    = new double[_nx * _ny * _nz];
      Xtmp = new double[_nx * _ny * _nz];
    }
    MPI_Gatherv(TempLoc, _NbElemtsLoc[_rank], MPI_DOUBLE,Xtmp, _NbElemtsLoc, _deplacement, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    delete[] TempLoc;

    if (_rank == 0) {
      for (int r = 0; r < _nbp; r++) {
        int rx = _NxLoc[r], ry = _NyLoc[r], rz = _NzLoc[r];
        for (int k = 0; k < rz; k++) {
          for (int j = 0; j < ry; j++) {
            for (int i = 0; i < rx; i++) {
              int I = _deplacement[r] + i + j*rx + k*rx*ry;
              int J = (_debutx[r]+i) + (_debuty[r]+j)*_nx + (_debutz[r]+k)*_nx*_ny;
              X[J] = Xtmp[I];
            }
          }
        }
      }

      std::stringstream mystream;
      mystream << "Data3D/Chal" << Ind << ".csv";

      std::ofstream File(mystream.str(), std::ofstream::out);
      File << "#Indice,#Point_x,#Point_y,#Point_z,#Solution_obtenue,#Solution_exacte\n";

      int m = 0;
      for (int k = 0; k <= _NbzElemtsGlob; k++) {
        for (int j = 0; j <= _NbyElemtsGlob; j++) {
          for (int i = 0; i <= _NbxElemtsGlob; i++) {
            File << m << "," << _Mesh0x[i] << "," << _Mesh0y[j] << "," << _Mesh0z[k] << "," << X[m] << "," << Sol->Exact(_Mesh0x[i], _Mesh0y[j], _Mesh0z[k], T) << "\n";
            m++;
          }
        }
      }

      File.close();
      delete[] X;
      delete[] Xtmp;
    }
  }
};
