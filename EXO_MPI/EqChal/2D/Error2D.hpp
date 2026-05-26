#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string>
#include <cmath>
#include <mpi.h>

class Error{  //Classe calculant les erreurs de notre méthode
private:
  double _L;             // Longueur du maillage 
  int _NbxElemtsGlob;    // Nombre d'éléments du maillageœ
  int _NbyElemtsGlob;    // Nombre d'éléments du maillage 
  int _k;                // Diffusivité thermique 
  double _Tfinal;        // Temps final 
  double* _Mesh0x;       // Maillage définiœ
  double* _Mesh0y;       // Maillage défini
  int _fin;              // Dernier indice de chaque partition
  int _rank;             // Indice de chaque proc
  int* _NyLoc;
  int* _NxLoc;
  int* _debutx;
  int* _debuty;

public:
  Error(const double& L,const int& Nx,const int& Ny, const int& k,const double& Tfinal,double* Initx ,double* Inity , int* debutx , int* debuty,int* NyLoc,int* NxLoc):_L(L), _NbxElemtsGlob(Nx), _NbyElemtsGlob(Ny), _k(k), _Tfinal(Tfinal), _Mesh0x(Initx),_Mesh0y(Inity),_debutx(debutx),_debuty(debuty),_NyLoc(NyLoc),_NxLoc(NxLoc){
    MPI_Comm_rank(MPI_COMM_WORLD,&_rank);   //Constructeur de la classe
  } 



  void ErrorInf(const double* Vect, Solution* Sol){  // Calcul de l'erreur en norme infinie
    double Errloc=0,Errglob=0; 
    for (int k =1; k<=_NyLoc[_rank];k++){
      for (int i=1; i<=_NxLoc[_rank];i++){
	int I = i+k*(_NxLoc[_rank]+2);
	Errloc=std::max(Errloc,std::abs(Vect[I]-Sol->Exact(_Mesh0x[_debutx[_rank]+i-1],_Mesh0y[_debuty[_rank]+k-1],_Tfinal)));
      }
    }
    MPI_Reduce(&Errloc,&Errglob,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD); //  Défnition du max entre tous les max de chaque proc 
    if(_rank==0) std::cout<<"ErreurInf="<<Errglob<<std::endl;
  }

  void ErrorL2(const double* Vect, Solution* Sol) {

    const double xi[2]  = { -1.0/sqrt(3.0), 1.0/sqrt(3.0) };
    const double w[2]   = { 1.0, 1.0 };

    const double dx = (_Mesh0x[1] - _Mesh0x[0]);
    const double dy = (_Mesh0y[1] - _Mesh0y[0]);

    double ErrlocL2 = 0.0, ErrglobL2 = 0.0;

    for (int k = 1; k < _NyLoc[_rank]; k++) {
        for (int i = 1; i < _NxLoc[_rank]; i++) {


            double x0 = _Mesh0x[_debutx[_rank] + i- 1];
            double x1 = _Mesh0x[_debutx[_rank] + i];
            double y0 = _Mesh0y[_debuty[_rank] + k-1];
            double y1 = _Mesh0y[_debuty[_rank] + k];

	    double cell_sum =0 ;
            for (int p = 0; p < 2; p++) {
                for (int q = 0; q < 2; q++) {
                    double xg = 0.5*(x0 + x1) + 0.5*dx*xi[p];
                    double yg = 0.5*(y0 + y1) + 0.5*dy*xi[q];

                    double s = (xg - x0) / dx;
                    double r = (yg - y0) / dy;

                    double uh = (1-s)*(1-r)*Vect[i+ k*(_NxLoc[_rank]+2)]+s*(1-r)*Vect[(i+1)+k*(_NxLoc[_rank]+2)]+(1-s)*r*Vect[i+(k+1)*(_NxLoc[_rank]+2)]+s*r*Vect[(i+1)+(k+1)*(_NxLoc[_rank]+2)];

                    double uex = Sol->Exact(xg, yg, _Tfinal);
		    
                    double diff = uh - uex;
                    cell_sum+= w[p] * w[q] * diff * diff;
                }
            }
            ErrlocL2 +=cell_sum* dx * dy / 4.0;
        }
    }

    MPI_Reduce(&ErrlocL2, &ErrglobL2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (_rank == 0)
        std::cout << "ErreurL2 = " << sqrt(ErrglobL2) << std::endl;
}
};
 
