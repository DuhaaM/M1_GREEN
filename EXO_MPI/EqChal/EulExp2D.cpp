#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cmath>
#include <vector>


typedef double (* Sol2DCondition)(const double& x,const double& y, const int& K,const double& T,const double& Lx,const double& Ly);
double SolSinus(const double& x,const double& y,const int& K,const double& T,const double& Lx,const double& Ly){

  return exp(-T*K*K*M_PI*M_PI/Lx/Lx)*exp(-T*K*K*M_PI*M_PI/Ly/Ly)*sin(K*M_PI*x/Lx)*sin(K*M_PI*y/Ly);
}



//Création de maillage 1D


std::vector<double> Mesh(const int& N,const double& L){

  std::vector<double> mesh;
  double dx=L/double(N);
  for (int i=0;i<=N;i++) mesh.push_back(i*dx);
  return mesh;
}



void File_writting(const int& k,const std::vector<double>& U, const std::vector<double>& x, const std::vector<double>& y,const int& N){
  std::stringstream mystream;
  mystream<<"Data2D/Chal"<<k<<".dat";
  
  std::ofstream File( mystream.str(), std::ofstream::out );
  for(int i=0;i<=N;i++){
    File <<U[i]<<'\n';
     }   
  File.close();

}



int main (int argc, char *argv[]){
  int Nx,Ny,K,kprint;
  double Tfinal,Lx,Ly;

  
  if (argc != 8){
     std::cerr << "usage is : " << argv[0] << " <Nx=mesh size> <Ny=mesh size> <t=final time> <K=integer> <kprint=print step> <Lx=Upperbound> <Ly=Upperbound> "
	      << std::endl;
     return 0;
    }

  else {
    Nx=std::atoi(argv[1]);
    Ny=std::atoi(argv[2]);
    Tfinal=std::atof(argv[3]);
    K=std::atoi(argv[4]);
    kprint=std::atoi(argv[5]);
    Lx=std::atof(argv[6]);
    Ly=std::atof(argv[7]);
  }
   
  //Définition des maillages
  std::vector<double> x=Mesh(Nx,Lx);
  std::vector<double> y=Mesh(Ny,Ly);

  int N=(Nx+1)*(Ny+1);
  std::vector<double> Ufinal,U,UTemp;
  double t=0;
  int kT=0,Ind=0;

  U.reserve(N);
  for (int i=0; i <= Nx; i++){
    for (int j=0; j <= Ny; j++){
      U.push_back(SolSinus(x[i],y[j],K,0,Lx,Ly));
    }
  }
  
  for(int i=0;i<N;i++){
    std::cout<<"U["<<i<<"]="<<U[i]<<std::endl;
    }

  
  double dy = Ly/double(Ny);
  double dx= Lx/double(Nx); // Pas d'espace
  
  double dtloc = 0.4* dx* dx; // Pas de temps

  UTemp.resize(N);
  
  while( t < Tfinal) {

    if ( t + dtloc > Tfinal ) dtloc = Tfinal - t ; // Assurance de ne pas dépasser Tfinal
    
    for (int i=1; i < Nx; i++){
      for (int j=1; j < Ny; j++){
	int I = i+j *(Nx+1);

	UTemp[I]= U[I]+ dtloc/dx/dx * ( U[I+1] - 2*U[I] + U[I-1] );
	UTemp[I]+= dtloc/dy/dy * ( U[I+Nx+1]- 2*U[I]+ U[I-Nx-1]);
      }
    }
    
    for(int i=1; i<N ;i++) U[i]=UTemp[i];	
    
    t+=dtloc;
    kT++;


    if( kT % kprint == 0){
      File_writting(Ind,U,x,y,N);
      Ind++;
     }

    if (t==Tfinal){
      std::stringstream mystream;
      mystream<<"Data2D/ChalFinal.dat";
      
      std::ofstream File( mystream.str(), std::ofstream::out );
      for(int i=0;i<N;i++){
	File <<U[i]<<'\n';
      }   
      File.close();
    }

  }


  return 0;
}
