#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>



int main (int argc, char *argv[]){
  int i,j,n,r,nb_procs,rang;
  int debut,fin;
  int N=100,etiquette=100;
  double dt=0.001,t=0,L=1,dx,T;
  double U[N+1],Uold[N+1],x[N+1];
  double temps_debut,temps_fin,temps;
 
 

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nb_procs);
  MPI_Comm_rank(MPI_COMM_WORLD,&rang);

  MPI_Request request[nb_procs];
  temps_debut=MPI_Wtime();
  
  dx=2*L/N;
  T=1/dt;

  n=N-2/nb_procs;
  r=N-2-nb_procs*n;
  
  if (rang<=r) {
    fin=(rang+1)*n+1;
    debut=rang*n;
  }
  else if (rang=0){
    debut=1;
      }
  else {
    fin=(rang+1)*n;
    debut=rang*n;
  }
  std::fstream File("Chal.dat");
  
  for(i=0;i<=N;i++){
    x[i]= -L+2*L*i*dx;
    U[i]=1;
    Uold[i]=1;
  }

  
  for(j=0;j<T;j++){
    t+=dt;
    U[0]=0;
    U[N]=0;
    Uold[0]=0;
    Uold[N]=0;

    
	
    for(i=debut;i<fin;i++){
      if (rang==0){
	MPI_Isend(&U[fin-1],1,MPI_DOUBLE,rang+1,etiquette,MPI_COMM_WORLD,request[rang]);
	MPI_Irecv(&U[fin],1,MPI_DOUBLE,rang+1,etiquette,MPI_COMM_WORLD,request[rang+1]);
      }
      else if (rang==nb_procs){
	MPI_Isend(&U[debut],1,MPI_DOUBLE,rang-1,etiquette,MPI_COMM_WORLD,request[rang]);
	MPI_Irecv(&U[debut-1],1,MPI_DOUBLE,rang-1,etiquette,MPI_COMM_WORLD,request[rang-1]);
      }
      else {
	MPI_Isend(&U[fin-1],1,MPI_DOUBLE,rang+1,etiquette,MPI_COMM_WORLD,request[rang]);
	MPI_Irecv(&U[fin],1,MPI_DOUBLE,rang+1,etiquette,MPI_COMM_WORLD,request[rang+1]);
	MPI_Isend(&U[debut],1,MPI_DOUBLE,rang-1,etiquette,MPI_COMM_WORLD,request[rang]);
	MPI_Irecv(&U[debut-1],1,MPI_DOUBLE,rang-1,etiquette,MPI_COMM_WORLD,request[rang-1]);
      }
	
      Uold[i]=U[i]+(dt/dx/dx)*(U[i+1]-2*U[i]+U[i-1]);
    }
    
    for(i=0;i<=N;i++){
      U[i]=Uold[i];
    }

    std::cout<<t<<std::endl;
    
    for(i=0;i<=N;i++){
      File<< x[i]<<"\t"<< U[i] <<"\n";
    }
  }
  temps_fin=MPI_Wtime();
  temps=temps_fin-temps_debut;
  std::cout<<temps<<"secondes"<<std::endl;

  
  return 0;
}
