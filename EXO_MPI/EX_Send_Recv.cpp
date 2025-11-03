//Ce programme a pour but d'envoyer m  entiers appartenant à [0,m] d'un processus à l'autre


#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>


int main (int argc,char *argv[]) {
  int rang,nb_procs;
  int m=1000,etiquette=11;
  double A[m];

  srand(time(NULL));

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rang);
  MPI_Comm_size(MPI_COMM_WORLD,&nb_procs);

  if( nb_procs != 2){
    if ( rang == 0 ) {
      std::cout <<"Ce programme est éxécutable uniquement pour 2 processus"
		<<std::endl;
    }
    MPI_Barrier( MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD,EPERM);
  }

  if( rang == 0 ){
    for(int n=0;n < m;n++) {
      A[n]= rand() %( m+1);
    }
  }
  else {
    memset(A,0,m*sizeof(double));
  }
  if(rang == 0) {
    MPI_Send(A,m,MPI_DOUBLE,1,etiquette,MPI_COMM_WORLD);}
  else {
    MPI_Recv(A,m,MPI_DOUBLE,0,etiquette,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    //printf("Moi,processus 1,ai recu %f du processus 0\n",A);
    std::cout << "Moi,processus 1,ai recu du processus 0"
	      << std::endl;
    for(int j = 0; j < 1000 ; j++) {
      std::cout << A[j] << "\t";
    }
    std::cout << std::endl;
  }
  MPI_Finalize();
}
