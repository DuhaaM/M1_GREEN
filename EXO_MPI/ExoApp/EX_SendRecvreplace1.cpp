//Ce programme a pour but d'envoyer m  entiers appartenant à [0,m] d'un processus à l'autre


#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>

int main (int argc,char *argv[]) {
  int rang,nb_procs;
  int m=1000,etiquette=11;
  double A[m];

  srand(time(NULL));
  std::fstream File0("Rec0.txt");
  std::fstream File1("Rec1.txt");

  MPI_Init(&argc,&argv);
  double temps_debut=MPI_Wtime();
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
  for(int n=0;n < m;n++) {
    A[n]= rand() %( m+1);
  }
    
  int peer=(rang+1)%2;
  MPI_Sendrecv_replace(A,m,MPI_DOUBLE,peer,etiquette,peer,etiquette,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  if (rang ==0) {
    for(int j = 0; j < 1000 ; j++) {
      File0<< A[j] << "\t";
    }
  }
  else {
    for(int j = 0; j < 1000 ; j++) {
      File1 << A[j] << "\t";
    }
  }
  double temps_fin=MPI_Wtime();
    printf("Le programme s'est réalise en %f secondes. \n",temps_fin-temps_debut);
  MPI_Finalize();
}
