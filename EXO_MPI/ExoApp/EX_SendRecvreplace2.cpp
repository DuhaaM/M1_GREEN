//Ce programme a pour but d'envoyer m  entiers appartenant à [0,m] d'un processus à l'autre


#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>

int main (int argc,char *argv[]) {
  int rang,nb_procs;
  int etiquette=11,nb_ech=9;

  srand(time(NULL));
  std::fstream File0("Rec0.txt");
  std::fstream File1("Rec1.txt");


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
  for (int j=0;j<nb_ech;j++){
    int m=rand() % (1001);
    int *p=&m;
    double A[*p];
    if ( rang == 0 ) {
    std::cout<<"Ce programme transmet un vecteur de taille"<<m<<"\n";
    }
    for(int n=0;n < m;n++) {
      A[n]= rand() %( 1001);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double temps_debut=MPI_Wtime();
    if (rang==0) {
      MPI_Sendrecv_replace(A,m,MPI_DOUBLE,1,etiquette,1,etiquette,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }

    if (rang==1){
      MPI_Sendrecv_replace(A,m,MPI_DOUBLE,0,etiquette,0,etiquette,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    double temps_fin=MPI_Wtime();
    MPI_Barrier( MPI_COMM_WORLD );
    if ( rang == 0) {
       for(int j = 0; j < 1000 ; j++) {
	 File0<< A[j] << "\t";
       }
    }
    MPI_Barrier( MPI_COMM_WORLD);
    if ( rang == 1) {
      for(int j = 0; j < 1000 ; j++) {
	File0<< A[j] << "\t";
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if( rang == 0 ) {
      printf("Le programme s'est réalise en %f secondes. \n",temps_fin-temps_debut);
    }
  }
  MPI_Finalize();
}
