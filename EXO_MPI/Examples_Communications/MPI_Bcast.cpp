//MPI_Bcast()git 

#include <mpi.h>
#include <stdio.h>

int main (int argc,char *argv[]){
  int rang,valeur;

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rang);

  if (rang == 2) valeur = rang+1000;

  MPI_Bcast(&valeur,1,MPI_INT,2,MPI_COMM_WORLD); //Le 2 processus diffuse sa valeur à tout le monde 

  printf("Moi,processus %d,ai recu %d du processus 2\n",rang,valeur);
  MPI_Finalize();
}
