#include <mpi.h>
#include <stdio.h>

int main(int argc,char *argv[]) {
  int rang,nb_procs,valeur,produit,i;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nb_procs);
  MPI_Comm_rank(MPI_COMM_WORLD,&rang);

  if (rang == 0)
    valeur = 10;
  else
    valeur = rang;

  MPI_Allreduce(&valeur,&produit,1,MPI_INT,MPI_PROD,MPI_COMM_WORLD);

  printf("Moi, processus %d, ai recu la valeur du produit globale %d\n",rang,produit);

  MPI_Finalize();
}
