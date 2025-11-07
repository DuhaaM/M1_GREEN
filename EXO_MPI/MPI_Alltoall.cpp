#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
int main(int argc,char *argv[]) {
  int nb_valeurs=8;
  int rang,nb_procs,longueur_tranche,i;
  float donnees[nb_valeurs],valeurs[nb_valeurs];
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nb_procs);
  MPI_Comm_rank(MPI_COMM_WORLD,&rang);

  for (i=0;i<nb_valeurs;i++) valeurs[i]=1001.+rang*nb_valeurs+i;
  longueur_tranche = nb_valeurs/nb_procs;
  printf("Moi, processus %d envoie mon tableau valeurs : ",rang);
  for (i=0; i<nb_valeurs;i++) printf("%f ",valeurs[i]);
  printf("\n");
  MPI_Alltoall(valeurs,longueur_tranche,MPI_FLOAT, donnees,longueur_tranche,MPI_FLOAT,MPI_COMM_WORLD);
 
  printf("Moi, processus %d, ai recu ",rang);
  for (i=0;i<nb_valeurs;i++) printf("%f ", donnees[i]);
  printf("\n");
  MPI_Finalize();
}
