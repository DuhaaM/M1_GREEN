#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc,char *argv[]) {
  int nb_valeurs=8,rang,nb_procs,longueur_tranche,i;
  float *valeurs, *donnees;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nb_procs);
  MPI_Comm_rank(MPI_COMM_WORLD,&rang);

  longueur_tranche = nb_valeurs/nb_procs;
  donnees = (float *) malloc(longueur_tranche*sizeof(float));

  if ( rang == 2 ) {
    valeurs = (float *) malloc(nb_valeurs*sizeof(float));

    for ( i=0 ; i < nb_valeurs ; i++) valeurs[i]=1001. +i;
    printf("Moi,processus %d envoie mon tableau valeurs:",rang);

    for ( i=0 ; i < nb_valeurs ; i++) {printf("%f ",valeurs[i]);} printf("\n");}

  MPI_Scatter(valeurs,longueur_tranche,MPI_FLOAT,donnees,longueur_tranche,MPI_FLOAT,2,MPI_COMM_WORLD);

  printf("Moi,processs %d, ai recu ",rang);

  for( i=0 ; i<longueur_tranche; i++) printf("%f\t",donnees[i]);

  printf("du processus 2\n");
  MPI_Finalize();
}
