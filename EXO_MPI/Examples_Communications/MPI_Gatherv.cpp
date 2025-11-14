#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc,char *argv[]) {
int nb_valeurs=10,rang,nb_procs,longueur_tranche,reste,i;
float donnees[nb_valeurs];
float *valeurs;
int *nb_elements_recus,*deplacements;

 MPI_Init(&argc,&argv);
 MPI_Comm_size(MPI_COMM_WORLD,&nb_procs);
 MPI_Comm_rank(MPI_COMM_WORLD,&rang);
 longueur_tranche = nb_valeurs/nb_procs;
 reste = nb_valeurs%nb_procs;
 if (rang < reste) longueur_tranche = longueur_tranche+1;
 valeurs = (float *) malloc(longueur_tranche*sizeof(float));
 for (i=0;i<longueur_tranche;i++)
   valeurs[i]=1001.+rang*(nb_valeurs/nb_procs)+(rang<reste?rang:reste)+i;
 printf("Moi, processus %d envoie mon tableau valeurs : ",rang);
 for (i=0; i<longueur_tranche;i++) {printf("%f ",valeurs[i]); }printf("\n");
 if (rang == 2) {
   nb_elements_recus = (int *) malloc(nb_procs*sizeof(int));
   deplacements = (int *) malloc(nb_procs*sizeof(int));
   nb_elements_recus[0] = nb_valeurs/nb_procs;
   if (reste > 0) nb_elements_recus[0] = nb_elements_recus[0]+1;
   deplacements[0] = 0;
   for (i=1;i<nb_procs;i++) {
     deplacements[i] = deplacements[i-1]+nb_elements_recus[i-1];
     nb_elements_recus[i] = nb_valeurs/nb_procs;
     if (i < reste) nb_elements_recus[i] = nb_elements_recus[i]+1;
   } }
 MPI_Gatherv(valeurs,longueur_tranche,MPI_FLOAT,
donnees,nb_elements_recus,deplacements,MPI_FLOAT,2,MPI_COMM_WORLD);
 if (rang==2) {
   printf("Moi, processus %d, ai recu ",rang);
   for (i=0;i<nb_valeurs;i++) printf("%f ", donnees[i]);
   printf("\n"); }
 MPI_Finalize();
}
