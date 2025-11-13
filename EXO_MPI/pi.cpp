#include <stdio.h>
#include <math.h>
#include <mpi.h>


int main(int argc, char *argv[]) {
  long long nbbloc,i,debut,fin;
  double n,largeur, somme, x,I;
  int ntmp,r,rang,nb_procs;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nb_procs);
  MPI_Comm_rank(MPI_COMM_WORLD,&rang);

  /* Nombre d'intervalles */
  nbbloc = 3*1000*1000LL*100;
  /* largeur des intervalles */
  largeur = 1.0/nbbloc;

  somme=0;

  ntmp=nbbloc/nb_procs;
  r=nbbloc - nb_procs*ntmp;
  n=ntmp;

  if (rang<=r) {
    fin=(rang+1)*n+1;
    debut=rang*n;
  }
  else {
    fin=nbbloc;
    debut=rang*n;
  }

  for (i=debut; i<fin; i++) {
    /* Point au milieu de l'intervalle */
    x = largeur*(i+0.5);
    /* Calcul de l'aire */
    somme = somme + largeur*(4.0 / (1.0 + x*x));
  }
 
  MPI_Reduce(&somme,&I,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(rang==0){
    printf("Pi = %.12lf\n", I);
    printf("Ecart = %g\n",somme-4.0*atan(1.0));
  }
  i=fin-debut;
  MPI_Finalize();
}
