#include <stdio.h>
#include <math.h>
#include <mpi.h>


int main(int argc, char *argv[]) {
  long long nbbloc;
  double  x,I;
  int rang,nb_procs;

 
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nb_procs);
  MPI_Comm_rank(MPI_COMM_WORLD,&rang);

  double temps_debut = MPI_Wtime();

  /* Nombre d'intervalles */
  nbbloc = 3*1000*1000LL*100;
  /* largeur des intervalles */
  double largeur = 1.0/nbbloc;

  double somme=0;

  int Q=nbbloc/nb_procs;
  int r=nbbloc % nb_procs;
  
  long long debut = rang * Q + std::min(rang,r);
  long long fin = debut+ Q ;
  
  if ( rang < r and r > 0) fin++;

  for (int i=debut; i<fin; i++) {
    /* Point au milieu de l'intervalle */
    x = largeur*(i+0.5);
    /* Calcul de l'aire */
    somme = somme + largeur*(4.0 / (1.0 + x*x));
  }
 
  MPI_Reduce(&somme,&I,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  double temps_fin  = MPI_Wtime();
  if(rang==0){
    printf("Pi = %.12lf\n", I);
    printf("Ecart = %g\n",I-4.0*atan(1.0));
    std::cout<<"Temps="<<temps_fin-temps_debut<<std::endl;
  }
  MPI_Finalize();
}






