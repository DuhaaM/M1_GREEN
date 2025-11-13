#include <mpi.h>
#include <iostream>
int main( int argc,char *argv[]){
  int rang,i,j,x;
  int n=4,m=5,etiquette = 100;
  int a[n][m];
  int longueurs_blocs[n],deplacement[n];
  MPI_Datatype type_colonne,type_transpose;

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rang);


  if (rang==0) {
    for ( i=0;i<n;i++){
      for ( j=0;j<m;j++){
	a[i][j]=1+i*m+j;
      }
      longueurs_blocs[i]=1;
      deplacement[i]=0;
    }
  }
  MPI_Type_vector(n,1,m,MPI_INT,&type_colonne);
  MPI_Type_commit(&type_colonne);
  MPI_Type_indexed(4,longueurs_blocs,deplacement,type_colonne,&type_transpose);
  MPI_Type_commit(&type_transpose);
  if(rang==0){
	MPI_Send(a,1,type_transpose,1,etiquette,MPI_COMM_WORLD);
  }
  if(rang==1){
    for ( i=0;i<n;i++){
      for ( j=0;j<m;j++){
	std::cout<<a[i][j]<<std::endl;
      }
    }
  }
  MPI_Type_free(&type_colonne);
  MPI_Type_free(&type_transpose);
  MPI_Finalize();
}

	
