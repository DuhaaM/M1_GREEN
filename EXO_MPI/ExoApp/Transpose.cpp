#include <mpi.h>
#include <iostream>
int main( int argc,char *argv[]){
  int rang,i,j,taille;
  int nb_lignes=4,nb_colonnes=5,etiquette = 100;
  int a[nb_lignes][nb_colonnes];
  int a_transpose[nb_colonnes][nb_lignes];
  MPI_Datatype type_colonne,type_transpose;
  MPI_Aint pas;
  MPI_Status statut;

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rang);
  


  if (rang==0) {
    std::cout<<"Matrice A :\n"<<std::endl;
    for ( i=0;i<nb_lignes;i++){
      for ( j=0;j<nb_colonnes;j++){
	a[i][j]=1+i*nb_colonnes+j;
	std::cout<<a[i][j]<<std::endl;
      }
    }
  }
  MPI_Type_vector(nb_colonnes,1,nb_lignes,MPI_INT,&type_colonne);
  MPI_Type_size(MPI_INT,&taille);
  pas=taille;
  MPI_Type_create_hvector(nb_lignes,1,pas,type_colonne,&type_transpose);
  MPI_Type_commit(&type_transpose);
  if(rang==0){
    MPI_Send(a,nb_colonnes*nb_lignes,MPI_INT,1,etiquette,MPI_COMM_WORLD);
  }
  if(rang==1){
    MPI_Recv(a_transpose,1,type_transpose,0,etiquette,MPI_COMM_WORLD,&statut);
    std::cout<<"Matrice A transposée:\n"<<std::endl;
    for ( i=0;i<nb_colonnes;i++){
      for ( j=0;j<nb_lignes;j++){
	std::cout<<a[i][j]<<std::endl;
      }
    }
  }
  MPI_Type_free(&type_colonne);
  MPI_Type_free(&type_transpose);
  MPI_Finalize();
}

	
