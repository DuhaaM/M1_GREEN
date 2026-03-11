#include <mpi.h>
#include <iostream>




int main (int argc,char *argv[]){
  int rang,A[3],B[3];
  int etiquette=100;
  int message_rec,message_env;

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rang);


  for (int i=0;i<3;i++){
    if (rang==0) A[i]=i+1;
    if (rang==1) B[i]=i+4;
    
    if (rang==0) std::cout<<"A["<<i<<"]="<<A[i]<<std::endl;
    if (rang==1) std::cout<<"B["<<i<<"]="<<B[i]<<std::endl;

  }
  if(rang==0) message_env=A[0];
  if(rang==1) message_env=B[2];

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Sendrecv(&message_env,1,MPI_INT,(rang+1)%2,etiquette,&message_rec,1,MPI_INT,(rang+1)%2,etiquette,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  if(rang==0) A[0]=message_rec;
  if(rang==1) B[2]=message_rec;

  
  for (int i=0 ;i<3;i++){
        
    if (rang==0) std::cout<<"A["<<i<<"]="<<A[i]<<std::endl;
    if (rang==1) std::cout<<"B["<<i<<"]="<<B[i]<<std::endl;
  }

  MPI_Finalize();
  return 0;
}
