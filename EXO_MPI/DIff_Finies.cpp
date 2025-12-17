#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cmath>
#include <mpi.h>

typedef double (* SolCondition)(const double& x, const int& K,const double& T);

double SolSinus(const double& x,const int& K,const double& T){

  if (T==0) return sin(K*M_PI*x);

  else return exp(-T*M_PI*M_PI)*sin(K*M_PI*x);
    

}



void File_writting(const int& k,const double* const U, const double* const x,const int& N){
  std::stringstream mystream;
  mystream<<"Chal"<<k<<".dat";
  
  std::ofstream File( mystream.str(), std::ofstream::out );
  for(int i=0;i<=N;i++){      
    File <<x[i]<< "\t" << U[i]<<'\n';
  }   
  File.close();

}

void Error(const int& N,const double* const U,const double* const Mesh,const int& Tfinal,const int& K,SolCondition Solution){
  double ERR=0;
  double xG1,xG2;
  double step = Mesh[1]-Mesh[0];
  for(int i = 0; i < N ; i++) {  

    xG1=step*(-1/sqrt(3)+1)/2+Mesh[i];
    xG2=step*(+1/sqrt(3)+1)/2+Mesh[i];
  
    ERR+=1/(4*double(N)*double(N)) * ( Solution( xG1,K,Tfinal ) - (U[i+1]-U[i])*(xG1-Mesh[i])/step + U[i]) * ( Solution( xG1,K,Tfinal ) -(U[i+1]-U[i])*(xG1-Mesh[i])/step + U[i]);
    ERR+=1/(4*double(N)*double(N)) * ( Solution( xG2,K,Tfinal ) - (U[i+1]-U[i])*(xG2-Mesh[i])/step + U[i]) * ( Solution( xG2,K,Tfinal ) -(U[i+1]-U[i])*(xG2-Mesh[i])/step + U[i]);
  }
  std::cout<<"Erreur="<<step<<"\t"<<sqrt(ERR)<<std::endl;

}
			       

void Error_ev(const int& N,const double* const U,const double* const Mesh,const int& K,const int& Tfinal,const int& dt){
  double step=Mesh[1]-Mesh[0];
  double Normesqr_U0=0,vp;
  for (int i=0; i<=N ;i++){
    Normesqr_U0+=step*U[i]*U[i];
  }
  vp=2/(step*step)*(cos(K*M_PI*step)-1);
  std::cout<<"Erreur_vp="<<step<<"\t"<<pow( 1+ dt*vp -exp(-Tfinal),N)*sqrt(Normesqr_U0)<<std::endl;
  
}


int main (int argc, char *argv[]){
  int i,j,N,K,kprint;
  double dt=0,Tfinal,dx,l;
  double Ufinal[N+1];

  if (argc != 5) {
    std::cerr << "usage is : " << argv[0] << " <N=mesh size> <t=final time> <K=integer> <kprint=print step>"
	      << std::endl;
    return 0;
  }
  else {
  
    N=std::atoi(argv[1]); 
    Tfinal=std::atof(argv[2]);
    K=std::atoi(argv[3]);
    kprint=std::atoi(argv[4]);
 }
  
  int nb_procs,rang;
  int etiquette=100;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nb_procs);
  MPI_Comm_rank(MPI_COMM_WORLD,&rang);
  
  
  dx=1./double(N);
  dt = 0.4 * dx * dx;

  double x[N+1];
  for (i=0;i<=N;i++) x[i]=i*dx;
  
  int Q,r,debut,fin;
  
  Q = (N+1) / nb_procs;
  r = (N+1) % nb_procs;
  
  if (r != 0){
    if (rang < r){
      debut = rang * (Q+1);
      fin = (rang+1)*Q+1;
     
    }
    else{
      debut = rang * Q + r;
      if (debut + Q >= N) fin = N+1;
      if (debut+Q < N) fin = r+ Q * (rang+1) ;
    }
  }
  else {
    debut=rang*(Q+1);
    fin=(rang+1)*Q;
  }




  int elem_recus[nb_procs],deplacement[nb_procs];
  
  if (rang==0){
    elem_recus[0]=Q;
    deplacement[0]=0;
    if(r > 0) elem_recus[0]++;
    for(j=1;j<nb_procs;j++){
      deplacement[j]=deplacement[j-1]+elem_recus[j-1];
      elem_recus[j]=Q;
      if (j < r) elem_recus[j]++;
    }
  }

  int Nproc=fin-debut;
  double U[Nproc],UTemp[Nproc];
  int s=0;

  for (i=debut;i<=fin;i++){
    U[s] = SolSinus(x[i],K,0);
    s++;
  }



  for( i = 0; i <Nproc; i++){
    std::cout << U[i] << "\t"<<rang<<std::endl;
  }
  std::cout << std::endl;

  double t = 0.;
  double dtloc = 0.;
  int kT = 0;
  double tmp;
  double Ind=0;
  double Send1,Send2,Recv1,Recv2;
  std::cout << "dt = " << dt << std::endl;

    
  while( t < Tfinal ) {
    dtloc = dt;
    if ( t + dtloc > Tfinal ) dtloc = Tfinal - t ;
    if (rang==0) UTemp[0]=0;
    if (rang == nb_procs-1)UTemp[Nproc-1]=0;


 
    for(i=1;i<Nproc-1;i++){
      UTemp[i]=U[i]+ dtloc / ( dx * dx ) * ( U[i+1] - 2*U[i] + U[i-1] );
    }
    Send1=U[Nproc-1];
    Send2=U[0];
      
    if (rang != nb_procs-1){
	
      MPI_Sendrecv(&Send1,1,MPI_DOUBLE,rang+1,etiquette,&Recv1,1,MPI_DOUBLE,rang+1,etiquette,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      UTemp[Nproc-1]=U[Nproc-1]+ dtloc / ( dx * dx ) * ( Recv1 - 2*U[Nproc-1] + U[Nproc-2] );
    }

    if ( rang != 0){
      
      MPI_Sendrecv(&Send2,1,MPI_DOUBLE,rang-1,etiquette,&Recv2,1,MPI_DOUBLE,rang-1,etiquette,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      UTemp[0]=U[0] +dtloc / (dx*dx) * ( U[1] -2*U[0] + Recv2);

    }
	  
    
    for(i=0;i<Nproc;i++){
      U[i]=UTemp[i];	
    }

    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gatherv(U,Nproc,MPI_DOUBLE,Ufinal,elem_recus,deplacement,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if (kT % kprint==0){
      File_writting(Ind,Ufinal,x,N);
      Ind++;
    }
    
    t+=dtloc;
    kT++;

  }


  if (rang == 0){
    Error(N,Ufinal,x,Tfinal,K,SolSinus);
    Error_ev(N,Ufinal,x,K,Tfinal,dt);
  
    for (i=0;i<=N;i++){
      std::cout<<Ufinal[i]<<std::endl;
    }
    std::ofstream File( "ChalExact.dat", std::ofstream::out );
    for(i=0;i<=N;i++){      
      File <<x[i]<< "\t" <<  SolSinus( x[i],K,Tfinal )<<'\n';
    }   
    File.close();
  }
  
  MPI_Finalize();
  return 0;
}

