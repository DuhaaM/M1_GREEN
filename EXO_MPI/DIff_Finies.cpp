#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cmath>
#include <mpi.h>
#include <vector>






typedef double (* SolCondition)(const double& x, const int& K,const double& T,const double& L);
double SolSinus(const double& x,const int& K,const double& T,const double& L){

   return exp(-T*K*K*M_PI*M_PI/L/L)*sin(K*M_PI*x/L);
    

}
class Cut{
private :


  int Q_;
  int r_;
  int debut_;
  int fin_;
  std::vector<int> elem_recus_;
  std::vector<int> deplacement_;
  int nb_procs;
  int rang;
  int Nproc_;
  std::vector<double> U_;
  

  
public :
  Cut(const int& N,const int& nbp,const int& rank) : nb_procs(nbp), rang(rank) {
    Q_ = (N + 1) / nbp;
    r_ = (N + 1) % nbp;
    elem_recus_.resize(nbp);
    deplacement_.resize(nbp);
  }

  void computeElemrecus(){
     for(int j=0;j<nb_procs;j++){
	elem_recus_[j]=Q_;
	if (j < r_ and r_ > 0) elem_recus_[j]++;
      }
  }

  void  computeDeplacement(){
    deplacement_[0]=0;
    for(int j=1;j<nb_procs;j++){
      deplacement_[j]=deplacement_[j-1]+elem_recus_[j-1];
    }
  }

    

  void setDebut(const int& rang){ debut_= rang*Q_ +std::min(rang,r_);}
 
  void setFin(const int& rang){ fin_ = debut_+ elem_recus_[rang]-1;}

  void computeNproc(const int& rang){
    Nproc_ = fin_-debut_+1;
  }

  
  void computeU(SolCondition Solution,const double& k, const double& L,const std::vector<double>&  mesh){
    U_.reserve(Nproc_);
      
     for (int i=debut_;i<=fin_;i++){
       U_.push_back( Solution(mesh[i],k,0,L));
     }
  }
    
  int Q() const {return Q_;}
  int r() const {return r_;}
  const std::vector<int>& elem_recus() const {return elem_recus_;}
  const std::vector<int>& deplacement() const { return deplacement_;}
  int debut() const { return debut_; }
  int fin()   const { return fin_; }
  int Nproc() const { return Nproc_; }
  const std::vector<double>& U() const {return U_;}
};





std::vector<double> Mesh(const int& N,const double& L){

  std::vector<double> mesh;
  double dx=L/double(N);
  for (int i=0;i<=N;i++) mesh.push_back(i*dx);
  return mesh;


}


  
  

void File_writting(const int& k,const double* const U, const std::vector<double>& x,const int& N){
  std::stringstream mystream;
  mystream<<"Chal"<<k<<".dat";
  
  std::ofstream File( mystream.str(), std::ofstream::out );
  for(int i=0;i<=N;i++){      
    File <<x[i]<< "\t" << U[i]<<'\n';
  }   
  File.close();

}





void Error(const int& N,const std::vector<double>& U,const std::vector<double>& Mesh,const double& Tfinal,const int& K,const double& L,SolCondition Solution){
  double ERR=0;
  double xG1,xG2;
  double step = Mesh[1]-Mesh[0];
  for(int i = 0; i < N ; i++) {  

    xG1=step*(-1/sqrt(3)+1)/2+Mesh[i];
    xG2=step*(+1/sqrt(3)+1)/2+Mesh[i];
  
    ERR+=1/(4*double(N)*double(N)) * ( Solution( xG1,K,Tfinal,L ) - (U[i+1]-U[i])*(xG1-Mesh[i])/step + U[i]) * ( Solution( xG1,K,Tfinal,L ) -(U[i+1]-U[i])*(xG1-Mesh[i])/step + U[i]);
    ERR+=1/(4*double(N)*double(N)) * ( Solution( xG2,K,Tfinal,L ) - (U[i+1]-U[i])*(xG2-Mesh[i])/step + U[i]) * ( Solution( xG2,K,Tfinal,L ) -(U[i+1]-U[i])*(xG2-Mesh[i])/step + U[i]);
  }
  std::cout<<"Erreur="<<sqrt(ERR)<<std::endl;

}





void Error_ev(const int& N,const std::vector<double>& U,const std::vector<double>& Mesh,const int& K,const double& Tfinal,const double& L, SolCondition Solution){
  double dt=0.4 * L * L / (N * N);
  double step=Mesh[1]-Mesh[0];
  double Normesqr_U0=0,vp,Err;
  for (int i=0; i<=N ;i++){
    Normesqr_U0+=step*Solution(Mesh[i],K,0,L)*Solution(Mesh[i],K,0,L);
  }
  vp=2*dt/(step*step)*(cos(K*M_PI*step)-1);
  Err=abs(pow(vp,N) -exp(-Tfinal*K*K*M_PI*M_PI/L/L));
  std::cout<<"Erreur_vp="<<sqrt(Normesqr_U0)*Err<<std::endl;
  
}











std::vector<double> EulExp_Sendrecv(SolCondition Solution, Cut cut, const double& Tfinal,const int& k,const std::vector<double>& mesh,const int& rang,const int& nb_procs,const int& N,const double& L,const double& dtloc,const std::vector<double>& U ){




  double dx=L/double(N);
  
  int Nproc=cut.Nproc();
  std::vector<double>  UTemp(Nproc);


  std::vector<int> deplacement(nb_procs), elem_recus(nb_procs);
  std::vector<double> Ufinal;
  double Send1,Send2,Recv1,Recv2;

        
    for(int i=1;i<Nproc-1;i++){
      UTemp[i]=U[i]+ dtloc / ( dx * dx ) * ( U[i+1] - 2*U[i] + U[i-1] );
    }

    int etiquette=100;
    
    Send1=U[Nproc-1];
    Send2=U[0];

    
    if (rang != nb_procs-1){
	
      MPI_Sendrecv(&Send1,1,MPI_DOUBLE,rang+1,etiquette,&Recv1,1,MPI_DOUBLE,rang+1,etiquette,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      UTemp[Nproc-1]=U[Nproc-1]+ dtloc / ( dx * dx ) * ( Recv1 - 2*U[Nproc-1] + U[Nproc-2] );
    }
    else UTemp[Nproc-1] = 0;
    
    if ( rang != 0){
      
      MPI_Sendrecv(&Send2,1,MPI_DOUBLE,rang-1,etiquette,&Recv2,1,MPI_DOUBLE,rang-1,etiquette,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      UTemp[0]=U[0] +dtloc / (dx*dx) * ( U[1] -2*U[0] + Recv2);
      
    }
    else UTemp[0]=0;

  return UTemp;
}














std::vector<double> EulExp_NoBlock(SolCondition Solution, Cut cut, const double& Tfinal,const int& k,const std::vector<double>& mesh,const int& rang,const int& nb_procs,const int& N,const double& L,const double& dtloc,const std::vector<double>& U ){




  double dx=L/double(N);
    
  int Nproc=cut.Nproc();
  std::vector<double> UTemp(Nproc);


  double Send1,Send2,Recv1,Recv2;
  MPI_Request request[4];
  int nb_req;


  int etiquette=100;
  
  Send1=U[Nproc-1];
  Send2=U[0];
  
  
  nb_req = 0;
  
  if (rang != nb_procs-1) MPI_Irecv(&Recv2,1,MPI_DOUBLE,rang+1,etiquette,MPI_COMM_WORLD,&request[nb_req++]);
  if (rang != 0)MPI_Irecv(&Recv1,1,MPI_DOUBLE,rang-1,etiquette,MPI_COMM_WORLD,&request[nb_req++]);
  
  Send1=U[Nproc-1];
  Send2=U[0];
  
        
  for(int i=1;i<Nproc-1;i++){
    UTemp[i]=U[i]+ dtloc / ( dx * dx ) * ( U[i+1] - 2*U[i] + U[i-1] );
  }
  
  if (rang != nb_procs-1) MPI_Isend(&Send1,1,MPI_DOUBLE,rang+1,etiquette,MPI_COMM_WORLD,&request[nb_req++]);
  if (rang != 0) MPI_Isend(&Send2,1,MPI_DOUBLE,rang-1,etiquette,MPI_COMM_WORLD,&request[nb_req++]);
  


  MPI_Waitall(nb_req,request,MPI_STATUS_IGNORE);
  
  if (rang != 0) UTemp[0]=U[0] +dtloc / (dx*dx) * ( U[1] -2*U[0] + Recv1);
  else UTemp[0]=0;
  
  
  if (rang != nb_procs-1)UTemp[Nproc-1]=U[Nproc-1]+ dtloc / ( dx * dx ) * (Recv2 - 2*U[Nproc-1] + U[Nproc-2] );
  else UTemp[Nproc-1]=0;



  return UTemp;
}










int main (int argc, char *argv[]){
  int i,j,N,K,kprint,method;
  double Tfinal,L;
  double temps_debut,temps_fin;

  if (argc != 7) {
    std::cerr << "usage is : " << argv[0] << " <N=mesh size> <t=final time> <K=integer> <kprint=print step> <L=Upperbound> <Comm_Method : Non-Blocking = 1 or Sendrecv = 2 >"
	      << std::endl;
    return 0;
  }
  else {
  
    N=std::atoi(argv[1]); 
    Tfinal=std::atof(argv[2]);
    K=std::atoi(argv[3]);
    kprint=std::atoi(argv[4]);
    L=std::atof(argv[5]);
    method = std::atoi(argv[6]);
 }

  int nb_procs,rang;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nb_procs);
  MPI_Comm_rank(MPI_COMM_WORLD,&rang);

  temps_debut=MPI_Wtime();
  

  std::vector<double> x=Mesh(N,L);

  Cut C1(N,nb_procs,rang);
  C1.computeElemrecus();
  C1.computeDeplacement();
  C1.setDebut(rang);
  C1.setFin(rang);
  C1.computeNproc(rang);
  C1.computeU(SolSinus,K,L,x);
  
  std::vector<double> Ufinal,U;

  double t=0;
  double dtloc;
  int kT=0,Ind=0;

  U=C1.U();
      
  while( t < Tfinal ) {
    dtloc = 0.4*L*L/(double(N)*double(N));
    if ( t + dtloc > Tfinal ) dtloc = Tfinal - t ;


    if(method == 1 ) U = EulExp_NoBlock(SolSinus,C1,Tfinal,K,x,rang,nb_procs,N,L,dtloc,U);
    
    
    if(method == 2) U = EulExp_Sendrecv(SolSinus,C1,Tfinal,K,x,rang,nb_procs,N,L,dtloc,U);

    if(kT % kprint == 0 or t==Tfinal){
      if (rang == 0) Ufinal.resize(N+1);
      MPI_Gatherv(U.data(),C1.Nproc(),MPI_DOUBLE,Ufinal.data(),C1.elem_recus().data(),C1.deplacement().data(),MPI_DOUBLE,0,MPI_COMM_WORLD);
    }

    if (rang == 0 and kT % kprint==0){
        File_writting(Ind,Ufinal.data(),x,N);
        Ind++;
    }

    t+=dtloc;
    kT++;
    
  } 
  if (rang == 0){
    Error(N,Ufinal,x,Tfinal,K,L,SolSinus);
    Error_ev(N,Ufinal,x,K,Tfinal,L,SolSinus);
    

    std::ofstream File( "ChalExact.dat", std::ofstream::out );
    for(i=0;i<=N;i++){      
      File <<x[i]<< "\t" <<  SolSinus( x[i],K,Tfinal,L )<<'\n';
    }   
    File.close();
  }

  if (rang ==0) {
  temps_fin = MPI_Wtime();
  std::cout<<"Temps="<<temps_fin-temps_debut<<std::endl;
  }
  
  MPI_Finalize();
  return 0;
}

