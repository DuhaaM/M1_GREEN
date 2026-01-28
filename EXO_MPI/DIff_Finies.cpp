#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cmath>
#include <mpi.h>
#include <vector>






typedef double (* SolCondition)(const double& x, const int& K,const double& T);
double SolSinus(const double& x,const int& K,const double& T){

   return exp(-T*K*M_PI*M_PI)*sin(K*M_PI*x);
    

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

  
public :
  Cut(const int& N,const int& nbp) : nb_procs(nbp) {
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
    
  int Q() const {return Q_;}
  int r() const {return r_;}
  const std::vector<int>& elem_recus() const {return elem_recus_;}
  const std::vector<int>& deplacement() const { return deplacement_;}
  int debut() const { return debut_; }
  int fin()   const { return fin_; }
  int Nproc() const { return Nproc_; }
};





std::vector<double> Mesh(const int& N,const int& Tfinal){

  std::vector<double> mesh;
  double dx=Tfinal/double(N);
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





void Error(const int& N,const std::vector<double>& U,const std::vector<double>& Mesh,const double& Tfinal,const int& K,SolCondition Solution){
  double ERR=0;
  double xG1,xG2;
  double step = Mesh[1]-Mesh[0];
  for(int i = 0; i < N ; i++) {  

    xG1=step*(-1/sqrt(3)+1)/2+Mesh[i];
    xG2=step*(+1/sqrt(3)+1)/2+Mesh[i];
  
    ERR+=1/(4*double(N)*double(N)) * ( Solution( xG1,K,Tfinal ) - (U[i+1]-U[i])*(xG1-Mesh[i])/step + U[i]) * ( Solution( xG1,K,Tfinal ) -(U[i+1]-U[i])*(xG1-Mesh[i])/step + U[i]);
    ERR+=1/(4*double(N)*double(N)) * ( Solution( xG2,K,Tfinal ) - (U[i+1]-U[i])*(xG2-Mesh[i])/step + U[i]) * ( Solution( xG2,K,Tfinal ) -(U[i+1]-U[i])*(xG2-Mesh[i])/step + U[i]);
  }
  std::cout<<"Erreur="<<sqrt(ERR)<<std::endl;

}





void Error_ev(const int& N,const std::vector<double>& U,const std::vector<double>& Mesh,const int& K,const double& Tfinal,const int& dt, SolCondition Solution){
  double step=Mesh[1]-Mesh[0];
  double Normesqr_U0=0,vp,Err;
  for (int i=0; i<=N ;i++){
    Normesqr_U0+=step*Solution(Mesh[i],K,0)*Solution(Mesh[i],K,0);
  }
  vp=2*dt/(step*step)*(cos(K*M_PI*step)-1);
  Err=abs(pow(vp,N) -exp(-Tfinal*K*M_PI*M_PI));
  std::cout<<"Erreur_vp="<<sqrt(Normesqr_U0)*Err<<std::endl;
  
}











std::vector<double> EulExp_Sendrecv(SolCondition Solution, Cut cut, const double& Tfinal,const int& k,const std::vector<double>& mesh,const int& rang,const int& nb_procs,const int& N ){




  double dx=Tfinal/double(N);
  double dt = 0.4 * dx * dx;
  
  int Nproc=cut.Nproc();
  std::vector<double> U, UTemp(Nproc);

  for (int i=cut.debut();i<=cut.fin();i++){
    U.push_back( Solution(mesh[i],k,0));
  }


  std::vector<int> deplacement(nb_procs), elem_recus(nb_procs);
  std::vector<double> Ufinal;
  double t = 0.;
  double dtloc;
  int kT = 0;
  double Ind=0;
  double Send1,Send2,Recv1,Recv2;

    
  while( t < Tfinal ) {
    dtloc = dt;
    if ( t + dtloc > Tfinal ) dtloc = Tfinal - t ;

        
    for(int i=1;i<Nproc-1;i++){
      UTemp[i]=U[i]+ k*dtloc / ( dx * dx ) * ( U[i+1] - 2*U[i] + U[i-1] );
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
	  
          
    for(int i=0;i<Nproc;i++){
      U[i]=UTemp[i];	
    }

   
    if (rang == 0) Ufinal.resize(N+1);

    MPI_Gatherv(U.data(),Nproc,MPI_DOUBLE,Ufinal.data(),cut.elem_recus().data(),cut.deplacement().data(),MPI_DOUBLE,0,MPI_COMM_WORLD);



    t+=dtloc;
    kT++;
  }

  return Ufinal;
}














std::vector<double> EulExp_NoBlock(SolCondition Solution, Cut cut, const double& Tfinal,const int& k,const std::vector<double>& mesh,const int& rang,const int& nb_procs,const int& N ){




  double dx=Tfinal/double(N);
  double dt = 0.4 * dx * dx;
    
  int Nproc=cut.Nproc();
  std::vector<double> U, UTemp(Nproc);

  for (int i=cut.debut();i<=cut.fin();i++){
    U.push_back( Solution(mesh[i],k,0));
  }


  std::vector<int> deplacement(nb_procs), elem_recus(nb_procs);
  std::vector<double> Ufinal;
  double t = 0.;
  double dtloc;
  int kT = 0;
  double Ind=0;
  double Send1,Send2,Recv1,Recv2;
  MPI_Request request[4];
  int nb_req;

    
  while( t < Tfinal ) {
    dtloc = dt;
    if ( t + dtloc > Tfinal ) dtloc = Tfinal - t ;


    int etiquette=100;
    
    Send1=U[Nproc-1];
    Send2=U[0];

    
    nb_req = 0;

    if (rang != nb_procs-1) MPI_Irecv(&Recv2,1,MPI_DOUBLE,rang+1,etiquette,MPI_COMM_WORLD,&request[nb_req++]);
    if (rang != 0)MPI_Irecv(&Recv1,1,MPI_DOUBLE,rang-1,etiquette,MPI_COMM_WORLD,&request[nb_req++]);

    Send1=U[Nproc-1];
    Send2=U[0];

        
    for(int i=1;i<Nproc-1;i++){
      UTemp[i]=U[i]+ k*dtloc / ( dx * dx ) * ( U[i+1] - 2*U[i] + U[i-1] );
    }
    
    if (rang != nb_procs-1) MPI_Isend(&Send1,1,MPI_DOUBLE,rang+1,etiquette,MPI_COMM_WORLD,&request[nb_req++]);
    if (rang != 0) MPI_Isend(&Send2,1,MPI_DOUBLE,rang-1,etiquette,MPI_COMM_WORLD,&request[nb_req++]);



    MPI_Waitall(nb_req,request,MPI_STATUS_IGNORE);
    
    if (rang != 0) UTemp[0]=U[0] +dtloc / (dx*dx) * ( U[1] -2*U[0] + Recv1);
    else UTemp[0]=0;

        
    if (rang != nb_procs-1)UTemp[Nproc-1]=U[Nproc-1]+ dtloc / ( dx * dx ) * (Recv2 - 2*U[Nproc-1] + U[Nproc-2] );
    else UTemp[Nproc-1]=0;

            
    for(int i=0;i<Nproc;i++){
      U[i]=UTemp[i];	
    }

    if (rang == 0) Ufinal.resize(N+1);

    MPI_Gatherv(U.data(),Nproc,MPI_DOUBLE,Ufinal.data(),cut.elem_recus().data(),cut.deplacement().data(),MPI_DOUBLE,0,MPI_COMM_WORLD);


    t+=dtloc;
    kT++;
  }

  return Ufinal;
}










int main (int argc, char *argv[]){
  int i,j,N,K,kprint,method;
  double dt=0,Tfinal,dx,l;
  double temps_debut,temps_fin;

  if (argc != 6) {
    std::cerr << "usage is : " << argv[0] << " <N=mesh size> <t=final time> <K=integer> <kprint=print step> <Comm_Method : Non-Blocking = 1 or Sendrecv = 2 >"
	      << std::endl;
    return 0;
  }
  else {
  
    N=std::atoi(argv[1]); 
    Tfinal=std::atof(argv[2]);
    K=std::atoi(argv[3]);
    kprint=std::atoi(argv[4]);
    method=std::atoi(argv[5]);
 }

  int nb_procs,rang;
  int etiquette=100;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nb_procs);
  MPI_Comm_rank(MPI_COMM_WORLD,&rang);

  temps_debut=MPI_Wtime();
  

  std::vector<double> x=Mesh(N,Tfinal);

  Cut C1(N,nb_procs);
  C1.computeElemrecus();
  C1.computeDeplacement();
  C1.setDebut(rang);
  C1.setFin(rang);
  C1.computeNproc(rang);
  
  std::vector<double> Ufinal;

  if(method == 1 ) Ufinal = EulExp_NoBlock(SolSinus,C1,Tfinal,K,x,rang,nb_procs,N);

  if(method == 2) Ufinal = EulExp_Sendrecv(SolSinus,C1,Tfinal,K,x,rang,nb_procs,N);
  
  if (rang == 0){
    Error(N,Ufinal,x,Tfinal,K,SolSinus);
    Error_ev(N,Ufinal,x,K,Tfinal,dt,SolSinus);

    std::ofstream File( "ChalExact.dat", std::ofstream::out );
    for(i=0;i<=N;i++){      
      File <<x[i]<< "\t" <<  SolSinus( x[i],K,Tfinal )<<'\n';
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

