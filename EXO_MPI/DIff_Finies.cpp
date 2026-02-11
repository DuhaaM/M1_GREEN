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

  //Premier et dernier indice de chaque partition
  
  void computeNproc(const int& rang){
    Nproc_ = fin_-debut_+1;    // NB d'éléments par partition
  }

  
  void computeU(SolCondition Solution,const double& k, const double& L,const std::vector<double>&  mesh){
    U_.reserve(Nproc_);
      
     for (int i=debut_;i<=fin_;i++){
       U_.push_back( Solution(mesh[i],k,0,L));  //Définition des partitions
     }
  }

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
  double step=L/double(N); 
  for(int i = 0; i < N ; i++) {  

    xG1=step*(-1/sqrt(3)+1) / 2 +Mesh[i];
    xG2=step*(+1/sqrt(3)+1) / 2 +Mesh[i];
  
    ERR+=L*L/(4*double(N)*double(N)) * ( -Solution( xG1,K,Tfinal,L ) + (U[i+1]-U[i])*(xG1-Mesh[i])/2 + U[i]) * ( -Solution( xG1,K,Tfinal,L ) +(U[i+1]-U[i])*(xG1-Mesh[i])/2 + U[i]);
    ERR+=L*L/(4*double(N)*double(N)) * ( -Solution( xG2,K,Tfinal,L ) + (U[i+1]-U[i])*(xG2-Mesh[i])/2 + U[i]) * ( -Solution( xG2,K,Tfinal,L ) +(U[i+1]-U[i])*(xG2-Mesh[i])/2 + U[i]);
  }
  std::cout<<"Erreur="<<sqrt(ERR)<<std::endl;

}

void Error2(const std::vector<double>&U,SolCondition Solution, const int& N,const int& K, const int& L, const double& Tfinal,std::vector<double>& x){
  double ERR=0;
  for(int i=0; i<=N;i++){
    ERR+=(U[i]-Solution(x[i],K,Tfinal,L))*(U[i]-Solution(x[i],K,Tfinal,L));
  }

  std::cout<<"Erreur="<<sqrt(ERR)<<std::endl;
}





void Error_ev(const int& N,const std::vector<double>& U,const std::vector<double>& Mesh,const int& K,const double& Tfinal,const double& L,const int& kT, SolCondition Solution){
  double dt=0.4 * L * L / (double(N) *double(N) );
  double dx=L/double(N);
  double Normesqr_U0=0,vp,Err;
  for (int i=0; i<=N ;i++){
    Normesqr_U0+=dx*Solution(Mesh[i],K,0,L)*Solution(Mesh[i],K,0,L);
  }
  vp=2*dt*(cos(K*M_PI*dx)-1)/(dx*dx) + 1;
  Err=abs(pow(vp,kT)-exp(-Tfinal*K*K*M_PI*M_PI/L/L));
  std::cout<<"Erreur_vp="<<sqrt(Normesqr_U0)*Err<<std::endl;
  
}












std::vector<double> EulExp_Sendrecv(Cut& cut,const int& rang,const int& nb_procs,const int& N,const double& L,const double& dtloc,const std::vector<double>& U ,std::vector<double>& UTemp,const int& Nproc,const double& dx){



  double Send1,Send2,Recv1,Recv2; //Stockage des communications
        
  for(int i=1;i<Nproc-1;i++){
    UTemp[i]=U[i]+ dtloc / ( dx * dx ) * ( U[i+1] - 2*U[i] + U[i-1] );
  } //Calcul au temps t+1
  
  int etiquette=100;
  
  Send1=U[Nproc-1];
  Send2=U[0];

  
  if (rang != nb_procs-1){
    
    MPI_Sendrecv(&Send1,1,MPI_DOUBLE,rang+1,etiquette,&Recv1,1,MPI_DOUBLE,rang+1,etiquette,MPI_COMM_WORLD,MPI_STATUS_IGNORE);  // Envoi du dernier élément et réception  du premier de la partition supèrieur
    UTemp[Nproc-1]=U[Nproc-1]+ dtloc / ( dx * dx ) * ( Recv1 - 2*U[Nproc-1] + U[Nproc-2] );  
  }
  else UTemp[Nproc-1] = U[Nproc-1];
  
  if ( rang != 0){
    
    MPI_Sendrecv(&Send2,1,MPI_DOUBLE,rang-1,etiquette,&Recv2,1,MPI_DOUBLE,rang-1,etiquette,MPI_COMM_WORLD,MPI_STATUS_IGNORE); // Inversion du premier échange
    UTemp[0]=U[0] +dtloc / (dx*dx) * ( U[1] -2*U[0] + Recv2); 
    
  }
  else UTemp[0]=U[0];
      
  return UTemp;
}








std::vector<double> EulExp_NoBlock(Cut& cut,const int& rang,const int& nb_procs,const int& N,const double& L,const double& dtloc,const std::vector<double>& U, std::vector<double>& UTemp,const int& Nproc, const double& dx ){

    
  double Send1,Send2,Recv1,Recv2;
  MPI_Request request[4];
  int nb_req;


  int etiquette=100;
  
  Send1=U[Nproc-1];
  Send2=U[0];
  
  nb_req = 0;
  
  if (rang != nb_procs-1) MPI_Irecv(&Recv2,1,MPI_DOUBLE,rang+1,etiquette,MPI_COMM_WORLD,&request[nb_req++]);
  if (rang != 0)MPI_Irecv(&Recv1,1,MPI_DOUBLE,rang-1,etiquette,MPI_COMM_WORLD,&request[nb_req++]);  //Création des zones de réception

  
  if (rang != nb_procs-1) MPI_Isend(&Send1,1,MPI_DOUBLE,rang+1,etiquette,MPI_COMM_WORLD,&request[nb_req++]);
  if (rang != 0) MPI_Isend(&Send2,1,MPI_DOUBLE,rang-1,etiquette,MPI_COMM_WORLD,&request[nb_req++]);


          
  for(int i=1;i<Nproc-1;i++){
    UTemp[i]=U[i]+ dtloc / dx/dx * ( U[i+1] - 2*U[i] + U[i-1] );
  }
  
  MPI_Waitall(nb_req,request,MPI_STATUS_IGNORE); //Attente de réalisation


  if (rang != 0) UTemp[0]=U[0] + dtloc / dx / dx * ( U[1] -2*U[0] + Recv1);
  else UTemp[0]=U[0];
  
  
  if (rang != nb_procs-1) UTemp[Nproc-1]=U[Nproc-1]+ dtloc /dx/ dx  * (Recv2 - 2*U[Nproc-1] + U[Nproc-2] );
  else UTemp[Nproc-1]=U[Nproc-1];


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


  std::vector<double> x=Mesh(N,L);

  Cut C1(N,nb_procs,rang);
  C1.computeElemrecus();
  C1.computeDeplacement();
  C1.setDebut(rang);
  C1.setFin(rang);
  C1.computeNproc(rang);
  C1.computeU(SolSinus,K,L,x);


  std::vector<double> Ufinal,U,UTemp;

  double t=0,time=0;
  int kT=0,Ind=0;

  U=C1.U();

  int Nproc=C1.Nproc();
  UTemp.resize(Nproc);

  
  double dx=L/double(N); // Pas d'espace
  double dtloc = 0.4* dx* dx; 

  while( t < Tfinal) {

    if ( t + dtloc > Tfinal ) dtloc = Tfinal - t ;

    if (rang==0) temps_debut=MPI_Wtime();
    if(method == 1 ) U = EulExp_NoBlock(C1,rang,nb_procs,N,L,dtloc,U,UTemp,Nproc,dx);
    if (rang == 0) temps_fin = MPI_Wtime();

    if (rang == 0) temps_fin = MPI_Wtime();

    if(method == 2) U = EulExp_Sendrecv(C1,rang,nb_procs,N,L,dtloc,U,UTemp,Nproc,dx);
    

    time+=temps_fin-temps_debut; 
    /* if(t+dtloc==Tfinal){

      if(rang == 0){
	File_writting(Ind,Ufinal.data(),x,N);
        Ind++;
	}	
    }*/

    t+=dtloc;
    kT++;
  }
  if (rang ==0) {
  std::cout<<"Temps="<<time<<std::endl;
  }
  
  if (rang == 0){
    Ufinal.resize(N+1);
  }
  MPI_Gatherv(U.data(),C1.Nproc(),MPI_DOUBLE,Ufinal.data(),C1.elem_recus().data(),C1.deplacement().data(),MPI_DOUBLE,0,MPI_COMM_WORLD);


  if (rang == 0){
    for (int i= 0; i<=N ; i++){
    std::cout<<"ucalc"<<i<<"\t"<<Ufinal[i]<<"\t"<<"uprevu"<<i<<"\t"<<SolSinus(x[i],K,0,L)*pow(0.8*(cos(K*M_PI*dx)-1)+1,kT)<<std::endl;
    }
    Error(N,Ufinal,x,Tfinal,K,L,SolSinus);
    //Error2(Ufinal,SolSinus,N,K,L,Tfinal,x);
    Error_ev(N,Ufinal,x,K,Tfinal,L,kT,SolSinus);
    /*std::ofstream File( "ChalExact.dat", std::ofstream::out );
    for(i=0;i<=N;i++){      
      File <<x[i]<< "\t" <<  SolSinus( x[i],K,Tfinal,L )<<'\n';
    }   
    File.close();*/
  }


  
  MPI_Finalize();
  return 0;
}

