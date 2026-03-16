#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string>
#include <cmath>
#include <mpi.h>

class Solution{
private:
  int _k;
  double _L;
  bool _Flag;
  int _ID;
  

public:
  Solution(const int& k, const double& L) : _k(k), _L(L){}

  double L() const {return _L;}

  int K() const {return _k;}
  
  double SinInit(const double& x){
    _ID = 1;
    _Flag = true;
    return sin(_k*M_PI*x/_L);
  }

  double ConstInit(){
    _ID=2;
    _Flag = true; // La condition initiale ne variera pas au cours du programme 
    return 18;
  }

  double LinInit(const double& x){
    _ID=3;
    _Flag = true; // La condition initiale ne variera pas au cours du programme
    return -150*x+sqrt(2);
  }


  bool SolExist(){return _Flag;}
  
  double SolExact(const double& t,const double& x){

    if(_ID==1)  return exp(-t*_k*_k*M_PI*M_PI/_L/_L)*SinInit(x);

    if(_ID==2)  return ConstInit();

    if(_ID==3)  return LinInit(x);

    else{
      std::cerr<<"Indeterminate function"<<std::endl;
      return 0;
    }
  }

};





class Mesh{
private:
  int _N;
  int _Tfinal;
  int _nbp;
  int _rank;
  int _Q;
  int _R;
  int _debut;
  int _fin;
  int* _elem_recus;
  int*  _deplacement;
  int _Nproc;
  double* _mesh;
  double*  _U; 
  Solution* Sol;

public:
  Mesh(const int& N,const int& nb_procs, const int& rang,const double& L,const double& K,const double& Tfinal ) {
    _N=N;
    _nbp=nb_procs;
    _rank=rang;
    _Tfinal= Tfinal;
    Sol = new Solution(K,L); 
  }
  
  void computeQ(){ _Q = (_N+1) / _nbp;}

  void computeR(){ _R = (_N+1) % _nbp;}


  void computeElem_recus(){
    _elem_recus= new int[_nbp];
    for(int j=0;j<_nbp;j++){
      _elem_recus[j]=_Q;
      if (j < _R and _R > 0) _elem_recus[j]++;
    }
  }

  void computeDeplacement(){
    _deplacement= new int[_nbp];
    _deplacement[0]=0;
    for(int j=1;j<_nbp;j++){
      _deplacement[j]=_deplacement[j-1]+_elem_recus[j-1];
    }
  }

  void setDebut(){ _debut= _rank*_Q +std::min(_rank,_R);}
 
  void setFin(){ _fin = _debut+ _elem_recus[_rank]-1;}

  void computeNproc(){  _Nproc = _fin-_debut+1; }    // NB d'éléments par partition


  double* mesh(){
    _mesh=new double[_N+1];
    for ( int i=0; i<=_N;i++) _mesh[i]=i*(Sol->L())/double(_N);
    return _mesh;
  }

  
  void computeU(){
    int j =1 ;
    _U= new double[_Nproc+2];
      for (int i=_debut; i<=_fin;i++){
	_U[j]=Sol->SinInit(_mesh[i]) ;
	j++;
      }
  }


  void AsynchroneComm( double* X,double* x, const double& dtloc){
    MPI_Request request[4];
    int nb_req=0;
    int etiquette=100;

    double dx= Sol->L()/double(_N);
    
    if (_rank != _nbp-1) MPI_Irecv(&x[_Nproc+1],1,MPI_DOUBLE,_rank+1,etiquette,MPI_COMM_WORLD,&request[nb_req++]);
    if (_rank != 0)MPI_Irecv(&x[0],1,MPI_DOUBLE,_rank-1,etiquette,MPI_COMM_WORLD,&request[nb_req++]);

    if (_rank != _nbp-1) MPI_Isend(&x[_Nproc],1,MPI_DOUBLE,_rank+1,etiquette,MPI_COMM_WORLD,&request[nb_req++]);
    if (_rank != 0) MPI_Isend(&x[1],1,MPI_DOUBLE,_rank-1,etiquette,MPI_COMM_WORLD,&request[nb_req++]);


    for (int i = 2; i<=_Nproc-1; i++){
      X[i]= x[i] + dtloc / dx / dx * ( x[i+1]- 2 * x[i]+ x[i-1] );
    }

    
    MPI_Waitall(nb_req,request,MPI_STATUS_IGNORE);

    if (_rank!=0) X[1]=x[1] + dtloc /dx /dx * ( x[2] - 2*x[1] + x[0]);
    else X[1]=0;

    if (_rank!=_nbp-1)X[_Nproc]=x[_Nproc] + dtloc /dx /dx * ( x[_Nproc+1] - 2*x[_Nproc] + x[_Nproc-1]);
    else X[_Nproc]=0;

  }

  void SynchroneComm(double* X, double* x, const double& dtloc){
    int etiquette=100;
    double dx= Sol->L()/double(_N);

    if (_rank != _nbp-1) MPI_Sendrecv(&_U[_Nproc],1,MPI_DOUBLE,_rank+1,etiquette,&_U[_Nproc+1],1,MPI_DOUBLE,_rank+1,etiquette,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    if ( _rank != 0)  MPI_Sendrecv(&_U[1],1,MPI_DOUBLE,_rank-1,etiquette,&_U[0],1,MPI_DOUBLE,_rank-1,etiquette,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    for (int i = 1; i<=_Nproc; i++){
      X[i]= x[i] + dtloc / dx / dx * ( x[i+1]- 2 * x[i]+ x[i-1] );
    }

    if (_rank==0) X[1]=0;
    
    if (_rank==_nbp-1) X[_Nproc]=0;
  }


  double* Gather(const double* const x){
    double* X = new double[_N+1];
    MPI_Gatherv(&x[1],_Nproc,MPI_DOUBLE,X,_elem_recus,_deplacement,MPI_DOUBLE,0,MPI_COMM_WORLD);
    return X;
  }


  void Verif(const int& kT){
    double dx= Sol->L() /double(_N);
    double dt= 0.4*dx*dx;
    double Normesqrt_U0=0,vp,Err;

    for (int i=0; i<=_N;i++){
      Normesqrt_U0 += dx*Sol->SinInit(_mesh[i])*Sol->SinInit(_mesh[i]);
    }
    
    vp= 2*dt* (cos (Sol->K()*M_PI*dx) -1) / (dx*dx) +1;

    Err=abs( pow(vp,kT) - exp(-_Tfinal* Sol->K()*Sol->K() *M_PI*M_PI /Sol->L()/Sol->L() ) );

    std::cout<<"Erreur_vp="<<sqrt(Normesqrt_U0)*Err<<std::endl;

  }

  void ErrorInf(const double* Vect){
    double Err=0;
    for (int i=1 ; i<_N ; i++) Err=std::max(Err,std::abs(Vect[i]-Sol->SolExact(_Tfinal,_mesh[i])));
    std::cout<<"ErreurInf="<<Err<<std::endl;
  }


  void ErrorL2(const double* Vect){
    double Err=0;
    double dx = Sol->L()/double(_N);
    double x_1=1 /sqrt(3), x_2=-1/sqrt(3);

    for (int i=0; i < _N ; i++){
      
      double xG1= (_mesh[i+1]+_mesh[i]) / 2 + (_mesh[i+1]-_mesh[i]) / 2 * x_1;
      double xG2= (_mesh[i+1]+_mesh[i]) / 2 + (_mesh[i+1]-_mesh[i]) / 2 * x_2;

      Err+= dx / 2 * (-Sol->SolExact(_Tfinal,xG1) + (Vect[i+1]+Vect[i])/2 + (Vect[i+1]-Vect[i])/2 *x_1) * (-Sol->SolExact(_Tfinal,xG1) + (Vect[i+1]+Vect[i])/2 + (Vect[i+1]-Vect[i])/2 *x_1);

      Err+= dx / 2 * (-Sol->SolExact(_Tfinal,xG2) + (Vect[i+1]+Vect[i])/2 + (Vect[i+1]-Vect[i])/2 *x_2) * (-Sol->SolExact(_Tfinal,xG2) + (Vect[i+1]+Vect[i])/2 + (Vect[i+1]-Vect[i])/2 *x_2);

    }
    std::cout<<"Erreur="<<sqrt(Err)<<std::endl;

  }
 

  void File_Writting(const int& Ind,const double* const X){
    std::stringstream mystream;
    mystream<<"Data/Chal"<<Ind<<".dat";
    
    std::ofstream File( mystream.str(), std::ofstream::out );
    for(int i=0;i<=_N;i++){      
      File << X[i]<<'\n';
    }   
    File.close();

}
    
    

  const int debut() const { return _debut; }
  const int fin()   const { return _fin; }
  const int Nproc() const { return _Nproc ; }
  const int* elem_recus() const {return _elem_recus ;}
  const int* deplacement() const { return _deplacement ;}
  double* U(){return _U ;}
};




int main (int argc, char *argv[]){
  int i,j,N,K,kprint,method;
  double Tfinal,L;
  double temps_debut,temps_fin;

  // Définiton manuelle des variables 
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

  Mesh Mesh(N,nb_procs,rang,L,K,Tfinal);
  //Défnition du maillage

  Mesh.computeQ();
  Mesh.computeR();
  Mesh.computeElem_recus();
  Mesh.computeDeplacement(); 
  Mesh.setDebut();
  Mesh.setFin();
  Mesh.computeNproc();
 

  double* x = Mesh.mesh();

  Mesh.computeU();

  double dx = L/double(N);
  double dtloc = 0.4* dx *dx;
  double t=0,time=0;
  int kT=0,Ind=0;

  double* U = Mesh.U();
  double* Utemp= new double[Mesh.Nproc()+2];
  
  while( t < Tfinal ) {

    if ( t + dtloc > Tfinal) dtloc = Tfinal - t;

    if (rang==0)temps_debut=MPI_Wtime();

    
    if ( method == 1) Mesh.AsynchroneComm(Utemp,U,dtloc);

    if ( method == 2) Mesh.SynchroneComm(Utemp,U,dtloc);

    if(rang==0){
      temps_fin=MPI_Wtime();
      time+=temps_fin-temps_debut;
    }



    for(int i = 1; i<=Mesh.Nproc(); i++){
      U[i]=Utemp[i];
      }

     if (kT % kprint == 0) {
       Mesh.File_Writting(Ind,Mesh.Gather(U));
       Ind++;
     }
     
    t+=dtloc;
    kT++;
  }


  double* Ufinal = Mesh.Gather(U);

  
  
  if (rang == 0 ){
    std::cout<<"Communication time  = "<<time<<"s"<<std::endl;
    Mesh.Verif(kT);
    Mesh.ErrorInf(Ufinal);
    Mesh.ErrorL2(Ufinal);
  }
  

  
  MPI_Finalize();
  return 0;

 
}
  


  

  
