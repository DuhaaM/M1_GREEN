#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cmath>

typedef double (* SolCondition)(const double& x, const int& K,const double& T);

double SolSinus(const double& x,const int& K,const double& T){

  if (T==0) return sin(K*M_PI*x);

  else return exp(-T)*sin(K*M_PI*x);
    

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
  

  double U[N+1],UTemp[N+1],x[N+1];
  
  dx=1./double(N);
  dt = 0.4 * dx * dx;
  
  for(i=0;i<=N;i++){
    x[i]= i*dx;
    U[ i ] = SolSinus(x[i],K,0); 
  }

  std::cout << "x = " ;
  for( i = 0; i <= N; i++)
    std::cout << x[i] << "\t";
  std::cout << std::endl;

  double t = 0.;
  double dtloc = 0.;
  int kT = 0;
  double ß;
  double Ind=0;
  std::cout << "dt = " << dt << std::endl;

  while( t < Tfinal ) {
    dtloc = dt;
    if ( t + dtloc > Tfinal ) dtloc = Tfinal - t ;
    U[0]=0;
    U[N]=0;


 
    for(i=1;i<N;i++){
      UTemp[i]=U[i]+ dtloc / ( dx * dx ) * ( U[i+1] - 2*U[i] + U[i-1] );
    }

    for(i=1;i<N;i++){
      U[i]=UTemp[i];	
    }


    if (kT % kprint==0){
      File_writting(Ind,U,x,N);
      Ind++;
    }

    t+=dtloc;
    kT++;
    
  }


  Error(N,U,x,Tfinal,K,SolSinus);
  Error_ev(N,U,x,K,Tfinal,dt);
  
   std::ofstream File( "ChalExact.dat", std::ofstream::out );
  for(i=0;i<=N;i++){      
    File <<x[i]<< "\t" <<  SolSinus( x[i],K,Tfinal )<<'\n';
  }   
  File.close();

  
  return 0;
}

