#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>

double g(double U,double t,double dt, double dx,double x){
  double G;

  G=((U-exp(-2*t)*sin(M_PI*(1/2+x/2)))*(U-exp(t*dt/dx/dx)*sin(M_PI*(1/2+x/2))))/2;

  return G;
}
    


int main (int argc, char *argv[]){
  int i,j,N;
  double dt=0,Tfinal,dx,l;

  if (argc !=3 ) {
    std::cerr << "usage is : " << argv[0] << " <N=mesh sise> <t=final time>"
	      << std::endl;
    return 0;
  }
  else {
    N=std::atoi(argv[1]); 
    Tfinal=std::atof(argv[2]); 
 }
  

  double U[N+1],Uold[N+1],x[N+1];
  
  dx=1./double(N);
  dt = 1./(dx * dx);
  //  T=1/dt;
  
  for(i=0;i<=N;i++){
    x[i]= i*dx;
    U[ i ] = sin( M_PI * x[i] ); 
    Uold[i]=1;
  }

  double t = 0.;
  double dtloc = 0.;
  while( t < Tfinal ) {
    dtloc = dt;
    if ( t + dtloc > Tfinal ) dtloc = Tfinal - t ;
    U[0]=0;
    U[N]=0;
    Uold[0]=0;
    Uold[N]=0;    

    
	
    for(i=1;i<N;i++){
	
      Uold[i]=U[i]+(dt/dx/dx)*(U[i+1]-2*U[i]+U[i-1]);
    }
    for(i=0;i<=N;i++){
      U[i]=Uold[i];	
    }
    t+=dt;
    std::stringstream mystream;
    mystream<<"Chal"<<i<<".dat";
   
    std::fstream File;
    File.open(mystream.str());
    for(i=0;i<=N;i++){      
      File <<x[i]<< "\t" << U[i]<<'\n';
    }   
    File.close();
  }

  double ERR = 0.;
  for( i = 0; i < N ; i++) {
    ERR += g(U[i],t,dt,dx,-1/pow(3,1/2))-g(U[i],t,dt,dx,1/pow(3,1/2));
  }
  
  std::cout<< dx << "\t" <<  sqrt(ERR )<<std::endl;

  
  return 0;
}
