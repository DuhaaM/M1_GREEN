#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>


int main (int argc, char *argv[]){
  int i,j,N;
  double dt=0.00001,t=0,T,dx,l;

  if (argc == 1 ) {
    return 0;
  }
  else {
    N=std::atoi(argv[1]);
  }
  

  double U[N+1],Uold[N+1],x[N+1],EX[N+1],ERR[N+1];
  
  dx=1./double(N);
  T=1/dt;
  
  for(i=0;i<=N;i++){
    x[i]= i*dx;
    U[i]=1;
    Uold[i]=1;
  }

  
  for(j=0;j<T;j++){
    t+=dt;
    U[0]=0;
    U[N]=0;
    Uold[0]=0;
    Uold[N]=0;    

    
	
    for(i=1;i<N;i++){
	
      Uold[i]=U[i]+(dt/dx/dx)*(U[i+1]-2*U[i]+U[i-1]);
      EX[i]=exp(double(-t))*sin(double(x[i])*M_PI);
    }
    for(i=0;i<=N;i++){
      U[i]=Uold[i];
      ERR[i]=abs(U[i]-EX[i]);
	
    }
    std::stringstream mystream;
    mystream<<"Chal"<<i<<".dat";
   
    std::fstream File(mystream.str());
    for(i=0;i<=N;i++){
      File <<x[i]<< "\t" << U[i]<<"\t"<<ERR[i]<<'\n';
    }   
    File.close();
  }
  
  std::cout<<std::max({ERR[i]});
  
  return 0;
}
