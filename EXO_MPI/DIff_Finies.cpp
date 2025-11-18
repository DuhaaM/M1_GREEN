#include <iostream>
#include <math.h>
#include <fstream>


int main (int argc, char *argv[]){
  int i,j;
  int N=10,K=1,L=1;
  double dt=1e-4,t=0,dx,T;
  double U[N+1],Uold[N+1],x[N+1];

  dx=2*L/N;
  T=1/dt;
  std::fstream File("Chal.dat");
  
  for(i=0;i<=N;i++){
    x[i]= -L+i*dx;
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
      Uold[i]=U[i]+(K*dt/dx/dx)*(U[i-1]-2*U[i]+U[i+1]);
    }
    
    for(i=0;i<=N;i++){
      U[i]=Uold[i];
    }

    std::cout<<t<<std::endl;
    
    for(i=0;i<=N;i++){
      File<<x[i],U[i],"\n";
    }
  }
  return 0;
}
