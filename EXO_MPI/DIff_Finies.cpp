#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>


double u_init(const double& x){

  return sin(M_PI*x);

}

double u_num(const double& l,const int& i, const double* const U,const double* const x){
  std::cout << "i = " << i << std::endl;
  std::cout << "U[i] = " << U[i] << "\t U[i+1] = " << U[i+1] << "\t x[i] = " << x[i]<< "\t x[i+1] = " << x[i+1] << std::endl;
  return (U[i+1]-U[i])*(l-x[i])/(x[i+1]-x[i]) + U[i];

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
  

  double U[N+1],UTemp[N+1],x[N+1];
  
  dx=1./double(N);
  dt = 0.4 * dx * dx;
  //  T=1/dt;
  
  for(i=0;i<=N;i++){
    x[i]= i*dx;
    U[ i ] = u_init(x[i]); 
    //Uold[i]=1;
  }

  std::cout << "x = " ;
  for( i = 0; i <= N; i++)
    std::cout << x[i] << "\t";
  std::cout << std::endl;

  double t = 0.;
  double dtloc = 0.;
  int kT = 0;
  std::cout << "dt = " << dt << std::endl;
  while( t < Tfinal ) {
    dtloc = dt;
    if ( t + dtloc > Tfinal ) dtloc = Tfinal - t ;
    U[0]=0;
    U[N]=0;
    //Uold[0]=0;
    //Uold[N]=0;    
	
    for(i=1;i<N;i++){
	
      UTemp[i]=U[i]+ dt / ( dx * dx ) * ( U[i+1] - 2*U[i] + U[i-1] );
    }
    for(i=1;i<N;i++){
      U[i]=UTemp[i];	
    }

    t+=dtloc;
    std::cout << "kT = " << kT << "\t U[N] = " << U[N] << std::endl;

    std::stringstream mystream;
    mystream<<"Chal"<<kT<<".dat";
   
    std::ofstream File( mystream.str(), std::ofstream::out );
    std::cout << "mystream.str() = " << mystream.str() << std::endl;
    //File.open(mystream.str());
    for(i=0;i<=N;i++){      
      File <<x[i]<< "\t" << U[i]<<'\n';
    }   
    File.close();
    kT++;
  }
  
  double ERR = 0.;

  std::cout << "U = " ;
  for( i = 0; i <= N; i++)
    std::cout << U[i] << "\t";
  std::cout << std::endl;
  
  for( i = 0; i < N ; i++) {
    double xG1,xG2;

    xG1=(x[i+1]-x[i])*(-1/sqrt(3)+1)/2+x[i];
    xG2=(x[i+1]-x[i])*(+1/sqrt(3)+1)/2+x[i];
    
	 
    //std::cout<<"x[i]:"<<x[i]<<"\t"<<"x[i+1]:"<<x[i+1]<<std::endl;
    //std::cout<<u_num(xG1,i,U,x)<<"\t"<<u_num(xG2,i,U,x)<<std::endl;
    std::cout << "exp(-Tfinal) * u_init( xG1 ) = " << exp(-Tfinal) * u_init( xG1 ) << "\t u_num(xG1,i,U,x) = " << u_num(xG1,i,U,x) << std::endl;
    ERR += 1/(2*double(N)) * (exp(-Tfinal) * u_init( xG1 ) - u_num(xG1,i,U,x)) * (exp(-Tfinal) * u_init( xG1 ) - u_num(xG1,i,U,x));
    ERR += 1/(2*double(N)) * (exp(-Tfinal) * u_init( xG2 ) - u_num(xG2,i,U,x)) * (exp(-Tfinal) * u_init( xG2 ) - u_num(xG2,i,U,x));
  }
  std::cout<< "Erreur = " << dx << "\t" <<  sqrt(ERR)<<std::endl;
   
  std::ofstream File( "ChakExact.dat", std::ofstream::out );
  for(i=0;i<=N;i++){      
    File <<x[i]<< "\t" << exp(-Tfinal) * u_init( x[i] )<<'\n';
  }   
  File.close();

  
  return 0;
}

