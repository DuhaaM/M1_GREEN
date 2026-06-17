#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <string>
#include <algorithm>
#include <cuda_runtime.h>
#include "Solution2D.hpp"

#define CUDA_CHECK(call) do { \
    cudaError_t err_ = (call); \
    if (err_ != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d : %s\n", __FILE__, __LINE__, \
                cudaGetErrorString(err_)); \
        exit(EXIT_FAILURE); \
    } \
} while(0)


__global__ void euler_exp_kernel(double* __restrict__ Utemp, const double* __restrict__ U, int NxLoc, int NyLoc, const double dtloc,const double dx, const double dy){

    int i = blockIdx.x * blockDim.x + threadIdx.x +2;
    int j = blockIdx.y * blockDim.y + threadIdx.y +2;

    if (i>= NxLoc || j >= NyLoc) return;

    const int sy = NxLoc+2;

    int I = i+ j*sy;

    Utemp[I]= U[I] + dtloc/dx/dx * (U[I+1]  - 2*U[I] + U[I-1])
    	      + dtloc/dy/dy * (U[I+sy] - 2*U[I] + U[I-sy]);
}

__global__ void copy_kernel ( double* __restrict__ U, const double* __restrict__ Utemp, int NxLoc, int NyLoc){

    int i = blockIdx.x * blockDim.x + threadIdx.x +2;
    int j = blockIdx.y * blockDim.y + threadIdx.y +2;
    const int sy = NxLoc+2;

    if (i>= NxLoc || j >= NyLoc) return;

    int I = i+ j*sy;

    U[I] = Utemp[I];
}


static void write(int Ind,const double* U, const double* x, const double* y, int NxLoc, int NyLoc, Solution* Sol , double T){

   
   std::stringstream mystream;
   mystream << "Data2D/Chal" << Ind << ".csv";
   
   std::ofstream File(mystream.str(), std::ofstream::out);
   File << "#Indice,#Point_x,#Point_y,#Solution_obtenue,#Solution_exacte\n";

   int m = 0;
   for (int j = 0; j <= NyLoc; j++) {
     for (int i = 0; i <= NxLoc; i++) {
       File << m << "," << x[i] << "," << y[j] << ","  << U[m] << "," << Sol->Exact(x[i], y[j], T) << "\n";
       m++;
     }
   }
   
   File.close();
}


static double ErrorInf (const double* U, const double* x, const double* y, int NxLoc, int NyLoc, Solution* Sol , double T){

   const int sy = NxLoc +2 ;
   double err=0.0;
   for (int j = 1; j <= NyLoc; j++){
     for (int i = 1; i <= NxLoc; i++) {
       int I = i + j*sy;
       err = std::max(err, std::abs(U[I]- Sol->Exact(x[i-1], y[j-1], T)));
     }
   }
   return err;
}


int main(int argc, char* argv[]){
    if (argc != 10) {
        std::cerr << "usage: " << argv[0]
                  << " <Nx> <Ny>  <Tfinal> <K> <kprint>"
                     " <Lx> <Ly> "
                     " <Function: 1=Sinus, 2=Linear, 3=Constant> <Blockdim>"
                  << std::endl;
        return 1;
    }

    int    Nx       = std::atoi(argv[1]);
    int    Ny       = std::atoi(argv[2]);
    double Tfinal   = std::atof(argv[3]);
    int    K        = std::atoi(argv[4]);
    int    kprint   = std::atoi(argv[5]);
    double Lx       = std::atof(argv[6]);
    double Ly       = std::atof(argv[7]);
    int    Function = std::atoi(argv[8]);
    int    Bdim     = std::atoi(argv[9]);
    
    Solution* Sol = nullptr;
    if (Function == 1) Sol = new Sin(K, Lx, Ly);
    if (Function == 2) Sol = new Lin();
    if (Function == 3) Sol = new Const();

    const double dx = Lx / double(Nx);
    const double dy = Ly / double(Ny);

    double dtloc= 0.2 * std::min(dx,dy) * std::min(dx,dy);      

    double* x = new double[Nx+1];
    double* y = new double[Ny+1];

    for (int i=0; i <= Nx; i++) x[i]= i * Lx / double(Nx);

    for (int j=0; j <= Ny; j++) y[j]= j * Ly / double(Ny);

    const int sy = Nx+1 +2;
    const int total = sy* ( Ny+1 +2);

    double* U_host = new double[total];
    for (int j=1; j<= Ny+1; j++){
      for (int i=1; i<=Nx+1; i++){
      	U_host[i+ j*sy ] = Sol ->Initial(x[i-1], y[j-1]);
      }
    }
    double *d_U, *d_Utemp;
    CUDA_CHECK(cudaMalloc(&d_U, total*sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_Utemp, total*sizeof(double)));
    CUDA_CHECK(cudaMemcpy(d_U,U_host,total*sizeof(double),cudaMemcpyHostToDevice));

    dim3 block(Bdim,Bdim);
    dim3 grid ((Nx - 1+ Bdim-1)/Bdim , (Ny - 1+Bdim-1)/Bdim);

    cudaEvent_t ev_start,ev_stop;

    CUDA_CHECK(cudaEventCreate(&ev_start));
    CUDA_CHECK(cudaEventCreate(&ev_stop));
    float elapsed_ms= 0.0f ,  gpu_ms=0.0f;

    double t =0.0;
    int kT=0, Ind=0;

    while( t<Tfinal ){

      if ( t+ dtloc > Tfinal) dtloc= Tfinal-t;

      CUDA_CHECK(cudaEventRecord(ev_start));

      euler_exp_kernel <<< grid,block >>> (d_Utemp,d_U,Nx+1,Ny+1,dtloc,dx,dy);

      CUDA_CHECK(cudaGetLastError());

      copy_kernel <<< grid,block >>>(d_U, d_Utemp,Nx+1,Ny+1);

      CUDA_CHECK(cudaGetLastError());

      CUDA_CHECK(cudaEventRecord(ev_stop));
      CUDA_CHECK(cudaEventSynchronize(ev_stop));
      CUDA_CHECK(cudaEventElapsedTime(&gpu_ms, ev_start, ev_stop));

      elapsed_ms += gpu_ms;

      t+= dtloc;
      kT++;

      if(kT % kprint == 0){

        CUDA_CHECK(cudaMemcpy(U_host, d_U, total*sizeof(double),cudaMemcpyDeviceToHost));
	write(Ind,U_host,x,y,Nx+1,Ny+1,Sol,t);
	Ind++;

      }
   }

   std::cout << "Computing time = " <<elapsed_ms/1000<<std::endl;

   CUDA_CHECK(cudaMemcpy(U_host, d_U, total*sizeof(double),cudaMemcpyDeviceToHost));
   	
   std::cout<< "ErreurInf = " << ErrorInf(U_host,x,y,Nx+1,Ny+1,Sol,Tfinal) << std::endl;

   CUDA_CHECK(cudaFree(d_U));                                 
   CUDA_CHECK(cudaFree(d_Utemp));                               
   CUDA_CHECK(cudaEventDestroy(ev_start));                      
   CUDA_CHECK(cudaEventDestroy(ev_stop));                      

   return 0;
}