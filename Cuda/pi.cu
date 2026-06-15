#include <iostream>
#include <cstdio>
#include <cmath>
#include <cuda_runtime.h>

#define CUDA_CHECK(call) do {\
	cudaError_t err_ = (call); \
	if (err_ != cudaSuccess) { \
	   fprintf(stderr,"CUDA error at %s:%d : %s\n",__FILE__,__LINE__, \
	   	   cudaGetErrorString(err_)); \
	   exit(EXIT_FAILURE);  \
	} \
} while(0)


__global__ void Somme_kernel(double* somme, double largeur, const double* __restrict__ x, long long N) {
    long long i = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    long long j = (long long)blockIdx.y * blockDim.y + threadIdx.y;
    long long k = (long long)blockIdx.z * blockDim.z + threadIdx.z;

    long long Nx = (long long)gridDim.x * blockDim.x;
    long long Ny = (long long)gridDim.y * blockDim.y;
    long long I  = i + j * Nx + k * Nx * Ny;

    if (I >= N) return;

    atomicAdd(somme, largeur * (4.0 / (1.0 + x[I] * x[I])));
}


int main(int argc, char* argv[]) {
    long long N = 1000LL * 1000 * 100;
    double largeur = 1.0 / N;
    double somme   = 0.0;

    double* x = new double[N];
    for (long long i = 0; i < N; i++) x[i] = largeur * (i + 0.5);

    double* d_x;
    CUDA_CHECK(cudaMalloc(&d_x, N * sizeof(double)));
    CUDA_CHECK(cudaMemcpy(d_x, x, N * sizeof(double), cudaMemcpyHostToDevice));

    double* d_somme;
    CUDA_CHECK(cudaMalloc(&d_somme, sizeof(double)));
    CUDA_CHECK(cudaMemcpy(d_somme, &somme, sizeof(double), cudaMemcpyHostToDevice));

    dim3 block(8, 8, 8);
    long long G = (long long)ceil(cbrt((double)N / 512.0));
    dim3 grid((unsigned int)G, (unsigned int)G, (unsigned int)G);

    cudaEvent_t ev_start, ev_stop;
    CUDA_CHECK(cudaEventCreate(&ev_start));
    CUDA_CHECK(cudaEventCreate(&ev_stop));
    float gpu_ms = 0.0f;

    CUDA_CHECK(cudaEventRecord(ev_start));

    Somme_kernel<<<grid, block>>>(d_somme, largeur, d_x, N);
    CUDA_CHECK(cudaGetLastError());

    CUDA_CHECK(cudaEventRecord(ev_stop));
    CUDA_CHECK(cudaEventSynchronize(ev_stop));
    CUDA_CHECK(cudaEventElapsedTime(&gpu_ms, ev_start, ev_stop));

    CUDA_CHECK(cudaMemcpy(&somme, d_somme, sizeof(double), cudaMemcpyDeviceToHost));

    std::cout << "Pi    = " << somme << std::endl;
    std::cout << "Temps = " << gpu_ms / 1000.0 << " s" << std::endl;

    CUDA_CHECK(cudaFree(d_x));
    CUDA_CHECK(cudaFree(d_somme));
    CUDA_CHECK(cudaEventDestroy(ev_start));
    CUDA_CHECK(cudaEventDestroy(ev_stop));
    delete[] x;

    return 0;
}
