#include <iostream>
#include <math.h>

using namespace std;

__global__ void MatrixIteration(double *d_A, double *d_B, int n) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < n * n) {
        int row = idx / n, col = idx % n;
        if (row == 0 || row == n - 1 || col == 0 || col == n - 1) {
            d_B[idx] = d_A[idx];
        }
        else {
            double around[4];
            around[0] = d_A[idx+n+1];
            around[1] = d_A[idx+n-1];
            around[2] = d_A[idx-n+1];
            around[3] = d_A[idx-n-1];
            if (around[1] < around[0]) {
                double tmp = around[1];
                around[1] = around[0];
                around[0] = tmp;
            }
            if (around[3] < around[2]) {
                double tmp = around[3];
                around[3] = around[2];
                around[2] = tmp;
            }
            double increment = (around[0] < around[2]) ? (min(around[1], around[2])) : (min(around[0], around[3]));
            d_B[idx] = d_A[idx] + increment;
        }
    }
}

__global__ void CalculateSum(double *d_C, int dist, int n) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < n * n && idx % (dist * 2) == 0 && idx + dist < n * n) {
        d_C[idx] += d_C[idx + dist];
    }
}

__global__ void VerifyResult(double *d_F, double *d_A, double *d_C, int n) {
    d_F[0] = d_C[0];
    d_F[1] = d_A[37 * n + 47];
}

int main(int argc, char** argv) {
    int n = atoi(argv[1]), t = atoi(argv[2]);
    int m_size = n * n * sizeof(double);
    double *A = (double*)malloc(m_size);
    for (int i = 0; i < n * n; ++i) {
        int row = i / n, col = i % n;
        A[i] = (1 + cos(2 * row) + sin(col)) * (1 + cos(2 * row) + sin(col));
    }
    int thread_per_block = 1024;
    int block_per_grid = n * n / thread_per_block + 1;
    dim3 gridDim(block_per_grid, 1, 1);
    dim3 blockDim(thread_per_block, 1, 1);
    double *d_A;
    double *d_B;
    cudaMalloc(&d_A, m_size);
    cudaMalloc(&d_B, m_size);
    cudaMemcpy(d_A, A, m_size, cudaMemcpyHostToDevice);
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
    for (int i = 0; i < t; ++i) {
        if (i % 2 == 0) {
            MatrixIteration<<<gridDim, blockDim>>>(d_A, d_B, n);
        }
        else {
            MatrixIteration<<<gridDim, blockDim>>>(d_B, d_A, n);
        }
        cudaDeviceSynchronize();
    }
    double *d_C = d_A;
    double *d_F;
    cudaMalloc(&d_F, 2 * sizeof(double));
    for (int dist = 1; dist <= n * n; dist *= 2) {
        CalculateSum<<<gridDim, blockDim>>>(d_C, dist, n);
        cudaDeviceSynchronize();
    }
    VerifyResult<<<1, 1>>>(d_F, d_A, d_C, n);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float total_time = 0.0;
    cudaEventElapsedTime(&total_time, start, stop);
    double *h_F = (double*)malloc(2 * sizeof(double));
    cudaMemcpy(h_F, d_F, 2 * sizeof(double), cudaMemcpyDeviceToHost);
    cout << "Results using " << n << '*' << n << " matrix are: " << endl;
    cout << endl;
    cout << "Sum of the matrix A is: " << h_F[0] << endl;
    cout << "Value of the matrix at A(37, 47): " << h_F[1] << endl;
    cout << "Elapsed time is: " << total_time << endl;
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
    cudaFree(d_F);
    free(A);
    free(h_F);
    return 0;
}
