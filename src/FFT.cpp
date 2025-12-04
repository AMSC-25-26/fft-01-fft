#include "FFT.hpp"
#include "utils.h"
#include <omp.h>

using namespace std::complex_literals;
using complex = std::complex<double>;

std::vector<std::complex<double>> FFT::recursive(const std::vector<std::complex<double>> A) {
    std::complex<double> N = A.size();
    int n = A.size(); 
    if (n == 1) return A;

    std::complex<double> Wn = std::exp((-2.0 * M_PI * 1i) / N);
    std::complex<double> W = 1;
    
    std::vector<std::complex<double>> A_even;
    std::vector<std::complex<double>> A_odd;
    
    for (int i = 0; i < n / 2; i++) {
        A_even.emplace_back(A[2 * i]);
        A_odd.emplace_back(A[2 * i + 1]);
    }

    std::vector<std::complex<double>> Y_even = recursive(A_even);
    std::vector<std::complex<double>> Y_odd = recursive(A_odd);

    std::vector<std::complex<double>> Y(n); 
    
    for (int j = 0; j < n / 2; j++) {
        std::complex<double> u = Y_even[j];
        std::complex<double> v = W * Y_odd[j];
        
        Y[j] = u + v;
        Y[j + n / 2] = u - v;
        W *= Wn;
    }

    return Y;
}

const std::complex<double> I(0.0, 1.0);

void FFT::iterative(std::vector<std::complex<double>>& A) {
    int N = A.size();
    if (N == 0 || (N & (N - 1)) != 0) {
        throw std::invalid_argument("Input vector size must be a power of two.");
    }
    applyBitReversalPermutation(A);

    int k = 2; 
    while (k <= N) { 
        int separation = k / 2; 
        int num_blocks = N / k; 

        for (int r = 0; r < num_blocks; ++r) {
            int block_start_index = r * k; 
            

            for (int j = 0; j < separation; ++j) {

            
                int top_idx    = block_start_index + j;
                int bottom_idx = block_start_index + j + separation;

               
                double exponent_val = -2.0 * M_PI * (double)j / (double)k;
 
                
                std::complex<double> X_top_old = A[top_idx];
                
             
                std::complex<double> T = W_k_j * A[bottom_idx];
                
               
                A[bottom_idx] = X_top_old - T;

                
                A[top_idx] = X_top_old + T; 
            }
        }
        
        k *= 2; 
    }
}

void FFT::inverse(std::vector<std::complex<double>>& A) {
    int N = A.size();
    if (N == 0 || (N & (N - 1)) != 0) {
        throw std::invalid_argument("Input vector size must be a power of two.");
    }

    applyBitReversalPermutation(A);     //reorders the input array in bit-reversed order to avoid recursion 

    int k = 2; 
    while (k <= N) { 
        int separation = k / 2; 
        int num_blocks = N / k; 

        for (int r = 0; r < num_blocks; ++r) {
            int block_start_index = r * k; 

            for (int j = 0; j < separation; ++j) {
            
                int top_idx    = block_start_index + j;
                int bottom_idx = block_start_index + j + separation;
               
                double exponent_val = 2.0 * M_PI * (double)j / (double)k;   // sign change
                std::complex<double> W_k_j = std::exp(I * exponent_val); 
                
                std::complex<double> X_top_old = A[top_idx];
                std::complex<double> T = W_k_j * A[bottom_idx];
               
                A[bottom_idx] = X_top_old - T;
                A[top_idx] = X_top_old + T; 
            }
        }
        
        k *= 2; 
    }
    
    // normalization 
    for (int i = 0; i < N; ++i) {
        A[i] /= N;

    }
}


void FFT::parallel_iterative(std::vector<std::complex<double>>& A) {
    int N = A.size();
    if (N == 0 || (N & (N - 1)) != 0) {
        throw std::invalid_argument("Input vector size must be a power of two.");
    }

    double wtime = omp_get_wtime();

    applyBitReversalPermutation(A);
    


    int k = 2; 
    while (k <= N) { 
        int separation = k / 2; 
        int num_blocks = N / k; 

        
        #pragma omp parallel for schedule(static)
        for (int r = 0; r < num_blocks; ++r) {
            int block_start_index = r * k; 
            

            for (int j = 0; j < separation; ++j) {

            
                int top_idx    = block_start_index + j;
                int bottom_idx = block_start_index + j + separation;

               
                double exponent_val = -2.0 * M_PI * (double)j / (double)k;
                std::complex<double> W_k_j = std::exp(I * exponent_val); //defining  

                
                std::complex<double> X_top_old = A[top_idx];
                
             
                std::complex<double> T = W_k_j * A[bottom_idx];
                
               
                A[bottom_idx] = X_top_old - T;

                
                A[top_idx] = X_top_old + T; 
            }
        }
        
        k *= 2; 
    }

    
    wtime = omp_get_wtime() - wtime;
    std::cout << "Parallel FFT Execution Time: " << wtime << " seconds" << std::endl;
}

void FFT::iterative(std::vector<std::complex<double>>& A) {
    int N = A.size();
    if (N == 0 || (N & (N - 1)) != 0) {
        throw std::invalid_argument("Input vector size must be a power of two.");
    }
    applyBitReversalPermutation(A);

    int k = 2; 
    while (k <= N) { 
        int separation = k / 2; 
        int num_blocks = N / k; 

        for (int r = 0; r < num_blocks; ++r) {
            int block_start_index = r * k; 
            

            for (int j = 0; j < separation; ++j) {

            
                int top_idx    = block_start_index + j;
                int bottom_idx = block_start_index + j + separation;

               
                double exponent_val = -2.0 * M_PI * (double)j / (double)k;
 
                
                std::complex<double> X_top_old = A[top_idx];
                
             
                std::complex<double> T = W_k_j * A[bottom_idx];
                
               
                A[bottom_idx] = X_top_old - T;

                
                A[top_idx] = X_top_old + T; 
            }
        }
        
        k *= 2; 
    }
}

void FFT::inverse_parallel(std::vector<std::complex<double>>& A) {
    int N = A.size();
    if (N == 0 || (N & (N - 1)) != 0) {
        throw std::invalid_argument("Input vector size must be a power of two.");
    }

    applyBitReversalPermutation(A);     //reorders the input array in bit-reversed order to avoid recursion 

    int k = 2; 
    while (k <= N) { 
        int separation = k / 2; 
        int num_blocks = N / k; 

        #pragma omp parallel for schedule(static)
        for (int r = 0; r < num_blocks; ++r) {
            int block_start_index = r * k; 

            for (int j = 0; j < separation; ++j) {
            
                int top_idx    = block_start_index + j;
                int bottom_idx = block_start_index + j + separation;
               
                double exponent_val = 2.0 * M_PI * (double)j / (double)k;   // sign change
                std::complex<double> W_k_j = std::exp(I * exponent_val); 
                
                std::complex<double> X_top_old = A[top_idx];
                std::complex<double> T = W_k_j * A[bottom_idx];
               
                A[bottom_idx] = X_top_old - T;
                A[top_idx] = X_top_old + T; 
            }
        }
        
        k *= 2; 
    }
    
    // normalization 
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; ++i) {
        A[i] /= N;

    }
}