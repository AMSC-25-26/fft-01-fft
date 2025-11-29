#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include "utils.h"

const std::complex<double> I(0.0, 1.0); // Define i

void fft_nonrecursive(std::vector<std::complex<double>>& A) {
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
                std::complex<double> W_k_j = std::exp(I * exponent_val); //defining  

                
                std::complex<double> X_top_old = A[top_idx];
                
             
                std::complex<double> T = W_k_j * A[bottom_idx];
                
               
                A[bottom_idx] = X_top_old - T;

                
                A[top_idx] = X_top_old + T; 
            }
        }
        
        k *= 2; 
    }
}


std::vector<std::complex<double>> DFT_Compute(const std::vector<std::complex<double>>& A) {
    int N = A.size();
    std::vector<std::complex<double>> X(N);

    for (int k = 0; k < N; k++) {
        X[k] = 0.0;
        for (int n = 0; n < N; n++) {
            double exponent = -2.0 * M_PI * (double)(k * n) / (double)N;
            std::complex<double> W = std::exp(I * exponent);
            X[k] += A[n] * W;
        }
    }
    return X;
}


int main() {
    // Set precision for printing complex numbers
    std::cout << std::fixed << std::setprecision(4);

    // TEST 1: Simple 8-point FFT (N=8)
    int N = 8;
    
    // Define the input signal x = [1, 2, 3, 4, 5, 6, 7, 8] + 0i
    std::vector<std::complex<double>> input_signal_fft = {
        {1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}, 
        {5.0, 0.0}, {6.0, 0.0}, {7.0, 0.0}, {8.0, 0.0}
    };
    // Copy for DFT verification
    std::vector<std::complex<double>> input_signal_dft = input_signal_fft;

    std::cout << "--- 8-Point FFT Test ---" << std::endl;
    std::cout << "Input Signal (x[n]):" << std::endl;
    for(size_t i = 0; i < N; ++i) {
        std::cout << "[" << i << "] " << input_signal_fft[i] << "\n";
    }
    std::cout << "\n";

    // 1. Compute FFT (Fast)
    fft_nonrecursive(input_signal_fft);
    
    // 2. Compute DFT (Slow, for verification)
    std::vector<std::complex<double>> result_dft = DFT_Compute(input_signal_dft);

    std::cout << "Results (FFT vs. DFT):" << std::endl;
    for (int i = 0; i < N; ++i) {
        // Output format: [k] FFT: (real + imag*i) | DFT: (real + imag*i)
        std::cout << "[" << i << "] FFT: " << input_signal_fft[i] 
                  << " | DFT: " << result_dft[i] << "\n";
    }
    
    std::cout << "\nVerification:\n";
    std::cout << "The FFT results should match the DFT results (within floating-point error)." << std::endl;

    return 0;
}