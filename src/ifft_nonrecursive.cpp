// the only differences from fft_nonrecursive are the sign of the exponent and the multiplication by 1/N

#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include "utils.h"

const std::complex<double> I(0.0, 1.0); // Define i

void ifft_nonrecursive(std::vector<std::complex<double>>& A) {
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




std::vector<std::complex<double>> IDFT_Compute(const std::vector<std::complex<double>>& A) {
    int N = A.size();
    std::vector<std::complex<double>> X(N);

    for (int k = 0; k < N; k++) {
        X[k] = 0.0;
        for (int n = 0; n < N; n++) {
            double exponent = 2.0 * M_PI * (double)(k * n) / (double)N;
            std::complex<double> W = std::exp(I * exponent);
            X[k] += A[n] * W;
        }
        X[k] /= N;
    }

    return X;
}


int main() {
    // ==========================================
    // TEST 2: Simple 8-point IFFT test (FFT^-1)
    // ==========================================

    std::cout << "\n--- 8-Point IFFT Test ---" << std::endl;

    // input spectrum X[k] = FFT{x[n]} where x[n] = [1,2,3,4,5,6,7,8]
    int N = 8;
    std::vector<std::complex<double>> input_signal_ifft = {
        {36.0, 0.0}, {-4.0, 9.6569}, {-4.0, 4.0}, {-4.0, 1.6569},
        {-4.0, 0.0}, {-4.0, -1.6569}, {-4.0, -4.0}, {-4.0, -9.6569}
    };

    // Copy for IDFT verification
    std::vector<std::complex<double>> input_signal_idft = input_signal_ifft;

    std::cout << "Input Spectrum (X[k]):" << std::endl;
    for (int i = 0; i < N; ++i) {
        std::cout << "[" << i << "] " << input_signal_ifft[i] << "\n";
    }
    std::cout << "\n";

    // 1. IFFT (Fast)
    ifft_nonrecursive(input_signal_ifft);

    // 2. IDFT (Slow)
    std::vector<std::complex<double>> result_idft = IDFT_Compute(input_signal_idft);

    std::cout << "Results (IFFT vs. IDFT):" << std::endl;
    for (int i = 0; i < N; ++i) {
        std::cout << "[" << i << "] IFFT: " << input_signal_ifft[i] 
                << " | IDFT: " << result_idft[i] << "\n";
    }

    std::cout << "\nVerification:\n";
    std::cout << "The IFFT results should match the IDFT results (within floating-point error)." << std::endl;
    return 0;
}
