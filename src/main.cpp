#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <iomanip>
#include "FFT.hpp"

std::vector<std::complex<double>> IDFT_Compute(const std::vector<std::complex<double>>& A);
std::vector<std::complex<double>> DFT_Compute(const std::vector<std::complex<double>>& A);

int main() {
    // Set precision for printing complex numbers
    std::cout << std::fixed << std::setprecision(4);

    // TEST 1: Simple 8-point FFT (N=8)
    int N = 8;
    
    // Define the input signal x = [1, 2, 3, 4, 5, 6, 7, 8] + 0i
    std::vector<std::complex<double>> input_signal = {
        {1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}, 
        {5.0, 0.0}, {6.0, 0.0}, {7.0, 0.0}, {8.0, 0.0}
    };
    // Copy for DFT verification
    std::vector<std::complex<double>> input_signal_dft = input_signal;
    std::vector<std::complex<double>> input_signal_fft = input_signal;

    

    std::cout << "--- 8-Point FFT Test ---" << std::endl;
    std::cout << "Input Signal (x[n]):" << std::endl;
    for(size_t i = 0; i < N; ++i) {
        std::cout << "[" << i << "] " << input_signal_fft[i] << "\n";
    }
    std::cout << "\n";

    /*
    // 1. Compute FFT
    FFT::iterative(input_signal_fft);

    // 2. Compute DFT
    std::vector<std::complex<double>> result_dft = DFT_Compute(input_signal_dft);

    std::cout << "Results (FFT vs. DFT):" << std::endl;
    for (int i = 0; i < N; ++i) {
        // Output format: [k] FFT: (real + imag*i) | DFT: (real + imag*i)
        std::cout << "[" << i << "] FFT: " << input_signal_fft[i] 
                  << " | DFT: " << result_dft[i] << "\n";
    }

    FFT::inverse(input_signal_fft);
    std::cout << "Results (IFFT vs. INPUT):" << std::endl;
    for (int i = 0; i < N; ++i) {
        // Output format: [k] FFT: (real + imag*i) | DFT: (real + imag*i)
        std::cout << "[" << i << "] IFFT: " << input_signal_fft[i] 
                  << " | INPUT: " << input_signal[i] << "\n";
    }

    */

    // 1. Compute Parallel FFT
    FFT::parallel_iterative(input_signal_fft);
    
    std::cout << "\nVerification:\n";
    std::cout << "The FFT results should match the DFT results (within floating-point error)." << std::endl;
    std::cout << "The IFFT results should match the input signal (within floating-point error)." << std::endl;

    return 0;
}
 
const std::complex<double> I(0.0, 1.0);

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