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
    const size_t N = 8;
    
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

    // 1. Compute FFT
    FFT::parallel_iterative(input_signal_fft);

    // 2. Compute DFT
    std::vector<std::complex<double>> result_dft = DFT_Compute(input_signal_dft);

    std::cout << "Results (FFT vs. DFT):" << std::endl;
    for (size_t i = 0; i < N; ++i) {
        // Output format: [k] FFT: (real + imag*i) | DFT: (real + imag*i)
        std::cout << "[" << i << "] FFT: " << input_signal_fft[i] 
                  << " | DFT: " << result_dft[i] << "\n";
    }

    FFT::parallel_inverse(input_signal_fft);
    std::cout << "Results (IFFT vs. INPUT):" << std::endl;
    for (size_t i = 0; i < N; ++i) {
        // Output format: [k] FFT: (real + imag*i) | DFT: (real + imag*i)
        std::cout << "[" << i << "] IFFT: " << input_signal_fft[i] 
                  << " | INPUT: " << input_signal[i] << "\n";
    }
    
    std::cout << "\nVerification:\n";
    std::cout << "The FFT results should match the DFT results (within floating-point error)." << std::endl;
    std::cout << "The IFFT results should match the input signal (within floating-point error)." << std::endl;

    //now we test our fft algorithm with a signal represented by a mathematical function 

    const size_t N2 = 1024;     // Number of samples (must be power of 2)
    double T_total = 2.0;    // Total duration of the signal (seconds)
    double dt = T_total / N2; // Sampling interval
    double a = 5.0;          


    
    // we choose a time domain function: f(t) = e^(-at)
    auto func = [a](double t) {
        return std::exp(-a * t);
    };

    // we know the frequency domain FT: |F(w)| = 1 / sqrt(a^2 + w^2)
    auto analytical_FT_magnitude = [a](double f) {
        double omega = 2 * M_PI * f;
        return 1.0 / std::sqrt(a * a + omega * omega);
    };

    // sample the function
    std::vector<std::complex<double>> sampled_signal(N2);
    for (size_t n = 0; n < N2; n++) {
        double t = n * dt;
        sampled_signal[n] = std::complex<double>(func(t), 0.0);
    }

    // compute fft
    std::vector<std::complex<double>> fft_result = sampled_signal;
    FFT::parallel_iterative(fft_result); 
    
    
    // compare results
    std::cout<<"\n FFT vs Analytical FT Comparison \n"
              << std::left << std::setw(10) << "Freq(Hz)" 
              << std::setw(15) << "FFT (Scaled)" 
              << std::setw(15) << "Analytical FT" 
              << std::setw(15) << "Error" << "\n";
    std::cout << "\n";

    for (size_t k = 0; k <= N2 / 2; k++) {
        // Calculate the frequency
        double f = (double)k / T_total; 

        // 1. Get FFT Magnitude 
        double fft_mag = std::abs(fft_result[k]) * dt;

        // 2. Get analytical Magnitude
        double analytical_mag = analytical_FT_magnitude(f);

        // 3. Compute Error
        double error = std::abs(fft_mag - analytical_mag);

        std::cout << std::left << std::setw(10) << std::fixed << std::setprecision(2) << f 
                  << std::setw(15) << std::setprecision(4) << fft_mag 
                  << std::setw(15) << std::setprecision(4) << analytical_mag 
                  << std::setw(15) << error << "\n";
    }

    return 0;
}
 
const std::complex<double> I(0.0, 1.0);

std::vector<std::complex<double>> IDFT_Compute(const std::vector<std::complex<double>>& A) {
    const size_t N = A.size();
    std::vector<std::complex<double>> X(N);

    for (size_t k = 0; k < N; k++) {
        X[k] = 0.0;
        for (size_t n = 0; n < N; n++) {
            double exponent = 2.0 * M_PI * (double)(k * n) / (double)N;
            std::complex<double> W = std::exp(I * exponent);
            X[k] += A[n] * W;
        }
        X[k] /= N;
    }

    return X;
}

std::vector<std::complex<double>> DFT_Compute(const std::vector<std::complex<double>>& A) {
    size_t N = A.size();
    std::vector<std::complex<double>> X(N);

    for (size_t k = 0; k < N; k++) {
        X[k] = 0.0;
        for (size_t n = 0; n < N; n++) {
            double exponent = -2.0 * M_PI * (double)(k * n) / (double)N;
            std::complex<double> W = std::exp(I * exponent);
            X[k] += A[n] * W;
        }
    }
    return X;
}