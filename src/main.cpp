#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <iomanip>
#include <random>
#include "FFT.hpp"

using complex = std::complex<double>; 
using complexVector = std::vector<complex>; 

complexVector IDFT_Compute(const complexVector& A);
complexVector DFT_Compute(const complexVector& A);
void test_polynomial_multiplication(size_t n1, size_t n2); 

int main(int argc, char* argv[]) {
    if(argc < 2){
        std::cerr << "Use: <int> <int>" << std::endl; 
        return -1; 
    }

    /**
     *  Set precision for printing complex numbers
     */
    std::cout << std::fixed << std::setprecision(4);

    /** 
     * TEST 1: Simple 8-point FFT (N=8) 
    */
    const size_t N = 8;
    
    /**
     *  Define the input signal x = [1, 2, 3, 4, 5, 6, 7, 8] + 0i
     */
    complexVector input_signal = {
        {1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}, 
        {5.0, 0.0}, {6.0, 0.0}, {7.0, 0.0}, {8.0, 0.0}
    };
    
    /**
     * Copy for DFT verification
     */
    complexVector input_signal_dft = input_signal;
    complexVector input_signal_fft = input_signal;

    std::cout << "=== 8-Point FFT Test ===" << std::endl;
    std::cout << "Input Signal (x[n]):" << std::endl;
    for(size_t i = 0; i < N; ++i) {
        std::cout << "[" << i << "] " << input_signal_fft[i] << "\n";
    }
    std::cout << "\n";

    // 1. Compute FFT
    FFT::parallel_iterative(input_signal_fft);

    // 2. Compute DFT
    complexVector result_dft = DFT_Compute(input_signal_dft);

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
    complexVector sampled_signal(N2);
    for (size_t n = 0; n < N2; n++) {
        double t = n * dt;
        sampled_signal[n] = std::complex<double>(func(t), 0.0);
    }

    // compute fft
    complexVector fft_result = sampled_signal;
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

    /**
     * Polynomial multiplication test
     */
    std::cout << "\n=== POLYNOMIAL MULTIPLICATION TEST ===" << std::endl; 
    const size_t n1 = atoi(argv[1]); 
    const size_t n2 = atoi(argv[2]); 
    test_polynomial_multiplication(n1, n2); 
    
    return 0;
}
 
const std::complex<double> I(0.0, 1.0);

complexVector generate_random_polynomial(size_t n)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-10.0, 10.0);
    
    std::vector<complex> poly(n);
    for(size_t i = 0; i < n; ++i) {
        poly[i] = {dis(gen), 0.0}; 
    }
    return poly;
}

size_t next_power_of_two(size_t n) {
    size_t power = 1;
    while (power < n) {
        power <<= 1;
    }
    return power;
}

complexVector naive_polynomial_multiplication(const complexVector& p1, const complexVector& p2)
{
    size_t n1 = p1.size();
    size_t n2 = p2.size(); 
    size_t nC = n1 + n2 - 1; 
    
    complexVector C(nC, {0.0, 0.0}); 

    for(size_t i = 0; i < n1; ++i) 
        for(size_t j = 0; j < n2; ++j) 
            C[i + j] += p1[i] * p2[j];
        
    return C; 
}

bool same_polynomial(complexVector p1, complexVector p2)
{
    size_t n1 = p1.size(), n2 = p2.size(); 
    if(n1 != n2) 
        return false;
    for(size_t i = 0; i < n1; ++i)
        if(double(abs(p1[i].real() - p2[i].real())) > 1.0e-12)
            return false; 

    return true; 
}

void print_polynomial(complexVector p)
{
    const size_t N = p.size(); 
    if(N > 0) 
    {
        std::cout << p[0].real() << " + "; 
        for(size_t i = 1; i < N-1; ++i)
            std::cout << p[i].real() << "x^" << i << " + "; 
        
        std::cout << p[N-1].real() << std::endl; 
    }
}

void test_polynomial_multiplication(size_t n1, size_t n2)
{
    complexVector poly1 = generate_random_polynomial(n1);    
    complexVector poly2 = generate_random_polynomial(n2);    
    complexVector naiveC = naive_polynomial_multiplication(poly1, poly2); 
    
    size_t final_size = n1 + n2 - 1; 
    size_t nC = next_power_of_two(final_size); 
    poly1.resize(nC, {0, 0}); 
    poly2.resize(nC, {0, 0}); 
    FFT::parallel_iterative(poly1); 
    FFT::parallel_iterative(poly2); 

    complexVector C(nC); 
    for(size_t i = 0; i < nC; ++i)
        C[i] = poly1[i] * poly2[i];     

    FFT::parallel_inverse(C); 

    C.resize(final_size); 
    std::cout << "FFT multiplication resulting polynomial: " << std::endl; 
    print_polynomial(C); 
    std::cout << "Naive multiplication resulting polynomial: " << std::endl; 
    print_polynomial(naiveC); 
    std::cout << (same_polynomial(C, naiveC) ? "Same polynomial - FFT Success" : "Not same polynomial - FFT Failed!") << std::endl; 
}


complexVector IDFT_Compute(const complexVector& A) {
    const size_t N = A.size();
    complexVector X(N);

    for (size_t k = 0; k < N; k++) {
        X[k] = 0.0;
        for (size_t n = 0; n < N; n++) {
            double exponent = 2.0 * M_PI * (double)(k * n) / (double)N;
            complex W = std::exp(I * exponent);
            X[k] += A[n] * W;
        }
        X[k] /= N;
    }

    return X;
}

complexVector DFT_Compute(const complexVector& A) {
    size_t N = A.size();
    complexVector X(N);

    for (size_t k = 0; k < N; k++) {
        X[k] = 0.0;
        for (size_t n = 0; n < N; n++) {
            double exponent = -2.0 * M_PI * (double)(k * n) / (double)N;
            complex W = std::exp(I * exponent);
            X[k] += A[n] * W;
        }
    }
    return X;
}

