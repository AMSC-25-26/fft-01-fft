#include "FFT.hpp"
#include "utils.h"

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
    size_t N = A.size();
    if (N == 0 || (N & (N - 1)) != 0) {
        throw std::invalid_argument("Input vector size must be a power of two.");
    }

    double theta = -2 * M_PI;
    int log2N = static_cast<int>(std::log2(N));

    applyBitReversalPermutation(A);

     // We iterate for log2(N) steps
    for (int k = 1; k <= log2N; ++k) { 
        int m = 1 << k; // 2^k So 
        int m2 = m >> 1; // m / 2 - equivalent to the separation used before

        // Principle root of nth complex
        // root of unity. We compute it once!
        // Then we build each twiddle using
        // the powers of it
        std::complex<double> wm = std::polar(1.0, theta / static_cast<double> (m)); // e^-(i*2*PI) / (2^k)

        // We iterate over the blocks
        for (size_t r = 0; r < N; r += m) {
            std::complex<double> w{1.0, 0};
            for (size_t j = 0; j < m2; ++j) {    
                std::complex<double> X_top_old = A[r + j];
                std::complex<double> T = w * A[r + j + m2];
               
                A[r + j + m2] = X_top_old - T;
                A[r + j] = X_top_old + T;
                w *= wm; // constructing e^(-2pi/k)^j
            }
        }         
    }
}

void FFT::inverse(std::vector<std::complex<double>>& A) {
    size_t N = A.size();
    if (N == 0 || (N & (N - 1)) != 0) {
        throw std::invalid_argument("Input vector size must be a power of two.");
    }

    double theta = 2 * M_PI;
    int log2N = static_cast<int>(std::log2(N));

    applyBitReversalPermutation(A);     //reorders the input array in bit-reversed order to avoid recursion 

    for (int k = 1; k <= log2N; ++k) { 
        int m = 1 << k; // 2^k So 
        int m2 = m >> 1; // m / 2 - equivalent to the separation used before

        // Principle root of nth complex
        // root of unity. We compute it once!
        // Then we build each twiddle using
        // the powers of it
        std::complex<double> wm = std::polar(1.0, theta / static_cast<double> (m)); // e^-(i*2*PI) / (2^k)

        // We iterate over the blocks
        for (size_t r = 0; r < N; r += m) {
            std::complex<double> w{1.0, 0};
            for (size_t j = 0; j < m2; ++j) {    
                std::complex<double> X_top_old = A[r + j];
                std::complex<double> T = w * A[r + j + m2];
               
                A[r + j + m2] = X_top_old - T;
                A[r + j] = X_top_old + T;
                w *= wm; // constructing e^(2pi/k)^j
            }
        }         
    }
    
    // normalization 
    for (size_t i = 0; i < N; ++i) {
        A[i] /= N;

    }
}