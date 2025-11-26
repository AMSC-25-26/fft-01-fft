#include "FFT.hpp"

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

