#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <iomanip> 

using namespace std::complex_literals;

std::vector<std::complex<double>> prova_fft(const std::vector<std::complex<double>> A) {
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

    std::vector<std::complex<double>> Y_even = prova_fft(A_even);
    std::vector<std::complex<double>> Y_odd = prova_fft(A_odd);

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

int main() {
    std::vector<std::complex<double>> input = {
        1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0
    };

    if ((input.size() & (input.size() - 1)) != 0) {
        std::cerr << "Error: N must be in the form 2^m" << std::endl;
        return -1;
    }

    std::cout << "Input Signal (time domain)" << std::endl;
    for (size_t i = 0; i < input.size(); ++i) {
        std::cout << "A[" << i << "] = " << input[i] << std::endl;
    }

    std::vector<std::complex<double>> output = prova_fft(input);

    std::cout << "\n Output FFT (frequency domain) ---" << std::endl;
    std::cout << std::fixed << std::setprecision(3);
    
    for (size_t k = 0; k < output.size(); ++k) {
        double real = std::abs(output[k].real()) < 1e-9 ? 0.0 : output[k].real();
        double imag = std::abs(output[k].imag()) < 1e-9 ? 0.0 : output[k].imag();
        
        std::cout << "Y[" << k << "] = " << real 
                  << (imag >= 0 ? " + " : " - ") 
                  << std::abs(imag) << "i" << std::endl;
    }

    return 0;
}