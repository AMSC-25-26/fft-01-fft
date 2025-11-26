#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <iomanip>
#include "FFT.cpp"


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

    std::vector<std::complex<double>> output = FFT::recursive(input);

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