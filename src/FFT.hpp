#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <iomanip> 

class FFT 
{
public:
    using complex = std::complex<double>;
    

    static std::vector<complex> recursive(const std::vector<complex> A);
    static std::vector<complex> iterative(const std::vector<complex> A);
    static std::vector<complex> inverse(const std::vector<complex> A);
};