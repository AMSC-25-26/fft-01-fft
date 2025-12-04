#include <iostream>
#include <vector>
#include <cmath>
#include <complex>


class FFT 
{
public:
    using complex = std::complex<double>;

    static std::vector<complex> recursive(const std::vector<complex> A);
    static void iterative(std::vector<complex> &A);
    static void inverse(std::vector<complex> &A);
    static void parallel_iterative(std::vector<complex> &A);
    static void parallel_inverse(std::vector<complex> &A);
};