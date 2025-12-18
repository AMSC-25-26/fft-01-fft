#include <iostream>
#include <vector>
#include <cmath>
#include <complex>

namespace FFT {
    using complex = std::complex<double>;

    void recursive(std::vector<complex> &A);
    void iterative(std::vector<complex> &A);
    void inverse(std::vector<complex> &A);
    void parallel_iterative(std::vector<complex> &A);
    void parallel_inverse(std::vector<complex> &A);
}