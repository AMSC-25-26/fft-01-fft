#include <iostream>
#include <complex>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using namespace std::complex_literals;
using namespace Eigen;
/*
 * A function to reverse the bits of a given index x.
 * The idea is to take the lsb of x and append it to
 * the msb of n
*/
unsigned int bitReverse(unsigned int x, int log2n) {
    int n = 0;
    for (int i = 0; i < log2n; i++) {
        n = n << 1;
        n = n | (x & 1);
        x = x >> 1;
    }
    return n;
}

MatrixXcd omega_n_half_builder(int N){
    MatrixXcd omega_n_half(N/2, N/2); 
    for(int j = 0; j < N/2; j++){      
            std::complex<double> exponent = -2.0 * M_PI * 1i * (double)(j) / (double)N;
            omega_n_half(j, j) = std::exp(exponent);
    }
    return omega_n_half;
}

int main(){
    int N = 8; 
    
    MatrixXcd omega_n_half = omega_n_half_builder(N);

    IOFormat CleanFmt(3, 0, "\t", "\n", "", "");

    std::cout << "--- Complete Fn ---" << std::endl;
    std::cout << omega_n_half.format(CleanFmt) << std::endl;

    return 0;
}