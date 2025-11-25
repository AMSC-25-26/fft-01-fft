#include <iostream>
#include <complex>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using namespace std::complex_literals;
using namespace Eigen;

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