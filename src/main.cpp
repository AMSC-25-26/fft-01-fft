#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <iomanip>
#include <random>
#include <fstream>
#include <ctime>
#include "FFT.hpp"

#define THR 9

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
                  << "\t| DFT: " << result_dft[i] << "\n";
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
    
    // === Polynomial multiplication test === 
    std::cout << "\n=== POLYNOMIAL MULTIPLICATION TEST ===" << std::endl; 
    const size_t n1 = atoi(argv[1]); 
    const size_t n2 = atoi(argv[2]); 
    test_polynomial_multiplication(n1, n2); 

    

    return 0;
}
 
const std::complex<double> I(0.0, 1.0);

/**
 * @brief Generates a vector of random polynomial coefficients.
 * * Creates a polynomial of size `n` where each coefficient is a complex number 
 * with a real part uniformly distributed in the range [-20.0, 20.0] 
 * and a zero imaginary part.
 * * @param n The number of coefficients of the polynomial to generate.
 * @return complexVector A vector containing the generated random coefficients.
 */
complexVector generate_random_polynomial(size_t n)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-20.0, 20.0);
    
    std::vector<complex> poly(n);
    for(size_t i = 0; i < n; ++i) {
        poly[i] = {dis(gen), 0.0}; 
    }
    return poly;
}

/**
 * @brief Computes the smallest power of two greater than or equal to the input.
 * * This function iterates to find the nearest upper power of two. 
 * If @p n is already a power of two, it returns @p n. 
 * Useful for padding data sizes for algorithms like FFT.
 * @param n The minimum required size.
 * @return size_t The calculated power of two (>= n).
 */
size_t next_power_of_two(size_t n) {
    size_t power = 1;
    while (power < n) {
        power <<= 1;
    }
    return power;
}

/**
 * @brief Computes the product of two polynomials using direct convolution.
 * This function performs a naive multiplication by iterating through all pairs 
 * of coefficients. The operation has a time complexity of O(N * M), where N and M 
 * are the sizes of the coefficient vectors. 
 * The resulting vector size is (size of p1 + size of p2 - 1).
 * @param p1 The vector of coefficients for the first polynomial.
 * @param p2 The vector of coefficients for the second polynomial.
 * @return complexVector The coefficients of the resulting product polynomial.
 */
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
/**
 * @brief Checks if two polynomials are approximately equal.
 * * Verifies equality by comparing the polynomial sizes and checking if the 
 * difference between the real parts of their coefficients is within a 
 * tolerance of 1.0e-12.
 * @note This function ignores the imaginary parts of the coefficients.
 * @param p1 The first polynomial vector.
 * @param p2 The second polynomial vector.
 * @return bool True if sizes match and real coefficients are within tolerance, false otherwise.
 */
bool same_polynomial(complexVector p1, complexVector p2)
{
    size_t n1 = p1.size(), n2 = p2.size(); 
    if(n1 != n2) 
        return false;
    for(size_t i = 0; i < n1; ++i)
        if(double(abs(p1[i].real() - p2[i].real())) > 1.e-12)
            return false; 

    return true; 
}

/**
 * @brief Prints the polynomial to standard output.
 * Displays the polynomial in the format "c0 + c1x^1 + ... + cNx^N".
 * Only the real parts of the coefficients are printed.
 * @param p The polynomial vector to print.
 */
void print_polynomial(complexVector p)
{
    const size_t N = p.size(); 

    // Conditional statement to avoid redundant prints 
    if(N >= THR) return; 
    if(N > 0) 
    {
        std::cout << p[0].real() << " + "; 
        size_t i; 
        for(i = 1; i < N-1; ++i)
            std::cout << p[i].real() << "x^" << i << " + "; 
        
        std::cout << p[N-1].real() << "x^" << i << std::endl; 
    }
}

/**
 * @brief Validates the FFT-based multiplication against the naive approach.
 * Generates two random polynomials, multiplies them using both the naive convolution 
 * method and the FFT-based convolution theorem, and asserts if the results are identical.
 * Prints both resulting polynomials and a success/failure message to stdout.
 * @param n1 Size of the first random polynomial.
 * @param n2 Size of the second random polynomial.
 */
void test_polynomial_multiplication(size_t n1, size_t n2)
{
    complexVector poly1 = generate_random_polynomial(n1);    
    complexVector poly2 = generate_random_polynomial(n2);  
    
    // Conditional statement to avoid redundant prints 
    if(n1 < THR && n2 < THR)
    {
        std::cout << "P1(x) = "; 
        print_polynomial(poly1); 
        std::cout << "P2(x) = "; 
        print_polynomial(poly2); 
    }

    clock_t initNaive = clock(); 
    complexVector naiveC = naive_polynomial_multiplication(poly1, poly2); 
    clock_t endNaive = clock() - initNaive; 

    size_t final_size = n1 + n2 - 1; 
    size_t nC = next_power_of_two(final_size); 
    poly1.resize(nC, {0, 0}); 
    poly2.resize(nC, {0, 0}); 

    clock_t initFFT = clock(); 
    FFT::parallel_iterative(poly1); 
    FFT::parallel_iterative(poly2); 

    complexVector C(nC); 
    for(size_t i = 0; i < nC; ++i)
        C[i] = poly1[i] * poly2[i];     

    FFT::parallel_inverse(C); 
    clock_t endFFT = clock() - initFFT; 

    C.resize(final_size); 
    
    // Conditional statement to avoid redundant prints 
    if(final_size < THR) 
    {
        std::cout << "\nFFT multiplication resulting polynomial: "; 
        print_polynomial(C); 
    
        std::cout << "Naive multiplication resulting polynomial: "; 
        print_polynomial(naiveC); 

        std::cout << "\nFFT Time: " << endNaive << std::endl; 
        std::cout << "Naive Multiplication Time: " << endFFT << std::endl << std::endl; 
    }
    std::cout << (same_polynomial(C, naiveC) ? "FFT Success" : "FFT Failed!") << std::endl; 
}


/**
 * @brief Computes the Inverse Discrete Fourier Transform (IDFT) directly.
 * Uses the direct summation definition to compute the inverse transform. 
 * Includes the 1/N normalization factor. 
 * @note This is an O(N^2) operation and is intended for reference or small N only.
 * @param A The input frequency-domain vector.
 * @return complexVector The time-domain signal vector.
 */
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

/**
 * @brief Computes the Discrete Fourier Transform (DFT) directly.
 * Uses the direct summation definition to compute the transform.
 * @note This is an O(N^2) operation. For large N, use an FFT algorithm.
 * @param A The input time-domain signal vector.
 * @return complexVector The frequency-domain vector.
 */
complexVector DFT_Compute(const complexVector& A) 
{
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