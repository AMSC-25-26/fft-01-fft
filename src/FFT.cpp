#include "FFT.hpp"
#include "utils.h"
#include <omp.h>

using namespace std::complex_literals;

namespace FFT
{

    /**
     * Recursive version of the Cooley-Tukey FFT Algorithm (Void / In-Place style)
     * @param A The input/output vector. 
     * Input: sampled values (time domain).
     * Output: DFT coefficients (frequency domain).
     */
    void recursive(std::vector<complex> &A) {
        int n = A.size();
        if (n <= 1) return;

        std::vector<complex> A_even;
        A_even.reserve(n/2);
        std::vector<complex> A_odd;
        A_odd.reserve(n/2);

        for (int i = 0; i < n / 2; i++) {
            A_even.emplace_back(A[2 * i]);
            A_odd.emplace_back(A[2 * i + 1]);
        }

        recursive(A_even);
        recursive(A_odd);

        complex Wn = std::exp((-2.0 * M_PI * 1i) / static_cast<double>(n));
        complex W = 1.0;

        for (int j = 0; j < n / 2; j++) {
            complex u = A_even[j];
            complex v = W * A_odd[j];
            
            A[j] = u + v;
            A[j + n / 2] = u - v;
            W *= Wn;
        }
    }

    /**
     * Iterative version of the Cooley-Tukey FFT Algorithm.
     * @param A The input vector which represents the sampled 
     *          values of an arbitrary signal (time domain). 
     * 
     * @return  The vector of coefficients yielded from the 
     *          Discrete Fourier Transform of the input signal (frequency domain). 
     */
    void iterative(std::vector<std::complex<double>>& A) {
        const size_t N = A.size();
        if (N == 0 || (N & (N - 1)) != 0) {
            throw std::invalid_argument("Input vector size must be a power of two.");
        }

        /**
         * Here we precompute the angle needed to
         * calculate the twiddle factors.
         */
        constexpr double theta = -2 * M_PI;

        /**
         * Cooley-Tukey is a divide and conquer algorithm
         * that splits the input vector in half at each
         * iteration resulting in log2(N) iterations.
         */
        const int log2N = static_cast<int>(std::log2(N));

        /**
         * Bit reversal permutation. Since the array
         * needs to holds the elements in the same 
         * order in wich they appear in the leaves of 
         * the tree.
         */
        applyBitReversalPermutation(A);

        /**
         * We iterate over the stages of the FFT computation.
         * In each stage we combine the results of smaller
         * DFTs to compute larger DFTs up to the N-point DFT.
         */
        for (int k = 1; k <= log2N; ++k) { 
            /**
             * The size of the subproblems in the 
             * current stage of the FFT. We use <<
             * (bitwise shift operator) to quickly
             * compute the power of 2. This operation
             * is faster than using pow(2, k).
             */
            const int m = 1 << k;

            /**
             * The midpoint value, that is half the size
             * of the subproblems in the current stage of 
             * the FFT. The value is essential for combining
             * the elements using the butterfly operation.
             * Again we use the bitwise shift operator to 
             * quickly compute the division by 2. 
             */
            const size_t mid = m >> 1;

            /**
             * Precompute the twiddle factor. This avoids
             * calling to std::exp repeatedly in the inner loop.
             */
            std::complex<double> wm = std::polar(1.0, theta / static_cast<double>(m));
            /**
             * The inner loop for computing the sub DFTs
             * in the current stage. Note that each sub DFTs
             * is computed indipendently from the others.
             * We iterate over the entire array with a step that
             * is the size of the subproblems.
             */
            for (size_t r = 0; r < N; r += m) {
                std::complex<double> w{1.0, 0};
                /**
                 * The loop performs the butterfly computation
                 * on two elements within the sub DFTs.
                 */
                for (size_t j = 0; j < mid; ++j) {    
                    /**
                     * Here we perform the so called phase
                     * adjustment, that is we rotate the phase
                     * of this element by multiplying it by 
                     * the twiddle factor.
                     */
                    std::complex<double> T = w * A[r + j + mid];
                    std::complex<double> X_top_old = A[r + j];
                    
                
                    /**
                     * The butterfly operation itself.
                     */
                    A[r + j + mid] = X_top_old - T;
                    A[r + j] = X_top_old + T;
                    /**
                     * Update of the twiddle factor.
                     * Constructing (e^-(i*2pi / m))^j
                     */
                    w *= wm;
                }
            }         
        }
    }

    /**
     * Iterative version of the Cooley-Tukey Inverse FFT Algorithm.
     * @param A The input vector which represents the DFT coefficients 
     *          of an arbitrary signal (frequency domain). 
     * 
     * @return  The vector of coefficients yielded from the 
     *          inverse FFT of the input signal (time domain). 
     */
    void inverse(std::vector<std::complex<double>>& A) {
        const size_t N = A.size();
        if (N == 0 || (N & (N - 1)) != 0) {
            throw std::invalid_argument("Input vector size must be a power of two.");
        }

        constexpr double theta = 2 * M_PI;
        const int log2N = static_cast<int>(std::log2(N));

        applyBitReversalPermutation(A);

        for (int k = 1; k <= log2N; ++k) { 
            const int m = 1 << k; 
            const size_t m2 = m >> 1; 

            std::complex<double> wm = std::polar(1.0, theta / static_cast<double> (m)); 

            // We iterate over the blocks
            for (size_t r = 0; r < N; r += m) {
                std::complex<double> w{1.0, 0};
                for (size_t j = 0; j < m2; ++j) {    
                    std::complex<double> X_top_old = A[r + j];
                    std::complex<double> T = w * A[r + j + m2];
                
                    A[r + j + m2] = X_top_old - T;
                    A[r + j] = X_top_old + T;
                    w *= wm;
                }
            }         
        }
        
        for (size_t i = 0; i < N; ++i) {
            A[i] /= N;
        }
    }


    /**
     * Parallel implementation of the Cooley-Tukey FFT Algorithm.
     * @param A The input vector which represents the sampled 
     *          values of an arbitrary signal (time domain). 
     */
    void parallel_iterative(std::vector<std::complex<double>>& A) {
        const size_t N = A.size();
        if (N == 0 || (N & (N - 1)) != 0) {
            throw std::invalid_argument("Input vector size must be a power of two.");
        }

        constexpr double theta = -2 * M_PI;
        const int log2N = static_cast<int>(std::log2(N));


        applyBitReversalPermutation(A);


        /**
        * Parallel region is defined before the outer loop to 
        * avoid the overhead of creating and destroying threads
        * at each iteration of the outer loop.
         */
        #pragma omp parallel
        {
        for (int k = 1; k <= log2N; ++k) { 
            const int m = 1 << k; 
            const size_t mid = m >> 1;

            std::complex<double> wm = std::polar(1.0, theta / static_cast<double> (m));

            /**
             * The array is divided into subproblems
             * of size m. Since each block is indipendent
             * and doesen't overlap with the others, different
             * threads can process them concurrently.
             */
            #pragma omp for schedule(static)
            for (size_t r = 0; r < N; r += m) {
                std::complex<double> w{1.0, 0};
                for (size_t j = 0; j < mid; ++j) {    
                    std::complex<double> X_top_old = A[r + j];
                    std::complex<double> T = w * A[r + j + mid];
                
                    A[r + j + mid] = X_top_old - T;
                    A[r + j] = X_top_old + T;
                    w *= wm;
                }
            }         
        }
    }
    }

    /**
     * Parallel implementation of the Cooley-Tukey Inverse FFT Algorithm.
     * @param A The input vector which represents the DFT coefficients 
     *          of an arbitrary signal (frequency domain). 
     */
    void parallel_inverse(std::vector<std::complex<double>>& A) {
        const size_t N = A.size();
        if (N == 0 || (N & (N - 1)) != 0) {
            throw std::invalid_argument("Input vector size must be a power of two.");
        }

        constexpr double theta = 2 * M_PI;
        const int log2N = static_cast<int>(std::log2(N));

        applyBitReversalPermutation(A);

        #pragma omp parallel
        {
        for (int k = 1; k <= log2N; ++k) { 
            const int m = 1 << k;
            const size_t mid = m >> 1;

            std::complex<double> wm = std::polar(1.0, theta / static_cast<double> (m));

            #pragma omp for schedule(static)
            for (size_t r = 0; r < N; r += m) {
                std::complex<double> w{1.0, 0};
                for (size_t j = 0; j < mid; ++j) {    
                    std::complex<double> X_top_old = A[r + j];
                    std::complex<double> T = w * A[r + j + mid];
                
                    A[r + j + mid] = X_top_old - T;
                    A[r + j] = X_top_old + T;
                    w *= wm;
                }
            }         
        }
        /**
         * Parallelization of the normalization step.
         * we  use nowait to avoid te implicit barrier
         * since this is the last operation in the parallel region.
         */
        #pragma omp for schedule(static) nowait
        for (size_t i = 0; i < N; ++i) {
            A[i] /= N;
        }
    }
}
}