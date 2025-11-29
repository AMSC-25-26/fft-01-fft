#ifndef UTILS_H
#define UTILS_H

unsigned int bitReverse(unsigned int x, int log2n) {
    int n = 0;
    for (int i = 0; i < log2n; i++) {
        n = n << 1;
        n = n | (x & 1);
        x = x >> 1;
    }
    return n;
}

template <typename T>
void applyBitReversalPermutation(std::vector<T>& input_vector) {
    int N = input_vector.size();
    
   
    int log2n = 0;
    if (N > 1) {
        log2n = std::log2(N);
    }
    

    // Iterate and Swap
    for (unsigned int i = 0; i < N; ++i) {
        // Calculate the bit-reversed index
        unsigned int j = bitReverse(i, log2n);
        
        // Only swap if the index j is greater than the current index i.
        
        if (i < j) {
            std::swap(input_vector[i], input_vector[j]);
        }
    }
}
#endif 