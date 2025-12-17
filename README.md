# Fast Fourier Transform – FFT

## Description

In this project, we implemented the **Fast Fourier Transform (FFT)** algorithm to efficiently compute the **Discrete Fourier Transform (DFT)** and its inverse, reducing the computational complexity from $O(n^2)$ to $O(n \log n)$.

We present:
- a **recursive FFT implementation**,
- an **iterative FFT implementation**,
- a **parallelized version** of the iterative algorithm, for both the forward and inverse transforms.

The project is developed in **C++**.

## Requirements 

To build and run the project, you need:
- a C++ compiler supporting **C++ 17**
- **CMake**
- a build system such as `make` or `ninja`

## Build Instructions

**1. Clone the repository** 

```bash
git clone https://github.com/AMSC-25-26/fft-01-fft.git
cd PATH/TO/fft-01-fft
```

**2. Create the build directory** 

```bash
mkdir build
cd build
```

**3. Configure and compile** 

```bash
cmake ..
make
```

After the compilation, the executable `fft` will be generated in the `build` directory.

## Execution

The executable requires **two integer arguments**:

```bash
./fft <n1> <n2>
```

The **input** parameters `n1` and `n2` represent the sizes of the two polynomials
used in the final test, which performs polynomial multiplication using FFT.

The **output** produced by the execution consists of the results of some tests performed by the program:
-  a fixed **8-Point FFT Test**
-  the **comparison** between FFT and DFT results
-  the **comparison** between the IFFT output signal and the original input points
-  a **FFT test** on a signal represented by a **mathematical function** compared with the **Analytical FT**
-  a **Polynomial multiplication test**, comparing the FFT-based multiplication with the naive one 

## Structure of the Repository 

```text
.
├── src/                  # Source files containing FFT implementations
│   ├── FFT.hpp           # Declaration of FFT algorithms
│   ├── FFT.cpp           # Implementation of recursive, iterative and parallel FFT
│   └── main.cpp          # Main program and test cases
├── Documentation/        # Project documentation
├── .github/              # GitHub configuration files
├── .vscode/              # VSCode configuration
├── CMakeLists.txt        # CMake build configuration
└── README.md             # Project description







