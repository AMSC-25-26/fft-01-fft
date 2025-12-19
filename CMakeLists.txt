cmake_minimum_required(VERSION 3.10)

project(FastFourierTransform-Handson)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(OpenMP REQUIRED)

add_executable(fft src/main.cpp src/FFT.cpp)

target_link_libraries(fft PUBLIC OpenMP::OpenMP_CXX)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(fft PRIVATE -Wall -Wextra -O3)
endif()