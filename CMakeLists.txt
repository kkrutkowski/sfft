cmake_minimum_required(VERSION 3.12)

project(GSS VERSION 1.0.0 LANGUAGES CXX)

# W tym przypadku icx/icpx dają sporo szybszy kod (~ 15% różnicy) od gcc/g++
#set(CMAKE_C_COMPILER "icx")
#set(CMAKE_CXX_COMPILER "icpx")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mfma -static -fopenmp -fno-sanitize=all") # -fsanitize=address -fnosanitize=all
# Add the source files for the executable

add_executable(sfft_demo src/main_b.cpp)

# Find FFTW package
target_link_libraries(sfft_demo PUBLIC fftw3f)

# Install the executable to /usr/local/bin
install(TARGETS sfft_demo DESTINATION /usr/local/bin/)
