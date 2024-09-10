cmake_minimum_required(VERSION 3.12)

project(GSS VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mfma -static -fopenmp") #-fno-sanitize=all -static -ipp=common -qopenmp  -fsanitize=address  -fopenmp

# Add the source files for the executable

add_executable(GLS_b src/main_b.cpp)
#target_link_libraries(GLS_b PUBLIC ${CMAKE_SOURCE_DIR}/lib/libx86simdsortcpp.a)

# Find 32-bit FFTW package
target_link_libraries(GLS_b PUBLIC fftw3f)

# Install the executable to /usr/local/bin
install(TARGETS GLS_b DESTINATION /usr/local/bin/)
