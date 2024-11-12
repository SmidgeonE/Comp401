cmake_minimum_required(VERSION 3.27.9)
project(Comp401)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)
include_directories(/software/local/languages/Intel-OneAPI)
include_directories(/software/spack/linux-rocky8-broadwell/gcc-12.3.0/gsl-2.7.1-w4io/include/gsl)
include_directories(/software/r22d_spack/linux-rocky8-broadwell/gcc-12.3.0/cuda-11.1.1-7pfi/targets/x86_64-linux/include/)
include_directories(/user/home/oz21652/glfw-3.4/include/GLFW)
add_executable(new_target BoidSim.cpp)
