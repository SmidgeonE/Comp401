cmake_minimum_required(VERSION 3.30)
project(Comp401)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(PUBLIC OpenMP::OpenMP_CXX)
endif()
