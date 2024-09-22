cmake_minimum_required(VERSION 3.30)
project(newProj)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(newProj
    helloTest.cpp
    helloTest.h)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(MyTarget PUBLIC OpenMP::OpenMP_CXX)
endif()
