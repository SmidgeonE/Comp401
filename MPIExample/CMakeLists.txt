cmake_minimum_required(VERSION 3.29)
project(Comp401)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)
include_directories(/software/local/languages/Intel-OneAPI)
add_executable(Comp401 multithreaded-pi-calc/picalc.cpp
        MPIAndQueueing.cpp
        multithreaded-pi-calc/picalc.h)