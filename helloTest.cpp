//
// Created by oz21652 on 9/22/24.
//

#include "helloTest.h"

#include <chrono>
#include <cmath>
#include <iostream>
#include <omp.h>

constexpr int iterations = 10000000;

double singleThreadedCalc() {
    double piByFour = 0;
    constexpr double dx = 1.0 / iterations;

    for (int i = 0; i < iterations; i++) {
        piByFour += sqrt(1-pow(i*dx, 2));
    }

    const double pi = piByFour * 4 * dx;

    return pi;
}

double multiThreadedCalc() {
    double piByFour = 0;
    constexpr double dx = 1.0 / iterations;
    int i;

    #pragma omp parallel private(i)
    {
        #pragma omp for reduction(+:piByFour)
        {
            for (i = 0; i < iterations; i++) {
                piByFour += sqrt(1-pow(i*dx, 2));
            }
        }
    }

    const double pi = piByFour * 4 * dx;

    return pi;
}

int main(const int argc, char* argv[]) {

    if (argc != 2) {
        std::cout << "Usage : <num threads>" << std::endl;
        return 0;
    }

    char* end;
    const long numThreads = strtol(argv[1], &end, 10);
    omp_set_num_threads(numThreads);

    std::cout << "This is the pi calculator. \n";

    auto startTime = omp_get_wtime();
    double pi = singleThreadedCalc();
    auto endTime = omp_get_wtime();

    std::cout << "The answer given by single threaded for pi was " << pi << std::endl;
    std::cout << "The single threaded application took " << endTime - startTime << "s \n" << std::endl;

    startTime = omp_get_wtime();
    pi = multiThreadedCalc();
    endTime = omp_get_wtime();

    std::cout << "The answer given by multi threaded for pi was " << pi << std::endl;
    std::cout << "The multi threaded application took " << endTime - startTime << "s \n" << std::endl;

    return 0;
}
