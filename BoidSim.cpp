//
// Created by oz21652 on 9/24/24.
//

#include "BoidSim.h"
#include <iostream>
#include <cmath>
#include <vector>
#include "mpi.h"
#include <gsl/gsl_vector.h>

std::vector<gsl_vector*> initializeBoidsToSquare(const int n, const double sideLength) {
    std::vector<gsl_vector*> boids; 
    double separation = sideLength / std::cbrt(n);
    int count = static_cast<int>(std::floor(sideLength / separation));
    std::cout << count << std::endl;
    
    for (int x = 0; x < count; ++x) {
        for (int y = 0; y < count; ++y) {
            for (int z = 0; z < count; ++z) {
                gsl_vector* boid = gsl_vector_alloc(3);
                gsl_vector_set(boid, 0, x * separation);
                gsl_vector_set(boid, 1, y * separation);
                gsl_vector_set(boid, 2, z * separation);
                boids.push_back(boid);
            }
        }
    }

    return boids;
}

int main(int argc, char* argv[]) {

    auto budgieBoids = initializeBoidsToSquare(105, 5);

    for (auto budgieBoid : budgieBoids) {
        std::cout << gsl_vector_get(budgieBoid, 0) << " " << gsl_vector_get(budgieBoid, 1) << " " << gsl_vector_get(budgieBoid, 2) << std::endl;
    }

    return 0;
}
