//
// Created by oz21652 on 9/24/24.
//

#include "BoidSim.h"
#include <iostream>
#include <cmath>
#include <vector>
#include "mpi.h"
#include <gsl/gsl_vector.h>
#include "array"

constexpr int SIZE_OF_SIMULATION = 1000;

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
    BoidSim* boidSim = new BoidSim();

    
    return 0;
}


class BoidSim {
private:
    VectorArray* BoidPositions;
    VectorArray* BoidDirections;

    std::array<double, SIZE_OF_SIMULATION>* BoidMasses;
    std::array<double, SIZE_OF_SIMULATION>* BoidSpeeds;

public:
    BoidSim() {
        BoidPositions = new VectorArray;
        BoidDirections = new VectorArray;
        BoidMasses = new std::array<double, SIZE_OF_SIMULATION>;
        BoidSpeeds = new std::array<double, SIZE_OF_SIMULATION>;
    }

    ~BoidSim() {
        delete BoidPositions;
        delete BoidDirections;
        delete[] BoidMasses;
        delete[] BoidSpeeds;
    }

    VectorArray* getBoidPositions() const {
        return BoidPositions;
    }

    VectorArray* getBoidDirections() const {
        return BoidDirections;
    }

    std::array<double, SIZE_OF_SIMULATION>* getBoidMasses() const {
        return BoidMasses;
    }

    std::array<double, SIZE_OF_SIMULATION>* getBoidSpeeds() const {
        return BoidSpeeds;
    }
};


class VectorArray {
private:
    std::array<double, SIZE_OF_SIMULATION>* arrayX;
    std::array<double, SIZE_OF_SIMULATION>* arrayY;
    std::array<double, SIZE_OF_SIMULATION>* arrayZ;

public:
    VectorArray() {
        arrayX = new std::array<double, SIZE_OF_SIMULATION>;
        arrayY = new std::array<double, SIZE_OF_SIMULATION>;
        arrayZ = new std::array<double, SIZE_OF_SIMULATION>;
    }

    ~VectorArray() {
        delete[] arrayX;
        delete[] arrayY;
        delete[] arrayZ;
    }

    std::array<double, SIZE_OF_SIMULATION>* getArrayX() const {
        return arrayX;
    }

    std::array<double, SIZE_OF_SIMULATION>* getArrayY() const {
        return arrayY;
    }

    std::array<double, SIZE_OF_SIMULATION>* getArrayZ() const {
        return arrayZ;
    }
};
    
