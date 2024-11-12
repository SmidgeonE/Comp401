//
// Created by oz21652 on 9/24/24.
//

#include <iostream>
#include <cmath>
#include "mpi.h"
#include <random>
#include <BoidSim.h>

int main(int argc, char* argv[]) {
    BoidSim* boidSim = new BoidSim();

    auto boidSpeeds = boidSim->getBoidSpeeds();
    
    return 0;
}



BoidSim::BoidSim() : BoidPositions(new VectorArray),
                     BoidDirections(new VectorArray),
                     BoidMasses(new std::array<double, SIZE_OF_SIMULATION>),
                     BoidSpeeds(new std::array<double, SIZE_OF_SIMULATION>) {}
            
BoidSim::~BoidSim() { delete BoidPositions; delete BoidDirections; delete BoidMasses; delete BoidSpeeds; }


VectorArray::VectorArray() : arrayX(new std::array<double, SIZE_OF_SIMULATION>),
                             arrayY(new std::array<double, SIZE_OF_SIMULATION>),
                             arrayZ(new std::array<double, SIZE_OF_SIMULATION>) {}
                            
VectorArray::~VectorArray() { delete arrayX; delete arrayY; delete arrayZ; }


void initializeRandomVectors(VectorArray* vectorArray, double lowerBound, double upperBound) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(lowerBound, upperBound);


    for (int i = 0; i < SIZE_OF_SIMULATION; ++i) {
        (*vectorArray->getArrayX())[i] = dis(gen);
        (*vectorArray->getArrayY())[i] = dis(gen);
        (*vectorArray->getArrayZ())[i] = dis(gen);
    }
}

    
