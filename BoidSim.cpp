//
// Created by oz21652 on 9/24/24.
//

#include <iostream>
#include <cmath>
#include "mpi.h"
#include <random>
#include "BoidSim.h"

int main(int argc, char* argv[]) {
    BoidSim* boidSim = new BoidSim();

    boidSim->SimView(5);

    return 0;
}



BoidSim::BoidSim() : BoidPositions(new VectorArray),
                     BoidDirections(new VectorArray),
                     BoidMasses(new std::array<double, SIZE_OF_SIMULATION>),
                     BoidSpeeds(new std::array<double, SIZE_OF_SIMULATION>) {

                        initialiseRandomVectors(BoidPositions, 0.0, 1.0, false);
                        initialiseRandomVectors(BoidDirections, -1.0, 1.0, true);
                        initialiseRandomScalars(BoidSpeeds, 0.0, 1.0);
                        BoidMasses->fill(1.0);
                     }
            
BoidSim::~BoidSim() { delete BoidPositions; delete BoidDirections; delete BoidMasses; delete BoidSpeeds; }

void BoidSim::SimView(const int viewNum) {
    auto boidSpeeds = this->BoidSpeeds;
    auto boidMasses = this->BoidMasses;

    auto BoidPositions = this->BoidPositions;
    auto boidDirections = this->BoidDirections;

    auto boidPositionsX = BoidPositions->getArrayX();
    auto boidPositionsY = BoidPositions->getArrayY();
    auto boidPositionsZ = BoidPositions->getArrayZ();

    auto boidDirectionsX = boidDirections->getArrayX();
    auto boidDirectionsY = boidDirections->getArrayY();
    auto boidDirectionsZ = boidDirections->getArrayZ();


    for (int i = 0; i < viewNum; ++i){
        std::cout << i << "th Boid:" 
        << "\n   -- position: " << (*boidPositionsX)[i] << ", " << (*boidPositionsY)[i] << ", " << (*boidPositionsZ)[i]
        << "\n   -- direction: " << (*boidDirectionsX)[i] << ", " << (*boidDirectionsY)[i] << ", " << (*boidDirectionsZ)[i]
        << "\n   -- direction magnitude: " << std::sqrt((*boidDirectionsX)[i]*(*boidDirectionsX)[i] + (*boidDirectionsY)[i]*(*boidDirectionsY)[i] + (*boidDirectionsZ)[i]*(*boidDirectionsZ)[i])
        << "\n   -- speed: " << (*boidSpeeds)[i]
        << "\n   -- mass: " << (*boidMasses)[i] << "\n";
    }
}


VectorArray::VectorArray() : arrayX(new std::array<double, SIZE_OF_SIMULATION>),
                             arrayY(new std::array<double, SIZE_OF_SIMULATION>),
                             arrayZ(new std::array<double, SIZE_OF_SIMULATION>) {}
                            
VectorArray::~VectorArray() { delete arrayX; delete arrayY; delete arrayZ; }


void initialiseRandomVectors(VectorArray* vectorArray, const double lowerBound, const double upperBound, const bool normalise) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(lowerBound, upperBound);

    for (int i = 0; i < SIZE_OF_SIMULATION; ++i) {
        double x = dis(gen);
        double y = dis(gen);
        double z = dis(gen);

        double normalisationFactor = 1;
    
        if (normalise) {
            normalisationFactor = 1 / std::sqrt((x*x + y*y + z*z));
        }

        (*vectorArray->getArrayX())[i] = x * normalisationFactor;
        (*vectorArray->getArrayY())[i] = y * normalisationFactor;
        (*vectorArray->getArrayZ())[i] = z * normalisationFactor;
    }
}


void initialiseRandomScalars(std::array<double, SIZE_OF_SIMULATION>* scalarArray, const double lowerBound, const double upperBound) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(lowerBound, upperBound);

    for (int i = 0; i < SIZE_OF_SIMULATION; ++i) {
        (*scalarArray)[i] = dis(gen);
    }
}

    
