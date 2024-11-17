#include "BoidSim.h"
#include <random>

void initialiseRandomScalars(std::array<double, SIZE_OF_SIMULATION>* scalarArray, const double lowerBound, const double upperBound) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(lowerBound, upperBound);

    for (int i = 0; i < SIZE_OF_SIMULATION; ++i) {
        (*scalarArray)[i] = dis(gen);
    }
}


std::array<double, 3> normaliseVector(const std::array<double, 3>& vector){
    double magnitude = std::sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
    return {vector[0]/magnitude, vector[1]/magnitude, vector[2]/magnitude};
}


double magSquared(const std::array<double, 3>& vector) {
    return vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2];
}


void initialiseRandomVectors(VectorArray* vectorArray, const double lowerBound, const double upperBound, const bool normalise) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(lowerBound, upperBound);

    for (int i = 0; i < SIZE_OF_SIMULATION; ++i) {
        const double x = dis(gen);
        const double y = dis(gen);
        const double z = dis(gen);

        auto vector = std::array<double, 3> {x, y, z};

        if (normalise) {
            vector = normaliseVector(vector);
        }


        (*vectorArray->getArrayX())[i] = vector[0];
        (*vectorArray->getArrayY())[i] = vector[1];
        (*vectorArray->getArrayZ())[i] = vector[2];
    }
}
    
