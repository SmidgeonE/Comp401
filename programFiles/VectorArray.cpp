#include "BoidSim.h"

VectorArray::VectorArray() : arrayX(new std::array<double, SIZE_OF_SIMULATION>),
                             arrayY(new std::array<double, SIZE_OF_SIMULATION>),
                             arrayZ(new std::array<double, SIZE_OF_SIMULATION>) {}
                            

VectorArray::~VectorArray() { delete arrayX; delete arrayY; delete arrayZ; }


std::array<double, 3> VectorArray::GetVectorAverage() {
    double averageX = 0;
    double averageY = 0;
    double averageZ = 0;

    for (int i = 0; i < SIZE_OF_SIMULATION; ++i){
        averageX += (*arrayX)[i];
        averageY += (*arrayY)[i];
        averageZ += (*arrayZ)[i];
    }

    averageX /= SIZE_OF_SIMULATION;
    averageY /= SIZE_OF_SIMULATION;
    averageZ /= SIZE_OF_SIMULATION;

    std::array<double, 3> averageDirection = {averageX, averageY, averageZ};

    return averageDirection;
}


void VectorArray::InitialiseRandomVectors(const double lowerBound, const double upperBound, const bool normalise) {
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


        GetArrayX()->at(i) = vector[0];
        GetArrayY()->at(i) = vector[1];
        GetArrayZ()->at(i) = vector[2];
    }
}

void VectorArray::InitialiseVectorsToLine(const int gridSize) {
    // Sets vectors to the line y=x=z

    auto spacing = gridSize / SIZE_OF_SIMULATION;

    for (double i = 0; i < gridSize; i += spacing) {
        GetArrayX()->at(i) = i;
        GetArrayY()->at(i) = i;
        GetArrayZ()->at(i) = i;
    }
}


