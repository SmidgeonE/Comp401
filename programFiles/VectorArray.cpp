#include "BoidSim.h"

VectorArray::VectorArray() : arrayX(new std::array<double, SIZE_OF_SIMULATION>),
                             arrayY(new std::array<double, SIZE_OF_SIMULATION>),
                             arrayZ(new std::array<double, SIZE_OF_SIMULATION>) {}
                            

VectorArray::~VectorArray() { delete arrayX; delete arrayY; delete arrayZ; }


std::array<double, 3> VectorArray::getVectorAverage() {
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

