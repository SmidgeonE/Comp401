#ifndef BOIDSIM_H
#define BOIDSIM_H

#include <vector>
#include "array"

constexpr long SIZE_OF_SIMULATION = 1000;

class VectorArray {
private:
    std::array<double, SIZE_OF_SIMULATION>* arrayX;
    std::array<double, SIZE_OF_SIMULATION>* arrayY;
    std::array<double, SIZE_OF_SIMULATION>* arrayZ;

public:
    VectorArray();
    ~VectorArray();

    std::array<double, SIZE_OF_SIMULATION>* getArrayX() const { return arrayX; };
    std::array<double, SIZE_OF_SIMULATION>* getArrayY() const { return arrayY; };
    std::array<double, SIZE_OF_SIMULATION>* getArrayZ() const { return arrayZ; };
};


class BoidSim {
private:
    VectorArray* BoidPositions;
    VectorArray* BoidDirections;
    std::array<double, SIZE_OF_SIMULATION>* BoidMasses;
    std::array<double, SIZE_OF_SIMULATION>* BoidSpeeds;

public:
    BoidSim();
    ~BoidSim();

    VectorArray* getBoidPositions() const { return BoidPositions; };
    VectorArray* getBoidDirections() const { return BoidDirections; };
    std::array<double, SIZE_OF_SIMULATION>* getBoidMasses() const { return BoidMasses; };
    std::array<double, SIZE_OF_SIMULATION>* getBoidSpeeds() const { return BoidSpeeds; };

    void SimView(const int viewNum);
};


void initialiseRandomVectors(VectorArray* vectorArray, const double lowerBound, const double upperBound, const bool normalise);

void initialiseRandomScalars(std::array<double, SIZE_OF_SIMULATION>* scalarArray, const double lowerBound, const double upperBound);

#endif // BOIDSIM_H