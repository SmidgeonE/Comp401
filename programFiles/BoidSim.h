#ifndef BOIDSIM_H
#define BOIDSIM_H

#include "array"
#include <iostream>
#include <fstream>

constexpr long SIZE_OF_SIMULATION = 20;

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

    bool writeToFile;
    std::ofstream outputStream;
    void WriteBoidSimulation();

    std::array<double, 3> getAverageFlockDirection();
    
    void applyAlignmentAlgorithm();
    void applyCohesionAlgorithm();
    void applySeparationAlgorithm();
    void applyTimeStep();

public:
    BoidSim();
    ~BoidSim();

    void setWriteToFile(const bool writeToFile);
    void SimView(int viewNum) const;
    void StartSimulation(long timeSteps);
};


void initialiseRandomVectors(VectorArray* vectorArray, double lowerBound, double upperBound, bool normalise);

void initialiseRandomScalars(std::array<double, SIZE_OF_SIMULATION>* scalarArray, double lowerBound, double upperBound);

std::array<double, 3> normaliseVector(const std::array<double, 3>& vector);

#endif // BOIDSIM_H