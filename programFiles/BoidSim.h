#ifndef BOIDSIM_H
#define BOIDSIM_H

#include "array"
#include <iostream>
#include <fstream>

#ifndef SIZE_OF_SIMULATION
#   define SIZE_OF_SIMULATION 5
#endif

constexpr double dt = 0.01f;

constexpr double SEPARATION_FORCE_CONSTANT = 0.00000001;
constexpr double ALIGNMENT_FORCE_CONSTANT = 0.00000001;
constexpr double COHESION_FORCE_CONSTANT = 0.00000001;

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

    std::array<double, 3> getVectorAverage();
};


class BoidSim {
private:
    VectorArray* BoidPositions;
    VectorArray* BoidDirections;
    VectorArray* BoidForces;

    std::array<double, SIZE_OF_SIMULATION>* BoidMasses;
    std::array<double, SIZE_OF_SIMULATION>* BoidSpeeds;

    bool writeToFile;
    std::ofstream outputStream;
    void WriteBoidSimulation();

    std::array<double, 3> getAverageFlockDirection();

    void calculateBoidVelocity();
    void resetForces();
    void addForce(const VectorArray& force);
    
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

double magSquared(const std::array<double, 3>& vector);

double runAndTimeSimulation(int timeSteps, bool writeToFile);

#endif // BOIDSIM_H