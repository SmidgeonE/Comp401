#ifndef BOIDSIM_H
#define BOIDSIM_H

#include "array"
#include <iostream>
#include <fstream>
#include <omp.h>
#include <cmath>
#include <random>


#ifndef SIZE_OF_SIMULATION
#   define SIZE_OF_SIMULATION 5
#endif

#ifndef NUM_THREADS
#   define NUM_THREADS -1
#endif

#ifndef DEBUG
#   define DEBUG false
#endif

constexpr double DT = 0.01f;

constexpr double SEPARATION_FORCE_CONSTANT = 0.00000001;
constexpr double ALIGNMENT_FORCE_CONSTANT = 0.00000001;
constexpr double COHESION_FORCE_CONSTANT = 0.00000001;

constexpr int NUM_SIMULATIONS = 6;

class VectorArray {
private:
    std::array<double, SIZE_OF_SIMULATION>* arrayX;
    std::array<double, SIZE_OF_SIMULATION>* arrayY;
    std::array<double, SIZE_OF_SIMULATION>* arrayZ;

public:
    VectorArray();
    ~VectorArray();

    std::array<double, SIZE_OF_SIMULATION>* GetArrayX() const { return arrayX; };
    std::array<double, SIZE_OF_SIMULATION>* GetArrayY() const { return arrayY; };
    std::array<double, SIZE_OF_SIMULATION>* GetArrayZ() const { return arrayZ; };

    std::array<double, 3> GetVectorAverage();
    void InitialiseRandomVectors(double lowerBound, double upperBound, bool normalise);
    void InitialiseVectorsToLine(int lineLength);
};


class BoidSim {
private:
    VectorArray* boidPositions;
    VectorArray* boidDirections;
    VectorArray* boidForces;

    std::array<double, SIZE_OF_SIMULATION>* boidMasses;
    std::array<double, SIZE_OF_SIMULATION>* boidSpeeds;

    bool writeToFile;
    std::ofstream outputStream;
    void writeBoidSimulation();

    std::array<double, 3> getAverageFlockDirection();

    void calculateBoidVelocity();
    void resetForces();
    void addForce(const VectorArray& force);
    
    void applyAlignmentForce();
    void applyCohesionForce();
    void applySeparationForce();
    void applyTimeStep();

public:
    BoidSim();
    ~BoidSim();

    void SetWriteToFile(const bool writeToFile);
    void SimView(int viewNum) const;
    void StartSimulation(long timeSteps);
};



void initialiseRandomScalars(std::array<double, SIZE_OF_SIMULATION>* scalarArray, double lowerBound, double upperBound);

std::array<double, 3> normaliseVector(const std::array<double, 3>& vector);

double magSquared(const std::array<double, 3>& vector);

double runAndTimeSimulation(int timeSteps, bool writeToFile);

#endif // BOIDSIM_H