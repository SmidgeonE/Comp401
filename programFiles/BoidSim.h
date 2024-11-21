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

#ifndef WRITE_SIM
#   define WRITE_SIM false
#endif

#define CELL_NUMBER 10

#ifndef DO_MULTIPLE_SIMS
#   define DO_MULTIPLE_SIMS false
#endif


constexpr double DT = 0.01f;

constexpr double SEPARATION_FORCE_CONSTANT = 1;
constexpr double ALIGNMENT_FORCE_CONSTANT = 1;
constexpr double COHESION_FORCE_CONSTANT = 1;


constexpr int NUM_SIMULATIONS = 6;

class VectorArray {
private:
    std::array<double, SIZE_OF_SIMULATION> arrayX;
    std::array<double, SIZE_OF_SIMULATION> arrayY;
    std::array<double, SIZE_OF_SIMULATION> arrayZ;

public:
    VectorArray();

    std::array<double, SIZE_OF_SIMULATION>& GetArrayX() { return arrayX; };
    std::array<double, SIZE_OF_SIMULATION>& GetArrayY() { return arrayY; };
    std::array<double, SIZE_OF_SIMULATION>& GetArrayZ() { return arrayZ; };

    std::array<double, 3> GetVectorAverage();

    void InitialiseRandomVectors(double lowerBound, double upperBound, bool normalise);
    void InitialiseVectorsToLine(const double lineLength);
    void InitaliseVectorsToZHat();

    void View(const int viewNum, const std::string& name);
};


class BoidSim {
private:
    VectorArray* boidPositions;
    VectorArray* boidDirections;
    VectorArray* boidForces;

    std::array<double, SIZE_OF_SIMULATION> boidMasses;
    std::array<double, SIZE_OF_SIMULATION> boidSpeeds;

    std::array<std::array<std::array<std::vector<int>, CELL_NUMBER>, CELL_NUMBER>, CELL_NUMBER> cellList;

    std::array<int, SIZE_OF_SIMULATION> filledCellsX;
    std::array<int, SIZE_OF_SIMULATION> filledCellsY;
    std::array<int, SIZE_OF_SIMULATION> filledCellsZ;

    std::array<double, 3> cellMinima;
    std::array<double, 3> cellMaxima;

    bool writeToFile;
    std::ofstream outputStream;
    void writeBoidSimulation();

    double startTime;
    double sepTime;
    double alignTime;
    double cohTime;
    double velTime;

    std::array<double, 3> getAverageFlockDirection();

    void calculateBoidVelocity();
    void resetForces();
    void addForce(VectorArray& force);
    
    void applyAlignmentForce();
    void applyCohesionForce();
    void applySeparationForce();
    void applySeparationForceCellList();
    void applyTimeStep();

    void constructCellList();
    void wipeCellList();

    std::vector<int> getAdjacentBoids(const int boidIndex);
    std::array<int, 3> getBoidCell(const int boidIndex);

public:
    BoidSim();
    ~BoidSim();

    void SetWriteToFile(const bool writeToFile);
    void SimView(int viewNum) const;
    void StartSimulation(long timeSteps);
};



void initialiseRandomScalars(std::array<double, SIZE_OF_SIMULATION>& scalarArray, double lowerBound, double upperBound);

std::array<double, 3> normaliseVector(const std::array<double, 3>& vector);

double magSquared(const std::array<double, 3>& vector);

double runAndTimeSimulation(int timeSteps, bool writeToFile);

std::array<double, 2> minMaxOfArray(std::array<double, SIZE_OF_SIMULATION>& array);

#endif // BOIDSIM_H