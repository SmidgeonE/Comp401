//
// Created by oz21652 on 9/24/24.
//

#include <cmath>
#include "mpi.h"
#include <random>
#include "BoidSim.h"

int main(int argc, char* argv[]) {
    auto* boidSim = new BoidSim();
    boidSim->setWriteToFile(true);
    // boidSim->SimView(5);

    boidSim->StartSimulation(20);

    delete boidSim;

    return 0;
}



BoidSim::BoidSim() : BoidPositions(new VectorArray),
                     BoidDirections(new VectorArray),
                     BoidMasses(new std::array<double, SIZE_OF_SIMULATION>),
                     BoidSpeeds(new std::array<double, SIZE_OF_SIMULATION>), writeToFile(false) {
    initialiseRandomVectors(BoidPositions, 0.0, 1.0, false);
    initialiseRandomVectors(BoidDirections, -1.0, 1.0, true);
    initialiseRandomScalars(BoidSpeeds, 0.0, 1.0);
    BoidMasses->fill(1.0);
}

void BoidSim::StartSimulation(const long timeSteps) {
    auto xPositions = BoidPositions->getArrayX();
    auto yPositions = BoidPositions->getArrayY();
    auto zPositions = BoidPositions->getArrayZ();

    auto xDirections = BoidDirections->getArrayX();
    auto yDirections = BoidDirections->getArrayY();
    auto zDirections = BoidDirections->getArrayZ();

    auto boidSpeeds = BoidSpeeds;

    if (writeToFile) {
        // We need to write the initial state of the simulation

        WriteBoidSimulation();
    }

    for (long i = 0; i < timeSteps; ++i) {
        // Calculate direction, change in speed for this time step here:

        // TO DO: Currently just leaving the direction and speed the same


        // Then we must apply the physical movements from this time step

        for (int j = 0; j < SIZE_OF_SIMULATION; ++j) {
            (*xPositions)[j] += (*xDirections)[j] * (*boidSpeeds)[j];
            (*yPositions)[j] += (*yDirections)[j] * (*boidSpeeds)[j];
            (*zPositions)[j] += (*zDirections)[j] * (*boidSpeeds)[j];
        }

        if (writeToFile) {
            WriteBoidSimulation();
        }
    }
}

void BoidSim::WriteBoidSimulation() {
    std::string outputStr = "";

    for (int i = 0; i < SIZE_OF_SIMULATION; ++i) {
        outputStr += "[[" + std::to_string((*BoidPositions->getArrayX())[i]) + "," 
            + std::to_string((*BoidPositions->getArrayY())[i]) + "," 
            + std::to_string((*BoidPositions->getArrayZ())[i]) + "],["
            + std::to_string((*BoidDirections->getArrayX())[i]) + ","
            + std::to_string((*BoidDirections->getArrayY())[i]) + ","
            + std::to_string((*BoidDirections->getArrayZ())[i]) + "],"
            + std::to_string((*BoidSpeeds)[i]) + ","
            + std::to_string((*BoidMasses)[i]) + "],";
    }

    outputStream << outputStr << "\n";
}

void BoidSim::setWriteToFile(const bool writeToFile) { 
    this->writeToFile = writeToFile;

    // Opening a stream to write the file

    outputStream.open("BoidSimulation.csv");
}

BoidSim::~BoidSim() { 
    delete BoidPositions; 
    delete BoidDirections; 
    delete BoidMasses; 
    delete BoidSpeeds; 

    if (writeToFile) {
        outputStream.close();
    }
}

void BoidSim::SimView(const int viewNum) const {
    const auto boidSpeeds = this->BoidSpeeds;
    const auto boidMasses = this->BoidMasses;

    const auto BoidPositions = this->BoidPositions;
    const auto boidDirections = this->BoidDirections;

    const auto boidPositionsX = BoidPositions->getArrayX();
    const auto boidPositionsY = BoidPositions->getArrayY();
    const auto boidPositionsZ = BoidPositions->getArrayZ();

    const auto boidDirectionsX = boidDirections->getArrayX();
    const auto boidDirectionsY = boidDirections->getArrayY();
    const auto boidDirectionsZ = boidDirections->getArrayZ();


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
        const double x = dis(gen);
        const double y = dis(gen);
        const double z = dis(gen);

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

    
