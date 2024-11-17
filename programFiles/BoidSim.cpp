//
// Created by oz21652 on 9/24/24.
//

#include <cmath>
#include "mpi.h"
#include "BoidSim.h"
#include <random>


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

    if (writeToFile) {
        // We need to write the initial state of the simulation

        WriteBoidSimulation();
    }

    // Now we must iterate through the time steps

    for (long i = 0; i < timeSteps; ++i) {
        applySeparationAlgorithm();
        applyAlignmentAlgorithm();
        applyCohesionAlgorithm();

        applyTimeStep();

        if (writeToFile) {
            WriteBoidSimulation();
        }
    }
}


void BoidSim::WriteBoidSimulation() {
    std::string outputStr = "";

    for (int i = 0; i < SIZE_OF_SIMULATION; ++i) {
        outputStr += "[" + std::to_string((*BoidPositions->getArrayX())[i]) + "," 
            + std::to_string((*BoidPositions->getArrayY())[i]) + "," 
            + std::to_string((*BoidPositions->getArrayZ())[i]) + ","
            + std::to_string((*BoidDirections->getArrayX())[i]) + ","
            + std::to_string((*BoidDirections->getArrayY())[i]) + ","
            + std::to_string((*BoidDirections->getArrayZ())[i]) + ","
            + std::to_string((*BoidSpeeds)[i]) + ","
            + std::to_string((*BoidMasses)[i]) + "],";
    }

    outputStream << outputStr << "\n";
}


void BoidSim::setWriteToFile(const bool writeToFile) { 
    this->writeToFile = writeToFile;

    // Opening a stream to write the file

    if (writeToFile){
        outputStream.open("BoidSimulation.csv");
    }
}


void BoidSim::SimView(const int viewNum) const {
    const auto boidPositionsX = BoidPositions->getArrayX();
    const auto boidPositionsY = BoidPositions->getArrayY();
    const auto boidPositionsZ = BoidPositions->getArrayZ();

    const auto boidDirectionsX = BoidDirections->getArrayX();
    const auto boidDirectionsY = BoidDirections->getArrayY();
    const auto boidDirectionsZ = BoidDirections->getArrayZ();


    for (int i = 0; i < viewNum; ++i){
        std::cout << i << "th Boid:" 
        << "\n   -- position: " << (*boidPositionsX)[i] << ", " << (*boidPositionsY)[i] << ", " << (*boidPositionsZ)[i]
        << "\n   -- direction: " << (*boidDirectionsX)[i] << ", " << (*boidDirectionsY)[i] << ", " << (*boidDirectionsZ)[i]
        << "\n   -- direction magnitude: " << std::sqrt((*boidDirectionsX)[i]*(*boidDirectionsX)[i] + (*boidDirectionsY)[i]*(*boidDirectionsY)[i] + (*boidDirectionsZ)[i]*(*boidDirectionsZ)[i])
        << "\n   -- speed: " << (*this->BoidSpeeds)[i]
        << "\n   -- mass: " << (*this->BoidMasses)[i] << "\n";
    }
}


void BoidSim::applyCohesionAlgorithm(){

}


void BoidSim::applySeparationAlgorithm(){

}


void BoidSim::applyAlignmentAlgorithm(){
    const auto directionsX = BoidDirections->getArrayX();
    const auto directionsY = BoidDirections->getArrayY();
    const auto directionsZ = BoidDirections->getArrayZ();

    const auto averageDirection = getAverageFlockDirection();

    for (int i = 0; i < SIZE_OF_SIMULATION; ++i){
        auto newDirectionX = ((*directionsX)[i] + 1 * averageDirection[0])/2;
        auto newDirectionY = ((*directionsY)[i] + 1 * averageDirection[1])/2;
        auto newDirectionZ = ((*directionsZ)[i] + 1 * averageDirection[2])/2;

        auto newDirection = normaliseVector({newDirectionX, newDirectionY, newDirectionZ});

        (*directionsX)[i] = newDirection[0];
        (*directionsY)[i] = newDirection[1];
        (*directionsZ)[i] = newDirection[2];


        if (i == 0){
            std::cout << "First boids Direction: "
             << (*BoidDirections->getArrayX())[0] 
             << ", " << (*BoidDirections->getArrayY())[0] 
             << ", " << (*BoidDirections->getArrayZ())[0] 
             << "\n";

            std::cout << "Average Direction: " 
            << averageDirection[0] 
            << ", " << averageDirection[1] 
            << ", " << averageDirection[2]
             << "\n";
        }
    }
}


void BoidSim::applyTimeStep(){
    auto xPositions = BoidPositions->getArrayX();
    auto yPositions = BoidPositions->getArrayY();
    auto zPositions = BoidPositions->getArrayZ();

    auto xDirections = BoidDirections->getArrayX();
    auto yDirections = BoidDirections->getArrayY();
    auto zDirections = BoidDirections->getArrayZ();

    auto boidSpeeds = BoidSpeeds;

    for (int j = 0; j < SIZE_OF_SIMULATION; ++j) {
        // Then we must apply the physical movements from this time step

        (*xPositions)[j] += (*xDirections)[j] * (*boidSpeeds)[j];
        (*yPositions)[j] += (*yDirections)[j] * (*boidSpeeds)[j];
        (*zPositions)[j] += (*zDirections)[j] * (*boidSpeeds)[j];
    }
}


std::array<double, 3> BoidSim::getAverageFlockDirection(){
    const auto directionsX = BoidPositions->getArrayX();
    const auto directionsY = BoidPositions->getArrayY();
    const auto directionsZ = BoidPositions->getArrayZ();

    double averageX = 0;
    double averageY = 0;
    double averageZ = 0;

    for (int i = 0; i < SIZE_OF_SIMULATION; ++i){
        averageX += (*directionsX)[i];
        averageY += (*directionsY)[i];
        averageZ += (*directionsZ)[i];
    }

    averageX /= SIZE_OF_SIMULATION;
    averageY /= SIZE_OF_SIMULATION;
    averageZ /= SIZE_OF_SIMULATION;

    std::cout << "Before normalisation, it is " << averageX << ", " << averageY << ", " << averageZ << "\n";

    std::array<double, 3> averageDirection = {averageX, averageY, averageZ};

    return normaliseVector(averageDirection);
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
