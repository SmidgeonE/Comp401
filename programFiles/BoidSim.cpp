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
                     BoidSpeeds(new std::array<double, SIZE_OF_SIMULATION>), 
                     BoidForces(new VectorArray),
                     writeToFile(false) {
    initialiseRandomVectors(BoidPositions, -10, 10, false);
    initialiseRandomVectors(BoidDirections, -1.0, 1.0, true);
    BoidSpeeds->fill(0.0);
    BoidMasses->fill(1.0);
}


void BoidSim::StartSimulation(const long timeSteps) {

    if (writeToFile) {
        // We need to write the initial state of the simulation

        WriteBoidSimulation();
    }

    // Now we must iterate through the time steps

    for (double i = 0.0; i < timeSteps * dt; i += dt) {
        resetForces();

        applySeparationAlgorithm();
        applyAlignmentAlgorithm();
        applyCohesionAlgorithm();

        calculateBoidVelocity();

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
        outputStream.open("/user/home/oz21652/Comp401/BoidSimulation.csv");
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
    const auto positionsX = BoidPositions->getArrayX();
    const auto positionsY = BoidPositions->getArrayY();
    const auto positionsZ = BoidPositions->getArrayZ();

    // Calculate COM of the flock

    std::array<double, 3> com = {0.0, 0.0, 0.0};
    double totalMass = 0.0;

    for (int i = 0; i < SIZE_OF_SIMULATION; ++i){
        com[0] += (*positionsX)[i] * (*BoidMasses)[i];
        com[1] += (*positionsY)[i] * (*BoidMasses)[i];
        com[2] += (*positionsZ)[i] * (*BoidMasses)[i];

        totalMass += (*BoidMasses)[i];
    }

    for (int i = 0; i < 3; ++i){
        com[i] /= totalMass;
    }

    // Now we have the COM we need to add a linear force based on the distance from the COM

    VectorArray* comForces = new VectorArray;

    for (int i = 0; i < SIZE_OF_SIMULATION; ++i){
        auto comDirectionX = com[0] - (*positionsX)[i];
        auto comDirectionY = com[1] - (*positionsY)[i];
        auto comDirectionZ = com[2] - (*positionsZ)[i];

        (*comForces->getArrayX())[i] = comDirectionX;
        (*comForces->getArrayY())[i] = comDirectionY;
        (*comForces->getArrayZ())[i] = comDirectionZ;
    }

    addForce(*comForces);

    delete comForces;
}


void BoidSim::applySeparationAlgorithm(){
    const auto positionsX = BoidPositions->getArrayX();
    const auto positionsY = BoidPositions->getArrayY();
    const auto positionsZ = BoidPositions->getArrayZ();

    VectorArray* repulsionForces = new VectorArray;


    for (int i = 0; i < SIZE_OF_SIMULATION; ++i){
        int numRejections = 0;

        for (int j = 0; j < SIZE_OF_SIMULATION; ++j){
            if (i == j) continue;
            
            auto separationX = (*positionsX)[j] - (*positionsX)[i];
            auto separationY = (*positionsY)[j] - (*positionsY)[i];
            auto separationZ = (*positionsZ)[j] - (*positionsZ)[i];
        
            auto magnitudeFactor = std::sqrt(magSquared({separationX, separationY, separationZ}));
            magnitudeFactor *= magnitudeFactor * magnitudeFactor;

            if (magnitudeFactor > 5) {
                ++numRejections;
                continue;
            }

            (*repulsionForces->getArrayX())[i] += -SEPARATION_FORCE_CONSTANT * separationX / magnitudeFactor;
            (*repulsionForces->getArrayY())[i] += -SEPARATION_FORCE_CONSTANT * separationY / magnitudeFactor;
            (*repulsionForces->getArrayZ())[i] += -SEPARATION_FORCE_CONSTANT * separationZ / magnitudeFactor;
        }

        if (numRejections == SIZE_OF_SIMULATION - 1) continue;

    }

    addForce(*repulsionForces);

    delete repulsionForces;
}



void BoidSim::applyAlignmentAlgorithm(){
    const auto averageDirection = getAverageFlockDirection();

    VectorArray* alignmentForces = new VectorArray;

    for (int i = 0; i < SIZE_OF_SIMULATION; ++i){
        auto alignmentForceX = ALIGNMENT_FORCE_CONSTANT * averageDirection[0];
        auto alignmentForceY = ALIGNMENT_FORCE_CONSTANT * averageDirection[1];
        auto alignmentForceZ = ALIGNMENT_FORCE_CONSTANT * averageDirection[2];


        (*alignmentForces->getArrayX())[i] = alignmentForceX;
        (*alignmentForces->getArrayY())[i] = alignmentForceY;
        (*alignmentForces->getArrayZ())[i] = alignmentForceZ;
    }

    addForce(*alignmentForces);

    delete alignmentForces;
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

void BoidSim::calculateBoidVelocity(){
    auto xDirections = BoidDirections->getArrayX();
    auto yDirections = BoidDirections->getArrayY();
    auto zDirections = BoidDirections->getArrayZ();


    for (int j = 0; j < SIZE_OF_SIMULATION; ++j) {
        // We have now calculated the total force for this time step, now we must calculate the new velocity because of that

        auto newVelocityX = (*xDirections)[j] * (*BoidMasses)[j] + (*BoidForces->getArrayX())[j] * dt / (*BoidMasses)[j];
        auto newVelocityY = (*yDirections)[j] * (*BoidMasses)[j] + (*BoidForces->getArrayY())[j] * dt / (*BoidMasses)[j];
        auto newVelocityZ = (*zDirections)[j] * (*BoidMasses)[j] + (*BoidForces->getArrayZ())[j] * dt / (*BoidMasses)[j];

        auto velocity = std::array<double, 3> {newVelocityX, newVelocityY, newVelocityZ};
        
        (*BoidSpeeds)[j] = std::sqrt(magSquared(velocity));
        
        auto velocityNormalised = normaliseVector(velocity);

        (*xDirections)[j] = velocityNormalised[0];
        (*yDirections)[j] = velocityNormalised[1];
        (*zDirections)[j] = velocityNormalised[2];
    }
}


void BoidSim::resetForces(){
    BoidForces->getArrayX()->fill(0.0);
    BoidForces->getArrayY()->fill(0.0);
    BoidForces->getArrayZ()->fill(0.0);
}


void BoidSim::addForce(const VectorArray& force){
    const auto forcesX = BoidForces->getArrayX();
    const auto forcesY = BoidForces->getArrayY();
    const auto forcesZ = BoidForces->getArrayZ();

    for (int i = 0; i < SIZE_OF_SIMULATION; ++i){
        (*forcesX)[i] += (*force.getArrayX())[i];
        (*forcesY)[i] += (*force.getArrayY())[i];
        (*forcesZ)[i] += (*force.getArrayZ())[i];
    }
}


std::array<double, 3> BoidSim::getAverageFlockDirection(){
    std::array<double, 3> averageDirection = BoidDirections->getVectorAverage();

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
