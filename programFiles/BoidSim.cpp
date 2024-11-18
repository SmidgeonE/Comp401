//
// Created by oz21652 on 9/24/24.
//

#include <cmath>
#include "mpi.h"
#include "BoidSim.h"
#include <random>


BoidSim::BoidSim() : boidPositions(new VectorArray),
                     boidDirections(new VectorArray),
                     boidMasses(new std::array<double, SIZE_OF_SIMULATION>),
                     boidSpeeds(new std::array<double, SIZE_OF_SIMULATION>), 
                     boidForces(new VectorArray),
                     writeToFile(false) {

    boidPositions->InitialiseVectorsToLine(10);
    boidDirections->InitialiseRandomVectors(-1.0, 1.0, true);
    boidSpeeds->fill(0.0);
    boidMasses->fill(1.0);
}


void BoidSim::StartSimulation(const long timeSteps) {
    if (DEBUG) {
        SimView(3);
    }

    if (writeToFile) {
        // We need to write the initial state of the simulation

        writeBoidSimulation();
    }

    // Now we must iterate through the time steps

    for (double i = 0.0; i < timeSteps * DT; i += DT) {
        if (DEBUG && i < DT) {
            std::cout << "Next Time Step: \n\n";
            SimView(3);
        }

        resetForces();

        applySeparationForce();
        // applyAlignmentForce();
        // applyCohesionForce();

        calculateBoidVelocity();

        applyTimeStep();

        if (writeToFile) {
            writeBoidSimulation();
        }
    }
}


void BoidSim::writeBoidSimulation() {
    std::string outputStr = "";

    for (int i = 0; i < SIZE_OF_SIMULATION; ++i) {
        outputStr += "[" + std::to_string((*boidPositions->GetArrayX())[i]) + "," 
            + std::to_string((*boidPositions->GetArrayY())[i]) + "," 
            + std::to_string((*boidPositions->GetArrayZ())[i]) + ","
            + std::to_string((*boidDirections->GetArrayX())[i]) + ","
            + std::to_string((*boidDirections->GetArrayY())[i]) + ","
            + std::to_string((*boidDirections->GetArrayZ())[i]) + ","
            + std::to_string((*boidSpeeds)[i]) + ","
            + std::to_string((*boidMasses)[i]) + "],";
    }

    outputStream << outputStr << "\n";
}


void BoidSim::SetWriteToFile(const bool writeToFile) { 
    this->writeToFile = writeToFile;

    // Opening a stream to write the file

    if (writeToFile){
        outputStream.open("/user/home/oz21652/Comp401/BoidSimulation.csv");
    }
}


void BoidSim::SimView(const int viewNum) const {
    const auto boidPositionsX = boidPositions->GetArrayX();
    const auto boidPositionsY = boidPositions->GetArrayY();
    const auto boidPositionsZ = boidPositions->GetArrayZ();

    const auto boidDirectionsX = boidDirections->GetArrayX();
    const auto boidDirectionsY = boidDirections->GetArrayY();
    const auto boidDirectionsZ = boidDirections->GetArrayZ();


    for (int i = 0; i < viewNum; ++i){
        std::cout << i << "th Boid:" 
        << "\n   -- position: " << (*boidPositionsX)[i] << ", " << (*boidPositionsY)[i] << ", " << (*boidPositionsZ)[i]
        << "\n   -- direction: " << (*boidDirectionsX)[i] << ", " << (*boidDirectionsY)[i] << ", " << (*boidDirectionsZ)[i]
        << "\n   -- direction magnitude: " << std::sqrt((*boidDirectionsX)[i]*(*boidDirectionsX)[i] + (*boidDirectionsY)[i]*(*boidDirectionsY)[i] + (*boidDirectionsZ)[i]*(*boidDirectionsZ)[i])
        << "\n   -- speed: " << (*this->boidSpeeds)[i]
        << "\n   -- mass: " << (*this->boidMasses)[i] << "\n";
    }
}


void BoidSim::applyCohesionForce(){
    const auto positionsX = boidPositions->GetArrayX();
    const auto positionsY = boidPositions->GetArrayY();
    const auto positionsZ = boidPositions->GetArrayZ();

    // Calculate COM of the flock

    std::array<double, 3> com = {0.0, 0.0, 0.0};
    double totalInvMass = 0.0;
    VectorArray* comForces = new VectorArray;

#pragma omp parallel private(totalInvMass, com)
{

#pragma omp for reduction(+:com, totalInvMass)
    for (int i = 0; i < SIZE_OF_SIMULATION; ++i){
        com[0] += (*positionsX)[i] * (*boidMasses)[i];
        com[1] += (*positionsY)[i] * (*boidMasses)[i];
        com[2] += (*positionsZ)[i] * (*boidMasses)[i];

        totalInvMass += (*boidMasses)[i];
    }

#pragma omp barrier

#pragma omp for reduction(*:com)
    for (int i = 0; i < 3; ++i){
        com[i] *= totalInvMass;
    }

#pragma omp barrier

    // Now we have the COM we need to add a linear force based on the distance from the COM

#pragma omp for
    for (int i = 0; i < SIZE_OF_SIMULATION; ++i){
        auto comDirectionX = com[0] - (*positionsX)[i];
        auto comDirectionY = com[1] - (*positionsY)[i];
        auto comDirectionZ = com[2] - (*positionsZ)[i];

        (*comForces->GetArrayX())[i] = comDirectionX;
        (*comForces->GetArrayY())[i] = comDirectionY;
        (*comForces->GetArrayZ())[i] = comDirectionZ;
    }
}

    addForce(*comForces);

    delete comForces;
}


void BoidSim::applySeparationForce(){
    const auto positionsX = boidPositions->GetArrayX();
    const auto positionsY = boidPositions->GetArrayY();
    const auto positionsZ = boidPositions->GetArrayZ();

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

            (*repulsionForces->GetArrayX())[i] += -SEPARATION_FORCE_CONSTANT * separationX / magnitudeFactor;
            (*repulsionForces->GetArrayY())[i] += -SEPARATION_FORCE_CONSTANT * separationY / magnitudeFactor;
            (*repulsionForces->GetArrayZ())[i] += -SEPARATION_FORCE_CONSTANT * separationZ / magnitudeFactor;
        }

        if (numRejections == SIZE_OF_SIMULATION - 1) continue;

    }

    addForce(*repulsionForces);

    delete repulsionForces;
}



void BoidSim::applyAlignmentForce(){
    const auto averageDirection = getAverageFlockDirection();

    VectorArray* alignmentForces = new VectorArray;

    for (int i = 0; i < SIZE_OF_SIMULATION; ++i){
        auto alignmentForceX = ALIGNMENT_FORCE_CONSTANT * averageDirection[0];
        auto alignmentForceY = ALIGNMENT_FORCE_CONSTANT * averageDirection[1];
        auto alignmentForceZ = ALIGNMENT_FORCE_CONSTANT * averageDirection[2];


        (*alignmentForces->GetArrayX())[i] = alignmentForceX;
        (*alignmentForces->GetArrayY())[i] = alignmentForceY;
        (*alignmentForces->GetArrayZ())[i] = alignmentForceZ;
    }

    addForce(*alignmentForces);

    delete alignmentForces;
}


void BoidSim::applyTimeStep(){
    auto xPositions = boidPositions->GetArrayX();
    auto yPositions = boidPositions->GetArrayY();
    auto zPositions = boidPositions->GetArrayZ();

    auto xDirections = boidDirections->GetArrayX();
    auto yDirections = boidDirections->GetArrayY();
    auto zDirections = boidDirections->GetArrayZ();

    for (int j = 0; j < SIZE_OF_SIMULATION; ++j) {
        // Then we must apply the physical movements from this time step

        (*xPositions)[j] += (*xDirections)[j] * (*boidSpeeds)[j];
        (*yPositions)[j] += (*yDirections)[j] * (*boidSpeeds)[j];
        (*zPositions)[j] += (*zDirections)[j] * (*boidSpeeds)[j];
    }
}

void BoidSim::calculateBoidVelocity(){
    auto xDirections = boidDirections->GetArrayX();
    auto yDirections = boidDirections->GetArrayY();
    auto zDirections = boidDirections->GetArrayZ();


    for (int j = 0; j < SIZE_OF_SIMULATION; ++j) {
        // We have now calculated the total force for this time step, now we must calculate the new velocity because of that

        auto newVelocityX = (*xDirections)[j] * (*boidMasses)[j] + (*boidForces->GetArrayX())[j] * DT / (*boidMasses)[j];
        auto newVelocityY = (*yDirections)[j] * (*boidMasses)[j] + (*boidForces->GetArrayY())[j] * DT / (*boidMasses)[j];
        auto newVelocityZ = (*zDirections)[j] * (*boidMasses)[j] + (*boidForces->GetArrayZ())[j] * DT / (*boidMasses)[j];

        auto velocity = std::array<double, 3> {newVelocityX, newVelocityY, newVelocityZ};
        
        (*boidSpeeds)[j] = std::sqrt(magSquared(velocity));
        
        auto velocityNormalised = normaliseVector(velocity);

        (*xDirections)[j] = velocityNormalised[0];
        (*yDirections)[j] = velocityNormalised[1];
        (*zDirections)[j] = velocityNormalised[2];
    }
}


void BoidSim::resetForces(){
    boidForces->GetArrayX()->fill(0.0);
    boidForces->GetArrayY()->fill(0.0);
    boidForces->GetArrayZ()->fill(0.0);
}


void BoidSim::addForce(const VectorArray& force){
    const auto forcesX = boidForces->GetArrayX();
    const auto forcesY = boidForces->GetArrayY();
    const auto forcesZ = boidForces->GetArrayZ();

    for (int i = 0; i < SIZE_OF_SIMULATION; ++i){
        (*forcesX)[i] += (*force.GetArrayX())[i];
        (*forcesY)[i] += (*force.GetArrayY())[i];
        (*forcesZ)[i] += (*force.GetArrayZ())[i];
    }
}


std::array<double, 3> BoidSim::getAverageFlockDirection(){
    std::array<double, 3> averageDirection = boidDirections->GetVectorAverage();

    return normaliseVector(averageDirection);
}


BoidSim::~BoidSim() { 
    delete boidPositions; 
    delete boidDirections; 
    delete boidMasses; 
    delete boidSpeeds; 

    if (writeToFile) {
        outputStream.close();
    }
}
