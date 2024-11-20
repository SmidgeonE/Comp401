//
// Created by oz21652 on 9/24/24.
//

#include <cmath>
#include "mpi.h"
#include "BoidSim.h"
#include <random>
#include <bits/algorithmfwd.h>


BoidSim::BoidSim(){
    boidPositions = new VectorArray;
    boidDirections = new VectorArray;
    boidForces = new VectorArray;

    boidPositions->InitialiseRandomVectors(-10.0, 10.0, false);
    boidDirections->InitialiseRandomVectors(-1.0, 1.0, true);

    boidSpeeds.fill(0.0);
    boidMasses.fill(1.0);
}


void BoidSim::StartSimulation(const long timeSteps) {
    if (writeToFile) {
        writeBoidSimulation();
    }


    for (int i = 0; i < DT; ++i) {
        if (DEBUG && i < 2 * DT) {
            std::cout << "\n\nNext Time Step: \n\n";
            SimView(3);
        }

        if (DEBUG && i > 2 * DT) {
            std::cout << "Debug has finished: 3 time steps have occured." << std::endl;
            return;
        }


        resetForces();

        if (DEBUG) startTime = omp_get_wtime();

        applySeparationForce();
        applyAlignmentForce();
        applyCohesionForce();

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
        outputStr += "[" + std::to_string(boidPositions->GetArrayX()[i]) + "," 
            + std::to_string(boidPositions->GetArrayY()[i]) + "," 
            + std::to_string(boidPositions->GetArrayZ()[i]) + ","
            + std::to_string(boidDirections->GetArrayX()[i]) + ","
            + std::to_string(boidDirections->GetArrayY()[i]) + ","
            + std::to_string(boidDirections->GetArrayZ()[i]) + ","
            + std::to_string(boidSpeeds[i]) + ","
            + std::to_string(boidMasses[i]) + "],";
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


void BoidSim::SimView(int viewNum) const {
    auto boidPositionsX = boidPositions->GetArrayX();
    auto boidPositionsY = boidPositions->GetArrayY();
    auto boidPositionsZ = boidPositions->GetArrayZ();

    auto boidDirectionsX = boidDirections->GetArrayX();
    auto boidDirectionsY = boidDirections->GetArrayY();
    auto boidDirectionsZ = boidDirections->GetArrayZ();

    if (viewNum >= SIZE_OF_SIMULATION) {
        std::cout << "View number is greater than the number of boids in the simulation\n";
        std::cout << "Defaulting to view all boids\n";
        return;
    }

    viewNum = std::min(viewNum, SIZE_OF_SIMULATION);


    for (int i = 0; i < viewNum; ++i){
        std::cout << i << "th Boid:" 
        << "\n   -- position: " << boidPositionsX[i] << ", " << boidPositionsY[i] << ", " << boidPositionsZ[i]
        << "\n   -- direction: " << boidDirectionsX[i] << ", " << boidDirectionsY[i] << ", " << boidDirectionsZ[i]
        << "\n   -- direction magnitude: " << std::sqrt(boidDirectionsX[i]*boidDirectionsX[i] + boidDirectionsY[i]*boidDirectionsY[i] + boidDirectionsZ[i]*boidDirectionsZ[i])
        << "\n   -- speed: " << this->boidSpeeds[i]
        << "\n   -- mass: " << this->boidMasses[i] << "\n";
    }
}


void BoidSim::applyCohesionForce(){    
    const auto positionsX = boidPositions->GetArrayX();
    const auto positionsY = boidPositions->GetArrayY();
    const auto positionsZ = boidPositions->GetArrayZ();
    
    // Calculate COM of the flock

    double comX = 0.0;
    double comY = 0.0;
    double comZ = 0.0;

    double totalMass = 0.0;
    VectorArray* comForces = new VectorArray;


#pragma omp parallel for reduction(+:comX, comY, comZ, totalMass)
    for (int i = 0; i < SIZE_OF_SIMULATION; ++i){
        comX += positionsX[i] * boidMasses[i];
        comY += positionsY[i] * boidMasses[i];
        comZ += positionsZ[i] * boidMasses[i];

        totalMass += boidMasses[i];
    }

    if (DEBUG) std::cout << "Total Mass: " << totalMass << std::endl;

    comX /= totalMass;
    comY /= totalMass;
    comZ /= totalMass;

    if (DEBUG) std::cout << "COM: " << comX << ", " << comY << ", " << comZ << std::endl;

    // Now we have the COM we need to add a linear force based on the distance from the COM

#pragma omp for
    for (int i = 0; i < SIZE_OF_SIMULATION; ++i){
        auto comDirectionX = comX - positionsX[i];
        auto comDirectionY = comY - positionsY[i];
        auto comDirectionZ = comZ - positionsZ[i];

        comForces->GetArrayX()[i] = COHESION_FORCE_CONSTANT * comDirectionX;
        comForces->GetArrayY()[i] = COHESION_FORCE_CONSTANT * comDirectionY;
        comForces->GetArrayZ()[i] = COHESION_FORCE_CONSTANT * comDirectionZ;
    }

    if (DEBUG) comForces->View(3, "COM Forces");

    addForce(*comForces);

    delete comForces;

    if (DEBUG) {
        cohTime += omp_get_wtime() - startTime;
        startTime = omp_get_wtime();
    }
}


// void BoidSim::applySeparationForceCellList(){
//     const auto positionsX = boidPositions->GetArrayX();
//     const auto positionsY = boidPositions->GetArrayY();
//     const auto positionsZ = boidPositions->GetArrayZ();

//     VectorArray* repulsionForces = new VectorArray;

//     constructCellList();
// }


void BoidSim::applySeparationForce(){
    const auto positionsX = boidPositions->GetArrayX();
    const auto positionsY = boidPositions->GetArrayY();
    const auto positionsZ = boidPositions->GetArrayZ();

    VectorArray* repulsionForces = new VectorArray;

#pragma omp parallel for
    for (int i = 0; i < SIZE_OF_SIMULATION; ++i){
        for (int j = 0; j < SIZE_OF_SIMULATION; ++j){
            if (i == j) continue;
            
            auto separationX = positionsX[j] - positionsX[i];
            auto separationY = positionsY[j] - positionsY[i];
            auto separationZ = positionsZ[j] - positionsZ[i];
        
            auto magnitudeFactor = std::sqrt(magSquared({separationX, separationY, separationZ}));
            magnitudeFactor *= magnitudeFactor * magnitudeFactor;

            auto repulsionX = std::min(-SEPARATION_FORCE_CONSTANT * separationX / magnitudeFactor, 10.0);
            auto repulsionY = std::min(-SEPARATION_FORCE_CONSTANT * separationY / magnitudeFactor, 10.0);
            auto repulsionZ = std::min(-SEPARATION_FORCE_CONSTANT * separationZ / magnitudeFactor, 10.0);

            repulsionForces->GetArrayX()[i] += repulsionX;
            repulsionForces->GetArrayY()[i] += repulsionY;
            repulsionForces->GetArrayZ()[i] += repulsionZ;
        }
    }

    if (DEBUG) repulsionForces->View(3, "Repulsion Forces");

    addForce(*repulsionForces);

    delete repulsionForces;

    if (DEBUG){ 
        sepTime += omp_get_wtime() - startTime;
        startTime = omp_get_wtime();
    }
}



void BoidSim::applyAlignmentForce(){
    const auto averageDirection = getAverageFlockDirection();

    VectorArray* alignmentForces = new VectorArray;

#pragma omp parallel for 
    for (int i = 0; i < SIZE_OF_SIMULATION; ++i){
        auto alignmentForceX = ALIGNMENT_FORCE_CONSTANT * (averageDirection[0] - boidDirections->GetArrayX()[i]);
        auto alignmentForceY = ALIGNMENT_FORCE_CONSTANT * (averageDirection[1] - boidDirections->GetArrayY()[i]);
        auto alignmentForceZ = ALIGNMENT_FORCE_CONSTANT * (averageDirection[2] - boidDirections->GetArrayZ()[i]);

        alignmentForces->GetArrayX()[i] = alignmentForceX;
        alignmentForces->GetArrayY()[i] = alignmentForceY;
        alignmentForces->GetArrayZ()[i] = alignmentForceZ;
    }

    if (DEBUG) alignmentForces->View(3, "Alignment Forces");

    addForce(*alignmentForces);

    delete alignmentForces;

    if (DEBUG) {
        alignTime += omp_get_wtime() - startTime;
        startTime = omp_get_wtime();
    }
}


void BoidSim::applyTimeStep(){
    auto xPositions = boidPositions->GetArrayX();
    auto yPositions = boidPositions->GetArrayY();
    auto zPositions = boidPositions->GetArrayZ();

    auto xDirections = boidDirections->GetArrayX();
    auto yDirections = boidDirections->GetArrayY();
    auto zDirections = boidDirections->GetArrayZ();


#pragma omp parallel for
    for (int j = 0; j < SIZE_OF_SIMULATION; ++j) {
        // Then we must apply the physical movements from this time step


        xPositions[j] += xDirections[j] * boidSpeeds[j];
        yPositions[j] += yDirections[j] * boidSpeeds[j];
        zPositions[j] += zDirections[j] * boidSpeeds[j];
    }

    if (DEBUG){
        velTime += omp_get_wtime() - startTime;
        startTime = omp_get_wtime();
    }
}

void BoidSim::calculateBoidVelocity(){
    auto xDirections = boidDirections->GetArrayX();
    auto yDirections = boidDirections->GetArrayY();
    auto zDirections = boidDirections->GetArrayZ();

#pragma omp parallel for
    for (int j = 0; j < SIZE_OF_SIMULATION; ++j) {
        // We have now calculated the total force for this time step, now we must calculate the new velocity because of that

        auto newVelocityX = xDirections[j] * boidMasses[j] + boidForces->GetArrayX()[j] * DT / boidMasses[j];
        auto newVelocityY = yDirections[j] * boidMasses[j] + boidForces->GetArrayY()[j] * DT / boidMasses[j];
        auto newVelocityZ = zDirections[j] * boidMasses[j] + boidForces->GetArrayZ()[j] * DT / boidMasses[j];

        auto velocity = std::array<double, 3> {newVelocityX, newVelocityY, newVelocityZ};
        
        boidSpeeds[j] = std::sqrt(magSquared(velocity));
        
        auto velocityNormalised = normaliseVector(velocity);

        xDirections[j] = velocityNormalised[0];
        yDirections[j] = velocityNormalised[1];
        zDirections[j] = velocityNormalised[2];
    }
}


void BoidSim::resetForces(){
    std::cout << "Resetting Forces" << std::endl;
    auto arrayX = boidForces->GetArrayX();
    std::cout << "Array: " << arrayX[0] << std::endl;
    boidForces->GetArrayX().fill(0.0);

    std::cout << "X Forces reset" << std::endl;
    boidForces->GetArrayY().fill(0.0);
    boidForces->GetArrayZ().fill(0.0);
}


void BoidSim::addForce(VectorArray& force){
    auto forcesX = boidForces->GetArrayX();
    auto forcesY = boidForces->GetArrayY();
    auto forcesZ = boidForces->GetArrayZ();

#pragma omp parallel for
    for (int i = 0; i < SIZE_OF_SIMULATION; ++i){
        forcesX[i] += force.GetArrayX()[i];
        forcesY[i] += force.GetArrayY()[i];
        forcesZ[i] += force.GetArrayZ()[i];
    }
}


std::array<double, 3> BoidSim::getAverageFlockDirection(){
    std::array<double, 3> averageDirection = boidDirections->GetVectorAverage();

    return normaliseVector(averageDirection);
}

// void BoidSim::constructCellList(){
//     auto positionsX = boidPositions->GetArrayX();
//     auto positionsY = boidPositions->GetArrayY();
//     auto positionsZ = boidPositions->GetArrayZ();

//     auto minMaxX = minMaxOfArray(positionsX);
//     auto minMaxY = minMaxOfArray(positionsY);
//     auto minMaxZ = minMaxOfArray(positionsZ);

//     (*cellMinima)[0] = minMaxX[0];
//     (*cellMinima)[1] = minMaxY[0];
//     (*cellMinima)[2] = minMaxZ[0];

//     (*cellMaxima)[0] = minMaxX[1];
//     (*cellMaxima)[1] = minMaxY[1];
//     (*cellMaxima)[2] = minMaxZ[1];

// #pragma omp parallel for
//     for (int i = 0; i < SIZE_OF_SIMULATION; ++i){
//         auto cellX = (int) ((((*positionsX)[i] - (*cellMinima)[0]) / ((*cellMaxima)[0] - (*cellMinima)[0])) * CELL_NUMBER);
//         auto cellY = (int) ((((*positionsY)[i] - (*cellMinima)[1]) / ((*cellMaxima)[1] - (*cellMinima)[1])) * CELL_NUMBER);
//         auto cellZ = (int) ((((*positionsZ)[i] - (*cellMinima)[2]) / ((*cellMaxima)[2] - (*cellMinima)[2])) * CELL_NUMBER);

//         (*filledCellsX)[i] = cellX;
//         (*filledCellsY)[i] = cellY;
//         (*filledCellsZ)[i] = cellZ;

// #pragma omp critical
//         {
//             (*cellList)[cellX][cellY][cellZ].push_back(i);
//         }
//     }
// }


// void BoidSim::wipeCellList(){
//     for (int i = 0; i < SIZE_OF_SIMULATION; ++i){
//         auto cellX = ((*filledCellsX)[i]);
//         auto cellY = ((*filledCellsY)[i]);
//         auto cellZ = ((*filledCellsZ)[i]);

//         (*cellList)[cellX][cellY][cellZ].clear();
//     }
// }


BoidSim::~BoidSim() { 
    delete boidPositions; 
    delete boidDirections; 
    delete boidForces;

    if (DEBUG) {
        std::cout << "Time taken for Separation Force: " << sepTime << " seconds" << std::endl;
        std::cout << "Time taken for Alignment Force: " << alignTime << " seconds" << std::endl;
        std::cout << "Time taken for Cohesion Force: " << cohTime << " seconds" << std::endl;
        std::cout << "Time taken for Velocity Calculation: " << velTime << " seconds" << std::endl;
    }

    if (writeToFile) {
        outputStream.close();
    }
}
