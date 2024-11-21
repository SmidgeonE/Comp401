//
// Created by oz21652 on 9/24/24.
//

#include <cmath>
#include "mpi.h"
#include "BoidSim.h"
#include <random>


BoidSim::BoidSim(int numProcesses, int thisProcess) {
    this->numProcesses = numProcesses;
    this->thisProcess = thisProcess;

    calculateProcessStartEndIndices();

    boidPositions = VectorArray();
    boidDirections = VectorArray();
    boidForces = VectorArray();

    filledCellsX.fill(0);
    filledCellsY.fill(0);
    filledCellsZ.fill(0);

    cellMinima.fill(0.0);
    cellMaxima.fill(0.0);

    simulationBoundaries[0] = BOX_SIZE/2;
    simulationBoundaries[1] = BOX_SIZE/2;
    simulationBoundaries[2] = BOX_SIZE/2;

    for (int i = 0; i < CELL_NUMBER; ++i){
        for (int j = 0; j < CELL_NUMBER; ++j){
            for (int k = 0; k < CELL_NUMBER; ++k){
                cellList[i][j][k] = {};
            }
        }
    }

    
    // Now we must either generate and broadcast state or receive state

    if (thisProcess == 0){
        generateInitialState();
    }

    // We broadcast the initial state to all other processes, note this doesnt need a MPI_Recv() call

    boidPositions.BroadcastVectorArray(0);
    boidDirections.BroadcastVectorArray(0);
}


void BoidSim::StartSimulation(const long timeSteps) {
    if (thisProcess == 0){
        std::cout << "Starting Simulation with on master process " << std::endl;
    }
    else {
        std::cout << "starting worker process " << thisProcess << std::endl;
        return;
    }

    if (writeToFile) {
        writeBoidSimulation();
    }


    for (int i = 0; i < timeSteps; ++i) {
        if (DEBUG && i <= 2) {
            std::cout << "\n\nNext Time Step: \n\n";
            SimView(3);
        }

        if (DEBUG && i > 2) {
            std::cout << "Debug has finished: 3 time steps have occured." << std::endl;
            return;
        }

        resetForces();

        if (DEBUG) startTime = omp_get_wtime();

        // applySeparationForceCellList();
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
        outputStr += "[" + std::to_string(boidPositions.GetArrayX()[i]) + "," 
            + std::to_string(boidPositions.GetArrayY()[i]) + "," 
            + std::to_string(boidPositions.GetArrayZ()[i]) + ","
            + std::to_string(boidDirections.GetArrayX()[i]) + ","
            + std::to_string(boidDirections.GetArrayY()[i]) + ","
            + std::to_string(boidDirections.GetArrayZ()[i]) + ","
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


void BoidSim::SimView(int viewNum) {
    auto& boidPositionsX = boidPositions.GetArrayX();
    auto& boidPositionsY = boidPositions.GetArrayY();
    auto& boidPositionsZ = boidPositions.GetArrayZ();

    auto& boidDirectionsX = boidDirections.GetArrayX();
    auto& boidDirectionsY = boidDirections.GetArrayY();
    auto& boidDirectionsZ = boidDirections.GetArrayZ();

    if (viewNum >= SIZE_OF_SIMULATION) {
        std::cout << "View number is greater than the number of boids in the simulation\n";
        std::cout << "Defaulting to view all boids\n";
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
    const auto& positionsX = boidPositions.GetArrayX();
    const auto& positionsY = boidPositions.GetArrayY();
    const auto& positionsZ = boidPositions.GetArrayZ();
    
    // Calculate COM of the flock

    double comX = 0.0;
    double comY = 0.0;
    double comZ = 0.0;

    double totalMass = 0.0;
    auto comForces = new VectorArray();


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

    if (DEBUG) {
        cohTime += omp_get_wtime() - startTime;
        startTime = omp_get_wtime();
    }

    delete comForces;
}


void BoidSim::applySeparationForceCellList(){
    const auto& positionsX = boidPositions.GetArrayX();
    const auto& positionsY = boidPositions.GetArrayY();
    const auto& positionsZ = boidPositions.GetArrayZ();

    auto repulsionForces = new VectorArray();
 
    constructCellList();

#pragma omp parallel for
    for (int i = 0; i < SIZE_OF_SIMULATION; ++i){
        auto adjacentBoids = getAdjacentBoids(i);

        for (auto adjacentBoid : adjacentBoids){
            auto separationX = positionsX[adjacentBoid] - positionsX[i];
            auto separationY = positionsY[adjacentBoid] - positionsY[i];
            auto separationZ = positionsZ[adjacentBoid] - positionsZ[i];

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

    if (DEBUG){ 
        sepTime += omp_get_wtime() - startTime;
        startTime = omp_get_wtime();
    }

    wipeCellList();

    delete repulsionForces;
}


void BoidSim::applySeparationForce(){
    const auto& positionsX = boidPositions.GetArrayX();
    const auto& positionsY = boidPositions.GetArrayY();
    const auto& positionsZ = boidPositions.GetArrayZ();

    auto repulsionForces = new VectorArray();

#pragma omp parallel for
    for (int i = 0; i < SIZE_OF_SIMULATION; ++i){
        for (int j = 0; j < SIZE_OF_SIMULATION; ++j){
            if (i == j) continue;
            
            auto separationX = positionsX[j] - positionsX[i];
            auto separationY = positionsY[j] - positionsY[i];
            auto separationZ = positionsZ[j] - positionsZ[i];
        
            auto distance = std::sqrt(magSquared({separationX, separationY, separationZ}));

            if (distance > 1 or std::abs(distance) < std::numeric_limits<double>::epsilon()) continue;

            auto magnitudeFactor = distance * distance * distance;

            auto repulsionX = std::min(-SEPARATION_FORCE_CONSTANT * separationX / magnitudeFactor, 1000.0);
            auto repulsionY = std::min(-SEPARATION_FORCE_CONSTANT * separationY / magnitudeFactor, 1000.0);
            auto repulsionZ = std::min(-SEPARATION_FORCE_CONSTANT * separationZ / magnitudeFactor, 1000.0);

            repulsionForces->GetArrayX()[i] += repulsionX;
            repulsionForces->GetArrayY()[i] += repulsionY;
            repulsionForces->GetArrayZ()[i] += repulsionZ;
        }
    }

    if (DEBUG) repulsionForces->View(3, "Repulsion Forces");

    addForce(*repulsionForces);

    if (DEBUG){ 
        sepTime += omp_get_wtime() - startTime;
        startTime = omp_get_wtime();
    }

    delete repulsionForces;
}



void BoidSim::applyAlignmentForce(){
    auto averageDirection = getAverageFlockDirection();

    auto alignmentForces = new VectorArray();

#pragma omp parallel for 
    for (int i = 0; i < SIZE_OF_SIMULATION; ++i){
        auto alignmentForceX = ALIGNMENT_FORCE_CONSTANT * (averageDirection[0] - boidDirections.GetArrayX()[i]);
        auto alignmentForceY = ALIGNMENT_FORCE_CONSTANT * (averageDirection[1] - boidDirections.GetArrayY()[i]);
        auto alignmentForceZ = ALIGNMENT_FORCE_CONSTANT * (averageDirection[2] - boidDirections.GetArrayZ()[i]);

        alignmentForces->GetArrayX()[i] = alignmentForceX;
        alignmentForces->GetArrayY()[i] = alignmentForceY;
        alignmentForces->GetArrayZ()[i] = alignmentForceZ;
    }

    if (DEBUG) alignmentForces->View(3, "Alignment Forces");

    addForce(*alignmentForces);

    if (DEBUG) {
        alignTime += omp_get_wtime() - startTime;
        startTime = omp_get_wtime();
    }

    delete alignmentForces;
}


void BoidSim::applyTimeStep(){
    auto& xPositions = boidPositions.GetArrayX();
    auto& yPositions = boidPositions.GetArrayY();
    auto& zPositions = boidPositions.GetArrayZ();

    auto& xDirections = boidDirections.GetArrayX();
    auto& yDirections = boidDirections.GetArrayY();
    auto& zDirections = boidDirections.GetArrayZ();


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
    auto& xDirections = boidDirections.GetArrayX();
    auto& yDirections = boidDirections.GetArrayY();
    auto& zDirections = boidDirections.GetArrayZ();

#pragma omp parallel for
    for (int j = 0; j < SIZE_OF_SIMULATION; ++j) {
        // We have now calculated the total force for this time step, now we must calculate the new velocity because of that

        auto addVelX = boidForces.GetArrayX()[j] / boidMasses[j];
        auto addVelY = boidForces.GetArrayY()[j] / boidMasses[j];
        auto addVelZ = boidForces.GetArrayZ()[j] / boidMasses[j];

        auto newVelocityX = xDirections[j] * boidSpeeds[j] + boidForces.GetArrayX()[j] / boidMasses[j];
        auto newVelocityY = yDirections[j] * boidSpeeds[j] + boidForces.GetArrayY()[j] / boidMasses[j];
        auto newVelocityZ = zDirections[j] * boidSpeeds[j] + boidForces.GetArrayZ()[j] / boidMasses[j];


        auto velocity = std::array<double, 3> {newVelocityX, newVelocityY, newVelocityZ};
        
        boidSpeeds[j] = std::sqrt(magSquared(velocity));

        auto velocityNormalised = normaliseVector(velocity);

        xDirections[j] = velocityNormalised[0];
        yDirections[j] = velocityNormalised[1];
        zDirections[j] = velocityNormalised[2];

        // Now we can apply the walls of the simulation

        if (std::abs(boidPositions.GetArrayX()[j]) > simulationBoundaries[0]){
            xDirections[j] *= -1.0;
        }

        if (std::abs(boidPositions.GetArrayY()[j]) > simulationBoundaries[1]){
            yDirections[j] *= -1.0;
        }

        if (std::abs(boidPositions.GetArrayZ()[j]) > simulationBoundaries[2]){
            zDirections[j] *= -1.0;
        }
    }
}


void BoidSim::resetForces(){
    auto& arrayX = boidForces.GetArrayX();
    auto& arrayY = boidForces.GetArrayY();
    auto& arrayZ = boidForces.GetArrayZ();

    arrayX.fill(0.0);
    arrayY.fill(0.0);
    arrayZ.fill(0.0);
}


void BoidSim::addForce(VectorArray& force){
    auto& forcesX = boidForces.GetArrayX();
    auto& forcesY = boidForces.GetArrayY();
    auto& forcesZ = boidForces.GetArrayZ();

#pragma omp parallel for
    for (int i = 0; i < SIZE_OF_SIMULATION; ++i){
        // std::cout << "adding force: " << force.GetArrayX()[i] << ", " << force.GetArrayY()[i] << ", " << force.GetArrayZ()[i] << std::endl;
        forcesX[i] += force.GetArrayX()[i];
        forcesY[i] += force.GetArrayY()[i];
        forcesZ[i] += force.GetArrayZ()[i];

        // std::cout << "forces after adding : " << forcesX[i] << ", " << forcesY[i] << ", " << forcesZ[i] << std::endl;
    }
}


std::array<double, 3> BoidSim::getAverageFlockDirection(){
    auto averageDirection = boidDirections.GetVectorAverage();

    return normaliseVector(averageDirection);
}


void BoidSim::constructCellList(){
    auto& positionsX = boidPositions.GetArrayX();
    auto& positionsY = boidPositions.GetArrayY();
    auto& positionsZ = boidPositions.GetArrayZ();

    auto minMaxX = minMaxOfArray(positionsX);
    auto minMaxY = minMaxOfArray(positionsY);
    auto minMaxZ = minMaxOfArray(positionsZ);

    cellMinima[0] = minMaxX[0];
    cellMinima[1] = minMaxY[0];
    cellMinima[2] = minMaxZ[0];

    cellMaxima[0] = minMaxX[1];
    cellMaxima[1] = minMaxY[1];
    cellMaxima[2] = minMaxZ[1];

#pragma omp parallel for
    for (int i = 0; i < SIZE_OF_SIMULATION; ++i){
        auto boidCell = getBoidCell(i);

        filledCellsX[i] = boidCell[0];
        filledCellsY[i] = boidCell[1];
        filledCellsZ[i] = boidCell[2];

#pragma omp critical
        {
            cellList[boidCell[0]][boidCell[1]][boidCell[2]].push_back(i);
        }
    }
}


std::array<int, 3> BoidSim::getBoidCell(const int boidIndex){
    auto positionsX = boidPositions.GetArrayX();
    auto positionsY = boidPositions.GetArrayY();
    auto positionsZ = boidPositions.GetArrayZ();

    auto cellX = (int) (((positionsX[boidIndex] - cellMinima[0]) / (cellMaxima[0] - cellMinima[0])) * CELL_NUMBER);
    auto cellY = (int) (((positionsY[boidIndex] - cellMinima[1]) / (cellMaxima[1] - cellMinima[1])) * CELL_NUMBER);
    auto cellZ = (int) (((positionsZ[boidIndex] - cellMinima[2]) / (cellMaxima[2] - cellMinima[2])) * CELL_NUMBER);

    // Special case: when the position is exactly the max or min, we must round down:

    cellX = std::min(cellX, CELL_NUMBER - 1);
    cellY = std::min(cellY, CELL_NUMBER - 1);
    cellZ = std::min(cellZ, CELL_NUMBER - 1);

    return {cellX, cellY, cellZ};
}


void BoidSim::wipeCellList(){
    for (int i = 0; i < SIZE_OF_SIMULATION; ++i){
        auto cellX = filledCellsX[i];
        auto cellY = filledCellsY[i];
        auto cellZ = filledCellsZ[i];

        cellList[cellX][cellY][cellZ].clear();
    }
}


std::vector<int> BoidSim::getAdjacentBoids(const int boidIndex){
    auto boidCell = getBoidCell(boidIndex);
    std::vector<int> adjacentBoids = {};


    for (int ix = boidCell[0]-1; ix < boidCell[0]+2; ++ix){
        for (int iy = boidCell[1]-1; iy < boidCell[1]+2; ++iy){
            for (int iz = boidCell[2]-1; iz < boidCell[2]+2; ++iz){
                if (ix < 0 || ix >= CELL_NUMBER || iy < 0 || iy >= CELL_NUMBER || iz < 0 || iz >= CELL_NUMBER) continue;

                auto currentAdjacentBoids = cellList[ix][iy][iz];

                if (currentAdjacentBoids.size() > 0){
                    for (int boid : currentAdjacentBoids){
                        if (boid != boidIndex){
                            adjacentBoids.push_back(boid);
                        }
                    }
                }
            }
        }
    }

    return adjacentBoids;
}


void BoidSim::generateInitialState() {
    if (DEBUG) { 
        boidPositions.InitialiseVectorsToLine(10);
        boidDirections.InitaliseVectorsToZHat();
        // boidSpeeds.fill(0.0);

        // boidPositions->InitialiseRandomVectors(-10.0, 10.0, false);
        // boidDirections->InitialiseRandomVectors(-1.0, 1.0, true);
    }
    else {
        boidPositions.InitialiseRandomVectors(-10.0, 10.0, false);
        boidDirections.InitialiseRandomVectors(-1.0, 1.0, true);
        initialiseRandomScalars(boidSpeeds, 0.0, 1.0);
    }

    boidMasses.fill(1.0);
}


BoidSim::~BoidSim() { 
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


void BoidSim::calculateProcessStartEndIndices(){
    thisProcessStartIndex = (int) ((thisProcess / (double)numProcesses) * SIZE_OF_SIMULATION);
    thisProcessEndIndex = (int) ((((thisProcess+1) / (double)numProcesses) * SIZE_OF_SIMULATION)-1);

    if (thisProcess == numProcesses - 1){
        thisProcessEndIndex = SIZE_OF_SIMULATION - 1;
    }

    if (DEBUG) std::cout << "Process " << thisProcess << " has indices " << thisProcessStartIndex << " to " << thisProcessEndIndex << std::endl;
}
