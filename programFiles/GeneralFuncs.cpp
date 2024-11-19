#include "BoidSim.h"

void initialiseRandomScalars(std::array<double, SIZE_OF_SIMULATION>* scalarArray, const double lowerBound, const double upperBound) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(lowerBound, upperBound);

    for (int i = 0; i < SIZE_OF_SIMULATION; ++i) {
        (*scalarArray)[i] = dis(gen);
    }
}


std::array<double, 3> normaliseVector(const std::array<double, 3>& vector){
    double magnitude = std::sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
    return {vector[0]/magnitude, vector[1]/magnitude, vector[2]/magnitude};
}


double magSquared(const std::array<double, 3>& vector) {
    return vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2];
}


double runAndTimeSimulation(int timeSteps, bool writeToFile){
    auto start = omp_get_wtime();
    auto boidSim = new BoidSim();
    
    if (writeToFile) {
        boidSim->SetWriteToFile(true);
    }

    boidSim->StartSimulation(timeSteps);

    delete boidSim;

    auto end = omp_get_wtime() - start;

    std::cout << "Time taken for time steps: " << timeSteps << " is " << end << " seconds" << std::endl;
    
    return end;
}
    
