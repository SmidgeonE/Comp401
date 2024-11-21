#include "BoidSim.h"

void initialiseRandomScalars(std::array<double, SIZE_OF_SIMULATION>* scalarArray, const double lowerBound, const double upperBound) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(lowerBound, upperBound);

#pragma omp parallel for
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


    auto end = omp_get_wtime() - start;

    std::cout << "---Time taken for time steps is: " << end << " seconds" << std::endl;

    delete boidSim;
    
    return end;
}


std::array<double, 2> minMaxOfArray(std::array<double, SIZE_OF_SIMULATION>& array){
    auto min = array[0];
    auto max = array[0];

    for (int i = 1; i < SIZE_OF_SIMULATION; ++i){
        auto value = array[i];

        if (value < min){
            min = value;
        }

        if (value > max){
            max = value;
        }
    }

    return {min, max};
}
    
