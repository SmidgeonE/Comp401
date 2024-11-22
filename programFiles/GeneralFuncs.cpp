#include "BoidSim.h"

void initialiseRandomScalars(std::array<double, SIZE_OF_SIMULATION>& scalarArray, const double lowerBound, const double upperBound) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(lowerBound, upperBound);

#pragma omp parallel for
    for (int i = 0; i < SIZE_OF_SIMULATION; ++i) {
        scalarArray[i] = dis(gen);
    }
}


std::array<double, 3> normaliseVector(const std::array<double, 3>& vector){
    double magnitude = std::sqrt(magSquared(vector));

    if (magnitude < std::numeric_limits<double>::epsilon()) {
        return {0.0, 0.0, 0.0};
    }

    return {vector[0]/magnitude, vector[1]/magnitude, vector[2]/magnitude};
}


double magSquared(const std::array<double, 3>& vector) {
    if (std::abs(vector[0]) < std::numeric_limits<double>::epsilon() &&
        std::abs(vector[1]) < std::numeric_limits<double>::epsilon() &&
        std::abs(vector[2]) < std::numeric_limits<double>::epsilon()) {
        return 0.0;
    }

    return vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2];
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
    
