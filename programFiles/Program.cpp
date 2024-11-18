#include "BoidSim.h"



int main(int argc, char* argv[]) {
    std::array<double, 6> timeTakenArray = {0, 0, 0, 0, 0, 0};

    std::cout << "Number of boids supplied : " << SIZE_OF_SIMULATION << std::endl;

    timeTakenArray[0] = runAndTimeSimulation(100, false);
    timeTakenArray[1] = runAndTimeSimulation(1000, false);
    timeTakenArray[2] = runAndTimeSimulation(10000, false);
    timeTakenArray[3] = runAndTimeSimulation(100000, false);
    timeTakenArray[4] = runAndTimeSimulation(1000000, false);
    timeTakenArray[5] = runAndTimeSimulation(10000000, false);

    std::cout << "Time taken in python array format:" << std::endl;
    
    std::cout << "[" << timeTakenArray[0] << ", " 
    << timeTakenArray[1] << ", " 
    << timeTakenArray[2] << ", " 
    << timeTakenArray[3] << ", " 
    << timeTakenArray[4] << ", " 
    << timeTakenArray[5] << "]" 
    << std::endl;

    return 0;
}