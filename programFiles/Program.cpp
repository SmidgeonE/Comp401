#include "BoidSim.h"


int main(int argc, char* argv[]) {
    int numThreads = 1;

    if (NUM_THREADS != -1) {
        numThreads = NUM_THREADS;
    }
    else{
        numThreads = omp_get_max_threads();
    }

    omp_set_num_threads(numThreads);
    std::cout << "Number of threads to use : " << numThreads << std::endl;
    std::cout << "Number of boids supplied : " << SIZE_OF_SIMULATION << std::endl;

    std::array<double, NUM_SIMULATIONS> timeTakenArray = {0, 0, 0, 0, 0, 0};

    for (int i = 0; i < NUM_SIMULATIONS; ++i){
        timeTakenArray[i] = runAndTimeSimulation(100 * std::pow(10, i), false);
    }


    std::cout << "Time taken in python array format:" << std::endl;
    std::cout << "[";

    for (int i = 0; i < NUM_SIMULATIONS; ++i){
        std::cout << timeTakenArray[i] << ", ";
    }
    
    std::cout << "]" << std::endl;

    return 0;
}