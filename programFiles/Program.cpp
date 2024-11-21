#include "BoidSim.h"


int main(int argc, char* argv[]) {
    int numThreads = NUM_THREADS == -1 ? omp_get_max_threads() : NUM_THREADS;

    int numProcesses, rank, nameLen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(processor_name, &nameLen);

    const std::string outputString = "Hello from process " + std::to_string(rank) + " out of " +
    std::to_string(numProcesses) + " process on node " + processor_name + "\n";
    std::cout << outputString;

    

    omp_set_num_threads(numThreads);
    std::cout << "Number of threads to use : " << numThreads << std::endl;
    std::cout << "Number of boids supplied : " << SIZE_OF_SIMULATION << std::endl;


    if (!DO_MULTIPLE_SIMS){
        std::cout << "Only doing one timeSteps = 1000" << std::endl;

        runAndTimeSimulation(1000, WRITE_SIM);

        return 0;
    }

    
    std::array<double, NUM_SIMULATIONS> timeTakenArray;
    timeTakenArray.fill(0.0);

    for (int i = 0; i < timeTakenArray.size(); ++i){
        if (DEBUG && i != 1) continue;

        std::cout << "Num of time Steps : " << 100 * std::pow(10, i) << std::endl;

        timeTakenArray[i] = runAndTimeSimulation(100 * std::pow(10, i), WRITE_SIM);

        if (WRITE_SIM) break;
    }


    std::cout << "Time taken in python array format:" << std::endl;
    std::cout << "[";

    for (int i = 0; i < NUM_SIMULATIONS; ++i){
        std::cout << timeTakenArray[i] << ", ";
    }
    
    std::cout << "]" << std::endl;

    
    MPI_Finalize();

   return 0;
}