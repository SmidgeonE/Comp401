#include "BoidSim.h"

int main(int argc, char* argv[]) {

    // Here we init OpenMP, MPI

    int numThreads = NUM_THREADS == -1 ? omp_get_max_threads() : NUM_THREADS;

    int nameLen, numProcesses, thisProcess;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &thisProcess);
    MPI_Get_processor_name(processor_name, &nameLen);

    omp_set_num_threads(numThreads);

    if (thisProcess == MASTER_PROCESS){
        std::cout << "Number of threads to use : " << numThreads << std::endl;
        std::cout << "Number of boids supplied : " << SIZE_OF_SIMULATION << std::endl;
    }

    auto startTime = MPI_Wtime();


    // Then we start the Sim using BoidSim
    if (!DO_MULTIPLE_SIMS){
        if (thisProcess == MASTER_PROCESS) std::cout << "Only doing one timeSteps = 1000" << std::endl;

        BoidSim newSim(numProcesses, thisProcess);

        newSim.SetWriteToFile(WRITE_SIM);

        newSim.StartSimulation(1000);
    }

    // std::array<double, NUM_SIMULATIONS> timeTakenArray;
    // timeTakenArray.fill(0.0);

    // for (int i = 0; i < timeTakenArray.size(); ++i){
    //     if (DEBUG && i != 1) continue;

    //     if (thisProcess == 0) std::cout << "Num of time Steps : " << 100 * std::pow(10, i) << std::endl;

    //     timeTakenArray[i] = runAndTimeSimulation(100 * std::pow(10, i), WRITE_SIM, numProcesses, thisProcess);

    //     if (WRITE_SIM) break;
    // }


    // std::cout << "Time taken in python array format:" << std::endl;
    // std::cout << "[";

    // for (int i = 0; i < NUM_SIMULATIONS; ++i){
    //     std::cout << timeTakenArray[i] << ", ";
    // }
    
    // std::cout << "]" << std::endl;


    if (thisProcess == MASTER_PROCESS) {
        std::cout << "Simulation Complete. Time Taken: " << MPI_Wtime() - startTime << " seconds" << std::endl;
    }
    
    MPI_Finalize();

    return 0;
}