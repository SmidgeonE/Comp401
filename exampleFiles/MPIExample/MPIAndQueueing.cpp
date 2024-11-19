//
// Created by oz21652 on 9/24/24.
//

#include <iostream>
#include "mpi.h"

int main(int argc, char* argv[]) {
    std::cout << "Multiprocess Application \n";

    int numProcesses, rank, nameLen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(processor_name, &nameLen);

    const std::string outputString = "Hello from process " + std::to_string(rank) + " out of " +
        std::to_string(numProcesses) + " process on node " + processor_name + "\n";
    std::cout << outputString;


    MPI_Finalize();

    return 0;
}
