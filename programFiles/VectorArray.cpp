#include "BoidSim.h"


VectorArray::VectorArray() {
    arrayX.fill(0.0);
    arrayY.fill(0.0);
    arrayZ.fill(0.0);
}
                            

std::array<double, 3> VectorArray::GetVectorAverage() {
    double averageX = 0;
    double averageY = 0;
    double averageZ = 0;

#pragma omp parallel reduction(+: averageX, averageY, averageZ) 
    for (int i = 0; i < SIZE_OF_SIMULATION; ++i){
        averageX += arrayX[i];
        averageY += arrayY[i];
        averageZ += arrayZ[i];
    }

    averageX /= SIZE_OF_SIMULATION;
    averageY /= SIZE_OF_SIMULATION;
    averageZ /= SIZE_OF_SIMULATION;

    std::array<double, 3> averageDirection = {averageX, averageY, averageZ};

    return averageDirection;
}


void VectorArray::InitialiseRandomVectors(const double lowerBound, const double upperBound, const bool normalise) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(lowerBound, upperBound);

    for (int i = 0; i < SIZE_OF_SIMULATION; ++i) {
        const double x = dis(gen);
        const double y = dis(gen);
        const double z = dis(gen);

        auto vector = std::array<double, 3> {x, y, z};

        if (normalise) {
            vector = normaliseVector(vector);
        }


        GetArrayX().at(i) = vector[0];
        GetArrayY().at(i) = vector[1];
        GetArrayZ().at(i) = vector[2];
    }
}


void VectorArray::InitialiseVectorsToLine(const double gridSize) {
    // Sets vectors to the line y=x=z

    auto spacing = gridSize / SIZE_OF_SIMULATION;

    for (int i = 0; i < SIZE_OF_SIMULATION; ++i) {
        GetArrayX()[i] = i * spacing;
        GetArrayY()[i] = i * spacing;
        GetArrayZ()[i] = i * spacing;
    }
}


void VectorArray::InitaliseVectorsToZHat() {
    for (int i = 0; i < SIZE_OF_SIMULATION; ++i) {
        GetArrayX().at(i) = 0;
        GetArrayY().at(i) = 0;
        GetArrayZ().at(i) = 1;
    }
}


void VectorArray::BroadcastVectorArray(const int rootProcess) {
    MPI_Bcast(GetArrayX().data(), SIZE_OF_SIMULATION, MPI_DOUBLE, rootProcess, MPI_COMM_WORLD);
    MPI_Bcast(GetArrayY().data(), SIZE_OF_SIMULATION, MPI_DOUBLE, rootProcess, MPI_COMM_WORLD);
    MPI_Bcast(GetArrayZ().data(), SIZE_OF_SIMULATION, MPI_DOUBLE, rootProcess, MPI_COMM_WORLD);
}


void VectorArray::SendVectorArray(const int destinationProcess, const int currentProcess) {
    // The tag schema is the following : 
    // currentProcess -> the process id
    // what is sent: "tag" + either "0", "1", or "2" for x, y, z respectively

    MPI_Send(GetArrayX().data(), SIZE_OF_SIMULATION, MPI_DOUBLE, destinationProcess, currentProcess*10, MPI_COMM_WORLD);
    MPI_Send(GetArrayY().data(), SIZE_OF_SIMULATION, MPI_DOUBLE, destinationProcess, currentProcess*10+1, MPI_COMM_WORLD);
    MPI_Send(GetArrayZ().data(), SIZE_OF_SIMULATION, MPI_DOUBLE, destinationProcess, currentProcess*10+2, MPI_COMM_WORLD);
}


VectorArray VectorArray::ReceiveVectorArray(const int sourceProcess, const int currentProcess) {
    // The tag schema is the following : 
    // sourceProcess -> the process id
    // what is received: "tag" + either "0", "1", or "2" for x, y, z respectively
    
    VectorArray receivedArray;

    MPI_Recv(receivedArray.GetArrayX().data(), SIZE_OF_SIMULATION, MPI_DOUBLE, sourceProcess, sourceProcess*10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(receivedArray.GetArrayY().data(), SIZE_OF_SIMULATION, MPI_DOUBLE, sourceProcess, sourceProcess*10+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(receivedArray.GetArrayZ().data(), SIZE_OF_SIMULATION, MPI_DOUBLE, sourceProcess, sourceProcess*10+2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    return receivedArray;
}
