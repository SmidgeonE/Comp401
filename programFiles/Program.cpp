#include "BoidSim.h"

int main(int argc, char* argv[]) {
    auto boidSim = new BoidSim();

    boidSim->setWriteToFile(true);
    boidSim->SimView(5);
    boidSim->StartSimulation(1000);

    delete boidSim;

    return 0;
}