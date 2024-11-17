#include "BoidSim.h"

VectorArray::VectorArray() : arrayX(new std::array<double, SIZE_OF_SIMULATION>),
                             arrayY(new std::array<double, SIZE_OF_SIMULATION>),
                             arrayZ(new std::array<double, SIZE_OF_SIMULATION>) {}
                            
VectorArray::~VectorArray() { delete arrayX; delete arrayY; delete arrayZ; }

