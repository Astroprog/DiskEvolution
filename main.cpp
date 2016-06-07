#include "Simulation.h"
#include <iostream>

int main(int argc, char** argv)
{
    if (argc == 2){
        Simulation::runOrdinarySimulation(argv[1]);
    } else {
        std::cout << "No parameter file specified" << std::endl;
    }


    return 0;
}