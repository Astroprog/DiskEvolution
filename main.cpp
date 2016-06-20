// main.cpp

#include "Simulation.h"
#include <iostream>
#include <mpi.h>

int main(int argc, char** argv)
{
    MPI::Init(argc, argv);

    if (argc == 2){
        Simulation::runDiskDispersalAnalysis(argv[1]);
    } else {
        std::cout << "No parameter file specified" << std::endl;
    }

    MPI::Finalize();


    return 0;
}
