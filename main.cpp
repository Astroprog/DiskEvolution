// main.cpp

#include "Simulation.h"
#include <iostream>
#include <mpi.h>
#include "ParameterParser.h"
#include "GridGeometry.h"

int main(int argc, char** argv)
{
    MPI::Init(argc, argv);

    if (argc == 2){
        std::map<std::string, double> pMap = ParameterParser::parseFile(argv[1]);
        int simulationType = (int)pMap["dispersalanalysis"];
        if (simulationType == 0) {
            Simulation::runOrdinarySimulation(argv[1]);
        } else {
            Simulation::runDiskDispersalAnalysis(argv[1]);
        }
    } else {
        std::cout << "No parameter file specified" << std::endl;
    }

    MPI::Finalize();

    return 0;
}
