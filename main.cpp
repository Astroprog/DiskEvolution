#include "DiskWind.h"
#include "ParameterParser.h"
#include <iostream>

int main(int argc, char** argv)
{
    if (argc == 2){
        std::map<std::string, double> pMap = ParameterParser::parseFile(argv[1]);

        bool restart = (bool)pMap["restart"];
        bool logscale = (bool)pMap["logscale"];
        int NGrid = (int)pMap["ngrid"];
        double rin = pMap["rin"];
        double rout = pMap["rout"];
        double a = pMap["a"];
        double mass = pMap["mass"];
        double density = pMap["density"];
        double cutoff = pMap["cutoff"];
        int frames = (int)pMap["frames"];
        int time = (int)pMap["time"];

        GridGeometry *g = new GridGeometry(rin, rout, NGrid, logscale);
        DiskWind *disk = new DiskWind(NGrid);
        disk->setParameters(a, mass);
        disk->setNumberOfFrames(frames);
        disk->setGeometry(g);

        if (restart) {
            int restartFrame = (int)pMap["restartframe"];
            disk->initWithRestartData(restartFrame);
            disk->computedx();
            disk->computedt();
            disk->restartSimulation(restartFrame, time);
        } else {
            disk->initWithDensityDistribution(density, cutoff);
            disk->computedx();
            disk->computedt();
            disk->runSimulation(time);
        }

    } else {
        std::cout << "No parameter file specified" << std::endl;
    }


    return 0;
}