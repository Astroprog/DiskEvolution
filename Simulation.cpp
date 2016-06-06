//
// Created by Peter Rodenkirch on 17.05.16.
//

#include <vector>
#include "Simulation.h"
#include "ParameterParser.h"
#include "DiskWind.h"


void Simulation::runDiskDispersalAnalysis(char* parseString)
{
    std::map<std::string, double> pMap = ParameterParser::parseFile(parseString);

    bool logscale = (bool)pMap["logscale"];
    int NGrid = (int)pMap["ngrid"];
    double rin = pMap["rin"];
    double rout = pMap["rout"];
    double a = pMap["a"];
    double mass = pMap["mass"];
    double initialDiskMassRatio = pMap["diskmass"];
    double radialScaleFactor = pMap["rscale"];
    double floorDensity = pMap["floor"];
    int frames = (int)pMap["frames"];
    int time = (int)pMap["time"];
    double luminosity = pMap["luminosity"];
    double rg = pMap["rg"];

    GridGeometry *g = new GridGeometry(rin, rout, NGrid, logscale);
    DiskWind *disk = new DiskWind(NGrid);
    disk->setParameters(a, mass, luminosity, rg, 1.0, frames, g);
    disk->initWithHCGADensityDistribution(initialDiskMassRatio * mass, radialScaleFactor, floorDensity);

    std::vector<double> *lambda = new std::vector<double>;
    for (double i = 2.4; i < 10.0; i = i + 0.1) {
        lambda->push_back(i);
    }

    disk->runDispersalAnalysis(10000000, lambda);
}

void Simulation::runOrdinarySimulation(char *parseString)
{
    std::map<std::string, double> pMap = ParameterParser::parseFile(parseString);

    bool restart = (bool)pMap["restart"];
    bool logscale = (bool)pMap["logscale"];
    int NGrid = (int)pMap["ngrid"];
    double rin = pMap["rin"];
    double rout = pMap["rout"];
    double a = pMap["a"];
    double mass = pMap["mass"];
    double initialDiskMassRatio = pMap["diskmass"];
    double radialScaleFactor = pMap["rscale"];
    double floorDensity = pMap["floor"];
    double luminosity = pMap["luminosity"];
    double rg = pMap["rg"];
    int frames = (int)pMap["frames"];
    int time = (int)pMap["time"];

    GridGeometry *g = new GridGeometry(rin, rout, NGrid, logscale);
    DiskWind *disk = new DiskWind(NGrid);
    disk->setParameters(a, mass, luminosity, rg, 4.0, frames, g);

    if (restart) {
        int restartFrame = (int)pMap["restartframe"];
        disk->initWithRestartData(restartFrame);
        disk->restartSimulation(restartFrame, time);
    } else {
        disk->initWithHCGADensityDistribution(initialDiskMassRatio * mass, radialScaleFactor, floorDensity);
        disk->runSimulation(time);
    }
}
