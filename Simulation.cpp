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
    double rg = pMap["rg"];
    double floorDensity = pMap["floor"];
    double luminosity = pMap["luminosity"];
    int frames = (int)pMap["frames"];
    int time = (int)pMap["time"];

    GridGeometry *g = new GridGeometry(rin, rout, NGrid, logscale, 1.49597871e13 * 1.49597871e13, 1.49597871e13, 3.15e7);
    DiskWind *disk = new DiskWind(NGrid);
    disk->setParameters(a, mass, luminosity, rg);
    disk->setNumberOfFrames(frames);
    disk->setGeometry(g);

    disk->initWithHCGADensityDistribution(initialDiskMassRatio * mass, radialScaleFactor, floorDensity);
    disk->computedx();
    disk->computedt();

    std::vector<double> *lambda = new std::vector<double>;
    lambda->push_back(1.3);
    lambda->push_back(1.4);
    lambda->push_back(1.5);

    disk->runDispersalAnalysis(10000000, lambda);
}

void Simulation::runOrdinarySimulation(char *parseString)
{
    std::map<std::string, double> pMap = ParameterParser::parseFile("parameter.par");

    bool restart = (bool)pMap["restart"];
    bool logscale = (bool)pMap["logscale"];
    int NGrid = (int)pMap["ngrid"];
    double rin = pMap["rin"];
    double rout = pMap["rout"];
    double a = pMap["a"];
    double mass = pMap["mass"];
    double initialDiskMassRatio = pMap["diskmass"];
    double radialScaleFactor = pMap["rscale"];
    double rg = pMap["rg"];
    double floorDensity = pMap["floor"];
    double luminosity = pMap["luminosity"];
    int frames = (int)pMap["frames"];
    int time = (int)pMap["time"];

    GridGeometry *g = new GridGeometry(rin, rout, NGrid, logscale, 1.49597871e13 * 1.49597871e13, 1.49597871e13, 3.15e7);
    DiskWind *disk = new DiskWind(NGrid);
    disk->setGeometry(g);
    disk->setParameters(a, mass, luminosity, rg);
    disk->setNumberOfFrames(frames);
    disk->setLeverArm(1.3);

    if (restart) {
        int restartFrame = (int)pMap["restartframe"];
        disk->initWithRestartData(restartFrame);
        disk->computedx();
        disk->computedt();
        disk->restartSimulation(restartFrame, time);
    } else {
        disk->initWithHCGADensityDistribution(initialDiskMassRatio * mass, radialScaleFactor, floorDensity);
        disk->computedx();
        disk->computedt();
        disk->runSimulation(time);
    }
}
