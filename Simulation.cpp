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

    bool logscale = (bool)pMap["logscale"]; // Logscale grid acitvated
    int NGrid = (int)pMap["ngrid"];         // Number of grid cells
    double rin = pMap["rin"];               // inner simulation boundary in AU
    double rout = pMap["rout"];             // outer simulation boundary in AU
    double a = pMap["a"];                   // viscosity alpha parameter
    double mass = pMap["mass"];             // mass of the central star in g
    double initialDiskMassRatio = pMap["diskmass"]; // initial disk mass in solar masses
    double radialScaleFactor = pMap["rscale"]; // scale factor for the initial surface density distribution
    double floorDensity = pMap["floor"];    // floor desnity as lower density boundary
    int frames = (int)pMap["frames"];       // Number of output frames
    int time = (int)pMap["time"];           // total duration of the Simulation in yr
    double luminosity = pMap["luminosity"]; // luminosity of the central star
    double rg = pMap["rg"];                 // gravitational radius - limit for photoevaporation
    int plasma = pMap["plasma"];            // plasma parameter as ration of gas pressure vs magnetic pressure

    GridGeometry *g = new GridGeometry(rin, rout, NGrid, logscale);
    DiskWind *disk = new DiskWind(NGrid);
    disk->setParameters(a, mass, luminosity, rg, 1.0, frames, g, plasma);
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
    int plasma = pMap["plasma"];
    int frames = (int)pMap["frames"];
    int time = (int)pMap["time"];

    GridGeometry *g = new GridGeometry(rin, rout, NGrid, logscale);
    DiskWind *disk = new DiskWind(NGrid);
    disk->setParameters(a, mass, luminosity, rg, 4.0, frames, g, plasma);

    if (restart) {
        int restartFrame = (int)pMap["restartframe"];
        disk->initWithRestartData(restartFrame);
        disk->restartSimulation(restartFrame, time);
    } else {
        disk->initWithHCGADensityDistribution(initialDiskMassRatio * mass, radialScaleFactor, floorDensity);
        disk->runSimulation(time);
    }
}
