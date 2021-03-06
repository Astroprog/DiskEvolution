//
// Simulation.cpp
// Created by Peter Rodenkirch on 17.05.16.
//

#include <vector>
#include <iostream>
#include <string>
#include <mpi.h>

#include "Simulation.h"
#include "ParameterParser.h"
#include "DiskWind.h"
#include "EvaporationModel.h"


void Simulation::runDiskDispersalAnalysis(char* parseString)
{
    std::map<std::string, double> pMap = ParameterParser::parseFile(parseString);

    bool logscale = (bool)pMap["logscale"]; // Logscale grid acitvated
    bool constLambda = (bool)pMap["constlambda"];
    bool constB = (bool)pMap["constb"];
    bool freezing = (bool)pMap["freezing"];
    bool pfreezing = (bool)pMap["pfreezing"];
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
    double plasma = pMap["plasma"];            // plasma parameter as ration of gas pressure vs magnetic pressure
    double lambda = pMap["lambda"];

    GridGeometry *g = new GridGeometry(rin, rout, NGrid, logscale);
    DiskWind *disk = new DiskWind(NGrid);
    disk->setParameters(a, mass, luminosity, rg, lambda, frames, g, plasma, constLambda, constB, freezing, pfreezing);
    disk->initWithHCGADensityDistribution(initialDiskMassRatio * mass, radialScaleFactor, floorDensity);

    std::vector<double> *parameters = new std::vector<double>;


    if (constLambda) {
        for (double i = 1.0; i < 20.0; i = i + 0.1)
        {
            parameters->push_back(i);
        }

        const std::string parameterType = "lambda";
        disk->runDispersalAnalysis(50000000, parameters, parameterType);
    } else {
        for (double i = 1e7; i > 100; i = i / 1.1) {
            parameters->push_back(i);
            std::cout << i << std::endl;
        }

        const std::string parameterType = "plasma";
        disk->runDispersalAnalysis(50000000, parameters, parameterType);
    }

}

void Simulation::runOrdinarySimulation(char *parseString)
{
    std::map<std::string, double> pMap = ParameterParser::parseFile(parseString);

    bool restart = (bool)pMap["restart"];
    bool logscale = (bool)pMap["logscale"];
    bool constLambda = (bool)pMap["constlambda"];
    bool constB = (bool)pMap["constb"];
    bool freezing = (bool)pMap["freezing"];
    bool pfreezing = (bool)pMap["pfreezing"];
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
    double plasma = pMap["plasma"];
    double lambda = pMap["lambda"];
    int frames = (int)pMap["frames"];
    int time = (int)pMap["time"];

    GridGeometry *g = new GridGeometry(rin, rout, NGrid, logscale);
    EvaporationModel* model = new EvaporationClarke(luminosity, rg);    
    //const std::string path = "surfaceLoss.dat";
    //EvaporationModel* model = new EvaporationIon100(path);
    DiskWind disk(NGrid);

    disk.setEvaporationModel(model);
    disk.setParameters(a, mass, luminosity, rg, lambda, frames, g, plasma, constLambda, constB, freezing, pfreezing);

    if (restart) {
        int restartFrame = (int)pMap["restartframe"];
        disk.initWithRestartData(restartFrame);
        disk.restartSimulation(restartFrame, time);
    } else {
        //disk.initWithHCGADensityDistribution(initialDiskMassRatio * mass, radialScaleFactor, floorDensity);
        disk.initWithCustomDensityDistribution(80.0, 5.2, 20.0, 0.9, floorDensity);
        disk.runSimulation(time);
    }

    delete g;
    delete model;
}
