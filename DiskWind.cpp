//
// Created by Peter Rodenkirch on 19.04.16.
//

#include "DiskWind.h"
#include <cstdlib>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <iterator>

DiskWind::DiskWind()
{
    NGrid = 2;
    frame = 0;
    frameStride = 1;
    outputFrame = 0;
    data = (Point *)malloc(NGrid * sizeof(Point));
}

DiskWind::DiskWind(int ncells)
{
    NGrid = ncells;
    frame = 0;
    frameStride = 1;
    outputFrame = 0;
    data = (Point *)malloc(NGrid * sizeof(Point));
}

DiskWind::~DiskWind()
{
    free(data);
    delete(g);
}

void DiskWind::writeFrame()
{
    std::cout << "Writing frame " << outputFrame << " (step " << frame << ", " << dt/year * frame << "yr)" << std::endl;
    std::vector<std::string> stringOutput;

    std::stringstream headerStream;
    headerStream << frame << " " << dt/year * frame;
    stringOutput.push_back(headerStream.str());

    for (int i = 0; i < NGrid; i++)
    {
        std::stringstream ss;
        ss << data[i].x << " " << data[i].y / data[i].x << " " << data[i].mdot;
        stringOutput.push_back(ss.str());
    }
    std::ostringstream tempStream;
    tempStream << "frame" << outputFrame << ".dat";
    std::ofstream outputFile(tempStream.str());
    std::ostream_iterator<std::string> output_iterator(outputFile, "\n");
    std::copy(stringOutput.begin(), stringOutput.end(), output_iterator);
    outputFrame++;
}

void DiskWind::setGeometry(GridGeometry *geometry)
{
    g = geometry;
}

void DiskWind::setNumberOfFrames(int NFrames)
{
    maxFrames = NFrames;
}

void DiskWind::computedx()
{
    dx = data[1].x - data[0].x;  //Determining the smallest dx in data
    if (NGrid > 2) {
        for (int i = 2; i < NGrid; i++) {
            double temp = data[i].x - data[i-1].x;
            if (temp < dx) {
                dx = temp;
            }
        }
    }
}

void DiskWind::computedt()
{
    double x = data[1].x - data[0].x;
    double c = 3 * alpha * kb * T0 / (sqrt(au) * 2.3 * mp * sqrt(G * M));
    dt = 0.25 * pow(x, 1.5)/c;

    for (int i = 2; i < NGrid; i++) {
        x = pow((data[i].x - data[i-1].x), 1.5);
        double temp = 0.25 * x / c;
        if (temp < dt) {
            dt = temp;
        }
    }

    std::cout << "Timestep: " << dt << std::endl;
}

double DiskWind::massLossAtRadius(double r, double rg)
{
    double rate = 1e25;
    if (r >= rg) {
        return rate * pow(r, -2.5);
    } else {
        return 0.0;
    }
}

double DiskWind::leverArmAtRadius(double r)
{
    return 1.0;
}

double DiskWind::constantLeverArm()
{
    return leverArm;
}

void DiskWind::setLeverArm(double arm)
{
    leverArm = arm;
}

double DiskWind::computeFluxDiff(int i)
{
    double rPlus = g->convertIndexToPosition(i + 1.0);
    double rPlusHalf = g->convertIndexToPosition(i + 0.5);
    double r = g->convertIndexToPosition(i);
    double rMinusHalf = g->convertIndexToPosition(i - 0.5);
    double rMinus = g->convertIndexToPosition(i - 1.0);

    double yPlus = data[i + 1].y;
    double y = data[i].y;
    double yMinus = data[i - 1].y;



    double Fleft, Fright;

    if (i == 0) {
        yMinus = y;
    } else if (i == NGrid - 1) {
        yPlus = 0.0;
    }

    Fright = (0.25 * (y + yPlus) + rPlusHalf * (yPlus - y) / (rPlus - r) + (constantLeverArm() - 1) * massLossAtRadius(r, 5.0) / (M_PI * au * year * (rPlus - r)) * r);
    Fleft = (0.25 * (y + yMinus) + rMinusHalf * (y - yMinus) / (r - rMinus) + (constantLeverArm() - 1) * massLossAtRadius(r, 5.0) / (M_PI * au * year * (r - rMinus)) * r);
    return Fright - Fleft;
}

void DiskWind::step()
{
    double *tempData = (double *)malloc(NGrid * sizeof(double));

    double c = 3 * alpha * kb * T0 / (sqrt(au) * 2.3 * mp * sqrt(G * M));

    for (int i = 0; i < NGrid; i++) {
        double fluxDiff = computeFluxDiff(i);
        double dr = g->convertIndexToPosition(i+0.5) - g->convertIndexToPosition(i-0.5);

        double densityLoss = massLossAtRadius(data[i].x, 5.0) / (2 * M_PI * au * au * dr) / year;
        tempData[i] = data[i].y + dt * (c / dr * fluxDiff - densityLoss);
    }

    for (int i = 0; i < NGrid; i++) {
        if (tempData[i] / data[i].x < floorDensity)
        {
            data[i].y = floorDensity * data[i].x;
        } else {
            data[i].y = tempData[i];
        }
    }

    free(tempData);

    frame++;
    if (frame % frameStride == 0) {
        writeFrame();
    }
}

void DiskWind::initWithRestartData(int lastFrame)
{
    std::cout << "Initializing with frame " << lastFrame << std::endl;
    std::stringstream fileStream;
    fileStream << "frame" << lastFrame << ".dat";
    std::ifstream input(fileStream.str());

    std::string line;
    int i = 0;
    while (std::getline(input, line))
    {
        std::istringstream iss(line);
        if (i == 0) {
            if(!(iss >> lastFrameDetail >> currentTime)) {
                std::cout << "Error while parsing restart file (first line)" << std::endl;
                break;
            }
            i++;
        } else {
            if (!(iss >> data[i-1].x >> data[i-1].y))
            {
                std::cout << "Error while parsing restart file" << std::endl;
                break;
            } else {
                data[i-1].y *= data[i-1].x;
                data[i-1].mdot = 0.0;
            }
            i++;
        }
    }
}


void DiskWind::setParameters(double a, double mass)
{
    alpha = a;
    M = mass;
}


void DiskWind::initWithHCGADensityDistribution(double initialDiskMass, double radialScaleFactor, double floor)
{
    std::cout << "Initializing grid of size " << NGrid << " with HCGA distribution" << std::endl;

    floorDensity = floor;
    radialScale = radialScaleFactor;
    diskMass = initialDiskMass;

    for (int i = 0; i < NGrid; i++) {
        data[i].x = g->convertIndexToPosition(i);
        data[i].y = initialDiskMass / (2 * M_PI * au * au * radialScaleFactor * data[i].x) * exp(-data[i].x / radialScaleFactor) * data[i].x;
        if (data[i].y / data[i].x < floorDensity) {
            data[i].y = floorDensity * data[i].x;
        }
        data[i].mdot = 0.0;
    }

    writeFrame();
}


void DiskWind::restartSimulation(int lastFrame, int years)
{
    double NSteps = (double)years * year / dt;
    frameStride = (int)(NSteps / (double)maxFrames);

    frame = lastFrameDetail;
    outputFrame = lastFrame + 1;

    for (int i = 0; dt/year * i < years; i++) {
        step();
    }
}

void DiskWind::runDispersalAnalysis(int timeLimit, std::vector<double>* leverArms)
{
    double NSteps = (double)timeLimit / 10 * year / dt;
    frameStride = (int)(NSteps / (double)maxFrames);

    for (unsigned int i = 0; i < leverArms->size(); i++) {
        leverArm = leverArms->at(i);
        frame = 0;
        outputFrame = 0;
        initWithHCGADensityDistribution(diskMass, radialScale, floorDensity);

        for (int k = 0; dt/year * k < timeLimit; k++) {
            step();
            bool dispersed = false;
            for (int j = 60; j < 66; j++) {
                if (data[j].y / data[j].x <= 2*floorDensity)
                {
                    std::cout << "For lambda = " << leverArms->at(i) << ", disk dispersal is reached after " << dt/year * k << " years." << std::endl;
                    dispersed = true;
                }
            }
            if (dispersed) {
                break;
            }
        }
    }
}

void DiskWind::runSimulation(int years)
{
    double NSteps = (double)years * year / dt;
    frameStride = (int)(NSteps / (double)maxFrames);

    for (int i = 0; dt/year * i < years; i++) {
        step();
    }
}