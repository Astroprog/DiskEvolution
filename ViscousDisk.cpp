//
// Created by Peter Rodenkirch on 13.04.16.
//

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iterator>
#include "ViscousDisk.h"


ViscousDisk::ViscousDisk()
{
    NGrid = 2;
    frame = 0;
    frameStride = 1;
    data = (Point *)malloc(NGrid * sizeof(Point));
}

ViscousDisk::ViscousDisk(int ncells)
{
    NGrid = ncells;
    frame = 0;
    frameStride = 1;
    data = (Point *)malloc(NGrid * sizeof(Point));
}

void ViscousDisk::writeFrame()
{
    std::cout << "Writing frame " << frame / frameStride << " (step " << frame << ", " << dt/year * frame << "yr)" << std::endl;
    std::vector<std::string> stringOutput;

    for(int i = 0; i < NGrid; i++)
    {
        std::stringstream ss;
        ss << data[i].x << " " << data[i].y / data[i].x;
        stringOutput.push_back(ss.str());
    }
    std::ostringstream tempStream;
    tempStream << "frame" << frame / frameStride << ".dat";
    std::ofstream outputFile(tempStream.str());
    std::ostream_iterator<std::string> output_iterator(outputFile, "\n");
    std::copy(stringOutput.begin(), stringOutput.end(), output_iterator);
}


void ViscousDisk::step()
{
    double *tempData = (double *)malloc(NGrid * sizeof(double));

    double c = 3 * alpha * kb * T0 / (sqrt(au) * 2.3 * mp * sqrt(G * M));

    for (int i = 0; i < NGrid; i++) {
        tempData[i] = data[i].y - c * dt / (g->convertIndexToPosition(i + 0.5) - g->convertIndexToPosition(i - 0.5)) * computeFluxDiff(i);
    }

    for (int i = 0; i < NGrid; i++) {
        data[i].y = tempData[i];
    }

    free(tempData);

    frame++;
    if (frame % frameStride == 0) {
        writeFrame();
    }
}

double ViscousDisk::computeFluxDiff(int i)
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

    Fright = -sqrt(rPlusHalf) * (0.5 * (y + yPlus) * (s+1)*pow(rPlusHalf, s) + pow(rPlusHalf, s+1) * (yPlus - y) / (rPlus - r));
    Fleft = -sqrt(rMinusHalf) * (0.5 * (y + yMinus) * (s+1)*pow(rMinusHalf, s) + pow(rMinusHalf, s+1) * (y - yMinus) / (r - rMinus));
    return Fright - Fleft;
}

void ViscousDisk::setParameters(double a, double mass)
{
    alpha = a;
    M = mass;
}

void ViscousDisk::initWithDensityDistribution(double densityAt1Au, double cutoff)
{
    std::cout << "Initializing grid with size " << NGrid << std::endl;

    for (int i = 0; i < NGrid; i++) {
        data[i].x = g->convertIndexToPosition(i);
        data[i].y = densityAt1Au * exp(-data[i].x / cutoff);
    }

    writeFrame();
}

void ViscousDisk::computedt()
{
    double D = pow(data[0].x, s+1);
    for (int i = 0; i < NGrid; i++) {
        double temp = pow(data[i].x, s+1);
        if (temp > D) {
            D = temp;
        }
    }

    dt = year / 3;
}


void ViscousDisk::runSimulation(int years)
{
    for (int i = 0; dt/year * i < years; i++) {
        step();
    }
}