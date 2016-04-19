//
// Created by Peter Rodenkirch on 12.04.16.
//

#include "Simulation.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <cmath>

void Simulation::writeFrame()
{
    std::cout << "Writing frame " << frame / frameStride << " (step " << frame << ", " << dt/year * frame << "yr)" << std::endl;
    std::vector<std::string> stringOutput;

    std::stringstream headerStream;
    headerStream << dt/year * frame;
    stringOutput.push_back(headerStream.str());

    for (int i = 0; i < NGrid; i++)
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

Simulation::Simulation()
{
    NGrid = 2;
    frame = 0;
    maxFrames = 100;
    data = (Point *)malloc(NGrid * sizeof(Point));
}

Simulation::Simulation(int ncells)
{
    NGrid = ncells;
    frame = 0;
    maxFrames = 100;
    data = (Point *)malloc(NGrid * sizeof(Point));
}

Simulation::~Simulation()
{
    free(data);
    delete(g);
}

void Simulation::setGeometry(GridGeometry *geometry)
{
    g = geometry;
}

void Simulation::setNumberOfFrames(int NFrames)
{
    maxFrames = NFrames;
}

void Simulation::computedx()
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