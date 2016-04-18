//
// Created by Peter Rodenkirch on 12.04.16.
//

#include "Simulation.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>

void Simulation::writeFrame()
{
    std::cout << "Writing frame " << frame / frameStride << " (step " << frame << ")" << std::endl;
    std::vector<std::string> stringOutput;
//    double sum = 0.0;
    for(int i = 0; i < NGrid; i++)
    {
        std::stringstream ss;
        ss << data[i].x << " " << data[i].y;
        stringOutput.push_back(ss.str());
//        if (i > 0) {
//            sum += (data[i].x - data[i - 1].x) * data[i].y;
//        }
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
    frameStride = 1;
    data = (Point *)malloc(NGrid * sizeof(Point));
}

Simulation::Simulation(int ncells)
{
    NGrid = ncells;
    frame = 0;
    frameStride = 1;
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

void Simulation::setFrameStride(int stride)
{
    frameStride = stride;
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

void Simulation::runSimulation(int steps)
{
    for (int i = 0; i < steps; i++) {
        step();
    }
}