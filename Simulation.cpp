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



