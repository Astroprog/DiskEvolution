//
// Created by Peter Rodenkirch on 12.04.16.
//

#ifndef DIFFUSION_SIMULATION_H
#define DIFFUSION_SIMULATION_H

#include "GridGeometry.h"

class Simulation {
public:

    struct Point {
        double x;
        double y;
    };

    int NGrid;
    int frame;
    int frameStride;
    double dx;
    double dt;
    GridGeometry *g;
    Point *data;

    Simulation();
    Simulation(int ncells);
    virtual ~Simulation();
    virtual void writeFrame();
    virtual void setGeometry(GridGeometry *geometry);
    virtual void setFrameStride(int stride);
    virtual void computedx();
    virtual void computedt() = 0;
    virtual void runSimulation(int years);
    virtual void step() = 0;



};


#endif //DIFFUSION_SIMULATION_H
