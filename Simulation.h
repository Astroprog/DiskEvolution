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
        double mdot;
    };

    double alpha;
    double au = 1.5e13;
    double G = 6.67e-8;
    double M;
    double kb = 1.38e-16;
    double mp = 1.67e-24;
    double T0 = 280;
    int NGrid;
    int frame;
    int maxFrames;
    int frameStride;
    int outputFrame;
    int lastFrameDetail;
    double currentTime;
    double dx;
    double dt;
    double year = 3.15e7;
    GridGeometry *g;
    Point *data;

    Simulation();
    Simulation(int ncells);
    virtual ~Simulation();
    virtual void writeFrame();
    virtual void setGeometry(GridGeometry *geometry);
    virtual void setNumberOfFrames(int NFrames);
    virtual void computedx();
    virtual void computedt() = 0;
    virtual void step() = 0;


};


#endif //DIFFUSION_SIMULATION_H
