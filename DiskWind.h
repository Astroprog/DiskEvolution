//
// Created by Peter Rodenkirch on 19.04.16.
//

#ifndef DISKEVOLUTION_DISKWIND_H
#define DISKEVOLUTION_DISKWIND_H

#include "GridGeometry.h"

class DiskWind {

private:
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

public:

    DiskWind();
    DiskWind(int ncells);
    ~DiskWind();

    void step();
    void computedt();
    void computedx();
    double computeFluxDiff(int i);

    double massLossAtRadius(double r, double rg);
    double leverArmAtRadius(double r);

    void setParameters(double a, double mass);
    void initWithDensityDistribution(double densityAt1Au, double cutoff);
    void setGeometry(GridGeometry *geometry);
    void setNumberOfFrames(int NFrames);

    void runSimulation(int years);
    void restartSimulation(int lastFrame, int years);
    void initWithRestartData(int lastFrame);

    void writeFrame();

};


#endif //DISKEVOLUTION_DISKWIND_H
