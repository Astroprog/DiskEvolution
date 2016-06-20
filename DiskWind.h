//
// DiskWind.h
// Created by Peter Rodenkirch on 19.04.16.
//

#ifndef DISKEVOLUTION_DISKWIND_H
#define DISKEVOLUTION_DISKWIND_H

#include <vector>
#include <string>
#include "GridGeometry.h"

class DiskWind {

private:

    struct Point {
        double x;
        double y;
        double mdot;
        double B2;
    };

    enum mpiTag {initSendLow, initSendHigh, finalRecvLow, finalRecvHigh, frameRecv, dispersalRecv, dispersalSend};

    double alpha;
    double au = 1.495978707e13;
    double G = 6.67408e-8;
    double M;
    double kb = 1.38064852e-16;
    double mp = 1.672621898e-24;
    double T0 = 280;
    double leverArm;
    double floorDensity;
    double diskMass;
    double radialScale;
    double luminosity;
    double normLuminosity = 1e41;
    double photoRadius;
    double viscousConstant;
    double plasma;
    double currentDiskExtent;
    int NGrid;
    int frame;
    int maxFrames;
    int frameStride;
    int outputFrame;
    int lastFrameDetail;
    double currentTime;
    double dx;
    double dt;
    double year = 3.1536e7;
    GridGeometry *g;
    Point *data;
    double *tempData;

public:

    DiskWind();
    DiskWind(int ncells);
    ~DiskWind();

    void step();
    void computedt();
    void computedx();
    void determineCurrentDiskExtent();
    double getUpdatedMagneticFluxDensityAtCell(int i);
    double computeFluxDiff(int i);

    double densityLossAtRadius(double r);
    double leverArmAtCell(double i);
    double constantLeverArm();

    void setParameters(double a, double mass, double lum, double rg, double lever, int NFrames, GridGeometry *geometry, double plasmaParameter);
    void initWithHCGADensityDistribution(double initialDiskMass, double radialScaleFactor, double floor);

    void runSimulation(int years);
    void runDispersalAnalysis(int timeLimit, std::vector<double>* parameters, const std::string parameterType);
    void restartSimulation(int lastFrame, int years);
    void initWithRestartData(int lastFrame);

    void writeFrame();

};


#endif //DISKEVOLUTION_DISKWIND_H
