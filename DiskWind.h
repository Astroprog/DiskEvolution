//
// DiskWind.h
// Created by Peter Rodenkirch on 19.04.16.
//

#ifndef DISKEVOLUTION_DISKWIND_H
#define DISKEVOLUTION_DISKWIND_H

#include <vector>
#include <string>
#include "GridGeometry.h"
#include "EvaporationModel.h"

class DiskWind {

private:

    struct Point {
        double x;
        double y;
        double mdot;
        double B2;
    };

    enum mpiTag {initSendLow, initSendHigh, finalRecvLow, finalRecvHigh, frameRecv, frameBRecv,
        frameMDotRecv, dispersalSend, massLossRight, totalWindloss};

    const double au = 1.495978707e13;
    const double G = 6.67408e-8;
    const double kb = 1.38064852e-16;
    const double mp = 1.672621898e-24;
    const double T0 = 280;
    const double normLuminosity = 1e41;
    const double year = 3.1536e7;

    EvaporationModel* photoevModel;

    bool constLambda;
    bool constB;
    bool fluxFreezing;
    bool perfectFluxFreezing;
    double alpha;
    double M;

    double leverArm;
    double floorDensity;
    double diskMass;
    double radialScale;
    double luminosity;
    double photoRadius;
    double viscousConstant;
    double plasma;
    double freezingPlasma;
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
    double accumulatedMassLossLeft;
    double accumulatedMassLossRight;
    double accumulatedWindLoss;
    GridGeometry *g;
    Point data[160];
    double *initialDensity;
    double *tempData;
    double *windloss;
    double *flux;

public:

    DiskWind();
    DiskWind(int ncells);
    ~DiskWind();

    void step();
    void computedt();
    void computedx();
    double getUpdatedMagneticFluxDensityAtCell(int i);
    void determineDiskExtent();
    double computeDiskMass();
    void computeFluxes(int minIndex, int maxIndex);


    double densityLossAtRadius(double r, int i);
    double leverArmAtCell(double i, double currentWindloss);
    double constantLeverArm();

    void setParameters(double a, double mass, double lum, double rg, double lever, int NFrames,
                       GridGeometry *geometry, double plasmaParameter, bool constlambda, bool constb, bool freezing, bool pfreezing);
    void setEvaporationModel(EvaporationModel* model);
    void initWithHCGADensityDistribution(double initialDiskMass, double radialScaleFactor, double floor);

    void runSimulation(int years);
    void runDispersalAnalysis(int timeLimit, std::vector<double>* parameters, const std::string parameterType);
    void restartSimulation(int lastFrame, int years);
    void initWithRestartData(int lastFrame);

    void writeFrame();

};


#endif //DISKEVOLUTION_DISKWIND_H
