//
// Created by Peter Rodenkirch on 13.04.16.
//

#ifndef DIFFUSION_VISCOUSDISK_H
#define DIFFUSION_VISCOUSDISK_H

#include "Diffusion.h"

class ViscousDisk : public Diffusion {
public:

    double alpha;
    double au = 1.5e13;
    double G = 6.67e-8;
    double M;
    double kb = 1.38e-16;
    double mp = 1.67e-24;
    double T0 = 280;
    double s = -0.5;

    ViscousDisk();
    ViscousDisk(int ncells);

    virtual void step();
    virtual double computeFluxDiff(int i);
    virtual void computedt();
    virtual void writeFrame();
    virtual void runSimulation(int years);
    virtual void setParameters(double a, double mass);
    virtual void initWithDensityDistribution(double densityAt1Au, double cutoff);
};


#endif //DIFFUSION_VISCOUSDISK_H
