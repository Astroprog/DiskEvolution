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
    double year = 3.15e7;

    ViscousDisk();
    ViscousDisk(int ncells);

    virtual void step();
    virtual double computeFluxDiff(int i);
    virtual void computedt();
    virtual void writeFrame();
    void setParameters(double a, double mass);
    void initWithDensityDistribution(double densityAt1Au, double cutoff);
    virtual void runSimulation(int years);
};


#endif //DIFFUSION_VISCOUSDISK_H
