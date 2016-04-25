//
// Created by Peter Rodenkirch on 13.04.16.
//

#ifndef DIFFUSION_VISCOUSDISK_H
#define DIFFUSION_VISCOUSDISK_H

#include "Diffusion.h"

class ViscousDisk : public Diffusion {
public:

    ViscousDisk();
    ViscousDisk(int ncells);

    virtual void step();
    virtual double computeFluxDiff(int i);
    virtual void computedt();
    virtual void setParameters(double a, double mass);
    virtual void runSimulation(int years);
    virtual void restartSimulation(int lastFrame, int years);
    virtual void initWithRestartData(int lastFrame);
    virtual void initWithDensityDistribution(double densityAt1Au, double cutoff);
};


#endif //DIFFUSION_VISCOUSDISK_H
