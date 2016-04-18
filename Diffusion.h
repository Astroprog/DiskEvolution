//
// Created by Peter Rodenkirch on 12.04.16.
//

#ifndef DIFFUSION_DIFFUSION_H
#define DIFFUSION_DIFFUSION_H

#include "Simulation.h"
#include "GridGeometry.h"

class Diffusion : public Simulation {

public:
    Diffusion();
    Diffusion(int ncells);

    double diffValue(double x);
    void initWithGaussian(double mu, double sigma);
    virtual void step();
    virtual void computedt();
    virtual double computeFluxDiff(int i);

};


#endif //DIFFUSION_DIFFUSION_H
