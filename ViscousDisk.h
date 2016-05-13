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

};


#endif //DIFFUSION_VISCOUSDISK_H
