//
// Created by Peter Rodenkirch on 19.04.16.
//

#ifndef DISKEVOLUTION_DISKWIND_H
#define DISKEVOLUTION_DISKWIND_H

#include "ViscousDisk.h"

class DiskWind : public ViscousDisk {

public:

    DiskWind();
    DiskWind(int ncells);

    virtual void step();
    virtual double computeFluxDiff(int i);
    double massLossAtRadius(double r, double rg);
    double leverArmAtRadius(double r);

};


#endif //DISKEVOLUTION_DISKWIND_H
