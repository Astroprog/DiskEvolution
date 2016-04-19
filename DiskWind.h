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
    virtual double massLossAtRadius(double r, double rg);

};


#endif //DISKEVOLUTION_DISKWIND_H
