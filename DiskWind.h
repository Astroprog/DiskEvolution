//
// Created by Peter Rodenkirch on 19.04.16.
//

#ifndef DISKEVOLUTION_DISKWIND_H
#define DISKEVOLUTION_DISKWIND_H

#include "ViscousDisk.h"

class DiskWind : public ViscousDisk {

    DiskWind();
    DiskWind(int ncells);

    virtual void step();
    virtual double massLossAtRadius(double r);

};


#endif //DISKEVOLUTION_DISKWIND_H
