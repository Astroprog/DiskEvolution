//
// Created by Peter Rodenkirch on 12.04.16.
//

#include "GridGeometry.h"
#include <math.h>
#include <iostream>

double GridGeometry::gaussian(double mu, double sigma, double x)
{
    return 1 / (sigma * sqrt(2 * M_PI)) * exp(-0.5 * pow((x - mu) / sigma, 2));
}

GridGeometry::GridGeometry()
{
    scaleFactor = 1.0;
    offset = 0.0;
}

GridGeometry::GridGeometry(double innerBound, double outerBound, int NCells, bool useLogscale)
{
    NGrid = NCells;
    logscale = useLogscale;

    if (useLogscale) {
        offset = log10(innerBound);
        scaleFactor = (log10(outerBound) - offset) / NGrid;
    } else {
        offset = innerBound;
        scaleFactor = outerBound / NGrid;
        std::cout << scaleFactor << ", " << offset << std::endl;
    }
}


