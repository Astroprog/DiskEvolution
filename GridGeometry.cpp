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

GridGeometry::GridGeometry(double innerBound, double outerBound, int NCells, bool useLogscale, double m0Value, double r0Value, double t0Value)
{
    NGrid = NCells;
    logscale = useLogscale;
    m0 = m0Value;
    r0 = r0Value;
    t0 = t0Value;

    if (useLogscale) {
        offset = log10(innerBound);
        scaleFactor = (log10(outerBound) - offset) / NGrid;
    } else {
        offset = innerBound;
        scaleFactor = outerBound / NGrid;
        std::cout << scaleFactor << ", " << offset << std::endl;
    }
}

double GridGeometry::mass(double cgsMass)
{
    return cgsMass / m0;
}

double GridGeometry::cgsMass(double mass)
{
    return mass * m0;
}

double GridGeometry::radius(double cgsRadius)
{
    return cgsRadius / r0;
}

double GridGeometry::cgsRadius(double radius)
{
    return radius * r0;
}

double GridGeometry::density(double cgsDensity)
{
    return cgsDensity * r0 * r0 / m0;
}

double GridGeometry::cgsDensity(double density)
{
    return density / (r0 * r0) * m0;
}

double GridGeometry::densityLoss(double cgsDensityLoss)
{
    return cgsDensityLoss * (r0 * r0 * t0) / m0;
}

double GridGeometry::cgsDensityLoss(double densityLoss)
{
    return densityLoss / (r0 * r0 * t0) * m0;
}

double GridGeometry::viscousConstant(double constant)
{
    return constant / sqrt(r0) * t0;
}

double GridGeometry::cgsViscousConstant(double constant)
{
    return constant * sqrt(r0) / t0;
}

double GridGeometry::viscosity(double cgsViscosity)
{
    return cgsViscosity * t0 / (r0 * r0);
}

double GridGeometry::cgsViscosity(double viscosity)
{
    return viscosity / t0 * r0 * r0;
}

double GridGeometry::simTime(double cgsTime)
{
    return cgsTime / t0;
}

double GridGeometry::cgsTime(double simTime)
{
    return simTime * t0;
}