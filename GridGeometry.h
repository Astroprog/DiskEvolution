//
// Created by Peter Rodenkirch on 12.04.16.
//

#ifndef DIFFUSION_GRIDGEOMETRY_H
#define DIFFUSION_GRIDGEOMETRY_H

#include <cmath>

class GridGeometry {
private:
    double scaleFactor;
    double offset;
    double NGrid;
    bool logscale;

    double au = 1.49597871e13;
    double G = 6.67e-8;
    double solarMass = 1.98855e33;
    double kb = 1.38e-16;
    double mp = 1.67e-24;
    double T0 = 280;
    double year = 3.15e7;
    double m0;
    double t0;
    double r0;

public:
    GridGeometry();
    GridGeometry(double innerBound, double outerBound, int NCells, bool useLogscale, double m0Value, double r0Value, double t0Value);

    double mass(double cgsMass);
    double cgsMass(double mass);
    double radius(double cgsRadius);
    double cgsRadius(double radius);
    double density(double cgsDensity);
    double cgsDensity(double density);
    double densityLoss(double cgsDensityLoss);
    double cgsDensityLoss(double densityLoss);
    double viscousConstant(double constant);
    double cgsViscousConstant(double constant);
    double viscosity(double cgsViscosity);
    double cgsViscosity(double viscosity);
    double simTime(double cgsTime);
    double cgsTime(double simTime);

    double convertIndexToPosition(double i)
    {
        if (logscale) {
            return pow(10, scaleFactor * i + offset);
        } else {
            return i * scaleFactor + offset;
        }
    }

    static double gaussian(double mu, double sigma, double x);
};


#endif //DIFFUSION_GRIDGEOMETRY_H
