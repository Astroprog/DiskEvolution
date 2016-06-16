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
public:
    GridGeometry();
    GridGeometry(double innerBound, double outerBound, int NCells, bool useLogscale);
    double convertIndexToPosition(double i)
    {
        if (logscale) {
            return pow(10, scaleFactor * i + offset);
        } else {
            return i * scaleFactor + offset;
        }
    }

    double convertPositionToIndex(double r)
    {
        if (logscale) {
            return (log10(r) - offset) / scaleFactor;
        } else {
            return (r - offset) / scaleFactor;
        }
    }

    static double gaussian(double mu, double sigma, double x);
};


#endif //DIFFUSION_GRIDGEOMETRY_H
