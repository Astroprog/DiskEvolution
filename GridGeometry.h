//
// Created by Peter Rodenkirch on 12.04.16.
//

#ifndef DIFFUSION_GRIDGEOMETRY_H
#define DIFFUSION_GRIDGEOMETRY_H


class GridGeometry {
private:
    double scaleFactor;
    double offset;
    double NGrid;
    bool logscale;
public:
    GridGeometry();
    GridGeometry(double innerBound, double outerBound, int NCells, bool useLogscale);
    double convertIndexToPosition(double i);
    static double gaussian(double mu, double sigma, double x);
};


#endif //DIFFUSION_GRIDGEOMETRY_H
