//
// Created by Peter Rodenkirch on 19.04.16.
//

#include "DiskWind.h"
#include <cstdlib>


DiskWind::DiskWind()
{
    NGrid = 2;
    frame = 0;
    frameStride = 1;
    outputFrame = 0;
    data = (Point *)malloc(NGrid * sizeof(Point));
}

DiskWind::DiskWind(int ncells)
{
    NGrid = ncells;
    frame = 0;
    frameStride = 1;
    outputFrame = 0;
    data = (Point *)malloc(NGrid * sizeof(Point));
}

double DiskWind::massLossAtRadius(double r, double rg)
{
    double rate = 1e25;
    if (r >= rg) {
        return rate * pow(r, -2.5);
    } else {
        return 0.0;
    }
}

double DiskWind::leverArmAtRadius(double r)
{
    return 1.2;
}

double DiskWind::computeFluxDiff(int i)
{
    double rPlus = g->convertIndexToPosition(i + 1.0);
    double rPlusHalf = g->convertIndexToPosition(i + 0.5);
    double r = g->convertIndexToPosition(i);
    double rMinusHalf = g->convertIndexToPosition(i - 0.5);
    double rMinus = g->convertIndexToPosition(i - 1.0);
    double dr = rPlusHalf - rMinusHalf;

    double yPlus = data[i + 1].y;
    double y = data[i].y;
    double yMinus = data[i - 1].y;



    double Fleft, Fright;

    if (i == 0) {
        yMinus = y;
    } else if (i == NGrid - 1) {
        yPlus = 0.0;
    }

    Fright = -(0.25 * (y + yPlus) + rPlusHalf * (yPlus - y) / (rPlus - r) + (leverArmAtRadius(r) - 1) * massLossAtRadius(r, 5.0) / (M_PI * au * year * (rPlus - r)) * r);
    Fleft = -(0.25 * (y + yMinus) + rMinusHalf * (y - yMinus) / (r - rMinus) + (leverArmAtRadius(r) - 1) * massLossAtRadius(r, 5.0) / (M_PI * au * year * (r - rMinus)) * r);
    return Fright - Fleft;
}

void DiskWind::step()
{
    double *tempData = (double *)malloc(NGrid * sizeof(double));

    double c = 3 * alpha * kb * T0 / (sqrt(au) * 2.3 * mp * sqrt(G * M));

    for (int i = 0; i < NGrid; i++) {
        double fluxDiff = computeFluxDiff(i);
        double dr = g->convertIndexToPosition(i+0.5) - g->convertIndexToPosition(i-0.5);
        data[i].mdot = 2 * M_PI * c * au * au * year / M * (0.5 * data[i].y + data[i].x * ((data[i+1].y - data[i].y)/2 - (data[i].y - data[i-1].y)/2) / dr);
        
        double densityLoss = massLossAtRadius(data[i].x, 5.0) / (2 * M_PI * au * au * dr) / year;
        tempData[i] = data[i].y - dt * (c / dr * fluxDiff + densityLoss);
    }

    for (int i = 0; i < NGrid; i++) {
        data[i].y = tempData[i];
    }

    free(tempData);

    frame++;
    if (frame % frameStride == 0) {
        writeFrame();
    }
}