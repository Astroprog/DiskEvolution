//
// Created by Peter Rodenkirch on 13.04.16.
//

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iterator>
#include "ViscousDisk.h"


ViscousDisk::ViscousDisk()
{
    NGrid = 2;
    frame = 0;
    frameStride = 1;
    currentTime = 0.0;
    data = (Point *)malloc(NGrid * sizeof(Point));
}

ViscousDisk::ViscousDisk(int ncells)
{
    NGrid = ncells;
    frame = 0;
    frameStride = 1;
    currentTime = 0.0;
    data = (Point *)malloc(NGrid * sizeof(Point));
}



void ViscousDisk::step()
{
    double *tempData = (double *)malloc(NGrid * sizeof(double));

    double c = 3 * alpha * kb * T0 / (sqrt(au) * 2.3 * mp * sqrt(G * M));

    for (int i = 0; i < NGrid; i++) {
        tempData[i] = data[i].y - c * dt / (g->convertIndexToPosition(i + 0.5) - g->convertIndexToPosition(i - 0.5)) * computeFluxDiff(i);
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

double ViscousDisk::computeFluxDiff(int i)
{
    double rPlus = g->convertIndexToPosition(i + 1.0);
    double rPlusHalf = g->convertIndexToPosition(i + 0.5);
    double r = g->convertIndexToPosition(i);
    double rMinusHalf = g->convertIndexToPosition(i - 0.5);
    double rMinus = g->convertIndexToPosition(i - 1.0);

    double yPlus = data[i + 1].y;
    double y = data[i].y;
    double yMinus = data[i - 1].y;



    double Fleft, Fright;

    if (i == 0) {
        yMinus = y;
    } else if (i == NGrid - 1) {
        yPlus = 0.0;
    }

    Fright = -0.25 * (y + yPlus) - rPlusHalf * (yPlus - y) / (rPlus - r);
    Fleft = -0.25 * (y + yMinus) - rMinusHalf * (y - yMinus) / (r - rMinus);
    return Fright - Fleft;
}



