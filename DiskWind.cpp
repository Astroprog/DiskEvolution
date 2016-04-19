//
// Created by Peter Rodenkirch on 19.04.16.
//

#include "DiskWind.h"
#include <cstdlib>
#include <cmath>



DiskWind::DiskWind()
{
    NGrid = 2;
    frame = 0;
    frameStride = 1;
    data = (Point *)malloc(NGrid * sizeof(Point));
}

DiskWind::DiskWind(int ncells)
{
    NGrid = ncells;
    frame = 0;
    frameStride = 1;
    data = (Point *)malloc(NGrid * sizeof(Point));
}

double DiskWind::massLossAtRadius(double r)
{

}


void DiskWind::step()
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