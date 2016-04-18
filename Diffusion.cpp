//
// Created by Peter Rodenkirch on 12.04.16.
//

#include <cstdlib>
#include <iostream>
#include <math.h>
#include "Diffusion.h"


Diffusion::Diffusion()
{
    NGrid = 2;
    frame = 0;
    maxFrames = 1;
    data = (Point *)malloc(NGrid * sizeof(Point));
}

Diffusion::Diffusion(int ncells)
{
    NGrid = ncells;
    frame = 0;
    maxFrames = 1;
    data = (Point *)malloc(NGrid * sizeof(Point));
}

void Diffusion::initWithGaussian(double mu, double sigma)
{
    std::cout << "Initializing grid with size " << NGrid << std::endl;

    for (int i = 0; i < NGrid; i++) {
        data[i].x = g->convertIndexToPosition(i);
        data[i].y = GridGeometry::gaussian(mu, sigma, data[i].x);
    }

    std::cout << data[100].y << std::endl;

    writeFrame();
}


void Diffusion::step()
{
    double *tempData = (double *)malloc(NGrid * sizeof(double));

    for (int i = 0; i < NGrid; i++) {
        tempData[i] = data[i].y - dt / (g->convertIndexToPosition(i + 0.5) - g->convertIndexToPosition(i - 0.5)) * computeFluxDiff(i);
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


double Diffusion::diffValue(double x)
{
    return 1.0;
}

double Diffusion::computeFluxDiff(int i)
{
    double Dleft = diffValue(g->convertIndexToPosition(i - 0.5));
    double Dright = diffValue(g->convertIndexToPosition(i + 0.5));


//Boundary Conditions (0)
    double Fleft = 0.0;
    double Fright = 0.0;
    if (i == 0) {
        Fleft = 0.0;
    } else if (i == NGrid - 1) {
        Fright = 0.0;
    } else {
        Fleft = -Dleft * (data[i].y - data[i - 1].y) / (data[i].x - data[i - 1].x);
        Fright = -Dright * (data[i + 1].y - data[i].y) / (data[i + 1].x - data[i].x);
    }

    return Fright - Fleft;
}

void Diffusion::computedt()
{
    dt = 0.1 * pow(dx, 2) / diffValue(0);
}