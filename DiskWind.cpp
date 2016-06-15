//
// Created by Peter Rodenkirch on 19.04.16.
//

#include "DiskWind.h"
#include <cstdlib>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <iterator>
#include <mpi.h>

DiskWind::DiskWind()
{
    NGrid = 2;
    frame = 0;
    frameStride = 1;
    outputFrame = 0;
    data = (double *)malloc(NGrid * sizeof(double));
}

DiskWind::DiskWind(int ncells)
{
    NGrid = ncells;
    frame = 0;
    frameStride = 1;
    outputFrame = 0;
    data = (double *)malloc(NGrid * sizeof(double));
}

DiskWind::~DiskWind()
{
    free(data);
    delete(g);
}

void DiskWind::writeFrame()
{
    std::cout << "Writing frame " << outputFrame << " (step " << frame << ", " << dt/year * frame << "yr)" << std::endl;
    std::vector<std::string> stringOutput;

    std::stringstream headerStream;
    headerStream << frame << " " << dt/year * frame;
    stringOutput.push_back(headerStream.str());

    for (int i = 0; i < NGrid; i++)
    {
        std::stringstream ss;
        double pos = g->convertIndexToPosition(i);
        ss << pos << " " << data[i] / pos;
        stringOutput.push_back(ss.str());
    }
    std::ostringstream tempStream;
    tempStream << "frame" << outputFrame << ".dat";
    std::ofstream outputFile(tempStream.str());
    std::ostream_iterator<std::string> output_iterator(outputFile, "\n");
    std::copy(stringOutput.begin(), stringOutput.end(), output_iterator);
    outputFrame++;
}

void DiskWind::computedx()
{
    dx = g->convertIndexToPosition(1) - g->convertIndexToPosition(0);  //Determining the smallest dx in data
    if (NGrid > 2) {
        for (int i = 2; i < NGrid; i++) {
            double temp = g->convertIndexToPosition(i) - g->convertIndexToPosition(i - 1);
            if (temp < dx) {
                dx = temp;
            }
        }
    }
}

void DiskWind::computedt()
{
    double x = g->convertIndexToPosition(1) - g->convertIndexToPosition(0);
    dt = 0.25 * pow(x, 1.5)/viscousConstant;

    for (int i = 2; i < NGrid; i++) {
        x = pow((g->convertIndexToPosition(i) - g->convertIndexToPosition(i - 1)), 1.5);
        double temp = 0.25 * x / viscousConstant;
        if (temp < dt) {
            dt = temp;
        }
    }

    std::cout << "Timestep: " << dt << std::endl;
}


double DiskWind::densityLossAtRadius(double r)
{
    if (r * au >= photoRadius) {
        double soundspeed = sqrt(kb * 1e4 / (2.3 * mp));
        double numberDensity = 5.7e4 * sqrt(luminosity / normLuminosity) * pow(photoRadius * 1e-14, -1.5) * pow(r * au / photoRadius, -2.5);
        double densityLoss = 2 * soundspeed * numberDensity * mp;
        return densityLoss;
    } else {
        return 0.0;
    }
}

double DiskWind::leverArmAtRadius(double r)
{
    return 1.0;
}

double DiskWind::constantLeverArm()
{
    return leverArm;
}


double DiskWind::computeFluxDiff(const int i)
{
    double rPlus = g->convertIndexToPosition(i + 1.0);
    double rPlusHalf = g->convertIndexToPosition(i + 0.5);
    double r = g->convertIndexToPosition(i);
    double rMinusHalf = g->convertIndexToPosition(i - 0.5);
    double rMinus = g->convertIndexToPosition(i - 1.0);

    double yPlus = data[i + 1];
    double y = data[i];
    double yMinus = data[i - 1];



    double Fleft, Fright;

    if (i == 0) {
        yMinus = y;
    } else if (i == NGrid - 1) {
        yPlus = 0.0;
    }

    Fright = viscousConstant * (0.25 * (y + yPlus) + rPlusHalf * (yPlus - y) / (rPlus - r)) - 2 * (constantLeverArm() - 1) * rPlusHalf * rPlusHalf * densityLossAtRadius(rPlusHalf);
    Fleft = viscousConstant * (0.25 * (y + yMinus) + rMinusHalf * (y - yMinus) / (r - rMinus)) - 2 * (constantLeverArm() - 1) * rMinusHalf * rMinusHalf * densityLossAtRadius(rMinusHalf);
    return Fright - Fleft;
}

void DiskWind::step()
{

    // Determine MPI Data
    const int root_process = 0;
    const int current_id = MPI::COMM_WORLD.Get_rank();
    const int processors = MPI::COMM_WORLD.Get_size();
    const int chunksize = NGrid / processors;

    // Root process controls the remaining cores
    if (current_id == root_process)
    {
        // temporary storage for all chunks
        double *tempData = (double *)malloc((chunksize + 1) * sizeof(double));

        // boundaries are exchanged
        MPI::COMM_WORLD.Send(&data[chunksize - 1], 1, MPI_DOUBLE, root_process + 1, 0);
        MPI::COMM_WORLD.Recv(&data[chunksize], 1, MPI_DOUBLE, root_process + 1, 1);

        // root process computes his chunk
        for (int i = 0; i < chunksize; i++) {
            double fluxDiff = computeFluxDiff(i);
            double dr = g->convertIndexToPosition(i+0.5) - g->convertIndexToPosition(i-0.5);
            double r  = g->convertIndexToPosition(i);
            tempData[i] = data[i] + dt * (fluxDiff / dr - densityLossAtRadius(r) * r);
        }

        for (int i = 0; i < chunksize; i++) {
            double r = g->convertIndexToPosition(i);
            if (tempData[i] / r < floorDensity)
            {
                data[i] = floorDensity * r;
            } else {
                data[i] = tempData[i];
            }
        }

        free(tempData);

        frame++;
        if (frame % frameStride == 0) {
            double *buffer = (double *)malloc((chunksize + 1) * sizeof(double));
            MPI::COMM_WORLD.Recv(buffer, chunksize + 1, MPI_DOUBLE, root_process + 1, 2);

            for (int i = chunksize; i < NGrid; i++) {
                data[i] = buffer[i - chunksize];
            }

            free(buffer);
            writeFrame();
        }
    } else {

        int minIndex = current_id * chunksize;
        int maxIndex = minIndex + chunksize;

        MPI::COMM_WORLD.Recv(&data[minIndex - 1], 1, MPI_DOUBLE, root_process, 0);
        MPI::COMM_WORLD.Send(&data[minIndex], 1, MPI_DOUBLE, root_process, 1);

        double *tempData = (double *)malloc((chunksize + 1) * sizeof(double));

        for (int i = minIndex; i < maxIndex; i++) {
            double fluxDiff = computeFluxDiff(i);
            double dr = g->convertIndexToPosition(i+0.5) - g->convertIndexToPosition(i-0.5);
            double r  = g->convertIndexToPosition(i);
            tempData[i - minIndex] = data[i] + dt * (fluxDiff / dr - densityLossAtRadius(r) * r);
        }

        for (int i = minIndex; i < maxIndex; i++) {
            double r = g->convertIndexToPosition(i);
            if (tempData[i - minIndex] / r < floorDensity)
            {
                data[i] = floorDensity * r;
                tempData[i - minIndex] = floorDensity * r;
            } else {
                data[i] = tempData[i - minIndex];
            }
        }

        frame++;
        if (frame % frameStride == 0) {
            MPI::COMM_WORLD.Send(tempData, chunksize + 1, MPI_DOUBLE, root_process, 2);
        }

        free(tempData);
    }
}

void DiskWind::initWithRestartData(int lastFrame)
{
    std::cout << "Initializing with frame " << lastFrame << std::endl;
    std::stringstream fileStream;
    fileStream << "frame" << lastFrame << ".dat";
    std::ifstream input(fileStream.str());

    std::string line;
    int i = 0;
    while (std::getline(input, line))
    {
        std::istringstream iss(line);
        if (i == 0) {
            if(!(iss >> lastFrameDetail >> currentTime)) {
                std::cout << "Error while parsing restart file (first line)" << std::endl;
                break;
            }
            i++;
        } else {
            double x;
            if (!(iss >> x >> data[i-1]))
            {
                std::cout << "Error while parsing restart file" << std::endl;
                break;
            } else {
                data[i-1] *= g->convertIndexToPosition(i - 1);
            }
            i++;
        }
    }
}


void DiskWind::setParameters(double a, double mass, double lum, double rg, double lever, int NFrames, GridGeometry *geometry)
{
    alpha = a;
    M = mass;
    luminosity = lum;
    photoRadius = rg;
    leverArm = lever;
    maxFrames = NFrames;
    g = geometry;

    viscousConstant = 3 * alpha * kb * T0 / (sqrt(au) * 2.3 * mp * sqrt(G * M));
}


void DiskWind::initWithHCGADensityDistribution(double initialDiskMass, double radialScaleFactor, double floor)
{
    std::cout << "Initializing grid of size " << NGrid << " with HCGA distribution" << std::endl;

    floorDensity = floor;
    radialScale = radialScaleFactor;
    diskMass = initialDiskMass;

    for (int i = 0; i < NGrid; i++) {
        double x = g->convertIndexToPosition(i);
        data[i] = initialDiskMass / (2 * M_PI * au * au * radialScaleFactor * x) * exp(-x / radialScaleFactor) * x;
        if (data[i] / x < floorDensity) {
            data[i] = floorDensity * x;
        }
    }

    writeFrame();
    computedx();
    computedt();
}


void DiskWind::restartSimulation(int lastFrame, int years)
{
    double NSteps = (double)years * year / dt;
    frameStride = (int)(NSteps / (double)maxFrames);

    frame = lastFrameDetail;
    outputFrame = lastFrame + 1;

    for (int i = 0; dt/year * i < years; i++) {
        step();
    }
}

void DiskWind::runDispersalAnalysis(int timeLimit, std::vector<double>* leverArms)
{
    double NSteps = (double)timeLimit * year / dt;
    frameStride = (int)(NSteps / (double)maxFrames);

    std::ofstream dispersalFile;
    dispersalFile.open("dispersal.dat");

    for (unsigned int i = 0; i < leverArms->size(); i++) {
        leverArm = leverArms->at(i);
        frame = 0;
        outputFrame = 0;
        initWithHCGADensityDistribution(diskMass, radialScale, floorDensity);

        for (int k = 0; dt/year * k < timeLimit; k++) {
            step();
            bool dispersed = false;
            for (int j = 60; j < 75; j++) {
                if (data[j] / g->convertIndexToPosition(j) <= 2*floorDensity)
                {
                    std::cout << "For lambda = " << leverArms->at(i) << ", disk dispersal is reached after " << dt/year * k << " years." << std::endl;
                    dispersalFile << leverArms->at(i) << " " << dt/year * k << std::endl;
                    dispersed = true;
                }
            }
            if (dispersed) {
                break;
            }
        }
    }

    dispersalFile.close();
}

void DiskWind::runSimulation(int years)
{
    double NSteps = (double)years * year / dt;
    frameStride = (int)(NSteps / (double)maxFrames);

    for (int i = 0; dt/year * i < years; i++) {
        step();
    }
}