//
// DiskWind.cpp
// Created by Peter Rodenkirch on 19.04.16.
//

#include "DiskWind.h"
#include <iostream>
#include <iomanip>
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
    accumulatedMassLossLeft = 0.0;
    accumulatedMassLossRight = 0.0;
    accumulatedWindLoss = 0.0;

    const int processors = MPI::COMM_WORLD.Get_size();
    const int chunksize = NGrid / processors;

    initialDensity = new double[NGrid];
    tempData = new double[chunksize + 1];
    windloss = new double[chunksize + 1];
    flux = new double[chunksize + 1];
}

DiskWind::DiskWind(int ncells)
{
    NGrid = ncells;
    frame = 0;
    frameStride = 1;
    outputFrame = 0;
    accumulatedMassLossLeft = 0.0;
    accumulatedMassLossRight = 0.0;
    accumulatedWindLoss = 0.0;

    const int processors = MPI::COMM_WORLD.Get_size();
    const int chunksize = NGrid / processors;

    initialDensity = new double[NGrid];
    tempData = new double[chunksize + 1];
    windloss = new double[NGrid + 1];
    flux = new double[NGrid + 1];
}

DiskWind::~DiskWind()
{
    delete[](initialDensity);
    delete[](tempData);
    delete[](windloss);
    delete[](flux);
}


// Writes the current state to a file
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
        ss << pos << " " << data[i].y / pos << " " << sqrt(data[i].B2) << " " << leverArmAtCell(i, densityLossAtRadius(g->convertIndexToPosition(i), i)) << " " << data[i].mdot;
        stringOutput.push_back(ss.str());
    }
    std::ostringstream tempStream;
    tempStream << "frame" << outputFrame << ".dat";
    std::ofstream outputFile(tempStream.str());
    std::ostream_iterator<std::string> output_iterator(outputFile, "\n");
    std::copy(stringOutput.begin(), stringOutput.end(), output_iterator);
    outputFrame++;
}


// Computes the smallest interval of space for the calculation of the time step
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


// Computes a stable time step
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


// Returns the density loss, computed by the photoevaporation model.
// First corrections, preventing the wind to carry more mass away than available.
double DiskWind::densityLossAtRadius(double r, int i)
{      
    double densityLoss = photoevModel->getLossAtRadius(r, data[i].y, dt);
    return densityLoss;
}


// Returns the magnetic lever arm, computed with the underlying magnetic field.
double DiskWind::leverArmAtCell(double i, double currentWindloss)
{
    if (currentWindloss > 0.0)
    {
        double B2Averaged = 0.5 * (data[(int)floor(i)].B2 + data[(int)ceil(i)].B2);
        if (ceil(i) == NGrid) {
            B2Averaged = data[NGrid - 1].B2;
        }

        double mu = 4.0 * M_PI * currentWindloss * sqrt(G * M / (g->convertIndexToPosition(i) * au)) / B2Averaged;
        return 1.5 * (1.0 + pow(mu, -2.0/3.0));
    } else {
        return 1.0;
    }
}

// Used for the simulation runs with a spacially constant magnetic lever arm.
double DiskWind::constantLeverArm()
{
    return leverArm;
}


// Computes the new fluxes for the current step, following the differential equations of the chosen model.
// Necessary corrections for low surface densities are made here.
void DiskWind::computeFluxes(int minIndex, int maxIndex)
{
    for (int i = minIndex; i < maxIndex; i++) {
        double r = g->convertIndexToPosition(i);
        double rMinusHalf = g->convertIndexToPosition(i - 0.5);
        double rMinus = g->convertIndexToPosition(i - 1.0);
        double dr = g->convertIndexToPosition(i + 0.5) - rMinusHalf;
        double drMinus = rMinusHalf - g->convertIndexToPosition(i - 1.5);

        double y;
        double yMinus;

        double currentWindloss;
        double currentWindlossMinus;

        if (i == 0) {
            y = data[i].y;
            yMinus = y;
            currentWindlossMinus = 0.0;
            currentWindloss = windloss[i];
        } else if (i == NGrid) {
            yMinus = data[i - 1].y;
            y = 0.0;
            currentWindlossMinus = windloss[i - 1];
            currentWindloss = 0.0;
        } else {
            y = data[i].y;
            yMinus = data[i - 1].y;
            currentWindloss = windloss[i];
            currentWindlossMinus = windloss[i - 1];
        }

        double yMinusHalf = 0.5 * (y + yMinus);
        double currentWindlossMinusHalf = 0.5 * (currentWindloss + currentWindlossMinus);

        double viscousTerm = viscousConstant * (0.5 * yMinusHalf + rMinusHalf * (y - yMinus) / (r - rMinus));

        double magneticTerm = 0.0;
        if (constLambda == true) {
            magneticTerm = 2 * (constantLeverArm() - 1) * rMinusHalf * rMinusHalf * currentWindlossMinusHalf;
        } else {
            magneticTerm = 2 * (leverArmAtCell(i - 0.5, currentWindlossMinusHalf) - 1) * rMinusHalf * rMinusHalf * currentWindlossMinusHalf;
        }

        double currentFlux = viscousTerm + magneticTerm;

        // Low surface density corrections
        if (currentFlux * dt / dr >= y - currentWindloss * r * dt) {
            currentFlux = (y - currentWindloss * r * dt) * dr / dt;
        } else if (-currentFlux * dt / drMinus >= yMinus - currentWindlossMinus * rMinus * dt) {
            currentFlux = -(yMinus - currentWindlossMinus * rMinus * dt) * drMinus / dt;
        }

        if (i == NGrid) {   // No inflow possible
            if (currentFlux > 0.0) {
                currentFlux = 0.0;
            }
        } else if (i == 0) {
            if (currentFlux < 0.0) {
                currentFlux = 0.0;
            }
        }

        flux[i] = currentFlux;
    }
}


// Advances the system by one time step
void DiskWind::step()
{
    // Determine MPI Data
    const int root_process = 0;
    const int current_id = MPI::COMM_WORLD.Get_rank();
    const int processors = MPI::COMM_WORLD.Get_size();
    const int chunksize = NGrid / processors;


    //Single processor

    if (processors == 1) {

        // Compute wind losses

        for (int i = 0; i < NGrid; i++) {
            windloss[i] = densityLossAtRadius(g->convertIndexToPosition(i), i);
        }

        // Compute fluxes

        computeFluxes(0, NGrid + 1);

        // Iterate over the array and generate the updated values.

        for (int i = 0; i < chunksize; i++) {
            double rPlusHalf = g->convertIndexToPosition(i + 0.5);
            double rMinusHalf = g->convertIndexToPosition(i - 0.5);
            double dr = rPlusHalf - rMinusHalf;
            double r  = g->convertIndexToPosition(i);
            double currentArea = M_PI * au * au * (pow(rPlusHalf, 2) - pow(rMinusHalf, 2));

            tempData[i] = data[i].y + dt * (flux[i + 1] - flux[i]) / dr;

            tempData[i] -= windloss[i] * r * dt;

            // Used for accuracy checks
            accumulatedWindLoss += windloss[i] * currentArea * dt;

            if (i == 0) {
                accumulatedMassLossLeft += flux[i] / (dr * data[i].x) * dt * currentArea;
            } else if (i == NGrid - 1) {
                accumulatedMassLossRight += -flux[i + 1] / (dr * data[i].x) * dt * currentArea;
            }

            // Accretion rate
            data[i].mdot = year * flux[i] / (dr * data[i].x) * currentArea;
        }

        // Updates with the new values

        for (int i = 0; i < chunksize; i++) {
            if (!constB) {
                data[i].B2 = getUpdatedMagneticFluxDensityAtCell(i);
            }

            data[i].y = tempData[i];
            if (data[i].y < 0.0) {
                data[i].y = 0.0;
            }
        }

        frame++;
        if (frame % frameStride == 0) {
            determineDiskExtent();
            writeFrame();
        }
    } else if (current_id == root_process)     // Root process controls the remaining cores when multiprocessing is used
    {

        // boundaries are exchanged
        for (int proc = root_process + 1; proc <= processors-1; proc++) {
            if (proc == processors - 1) {
                MPI::COMM_WORLD.Send(&data[proc * chunksize - 1].y, 1, MPI_DOUBLE, proc, initSendLow); // sending lower boundary
                MPI::COMM_WORLD.Send(&data[proc * chunksize - 1].B2, 1, MPI_DOUBLE, proc, initSendLow); // sending lower boundary
            } else {
                MPI::COMM_WORLD.Send(&data[proc * chunksize - 1].y, 1, MPI_DOUBLE, proc, initSendLow); // sending lower boundary
                MPI::COMM_WORLD.Send(&data[(proc + 1) * chunksize].y, 1, MPI_DOUBLE, proc, initSendHigh); // sending upper boundary
                MPI::COMM_WORLD.Send(&data[proc * chunksize - 1].B2, 1, MPI_DOUBLE, proc, initSendLow); // sending lower boundary
                MPI::COMM_WORLD.Send(&data[(proc + 1) * chunksize].B2, 1, MPI_DOUBLE, proc, initSendHigh); // sending upper boundary
            }
        }


        // root process computes his chunk

        for (int i = 0; i < chunksize + 1; i++) {
            windloss[i] = densityLossAtRadius(g->convertIndexToPosition(i), i);
        }

        computeFluxes(0, chunksize + 1);

        for (int i = 0; i < chunksize; i++) {
            double rPlusHalf = g->convertIndexToPosition(i + 0.5);
            double rMinusHalf = g->convertIndexToPosition(i - 0.5);
            double dr = rPlusHalf - rMinusHalf;
            double r  = g->convertIndexToPosition(i);
            double currentArea = M_PI * au * au * (pow(rPlusHalf, 2) - pow(rMinusHalf, 2));

            tempData[i] = data[i].y + dt * (flux[i + 1] - flux[i]) / dr;

            tempData[i] -= windloss[i] * r * dt;
            accumulatedWindLoss += windloss[i] * currentArea * dt;
            data[i].mdot = year * flux[i] / (dr * data[i].x) * currentArea;

            if (i == 0) {
                accumulatedMassLossLeft += flux[i] / (dr * data[i].x) * dt * currentArea;
            } else if (i == NGrid - 1) {
                accumulatedMassLossRight += -flux[i + 1] / (dr * data[i].x) * dt * currentArea;
            }
        }

        for (int i = 0; i < chunksize; i++) {
            if (!constB) {
                data[i].B2 = getUpdatedMagneticFluxDensityAtCell(i);
            }

            data[i].y = tempData[i];
            if (data[i].y < 0.0) {
                data[i].y = 0.0;
            }
        }



        for (int proc = root_process + 1; proc <= processors-1; proc++) {
            MPI::COMM_WORLD.Recv(&data[proc * chunksize].y, 1, MPI_DOUBLE, proc, finalRecvLow);
            MPI::COMM_WORLD.Recv(&data[proc * chunksize].B2, 1, MPI_DOUBLE, proc, finalRecvLow);
            if (proc < processors - 1) {
                MPI::COMM_WORLD.Recv(&data[(proc + 1) * chunksize - 1].y, 1, MPI_DOUBLE, proc, finalRecvHigh);
                MPI::COMM_WORLD.Recv(&data[(proc + 1) * chunksize - 1].B2, 1, MPI_DOUBLE, proc, finalRecvHigh);
            }
        }


        frame++;
        if (frame % frameStride == 0) {
            for (int proc = root_process + 1; proc <= processors-1; proc++) {

                double *buffer = new double[chunksize];
                double *magneticBuffer = new double[chunksize];
                double *accretionBuffer = new double[chunksize];

                MPI::COMM_WORLD.Recv(buffer, chunksize, MPI_DOUBLE, proc, frameRecv);
                MPI::COMM_WORLD.Recv(magneticBuffer, chunksize, MPI_DOUBLE, proc, frameBRecv);
                MPI::COMM_WORLD.Recv(accretionBuffer, chunksize, MPI_DOUBLE, proc, frameMDotRecv);

                for (int i = proc * chunksize; i < (proc + 1) * chunksize; i++) {
                    data[i].y = buffer[i - proc * chunksize];
                    data[i].B2 = magneticBuffer[i - proc * chunksize];
                    data[i].mdot = accretionBuffer[i - proc * chunksize];
                }

                delete[](buffer);
                delete[](magneticBuffer);
                delete[](accretionBuffer);
            }
            determineDiskExtent();
            writeFrame();
        }



    } else {  // additional processors
        int minIndex = current_id * chunksize;
        int maxIndex = minIndex + chunksize;

        if (current_id == processors - 1) {
            MPI::COMM_WORLD.Recv(&data[minIndex - 1].y, 1, MPI_DOUBLE, root_process, initSendLow); // receiving lower boundary
            MPI::COMM_WORLD.Recv(&data[minIndex - 1].B2, 1, MPI_DOUBLE, root_process, initSendLow); // receiving lower boundary
        } else {
            MPI::COMM_WORLD.Recv(&data[minIndex - 1].y, 1, MPI_DOUBLE, root_process, initSendLow); // receiving lower boundary
            MPI::COMM_WORLD.Recv(&data[maxIndex].y, 1, MPI_DOUBLE, root_process, initSendHigh); // receiving upper boundary
            MPI::COMM_WORLD.Recv(&data[minIndex - 1].B2, 1, MPI_DOUBLE, root_process, initSendLow); // receiving lower boundary
            MPI::COMM_WORLD.Recv(&data[maxIndex].B2, 1, MPI_DOUBLE, root_process, initSendHigh); // receiving upper boundary
        }

        for (int i = minIndex - 1; i < maxIndex + 1; i++) {
            windloss[i] = densityLossAtRadius(g->convertIndexToPosition(i), i);
        }


        computeFluxes(minIndex, maxIndex + 1);

        for (int i = minIndex; i < maxIndex; i++) {
            double rPlusHalf = g->convertIndexToPosition(i + 0.5);
            double rMinusHalf = g->convertIndexToPosition(i - 0.5);
            double dr = rPlusHalf - rMinusHalf;
            double r  = g->convertIndexToPosition(i);
            double currentArea = M_PI * au * au * (pow(rPlusHalf, 2) - pow(rMinusHalf, 2));

            tempData[i - minIndex] = data[i].y + dt * (flux[i + 1] - flux[i]) / dr;

            tempData[i - minIndex] -= windloss[i] * r * dt;
            accumulatedWindLoss += windloss[i] * currentArea * dt;

            if (i == 0) {
                accumulatedMassLossLeft += flux[i] / (dr * data[i].x) * dt * currentArea;
            } else if (i == NGrid - 1) {
                accumulatedMassLossRight += -flux[i + 1] / (dr * data[i].x) * dt * currentArea;
            }

            data[i].mdot = year * flux[i] / (dr * data[i].x) * currentArea;
        }

        for (int i = minIndex; i < maxIndex; i++) {

            if (!constB) {
                data[i].B2 = getUpdatedMagneticFluxDensityAtCell(i);
            }

            data[i].y = tempData[i - minIndex];
            if (data[i].y < 0.0) {
                data[i].y = 0.0;
            }
        }



        MPI::COMM_WORLD.Send(&data[minIndex].y, 1, MPI_DOUBLE, root_process, finalRecvLow);
        MPI::COMM_WORLD.Send(&data[minIndex].B2, 1, MPI_DOUBLE, root_process, finalRecvLow);

        if (current_id < processors - 1) {
            MPI::COMM_WORLD.Send(&data[maxIndex - 1].y, 1, MPI_DOUBLE, root_process, finalRecvHigh);
            MPI::COMM_WORLD.Send(&data[maxIndex - 1].B2, 1, MPI_DOUBLE, root_process, finalRecvHigh);
        }

        frame++;
        if (frame % frameStride == 0) {
            double *bfield = new double[chunksize];
            double *accretionRate = new double[chunksize];
            for (int i = minIndex; i < maxIndex; i++) {
                bfield[i - minIndex] = data[i].B2;
                accretionRate[i - minIndex] = data[i].mdot;
            }
            MPI::COMM_WORLD.Send(tempData, chunksize, MPI_DOUBLE, root_process, frameRecv);
            MPI::COMM_WORLD.Send(bfield, chunksize, MPI_DOUBLE, root_process, frameBRecv);
            MPI::COMM_WORLD.Send(accretionRate, chunksize, MPI_DOUBLE, root_process, frameMDotRecv);

            delete[](bfield);
            delete[](accretionRate);
        }
    }
}

// Returns the updated magnetic flux density for a given cell i.

double DiskWind::getUpdatedMagneticFluxDensityAtCell(int i)
{
    double soundSpeed = sqrt(kb * T0 / (2.3 * mp * sqrt(data[i].x)));
    double omega = sqrt(G*M/pow(data[i].x * au, 3));
    double scaleHeight = soundSpeed / omega;
    double currentDensity = data[i].y / data[i].x;
    double midplaneDensity = currentDensity / (sqrt(2 * M_PI) * scaleHeight);
    if (fluxFreezing) {
        return 8 * M_PI * midplaneDensity * currentDensity / initialDensity[i] * pow(soundSpeed, 2) / plasma;
    } else if (perfectFluxFreezing) {
        return midplaneDensity * midplaneDensity / freezingPlasma;
    } else {
        return 8 * M_PI * midplaneDensity * soundSpeed * soundSpeed / plasma;
    }
}

// Computes the extent of the disk. Used for density bump detection
void DiskWind::determineDiskExtent()
{
    for (int i = NGrid - 1; i >= 0; i--) {
        if (data[i].y / data[i].x > 2 * floorDensity) {
            currentDiskExtent = g->convertIndexToPosition(i);
            break;
        }
    }
}


// Computes the disk mass
double DiskWind::computeDiskMass()
{
    double mass = 0.0;
    for (int i = 0; i < NGrid; i++) {
        mass += data[i].y / data[i].x * M_PI * au * au * (pow(g->convertIndexToPosition(i + 0.5), 2) - pow(g->convertIndexToPosition(i - 0.5), 2));
    }
    return mass;
}


// Currently not working
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
            if (!(iss >> x >> data[i-1].y))
            {
                std::cout << "Error while parsing restart file" << std::endl;
                break;
            } else {
                data[i-1].y *= g->convertIndexToPosition(i - 1);
            }
            i++;
        }
    }
}


// Initializes the system with the chosen parameters from the Simulation class
void DiskWind::setParameters(double a, double mass, double lum, double rg, double lever, int NFrames,
                             GridGeometry *geometry, double plasmaParameter, bool constlambda,
                             bool constb, bool freezing, bool pfreezing)
{
    alpha = a;
    M = mass;
    luminosity = lum;
    photoRadius = rg;
    leverArm = lever;
    maxFrames = NFrames;
    g = geometry;
    plasma = plasmaParameter;
    constLambda = constlambda;
    constB = constb;
    fluxFreezing = freezing;
    perfectFluxFreezing = pfreezing;

    viscousConstant = 3 * alpha * kb * T0 / (sqrt(au) * 2.3 * mp * sqrt(G * M));
}

void DiskWind::setEvaporationModel(EvaporationModel* model) {
    photoevModel = model;
}


void DiskWind::initWithCustomDensityDistribution(double sigma0, double r0, double r_char, double slope, double floor)
{
    std::cout << "Initializing grid of size " << NGrid << " with custom distribution" << std::endl;

    computedx();
    computedt();

    floorDensity = floor;

    for (int i = 0; i < NGrid; i++) {
        data[i].x = g->convertIndexToPosition(i);
        data[i].y = sigma0 * pow(data[i].x / r0, -slope) * exp(-pow(data[i].x / r_char, 2.0 - slope)) * data[i].x;
        if (data[i].y / data[i].x < floorDensity) {
            data[i].y = floorDensity * data[i].x;
        }
        data[i].mdot = 0.0;

        double soundSpeed = sqrt(kb * T0 / (2.3 * mp * sqrt(data[i].x)));
        double scaleHeight = soundSpeed / sqrt(G*M/pow(data[i].x * au, 3));
        double midplaneDensity = data[i].y / data[i].x / (sqrt(2 * M_PI) * scaleHeight);
        double B2 = 8 * M_PI * midplaneDensity * soundSpeed * soundSpeed / plasma;
        data[i].B2 = B2;

        if (perfectFluxFreezing)
        {
            // At R_g AU the value of reference is set
            if (i == 68) {
                freezingPlasma = midplaneDensity * midplaneDensity / B2;
                std::cout << sqrt(B2) << std::endl;
                std::cout << "freezing Plasma: " << midplaneDensity / freezingPlasma << std::endl;
            }
        }

        if (fluxFreezing) {
            initialDensity[i] = data[i].y / data[i].x;
        }
    }

    if (perfectFluxFreezing)
    {
        for (int i = 0; i < NGrid; i++) {
            data[i].B2 = getUpdatedMagneticFluxDensityAtCell(i);
        }

    }

    determineDiskExtent();
    writeFrame();
}

// Initializes the density distribution
void DiskWind::initWithHCGADensityDistribution(double initialDiskMass, double radialScaleFactor, double floor)
{
    std::cout << "Initializing grid of size " << NGrid << " with HCGA distribution" << std::endl;

    computedx();
    computedt();

    floorDensity = floor;
    radialScale = radialScaleFactor;
    diskMass = initialDiskMass;

    for (int i = 0; i < NGrid; i++) {
        data[i].x = g->convertIndexToPosition(i);
        data[i].y = initialDiskMass / (2 * M_PI * au * au * radialScaleFactor * data[i].x) * exp(-data[i].x / radialScaleFactor) * data[i].x;
        if (data[i].y / data[i].x < floorDensity) {
            data[i].y = floorDensity * data[i].x;
        }
        data[i].mdot = 0.0;

        double soundSpeed = sqrt(kb * T0 / (2.3 * mp * sqrt(data[i].x)));
        double scaleHeight = soundSpeed / sqrt(G*M/pow(data[i].x * au, 3));
        double midplaneDensity = data[i].y / data[i].x / (sqrt(2 * M_PI) * scaleHeight);
        double B2 = 8 * M_PI * midplaneDensity * soundSpeed * soundSpeed / plasma;
        data[i].B2 = B2;

        if (perfectFluxFreezing)
        {
            // At R_g AU the value of reference is set
            if (i == 68) {
                freezingPlasma = midplaneDensity * midplaneDensity / B2;
                std::cout << sqrt(B2) << std::endl;
                std::cout << "freezing Plasma: " << midplaneDensity / freezingPlasma << std::endl;
            }
        }

        if (fluxFreezing) {
            initialDensity[i] = data[i].y / data[i].x;
        }
    }

    if (perfectFluxFreezing)
    {
        for (int i = 0; i < NGrid; i++) {
            data[i].B2 = getUpdatedMagneticFluxDensityAtCell(i);
        }

    }

    determineDiskExtent();
    writeFrame();
}


// Currently not working
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


// Triggers a dispersal analysis with a list of parameter values, given by the Simulation class
void DiskWind::runDispersalAnalysis(int timeLimit, std::vector<double>* parameters, const std::string parameterType)
{
    double NSteps = (double)timeLimit * year / dt;
    frameStride = (int)(NSteps / (double)maxFrames);

    const int root_process = 0;
    const int current_id = MPI::COMM_WORLD.Get_rank();
    const int processors = MPI::COMM_WORLD.Get_size();
    const int chunksize = NGrid / processors;

    std::ofstream dispersalFile;
    dispersalFile.open("dispersal.dat");

    const int minIndex = (int)g->convertPositionToIndex(5.0);
    determineDiskExtent();
    std::cout << "Extent: " << currentDiskExtent << std::endl;

    for (unsigned int i = 0; i < parameters->size(); i++) {

        if (parameterType == "lambda") {
            leverArm = parameters->at(i);
        } else if (parameterType == "plasma") {
            plasma = parameters->at(i);
        }

        frame = 0;
        outputFrame = 0;
        initWithHCGADensityDistribution(diskMass, radialScale, floorDensity);

        int maxIndex = (int)g->convertPositionToIndex(currentDiskExtent - 50.0);

        if (processors == 1) {

            for (int k = 0; dt / year * k < timeLimit; k++) {
                step();
                bool dispersed = false;
                for (int j = minIndex; j < maxIndex; j++) {
                    if (data[j].y / g->convertIndexToPosition(j) <= floorDensity) {
                        std::cout << "For " << parameterType << " = " << parameters->at(i) <<
                        ", disk dispersal is reached after " << dt / year * k << " years, at radius " << g->convertIndexToPosition(j) << std::endl;
                        dispersalFile << parameters->at(i) << " " << dt / year * k << " " << g->convertIndexToPosition(j) << std::endl;
                        dispersed = true;
                    }
                }

                if (dispersed) {
                    break;
                }
            }
        } else if (current_id == root_process) {

            for (int k = 0; dt / year * k < timeLimit; k++) {
                step();

                for (int proc = root_process + 1; proc <= processors-1; proc++) {
                    double *buffer = new double[chunksize];
                    MPI::COMM_WORLD.Recv(buffer, chunksize, MPI_DOUBLE, proc, frameRecv);

                    for (int j = proc * chunksize; j < (proc + 1) * chunksize; j++) {
                        data[j].y = buffer[j - proc * chunksize];
                    }
                    delete[](buffer);
                }

                bool dispersed = false;
                for (int j = minIndex; j < maxIndex; j++) {
                    if (data[j].y / g->convertIndexToPosition(j) <= 0.01 * data[j + 5].y / g->convertIndexToPosition(j + 5)) {
                        std::cout << "For " << parameterType << " = " << parameters->at(i) <<
                        ", disk dispersal is reached after " << dt / year * k << " years, at radius " << g->convertIndexToPosition(j) << std::endl;
                        dispersalFile << parameters->at(i) << " " << dt / year * k << " " << g->convertIndexToPosition(j) << std::endl;
                        dispersed = true;
                    }
                }

                for (int proc = root_process + 1; proc <= processors-1; proc++)
                {
                    MPI::COMM_WORLD.Send(&dispersed, 1, MPI_C_BOOL, proc, dispersalSend);
                }

                if (dispersed) {
                    break;
                }


            }
        } else {

            for (int k = 0; dt / year * k < timeLimit; k++) {
                step();
                MPI::COMM_WORLD.Send(tempData, chunksize, MPI_DOUBLE, root_process, frameRecv);
                bool dispersed = false;
                MPI::COMM_WORLD.Recv(&dispersed, 1, MPI_C_BOOL, root_process, dispersalSend);
                if (dispersed) {
                    break;
                }
            }
        }
    }

    dispersalFile.close();

}


// Runs the simulation after it has been initialized
void DiskWind::runSimulation(int years)
{
    double NSteps = (double)years * year / dt;
    frameStride = (int)(NSteps / (double)maxFrames);

    const int root_process = 0;
    const int current_id = MPI::COMM_WORLD.Get_rank();
    const int processors = MPI::COMM_WORLD.Get_size();
    const int chunksize = NGrid / processors;

    std::ofstream massFile;
    massFile.open("massloss.dat");

    std::ofstream checkFile;
    checkFile.open("check.dat");


    if (processors == 1) {

        for (int i = 0; dt/year * i < years; i++) {
            step();
            double currentMass = computeDiskMass();
            if (frame % frameStride == 0) {
                //std::cout << std::setprecision(16) << frame << ": " << currentMass << ", " << accumulatedMassLossLeft
                //          << ", " << accumulatedMassLossRight << ", " << accumulatedWindLoss << std::endl;
                //std::cout << std::setprecision(16) << "Total: " << currentMass + accumulatedMassLossLeft + accumulatedMassLossRight
                                                              //     + accumulatedWindLoss << std::endl << std::endl;
                massFile << std::setprecision(16) << dt/year * frame << "\t" << currentMass << "\t"
                         << accumulatedMassLossLeft << "\t" << accumulatedMassLossRight << "\t"
                         << accumulatedWindLoss << "\t" << currentMass + accumulatedMassLossLeft +
                        accumulatedMassLossRight + accumulatedWindLoss << std::endl;
                checkFile << dt/year * i << "\t" << data[40].y / data[40].x << "\t" << data[40].x << std::endl;

            }
        }
    } else if (current_id == root_process) {

        for (int i = 0; dt/year * i < years; i++) {
            step();

            double totalAccumulatedWindLoss = accumulatedWindLoss;

            for (int proc = root_process + 1; proc <= processors-1; proc++) {
                double *buffer = new double[chunksize];
                double windlossBuffer = 0.0;

                MPI::COMM_WORLD.Recv(buffer, chunksize, MPI_DOUBLE, proc, frameRecv);
                MPI::COMM_WORLD.Recv(&windlossBuffer, 1, MPI_DOUBLE, proc, totalWindloss);

                totalAccumulatedWindLoss += windlossBuffer;

                for (int j = proc * chunksize; j < (proc + 1) * chunksize; j++) {
                    data[j].y = buffer[j - proc * chunksize];
                }
                delete[](buffer);
            }

            MPI::COMM_WORLD.Recv(&accumulatedMassLossRight, 1, MPI_DOUBLE, processors - 1, massLossRight);

            double currentMass = computeDiskMass();
            if (frame % frameStride == 0) {
                //std::cout << std::setprecision(16) << frame << ": " << currentMass << ", " << accumulatedMassLossLeft
                //          << ", " << accumulatedMassLossRight << ", " << totalAccumulatedWindLoss << std::endl;
                //std::cout << std::setprecision(16) << "Total: "
                //          << currentMass + accumulatedMassLossLeft + accumulatedMassLossRight + totalAccumulatedWindLoss
                //          << std::endl << std::endl;
                massFile << std::setprecision(16) << dt/year * frame << "\t" << currentMass
                         << "\t" << accumulatedMassLossLeft << "\t"
                         << accumulatedMassLossRight << "\t" << totalAccumulatedWindLoss << "\t"
                         << currentMass + accumulatedMassLossLeft + accumulatedMassLossRight + totalAccumulatedWindLoss
                         << "\t" << data[0].mdot << std::endl;
                checkFile << dt/year * i << "\t" << data[40].y / data[40].x << "\t" << data[40].x << std::endl;
            }

        }
    } else {
        for (int i = 0; dt/year * i < years; i++) {
            step();
            MPI::COMM_WORLD.Send(tempData, chunksize, MPI_DOUBLE, root_process, frameRecv);
            MPI::COMM_WORLD.Send(&accumulatedWindLoss, 1, MPI_DOUBLE, root_process, totalWindloss);

            if (current_id == processors - 1)
            {
                MPI::COMM_WORLD.Send(&accumulatedMassLossRight, 1, MPI_DOUBLE, root_process, massLossRight);
            }
        }
    }

    massFile.close();
    checkFile.close();
}
