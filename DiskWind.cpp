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

    data = new Point[NGrid];
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

    data = new Point[NGrid];
    tempData = new double[chunksize + 1];
    windloss = new double[NGrid + 1];
    flux = new double[NGrid + 1];
}

DiskWind::~DiskWind()
{
    delete[](data);
    delete[](tempData);
    delete[](windloss);
    delete[](flux);
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
        ss << pos << " " << data[i].y / pos << " " << sqrt(data[i].B2) << " " << leverArmAtCell(i, densityLossAtRadius(g->convertIndexToPosition(i), i));
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

double DiskWind::densityLossAtRadius(double r, int i)
{
    if (r * au >= photoRadius) {
        double soundspeed = sqrt(kb * 1e4 / (2.3 * mp));
        double numberDensity = 5.7e4 * sqrt(luminosity / normLuminosity) * pow(photoRadius * 1e-14, -1.5) * pow(r * au / photoRadius, -2.5);
        double densityLoss = 2 * soundspeed * numberDensity * mp;
        if (densityLoss * r * dt >= data[i].y) {
            densityLoss = data[i].y / (r * dt);
        }
        return densityLoss;
    } else {
        return 0.0;
    }
}

double DiskWind::leverArmAtCell(double i, double currentWindloss)
{
    if (currentWindloss != 0.0)
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

double DiskWind::constantLeverArm()
{
    return leverArm;
}



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
        double magneticTerm = 2 * (leverArmAtCell(i - 0.5, currentWindlossMinusHalf) - 1) * rMinusHalf * rMinusHalf * currentWindlossMinusHalf;

        double currentFlux = viscousTerm + magneticTerm;


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

void DiskWind::step()
{
    // Determine MPI Data
    const int root_process = 0;
    const int current_id = MPI::COMM_WORLD.Get_rank();
    const int processors = MPI::COMM_WORLD.Get_size();
    const int chunksize = NGrid / processors;


    //Single processor

    if (processors == 1) {

        for (int i = 0; i < NGrid; i++) {
            windloss[i] = densityLossAtRadius(g->convertIndexToPosition(i), i);
        }

        computeFluxes(0, NGrid + 1);

        for (int i = 0; i < chunksize; i++) {
            double rPlusHalf = g->convertIndexToPosition(i + 0.5);
            double rMinusHalf = g->convertIndexToPosition(i - 0.5);
            double dr = rPlusHalf - rMinusHalf;
            double r  = g->convertIndexToPosition(i);
            double currentArea = M_PI * au * au * (pow(rPlusHalf, 2) - pow(rMinusHalf, 2));

            tempData[i] = data[i].y + dt * (flux[i + 1] - flux[i]) / dr;

            tempData[i] -= windloss[i] * r * dt;
            accumulatedWindLoss += windloss[i] * currentArea * dt;

            if (i == 0) {
                accumulatedMassLossLeft += flux[i] / (dr * data[i].x) * dt * currentArea;
            } else if (i == NGrid - 1) {
                accumulatedMassLossRight += -flux[i + 1] / (dr * data[i].x) * dt * currentArea;
            }
        }

        for (int i = 0; i < chunksize; i++) {
            data[i].y = tempData[i];
            if (data[i].y < 0.0) {
                data[i].y = 0.0;
            }
            data[i].B2 = getUpdatedMagneticFluxDensityAtCell(i);
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

            if (i == 0) {
                accumulatedMassLossLeft += flux[i] / (dr * data[i].x) * dt * currentArea;
            } else if (i == NGrid - 1) {
                accumulatedMassLossRight += -flux[i + 1] / (dr * data[i].x) * dt * currentArea;
            }
        }

        for (int i = 0; i < chunksize; i++) {
            data[i].y = tempData[i];
            if (data[i].y < 0.0) {
                data[i].y = 0.0;
            }
            data[i].B2 = getUpdatedMagneticFluxDensityAtCell(i);
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
                MPI::COMM_WORLD.Recv(buffer, chunksize, MPI_DOUBLE, proc, frameRecv);
                MPI::COMM_WORLD.Recv(magneticBuffer, chunksize, MPI_DOUBLE, proc, frameBRecv);

                for (int i = proc * chunksize; i < (proc + 1) * chunksize; i++) {
                    data[i].y = buffer[i - proc * chunksize];
                    data[i].B2 = magneticBuffer[i - proc * chunksize];
                }

                delete[](buffer);
                delete[](magneticBuffer);
            }
            determineDiskExtent();
            writeFrame();
        }



    } else {
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
        }

        for (int i = minIndex; i < maxIndex; i++) {
            data[i].y = tempData[i - minIndex];
            if (data[i].y < 0.0) {
                data[i].y = 0.0;
            }
            data[i].B2 = getUpdatedMagneticFluxDensityAtCell(i);
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
            for (int i = minIndex; i < maxIndex; i++) {
                bfield[i - minIndex] = data[i].B2;
            }
            MPI::COMM_WORLD.Send(tempData, chunksize, MPI_DOUBLE, root_process, frameRecv);
            MPI::COMM_WORLD.Send(bfield, chunksize, MPI_DOUBLE, root_process, frameBRecv);

            delete[](bfield);
        }
    }
}

double DiskWind::getUpdatedMagneticFluxDensityAtCell(int i)
{
    double soundSpeed = sqrt(kb * T0 / (2.3 * mp * sqrt(data[i].x)));
    double scaleHeight = soundSpeed / sqrt(G*M/pow(data[i].x * au, 3));
    double midplaneDensity = data[i].y / data[i].x / (sqrt(2 * M_PI) * scaleHeight);
    return 8 * M_PI * midplaneDensity * soundSpeed * soundSpeed / plasma;
}

void DiskWind::determineDiskExtent()
{
    for (int i = NGrid - 1; i >= 0; i--) {
        if (data[i].y / data[i].x > floorDensity) {
            currentDiskExtent = g->convertIndexToPosition(i);
            std::cout << currentDiskExtent << std::endl;
            break;
        }
    }
}

double DiskWind::computeDiskMass()
{
    double mass = 0.0;
    for (int i = 0; i < NGrid; i++) {
        mass += data[i].y / data[i].x * M_PI * au * au * (pow(g->convertIndexToPosition(i + 0.5), 2) - pow(g->convertIndexToPosition(i - 0.5), 2));
    }
    return mass;
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


void DiskWind::setParameters(double a, double mass, double lum, double rg, double lever, int NFrames, GridGeometry *geometry, double plasmaParameter)
{
    alpha = a;
    M = mass;
    luminosity = lum;
    photoRadius = rg;
    leverArm = lever;
    maxFrames = NFrames;
    g = geometry;
    plasma = plasmaParameter;

    viscousConstant = 3 * alpha * kb * T0 / (sqrt(au) * 2.3 * mp * sqrt(G * M));
}


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
    }

    determineDiskExtent();
    writeFrame();
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

    for (unsigned int i = 0; i < parameters->size(); i++) {

        if (parameterType == "lambda") {
            leverArm = parameters->at(i);
        } else if (parameterType == "plasma") {
            plasma = parameters->at(i);
        }

        frame = 0;
        outputFrame = 0;
        initWithHCGADensityDistribution(diskMass, radialScale, floorDensity);

        int maxIndex = (int) g->convertPositionToIndex(currentDiskExtent - 100.0);

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
                    if (data[j].y / g->convertIndexToPosition(j) <= floorDensity) {
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

void DiskWind::runSimulation(int years)
{
    double NSteps = (double)years * year / dt;
    frameStride = (int)(NSteps / (double)maxFrames);

    const int root_process = 0;
    const int current_id = MPI::COMM_WORLD.Get_rank();
    const int processors = MPI::COMM_WORLD.Get_size();
    const int chunksize = NGrid / processors;


    if (processors == 1) {

        for (int i = 0; dt/year * i < years; i++) {
            step();
            double currentMass = computeDiskMass();
            if (frame % frameStride == 0) {
                std::cout << std::setprecision(16) << frame << ": " << currentMass << ", " << accumulatedMassLossLeft << ", " << accumulatedMassLossRight << ", " << accumulatedWindLoss << std::endl;
                std::cout << std::setprecision(16) << "Total: " << currentMass + accumulatedMassLossLeft + accumulatedMassLossRight + accumulatedWindLoss << std::endl << std::endl;
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
                std::cout << std::setprecision(16) << frame << ": " << currentMass << ", " << accumulatedMassLossLeft << ", " << accumulatedMassLossRight << ", " << totalAccumulatedWindLoss << std::endl;
                std::cout << std::setprecision(16) << "Total: " << currentMass + accumulatedMassLossLeft + accumulatedMassLossRight + totalAccumulatedWindLoss << std::endl << std::endl;
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
}
