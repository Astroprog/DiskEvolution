#include "DiskWind.h"


int main()
{
    int NGrid = 200;
    bool restart = false;

    GridGeometry *g = new GridGeometry(0.1, 1000, 200, true);
    DiskWind *disk = new DiskWind(NGrid);
    disk->setParameters(0.01, 2e33);
    disk->setNumberOfFrames(100);
    disk->setGeometry(g);

    if (restart) {
        disk->initWithRestartData(455);
    } else {
        disk->initWithDensityDistribution(500.0, 100);
    }

    disk->computedx();
    disk->computedt();

    if (restart) {
        disk->restartSimulation(455, 100000);
    } else {
        disk->runSimulation(1000000);
    }

    return 0;
}