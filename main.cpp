#include "DiskWind.h"


int main()
{
    int NGrid = 200;
    bool restart = true;

    GridGeometry *g = new GridGeometry(0.1, 1000, 200, true);
    DiskWind *disk = new DiskWind(NGrid);
    disk->setParameters(0.01, 2e33);
    disk->setNumberOfFrames(100);
    disk->setGeometry(g);

    if (restart) {
        disk->initWithRestartData(10);
    } else {
        disk->initWithDensityDistribution(500.0, 100);
    }

    disk->computedx();
    disk->computedt();

    if (restart) {
        disk->restartSimulation(10, 10000);
    } else {
        disk->runSimulation(10000);
    }

    return 0;
}