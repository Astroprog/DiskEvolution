#include "DiskWind.h"


int main()
{
    int NGrid = 150;
    bool restart = false;

    GridGeometry *g = new GridGeometry(0.1, 2000, NGrid, true);
    DiskWind *disk = new DiskWind(NGrid);
    disk->setParameters(0.01, 2e33);
    disk->setNumberOfFrames(500);
    disk->setGeometry(g);

    if (restart) {
        disk->initWithRestartData(700);
    } else {
        disk->initWithDensityDistribution(500.0, 100);
    }

    disk->computedx();
    disk->computedt();

    if (restart) {
        disk->restartSimulation(700, 2000000);
    } else {
        disk->runSimulation(5000000);
    }

    return 0;
}