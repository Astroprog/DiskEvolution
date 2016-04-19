#include "DiskWind.h"


int main()
{
    int NGrid = 200;

    GridGeometry *g = new GridGeometry(0.1, 1000, 200, true);
    DiskWind *disk = new DiskWind(NGrid);
    disk->setParameters(0.01, 2e33);
    disk->setNumberOfFrames(100);
    disk->setGeometry(g);
//    disk->initWithDensityDistribution(500.0, 100);
    disk->initWithRestartData(198);
    disk->computedx();
    disk->computedt();
//    disk->runSimulation(100000);
    disk->restartSimulation(198, 100000);


    return 0;
}