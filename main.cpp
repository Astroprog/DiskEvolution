#include "DiskWind.h"


int main()
{
    int NGrid = 200;
//    GridGeometry *g = new GridGeometry(1, -250.0);
//    Diffusion *diffusion = new Diffusion(NGrid);
//    diffusion->setFrameStride(500);
//    diffusion->setGeometry(g);
//    diffusion->initWithGaussian(0.0, 0.5);
//    diffusion->computedx();
//    diffusion->computedt();
//    diffusion->runSimulation(500000);

    GridGeometry *g = new GridGeometry(0.1, 100, 200, true);
    DiskWind *disk = new DiskWind(NGrid);
    disk->setParameters(0.01, 2e33);
    disk->setNumberOfFrames(1000);
    disk->setGeometry(g);
    disk->initWithDensityDistribution(500.0, 20);
    disk->computedx();
    disk->computedt();
    disk->runSimulation(1000000);


    return 0;
}