#include "ViscousDisk.h"


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

    GridGeometry *g = new GridGeometry(0.1, 1000, 200, true);
    ViscousDisk *disk = new ViscousDisk(NGrid);
    disk->setParameters(0.01, 1e33);
    disk->setNumberOfFrames(1000);
    disk->setGeometry(g);
    disk->initWithDensityDistribution(500.0, 100);
    disk->computedx();
    disk->computedt();
    disk->runSimulation(1000000);


    return 0;
}