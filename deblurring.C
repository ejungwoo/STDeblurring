#include "STDeblurring.h"
#include "STDeblurring.cpp"

void deblurring()
{
    auto db = new STDeblurring();
    //double* CreateUniformExample(int numEvents, int numParticles, double xRange1, double xRange2, int numBins, double histRange1, double histRange2);
    db -> CreateUniformExample(1000, 10, -10, 10, 40, -20, 20);
}
