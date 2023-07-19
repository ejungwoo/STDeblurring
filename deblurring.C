#include "STDeblurring.h"
#include "STDeblurring.cpp"

void deblurring()
{
    auto bl = new STDeblurring();

    bl -> SetBins(41,-20.5,20.5);
    //bl -> SetBins(40,-20,20);
    bl -> SetLambda(0.01);
    bl -> SetAccelerationPower(1.99);
    bl -> SetNumIterations(1000);

    const int numParticles = 8;
    const int numEvents = 10000;
    double* array = new double[numParticles];

    for (auto iEvent=0; iEvent<numEvents; ++iEvent) {
        bl -> CreateUniformExample(array, numParticles, -10, 10);
        bl -> SetData(STDeblurring::EDataType::kAverage, array,numParticles);
        if (iEvent<10) {
            auto cvs = new TCanvas(Form("cvs%d",iEvent),Form("cvs%d",iEvent),1000,450);
            cvs -> Divide(2,1);
            auto histDataSets = bl -> GetHistDataSets();
            auto histMeasured = bl -> GetHistMeasured();
            cvs -> cd(1); histDataSets -> Clone() -> Draw();
            cvs -> cd(2); histMeasured -> Clone() -> Draw();
        }
    }
    return;

    bl -> Run();
    bl -> Draw();
}
