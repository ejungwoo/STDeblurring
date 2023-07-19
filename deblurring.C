#include "STDeblurring.h"
#include "STDeblurring.cpp"

void deblurring()
{
    bool useInputData = true;

    auto bl = new STDeblurring();

    bl -> SetBins(41,-20.5,20.5);
    //bl -> SetBins(40,-20,20);
    bl -> SetLambda(0.01);
    bl -> SetAccelerationPower(1.99);
    bl -> SetNumIterations(1000);

    const int numParticles = 10;
    const int numEvents = 10000;
    double* array = new double[numParticles];

    ifstream fileGenIn;
    ofstream fileGenOut;
    if (useInputData) fileGenIn.open("data/generated_data.dat");
    else fileGenOut.open("data/generated_data.dat");

    for (auto iEvent=0; iEvent<numEvents; ++iEvent) {
        if (useInputData) {
            double value;
            for (auto iParticle=0; iParticle<numParticles; ++iParticle) {
                fileGenIn >> value;
                array[iParticle]  = value;
            }
        }
        else {
            bl -> CreateUniformExample(array, numParticles, -10, 10);
            for (auto iParticle=0; iParticle<numParticles; ++iParticle)
                fileGenOut << array[iParticle] << " ";
            fileGenOut << endl;
        }
        bl -> SetData(STDeblurring::EDataType::kAverage, array, numParticles);
        if (iEvent<5) {
            auto cvs = new TCanvas(Form("cvs_check_%d",iEvent),Form("cvs_check_%d",iEvent),1000,450);
            cvs -> Divide(2,1);
            auto histDataSets = bl -> GetHistDataSets();
            auto histMeasured = bl -> GetHistMeasured();
            cvs -> cd(1); histDataSets -> Clone(Form("c_hist_data_%d",iEvent)) -> Draw();
            cvs -> cd(2); histMeasured -> Clone(Form("c_hist_meas_%d",iEvent)) -> Draw();
        }
    }

    bl -> Run();
    bl -> Draw();
    bl -> SaveFigures();
}
