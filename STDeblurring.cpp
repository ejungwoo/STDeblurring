#include "STDeblurring.h"

#include <iostream>
using namespace std;

//ClassImp(STDeblurring);

STDeblurring::STDeblurring(const char* name, const char* title)
: TNamed(name,title)
{
    Init();
}

bool STDeblurring::Init()
{
    // Put intialization todos here which are not iterative job though event
    cout << "Initializing STDeblurring" << std::endl;

    fHistMeasured = new TH1D(fName+"_measured",fTitle+" initial distribution;x",fNumBins,fXMin,fXMax);
    fHistSmearing = new TH1D(fName+"_smearing",fTitle+" smearing;x",fNumBins,fXMin,fXMax);

    fArrayMeasured = new double[fNumBins];
    fArraySmearing = new double[fNumBins];
    fArrayProcess2 = new double[fNumBins];
    fArrayProcess3 = new double[fNumBins];
    fArrayRestored = new double[fNumBins];

    return true;
}

void STDeblurring::Clear(Option_t *option)
{
    TObject::Clear(option);
    fNumIterations = 100;
    fNumBins = 100;
    fXMin = 10;
    fXMax = 0;
    fLambda = 0.01;
    fAccelerationPower = 1.99;

    memset(fArrayMeasured, 0, sizeof(double)*fNumBins);
    memset(fArraySmearing, 0, sizeof(double)*fNumBins);
    memset(fArrayProcess2, 0, sizeof(double)*fNumBins);
    memset(fArrayProcess3, 0, sizeof(double)*fNumBins);
    memset(fArrayRestored, 0, sizeof(double)*fNumBins);

    fMeanAverageErrorSquared = 0;

    fHistMeasured -> Reset();
    fHistSmearing -> Reset();
}

void STDeblurring::Print(Option_t *option) const
{
    // You will probability need to modify here
    cout << "STDeblurring" << std::endl;
    cout << "fNumIterations : " << fNumIterations << std::endl;
    cout << "fNumBins : " << fNumBins << std::endl;
    cout << "fXMin : " << fXMin << std::endl;
    cout << "fXMax : " << fXMax << std::endl;
    cout << "fLambda : " << fLambda << std::endl;
    cout << "fAccelerationPower : " << fAccelerationPower << std::endl;
    for (auto bin=0; bin<fNumBins; ++bin) cout << fArrayMeasured[bin] << " "; cout << endl;
    for (auto bin=0; bin<fNumBins; ++bin) cout << fArraySmearing[bin] << " "; cout << endl;
    for (auto bin=0; bin<fNumBins; ++bin) cout << fArrayProcess2[bin] << " "; cout << endl;
    for (auto bin=0; bin<fNumBins; ++bin) cout << fArrayProcess3[bin] << " "; cout << endl;
    for (auto bin=0; bin<fNumBins; ++bin) cout << fArrayRestored[bin] << " "; cout << endl;
}

void STDeblurring::Draw(Option_t *option)
{
    if (fHistMeasured==nullptr) fHistMeasured = new TH1D(fName+"_measured",fTitle+" initial distribution;x",fNumBins,fXMin,fXMax);
    if (fHistSmearing==nullptr) fHistSmearing = new TH1D(fName+"_smearing",fTitle+" smearing;x",fNumBins,fXMin,fXMax);
    if (fHistProcess2==nullptr) fHistProcess2 = new TH1D(fName+"_process2",fTitle+";x",fNumBins,fXMin,fXMax);
    if (fHistProcess3==nullptr) fHistProcess3 = new TH1D(fName+"_process3",fTitle+";x",fNumBins,fXMin,fXMax);
    if (fHistRestored==nullptr) fHistRestored = new TH1D(fName+"_restored",fTitle+" restored distribution;x",fNumBins,fXMin,fXMax);

    for (auto bin=0; bin<fNumBins; ++bin) fHistMeasured -> SetBinContent(bin+1,fArrayMeasured[bin]);
    for (auto bin=0; bin<fNumBins; ++bin) fHistSmearing -> SetBinContent(bin+1,fArraySmearing[bin]);
    for (auto bin=0; bin<fNumBins; ++bin) fHistProcess2 -> SetBinContent(bin+1,fArrayProcess2[bin]);
    for (auto bin=0; bin<fNumBins; ++bin) fHistProcess3 -> SetBinContent(bin+1,fArrayProcess3[bin]);
    for (auto bin=0; bin<fNumBins; ++bin) fHistRestored -> SetBinContent(bin+1,fArrayRestored[bin]);

    auto cvs = new TCanvas("cvs","",1000,800);
    cvs -> Divide(3,2);
    cvs -> cd(1); fHistMeasured -> Draw();
    cvs -> cd(2); fHistRestored -> Draw();
    cvs -> cd(3); fHistSmearing -> Draw();
    cvs -> cd(4); fHistProcess2 -> Draw();
    cvs -> cd(5); fHistProcess3 -> Draw();
}

double* STDeblurring::CreateUniformExample(int numEvents, int numParticles, double xRange1, double xRange2, int numBins, double histRange1, double histRange2)
{
    cout << "Creating uniform example " << endl;
    cout << "numEvents    : " << numEvents << endl;
    cout << "numParticles : " << numParticles << endl;
    cout << "xRange1      : " << xRange1 << endl;
    cout << "xRange2      : " << xRange2 << endl;
    cout << "numBins      : " << numBins << endl;
    cout << "histRange1   : " << histRange1 << endl;
    cout << "histRange2   : " << histRange2 << endl;

    if (fHistReference==nullptr)
        fHistReference = new TH1D(fName+"_reference",fTitle+"initial distribution;x",numBins,histRange1,histRange2);

    double sum_variance = 0.;
    for (auto iEvent=0; iEvent<numEvents; ++iEvent)
    {
        double xMean=0.;
        double xSquaredMean=0.;
        double xInitArray[numParticles]; 
        for (auto iParticle=0; iParticle<numParticles; ++iParticle)
        {
            double xInit = gRandom -> Uniform(xRange1, xRange2);
            xInitArray[iParticle] = xInit;
            xMean         += xInit;
            xSquaredMean  += xInit*xInit;
        }
        xMean             = xMean        / numParticles;
        xSquaredMean      = xSquaredMean / numParticles;
        double x_variance = (xSquaredMean-xMean*xMean) / (numParticles-1);
        sum_variance += x_variance;
        for (auto iParticle=0; iParticle<numParticles; ++iParticle)
        {
            double xInit = xInitArray[iParticle];
            double xNormalized = (numParticles*xMean - xInit) / (numParticles-1);
            double xRenormalized = xInit - xNormalized;
            fHistReference -> Fill(xRenormalized);
            fHistReference -> Fill(xRenormalized);
        }
    }

    double* array = new double[fNumBins];
    for (auto bin=0; bin<numBins; ++bin)
        array[bin] = fHistReference -> GetBinContent(bin+1);

    return array;
}

void STDeblurring::SetData(double *array, int numData)
{
    double x_mean=0.;
    double x_squared_mean=0.;
    for (int iData=0; iData<numData; ++iData) {
        auto value = array[iData];
        x_mean         += value;
        x_squared_mean += value*value;
    }
    x_mean            = x_mean         / numData;
    x_squared_mean    = x_squared_mean / numData;
    double x_variance = (x_squared_mean-x_mean*x_mean);
    double x_average_error_squared = (x_squared_mean-x_mean*x_mean) / (numData-1);

    fMeanAverageErrorSquared += x_average_error_squared;
    for (auto i_particle=0; i_particle<numData; ++i_particle)
    {
        double x_init = array[i_particle];
        double x_normalized = (numData*x_mean - x_init) / (numData-1);
        double x_renormalized = x_init - x_normalized;
        fHistMeasured -> Fill(x_renormalized);
        fHistRestored -> Fill(x_renormalized);
    }
}

void STDeblurring::Run()
{
    //
    int num_off = 0;
    if (fNumBins%2==0) num_off = fNumBins/2;
    else num_off = (fNumBins+1)/2;
    //

    //
    double varianceOfUniformDist = (fXMax-fXMin)*(fXMax-fXMin)/12.; // variance;
    fMeanAverageErrorSquared = varianceOfUniformDist / 10;
    //
    //

    double sumRenormalized=0.;
    for (auto bin=0; bin<fNumBins; ++bin) {
        fArrayMeasured[bin] = fHistMeasured -> GetBinContent(bin+1);
        fArrayRestored[bin] = fHistRestored -> GetBinContent(bin+1);
        sumRenormalized += fArrayMeasured[bin];
    }

    double sigmaSmear = TMath::Sqrt(fMeanAverageErrorSquared);
    double normSmear = 1. / (sigmaSmear * TMath::Sqrt(2*TMath::Pi()));

    for (auto bin=0; bin<fNumBins; ++bin)
    {
        double xBin = fHistReference -> GetBinCenter(bin);
        double smearing = normSmear * TMath::Exp(-.5*xBin*xBin/sigmaSmear)*fBinWidth;
        fArraySmearing[bin] = smearing;
    }

    for (auto iIterate=0; iIterate<fNumIterations; ++iIterate)
    {
        for (auto bin2=0; bin2<fNumBins; ++bin2)
        {
            double xSmearRestore = 0.;
            for (auto binS=0; binS<fNumBins; ++binS)
            {
                int binR = bin2 - binS + num_off;
                if (binR>0 && binR<=fNumBins)
                {
                    double xSmeared = fArraySmearing[binS];
                    double xRestore = fArrayRestored[binR];
                    xSmearRestore += xSmeared * xRestore;
                }
            }
            /// convoluted current restore w/blurring function
            fArrayProcess2[bin2] = xSmearRestore;
        }

        for (auto binR=0; binR<fNumBins; ++binR)
        {
            double x_restore1=0.;
            for (auto binS=0; binS<fNumBins; ++binS)
            {
                int binM = binR + binS - num_off;
                if (binM>0 && binM<=fNumBins)
                {
                    double xSmeared      = fArraySmearing[binS];
                    double xReNormalized = fArrayMeasured[binM];
                    double xSmearRestore = fArrayProcess2[binM];

                    if (xSmearRestore>0.)
                        x_restore1 += xSmeared * xReNormalized / xSmearRestore;
                    else
                        x_restore1 += xSmeared;
                }
            }


            double regularization_factor = 1.;
            if (binR>1 && binR<=fNumBins-1)
            {
                double x_reno_0 = fArrayRestored[binR];
                double x_reno_b = fArrayRestored[binR-1];
                double x_reno_a = fArrayRestored[binR+1];
                if      (x_reno_0 > x_reno_b && x_reno_0 > x_reno_a) regularization_factor = 1./(1. + fLambda); // peak
                else if (x_reno_0 < x_reno_b && x_reno_0 < x_reno_a) regularization_factor = 1./(1. - fLambda); // bump
            }

            double x_restore2 = fArrayRestored[binR];
            double x_restore3 = x_restore2 * regularization_factor * TMath::Power(x_restore1, fAccelerationPower);
            fArrayProcess3[binR] = x_restore3;
        }

        double sum_of_xsr=0.;
        for (auto bin=0; bin<fNumBins; ++bin)
            sum_of_xsr += fArrayProcess3[bin];

        for (auto bin=0; bin<fNumBins; ++bin) {
            double content = fArrayProcess3[bin] * sumRenormalized / sum_of_xsr;
            fArrayRestored[bin] = content;
        }

        /*
        if (draw_process1)
        {
            if (
                    (fNumIterations<=15&&iIterate<15) ||
                    (iIterate==0 ||
                     iIterate==1 ||
                     iIterate==3 ||
                     iIterate==5 ||
                     iIterate==10 ||
                     iIterate==50 ||
                     iIterate==100 ||
                     iIterate==200 ||
                     iIterate==500 ||
                     iIterate==1000 ||
                     iIterate==1500 ||
                     iIterate==2000 ||
                     iIterate==3000 ||
                     iIterate==4000 ||
                     iIterate==5000)
               )
            {
                ++i_pad_restore;
                cvs_restore -> cd(i_pad_restore);
                auto hist_restored_clone = (TH1D *) fArrayRestored -> Clone(Form("restored_it%d_%d",iIterate,iseed));
                hist_restored_clone -> SetTitle(Form("restored with iteration=%d",iIterate));
                hist_restored_clone -> Draw();
            }
        }
        */
    }
}
