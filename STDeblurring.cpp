#include "STDeblurring.h"

#include <iostream>
using namespace std;

//ClassImp(STDeblurring);

STDeblurring::STDeblurring(const char* name, const char* title)
: TNamed(name,title)
{
    Init();
    Clear();
}

bool STDeblurring::Init()
{
    // Put intialization todos here which are not iterative job though event
    cout << "Initializing STDeblurring" << std::endl;

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

    if (fHistMeasured!=nullptr) fHistMeasured -> Reset();
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
    if (fHistDataSets==nullptr) fHistDataSets = new TH1D(fName+"_datasets",fTitle+" example data set;x",fNumBins,fXMin,fXMax);
    if (fHistMeasured==nullptr) fHistMeasured = new TH1D(fName+"_measured",fTitle+" measured distribution;x",fNumBins,fXMin,fXMax);
    if (fHistSmearing==nullptr) fHistSmearing = new TH1D(fName+"_smearing",fTitle+" smearing;x",fNumBins,fXMin,fXMax);
    if (fHistProcess2==nullptr) fHistProcess2 = new TH1D(fName+"_process2",fTitle+" dummy 2;x",fNumBins,fXMin,fXMax);
    if (fHistProcess3==nullptr) fHistProcess3 = new TH1D(fName+"_process3",fTitle+" dummy 3;x",fNumBins,fXMin,fXMax);
    if (fHistRestored==nullptr) fHistRestored = new TH1D(fName+"_restored",fTitle+" restored distribution;x",fNumBins,fXMin,fXMax);

    //for (auto bin=0; bin<fNumBins; ++bin) fHistMeasured -> SetBinContent(bin+1,fArrayMeasured[bin]);
    for (auto bin=0; bin<fNumBins; ++bin) fHistSmearing -> SetBinContent(bin+1,fArraySmearing[bin]);
    for (auto bin=0; bin<fNumBins; ++bin) fHistProcess2 -> SetBinContent(bin+1,fArrayProcess2[bin]);
    for (auto bin=0; bin<fNumBins; ++bin) fHistProcess3 -> SetBinContent(bin+1,fArrayProcess3[bin]);
    for (auto bin=0; bin<fNumBins; ++bin) fHistRestored -> SetBinContent(bin+1,fArrayRestored[bin]);

    if (fCanvas==nullptr) {
        fCanvas = new TCanvas("cvs","",1250,700);
        fCanvas -> Divide(3,2);
    }
    fCanvas -> cd(1); fHistDataSets -> Draw();
    fCanvas -> cd(2); fHistMeasured -> Draw();
    fCanvas -> cd(3);
    fHistRestored -> Draw();
    auto cHistMeasured = (TH1D *) fHistMeasured -> Clone(fName+"_measured2");
    cHistMeasured -> SetLineStyle(2);
    cHistMeasured -> SetLineColor(kRed);
    cHistMeasured -> Draw("same");
    fCanvas -> cd(4); fHistSmearing -> Draw();
    fCanvas -> cd(5); fHistProcess2 -> Draw();
    fCanvas -> cd(6); fHistProcess3 -> Draw();
}

void STDeblurring::SaveFigures() {
    fCanvas -> SaveAs(Form("figures/%s_%s.pdf",fName.Data(),fCanvas->GetName()));
}

void STDeblurring::CreateUniformExample(double *array, int numPoints, double xRange1, double xRange2)
{
    //cout << "Creating uniform example " << endl;
    //cout << "numPoints : " << numPoints << endl;
    //cout << "xRange1 : " << xRange1 << endl;
    //cout << "xRange2 : " << xRange2 << endl;

    for (auto iPoint=0; iPoint<numPoints; ++iPoint)
        array[iPoint] = gRandom -> Uniform(xRange1, xRange2);
}

void STDeblurring::SetData(EDataType dataType, double *array, int numPoints)
{
    bool thisIsFirstDataSet = false;
    if (fHistMeasured==nullptr) {
        thisIsFirstDataSet = true;
        fHistMeasured = new TH1D(fName+"_measured",fTitle+" measured distribution;x",fNumBins,fXMin,fXMax);
        fHistDataSets = new TH1D(fName+"_datasets",fTitle+" example data set;x",fNumBins,fXMin,fXMax);
    }
    //fHistDataSets -> Reset("ICES");

    double xMean=0.;
    double xSquaredMean=0.;
    for (int iPoint=0; iPoint<numPoints; ++iPoint) {
        auto value = array[iPoint];
        xMean += value;
        xSquaredMean += value*value;
        //if (thisIsFirstDataSet)
        fHistDataSets -> Fill(value);
    }
    xMean = xMean / numPoints;
    xSquaredMean = xSquaredMean / numPoints;
    double xAverageErrorSquared = (xSquaredMean-xMean*xMean) / (numPoints-1);

    //
    fMeanAverageErrorSquared += xAverageErrorSquared / numPoints;
    //
    for (auto iPoint=0; iPoint<numPoints; ++iPoint)
    {
        double value = array[iPoint];
        double xNormalized = (numPoints*xMean - value) / (numPoints-1);
        double xMeasured = value - xNormalized;
        fHistMeasured -> Fill(xMeasured);
    }
}

void STDeblurring::Run()
{
    double sumRenormalized=0.;
    for (auto bin=0; bin<fNumBins; ++bin) {
        fArrayMeasured[bin] = fHistMeasured -> GetBinContent(bin+1);
        fArrayRestored[bin] = fHistMeasured -> GetBinContent(bin+1);
        sumRenormalized += fArrayMeasured[bin];
    }

    //
    int num_off = 0;
    if (fNumBins%2==0) num_off = fNumBins/2;
    else num_off = (fNumBins+1)/2;
    //

    //
    double varianceOfUniformDist = (fXMax-fXMin)*(fXMax-fXMin)/12.; // variance;
    fMeanAverageErrorSquared = varianceOfUniformDist / 10;
    //

    double sigmaSmear = TMath::Sqrt(fMeanAverageErrorSquared);
    double normSmear = 1. / (sigmaSmear * TMath::Sqrt(2*TMath::Pi()));

    for (auto bin=0; bin<fNumBins; ++bin)
    {
        double xBin = fHistMeasured -> GetBinCenter(bin);
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
    }
}
