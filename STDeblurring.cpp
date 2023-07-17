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

    fArrayInitData = new double[fNumBins];
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
    memset(fArrayInitData, 0, sizeof(double)*fNumBins);
    memset(fArraySmearing, 0, sizeof(double)*fNumBins);
    memset(fArrayProcess2, 0, sizeof(double)*fNumBins);
    memset(fArrayProcess3, 0, sizeof(double)*fNumBins);
    memset(fArrayRestored, 0, sizeof(double)*fNumBins);
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
    for (auto iBin=0; iBin<fNumBins; ++iBin) cout << fArrayInitData[iBin] << " "; cout << endl;
    for (auto iBin=0; iBin<fNumBins; ++iBin) cout << fArraySmearing[iBin] << " "; cout << endl;
    for (auto iBin=0; iBin<fNumBins; ++iBin) cout << fArrayProcess2[iBin] << " "; cout << endl;
    for (auto iBin=0; iBin<fNumBins; ++iBin) cout << fArrayProcess3[iBin] << " "; cout << endl;
    for (auto iBin=0; iBin<fNumBins; ++iBin) cout << fArrayRestored[iBin] << " "; cout << endl;
}

void STDeblurring::Draw(Option_t *option)
{
    if (fHistInitData==nullptr) fHistInitData = new TH1D(fName+"_initData",fTitle+" initial distribution;x",fNumBins,fXMin,fXMax);
    if (fHistSmearing==nullptr) fHistSmearing = new TH1D(fName+"_smearing",fTitle+" smearing;x",fNumBins,fXMin,fXMax);
    if (fHistProcess2==nullptr) fHistProcess2 = new TH1D(fName+"_process2",fTitle+";x",fNumBins,fXMin,fXMax);
    if (fHistProcess3==nullptr) fHistProcess3 = new TH1D(fName+"_process3",fTitle+";x",fNumBins,fXMin,fXMax);
    if (fHistRestored==nullptr) fHistRestored = new TH1D(fName+"_restored",fTitle+" restored distribution;x",fNumBins,fXMin,fXMax);

    for (auto iBin=0; iBin<fNumBins; ++iBin) fHistInitData -> SetBinContent(iBin+1,fArrayInitData[iBin]);
    for (auto iBin=0; iBin<fNumBins; ++iBin) fHistSmearing -> SetBinContent(iBin+1,fArraySmearing[iBin]);
    for (auto iBin=0; iBin<fNumBins; ++iBin) fHistProcess2 -> SetBinContent(iBin+1,fArrayProcess2[iBin]);
    for (auto iBin=0; iBin<fNumBins; ++iBin) fHistProcess3 -> SetBinContent(iBin+1,fArrayProcess3[iBin]);
    for (auto iBin=0; iBin<fNumBins; ++iBin) fHistRestored -> SetBinContent(iBin+1,fArrayRestored[iBin]);

    auto cvs = new TCanvas("cvs","",1000,800);
    cvs -> Divide(3,2);
    cvs -> cd(1); fHistInitData -> Draw();
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
    for (auto iBin=0; iBin<numBins; ++iBin)
        array[iBin] = fHistReference -> GetBinContent(iBin+1);

    //const double *array = fHistReference -> GetArray();
    //for (auto iBin=0; iBin<numBins; ++iBin)
    //    cout << fHistReference -> GetBinContent(iBin) << " " << array[iBin] << endl;

    /*
    double x_mean_of_variance = sum_variance / numEvents;
    double variance_uniform_dist = x_max/3.;
    double variance_div_num_particles = variance_uniform_dist / numParticles;
    double sigma_new = TMath::Sqrt(variance_div_num_particles);
    double norm_factor = 1. / (sigma_new * TMath::Sqrt(2*TMath::Pi()));

    /// renormalized width of average position: ptcle vs rest
    double variance_norm = x_mean_of_variance * numParticles/(numParticles-1);
    double sqrt_vn = TMath::Sqrt(variance_norm);
    double norm_factor2 = 1. / (sqrt_vn * TMath::Sqrt(2*TMath::Pi()));

    for (auto iBin=1; iBin<=num_bins; ++iBin) {
        double bin_x = hist_x_smearing -> GetBinCenter(iBin);
        double smearing = norm_factor2 * TMath::Exp(-.5*bin_x*bin_x/variance_norm)*bin_width;
        hist_x_smearing -> SetBinContent(iBin,smearing);
    }

    double sum_x_renormalized=0.;
    for (auto iBin=1; iBin<=num_bins; ++iBin)
        sum_x_renormalized += fHistReference -> GetBinContent(iBin);
        */

    return array;
}


void Run()
{
    double xMeanOfVariance = sum_variance / numEvents;
    double varianceOfUniformDist = x_max/3.;
    double variance_div_num_particles = varianceOfUniformDist / numParticles;
    double sigma_new = TMath::Sqrt(variance_div_num_particles);
    double norm_factor = 1. / (sigma_new * TMath::Sqrt(2*TMath::Pi()));

    /// renormalized width of average position: ptcle vs rest
    double variance_norm = xMeanOfVariance * numParticles/(numParticles-1);
    double sqrt_vn = TMath::Sqrt(variance_norm);
    double norm_factor2 = 1. / (sqrt_vn * TMath::Sqrt(2*TMath::Pi()));

    for (auto iBin=0; iBin<fNumBins; ++iBin)
    {
        double xBin = fHistReference -> GetBinCenter(iBin);
        double smearing = norm_factor2 * TMath::Exp(-.5*xBin*xBin/variance_norm)*bin_width;
        fArraySmearing -> SetBinContent(iBin,smearing);
    }

    for (auto iIterate=0; iIterate<fNumIterations; ++iIterate)
    {
        for (auto iBin2=0; iBin2<num_bins; ++iBin2)
        {
            double xSmearRestore = 0.;
            for (auto hBinS=0; hBinS<num_bins; ++hBinS)
            {
                int iBinR = iBin2 - hBinS + num_off;
                if (iBinR>0 && iBinR<=num_bins)
                {
                    double xSmeared = fArraySmearing -> GetBinContent(hBinS);
                    double xRestore = hist_x_restored -> GetBinContent(iBinR);
                    xSmearRestore += xSmeared * xRestore;
                }
            }
            /// convoluted current restore w/blurring function
            fArrayProcess2[iBin2] = xSmearRestore;
        }
        //break;

        for (auto iBinR=0; iBinR<num_bins; ++iBinR)
        {
            double x_restore1=0.;
            for (auto hBinS=0; hBinS<num_bins; ++hBinS)
            {
                int iBin_m = iBinR + hBinS - num_off;
                if (iBin_m>0 && iBin_m<=num_bins)
                {
                    double xSmeared        = hist_x_smearing -> GetBinContent(hBinS);
                    double x_renormalized   = hist_x_measured -> GetBinContent(iBin_m);
                    double xSmearRestore  = fArrayProcess2 -> GetBinContent(iBin_m);
                    if (xSmearRestore>0.)
                        x_restore1 += xSmeared * x_renormalized / xSmearRestore;
                    else
                        x_restore1 += xSmeared;
                }
            }


            double regularization_factor = 1.;
            if (iBinR>1 && iBinR<=num_bins-1)
            {
                double x_reno_0 = hist_x_restored -> GetBinContent(iBinR);
                double x_reno_b = hist_x_restored -> GetBinContent(iBinR-1);
                double x_reno_a = hist_x_restored -> GetBinContent(iBinR+1);
                if      (x_reno_0 > x_reno_b && x_reno_0 > x_reno_a) regularization_factor = 1./(1. + lambda); // peak
                else if (x_reno_0 < x_reno_b && x_reno_0 < x_reno_a) regularization_factor = 1./(1. - lambda); // bump
            }

            double x_restore2 = hist_x_restored -> GetBinContent(iBinR);
            double x_restore3 = x_restore2 * regularization_factor * TMath::Power(x_restore1, acc_power);
            hist_x_in_process_3 -> SetBinContent(iBinR, x_restore3);
        }

        double sum_of_xsr=0.;
        for (auto iBin=0; iBin<num_bins; ++iBin)
            sum_of_xsr += hist_x_in_process_3 -> GetBinContent(iBin);

        for (auto iBin=0; iBin<num_bins; ++iBin) {
            double content = hist_x_in_process_3 -> GetBinContent(iBin) * sum_x_renormalized / sum_of_xsr;
            hist_x_restored -> SetBinContent(iBin,content);
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
                auto hist_restored_clone = (TH1D *) hist_x_restored -> Clone(Form("restored_it%d_%d",iIterate,iseed));
                hist_restored_clone -> SetTitle(Form("restored with iteration=%d",iIterate));
                hist_restored_clone -> Draw();
            }
        }
        */
    }
}
