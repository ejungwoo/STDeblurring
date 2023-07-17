#ifndef STDEBLURRING_HH
#define STDEBLURRING_HH

#include "TNamed.h"

class STDeblurring : public TNamed
{
    public:
        STDeblurring(const char* name="Deblurring", const char* title="Deblurring");
        virtual ~STDeblurring() { ; }

        bool Init();
        void Clear(Option_t *option="");
        void Print(Option_t *option="") const;
        void Draw(Option_t *option="");


        int GetNumIterations() const  { return fNumIterations; }
        int GetNumBins() const  { return fNumBins; }
        int GetXMin() const  { return fXMin; }
        int GetXMax() const  { return fXMax; }
        int GetLambda() const  { return fLambda; }
        int GetAccelerationPower() const  { return fAccelerationPower; }
        void SetNumIterations(int numIterations) { fNumIterations = numIterations; }
        //void SetNumBins(int numBins) { fNumBins = numBins; }
        //void SetXMin(int xMin) { fXMin = xMin; }
        //void SetXMax(int xMax) { fXMax = xMax; }
        void SetBins(int numBins, double xMin, double xMax) { fNumBins = numBins; fXMin = xMin; fXMax = xMax; ConfigureBinWidth(); }
        void SetLambda(int lambda) { fLambda = lambda; }
        void SetAccelerationPower(int accelerationPower) { fAccelerationPower = accelerationPower; }
        void ConfigureBinWidth() { fBinWidth = (fXMax-fXMin)/fNumBins; }


        double GetMeanAverageErrorSquared() const { return fMeanAverageErrorSquared; }
        void SetMeanAverageErrorSquared(double val) { fMeanAverageErrorSquared = val; }


        void SetArrayMeasured(double *inputArray) { memcpy(fArrayMeasured, inputArray, fNumBins * sizeof(double)); }
        void SetArraySmearing(double *inputArray) { memcpy(fArraySmearing, inputArray, fNumBins * sizeof(double)); }
        void SetArrayProcess2(double *inputArray) { memcpy(fArrayProcess2, inputArray, fNumBins * sizeof(double)); }
        void SetArrayProcess3(double *inputArray) { memcpy(fArrayProcess3, inputArray, fNumBins * sizeof(double)); }
        void SetArrayRestored(double *inputArray) { memcpy(fArrayRestored, inputArray, fNumBins * sizeof(double)); }
        double* GetArrayMeasured() const  { return fArrayMeasured; }
        double* GetArraySmearing() const  { return fArraySmearing; }
        double* GetArrayProcess2() const  { return fArrayProcess2; }
        double* GetArrayProcess3() const  { return fArrayProcess3; }
        double* GetArrayRestored() const  { return fArrayRestored; }


        double* CreateUniformExample(int numEvents, int numParticles, double xRange1, double xRange2, int numBins, double histRange1, double histRange2);

        void SetData(double *array, int numData);
        void Run();

    private:
        int          fNumIterations = 100;
        int          fNumBins = 100;
        double       fBinWidth = 0;
        double       fXMin = 10;
        double       fXMax = 0;
        double       fLambda = 0.01;
        double       fAccelerationPower = 1.99;

        double       fMeanAverageErrorSquared = 0;

        double*      fArrayMeasured;
        double*      fArraySmearing;
        double*      fArrayProcess2;
        double*      fArrayProcess3;
        double*      fArrayRestored;

        TH1D*        fHistMeasured = nullptr;
        TH1D*        fHistSmearing = nullptr;
        TH1D*        fHistProcess2 = nullptr;
        TH1D*        fHistProcess3 = nullptr;
        TH1D*        fHistRestored = nullptr;

        TH1D*        fHistReference = nullptr;

    //ClassDef(STDeblurring,1);
};

#endif
