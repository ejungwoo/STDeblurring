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

        //void SetNumBins(int numBins) { fNumBins = numBins; }
        //void SetXMin(int xMin) { fXMin = xMin; }
        //void SetXMax(int xMax) { fXMax = xMax; }
        void SetNumIterations(int numIterations) { fNumIterations = numIterations; }
        void SetBins(int numBins, double xMin, double xMax) { fNumBins = numBins; fXMin = xMin; fXMax = xMax; ConfigureBinWidth(); }
        void SetLambda(double lambda) { fLambda = lambda; }
        void SetAccelerationPower(double accelerationPower) { fAccelerationPower = accelerationPower; }
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

        TH1D* GetHistDataSets() { return fHistDataSets; }
        TH1D* GetHistMeasured() { return fHistMeasured; }
        TH1D* GetHistSmearing() { return fHistSmearing; }
        TH1D* GetHistProcess2() { return fHistProcess2; }
        TH1D* GetHistProcess3() { return fHistProcess3; }
        TH1D* GetHistRestored() { return fHistRestored; }


        void CreateUniformExample(double *array, int numParticles, double xRange1, double xRange2);

        enum class EDataType {
            kAverage,
            kReactionPlane,
        };

        void SetData(EDataType dataType, double *array, int numData);
        void Run();

    private:
        int          fNumIterations = 1000;
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

        TH1D*        fHistDataSets = nullptr;

        TH1D*        fHistMeasured = nullptr;
        TH1D*        fHistSmearing = nullptr;
        TH1D*        fHistProcess2 = nullptr;
        TH1D*        fHistProcess3 = nullptr;
        TH1D*        fHistRestored = nullptr;

    //ClassDef(STDeblurring,1);
};

#endif
