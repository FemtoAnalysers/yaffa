#ifndef SUPERFITTER_H
#define SUPERFITTER_H

#include <cmath>
#include <map>

#include "Observable.h"
#include "Riostream.h"
#include "TF1.h"
#include "TFormula.h"
#include "TH1.h"
#include "TObject.h"

double Gaus(double* x, double* p) {
    double xx = x[0];
    double norm = p[0];
    double mean = p[1];
    double sigma = p[2];

    double normFactor = norm / (std::sqrt(2 * M_PI) * sigma);
    double exponent = -0.5 * std::pow((xx - mean) / sigma, 2);
    return normFactor * std::exp(exponent);
}

double Pol0(double* x, double* p) { return p[0]; }
double Pol1(double* x, double* p) { return Pol0(x, p) + p[1] * pow(x[0], 1); }
double Pol2(double* x, double* p) { return Pol1(x, p) + p[2] * pow(x[0], 2); }
double Pol3(double* x, double* p) { return Pol2(x, p) + p[3] * pow(x[0], 3); }
double Pol4(double* x, double* p) { return Pol3(x, p) + p[4] * pow(x[0], 4); }
double Pol5(double* x, double* p) { return Pol4(x, p) + p[5] * pow(x[0], 5); }
double Pol6(double* x, double* p) { return Pol5(x, p) + p[6] * pow(x[0], 6); }
double Pol7(double* x, double* p) { return Pol6(x, p) + p[7] * pow(x[0], 7); }
double Pol8(double* x, double* p) { return Pol7(x, p) + p[8] * pow(x[0], 8); }
double Pol9(double* x, double* p) { return Pol8(x, p) + p[9] * pow(x[0], 9); }

class SuperFitter : public TObject {
   private:
    Observable* fObs;
    TF1* fFit;
    double fMin;
    double fMax;

   public:
    // Empty Contructor
    SuperFitter() : fObs(nullptr), fFit(nullptr) {};

    // Standard Contructor
    SuperFitter(Observable* hObs, double xMin, double xMax);

    // Destructor
    ~SuperFitter();

    // Fit
    void Fit(const char* opt = "");

    // Add fit component
    void Add(const char* name = "");

    // Draw
    void Draw(const char* opt = "") const;

    ClassDef(SuperFitter, 1)
};

ClassImp(SuperFitter);

// Empty Constructor
SuperFitter::~SuperFitter() { delete fObs; }

// Standard Constructor
SuperFitter::SuperFitter(Observable* hObs, double xMin, double xMax) {
    this->fObs = hObs;
    this->fMin = xMin;
    this->fMax = xMax;
}

// Fit
void SuperFitter::Fit(const char* opt) {
    this->fFit->SetParameter(0, 1);
    this->fFit->FixParameter(1, 0.2);
    this->fFit->FixParameter(2, 0.1);
    this->fObs->Fit(this->fFit, opt);
}

// Add fit component
void SuperFitter::Add(const char* name) {
    if (strcmp(name, "pol0") == 0) {
        this->fFit = new TF1("fFit", Pol0, this->fMin, this->fMax, 1);
    } else if (strcmp(name, "pol1") == 0) {
        this->fFit = new TF1("fFit", Pol1, this->fMin, this->fMax, 2);
    } else if (strcmp(name, "pol2") == 0) {
        this->fFit = new TF1("fFit", Pol2, this->fMin, this->fMax, 3);
    } else if (strcmp(name, "pol3") == 0) {
        this->fFit = new TF1("fFit", Pol3, this->fMin, this->fMax, 4);
    } else if (strcmp(name, "pol4") == 0) {
        this->fFit = new TF1("fFit", Pol4, this->fMin, this->fMax, 5);
    } else if (strcmp(name, "pol5") == 0) {
        this->fFit = new TF1("fFit", Pol5, this->fMin, this->fMax, 6);
    } else if (strcmp(name, "pol6") == 0) {
        this->fFit = new TF1("fFit", Pol6, this->fMin, this->fMax, 7);
    } else if (strcmp(name, "pol7") == 0) {
        this->fFit = new TF1("fFit", Pol7, this->fMin, this->fMax, 8);
    } else if (strcmp(name, "pol8") == 0) {
        this->fFit = new TF1("fFit", Pol8, this->fMin, this->fMax, 9);
    } else if (strcmp(name, "pol9") == 0) {
        this->fFit = new TF1("fFit", Pol9, this->fMin, this->fMax, 10);
    } else if (strcmp(name, "gaus") == 0) {
        this->fFit = new TF1("fFit", Gaus, this->fMin, this->fMax, 3);
    } else {
        printf("Not implemented. Exit!\n");
        exit(1);
    }
}

// Draw
void SuperFitter::Draw(const char* opt) const { this->fObs->Draw(opt); }

#endif
