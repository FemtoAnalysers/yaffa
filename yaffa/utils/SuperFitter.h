#ifndef SUPERFITTER_H
#define SUPERFITTER_H

#include <cmath>

#include "Observable.h"
#include "TFormula.h"
#include "TH1.h"
#include "TObject.h"


double Pol0(double* x, double* p) {
    return p[0];
}

std::map<const char*, double (*)(double*, double*)> functions {
    {"mypol0", Pol0},
};

class SuperFitter : public TObject {
   private:
    Observable* fObs;
    TF1* fFit;

   public:
    // Empty Contructor
    SuperFitter() : fObs(nullptr) {};

    // Standard Contructor
    SuperFitter(Observable* hObs);

    // Destructor
    ~SuperFitter();

    // Add function
    // void Add(const char* opt);

    // Fit
    void Fit(const char* opt ="" );

    // Draw
    void Draw(const char* opt = "") const;

    ClassDef(SuperFitter, 1)
};

ClassImp(SuperFitter);

// Empty Constructor
SuperFitter::~SuperFitter() { delete fObs; }

// Standard Constructor
SuperFitter::SuperFitter(Observable* hObs) { this->fObs = hObs; }

// Add
// void SuperFitter::Add(const char* name) {
//     TFormula::DefineFunction(name, functions[name]);
// }

// // Fit
void SuperFitter::Fit(const char* opt) { this->fObs->Fit(opt); }

// Draw
void SuperFitter::Draw(const char* opt) const { this->fObs->Draw(opt); }

#endif
