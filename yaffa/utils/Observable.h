#ifndef OBSERVABLE_H
#define OBSERVABLE_H

#include "TH1.h"
#include "TF1.h"
#include "TObject.h"

class Observable : public TObject {
   private:
    TH1* fHObs;

   public:
    // Empty Contructor
    Observable() : fHObs(nullptr) {};

    // Standard Contructor
    Observable(TH1* hObs);

    // Destructor
    ~Observable();

    // Draw
    void Draw(const char* opt = "") const;

    // Fit
    void Fit(TF1 *fFit, const char* opt = "", double xMin=std::nan(""), double xMax=std::nan("")) const;

    ClassDef(Observable, 1)
};

ClassImp(Observable);

// Empty Constructor
Observable::~Observable() { delete fHObs; }

// Standard Constructor
Observable::Observable(TH1* hObs) { this->fHObs = hObs; }

// Draw
void Observable::Draw(const char* opt) const { this->fHObs->Draw(opt); }

// Fit
void Observable::Fit(TF1 *fFit, const char* opt, double xMin, double xMax) const {
    if (xMin == xMin && xMax == xMax) {
        this->fHObs->Fit(fFit, opt, "", xMin, xMax);    
    } else if (xMin != xMin && xMax != xMax) {
        this->fHObs->Fit(fFit, opt);
    } else {
        throw std::invalid_argument("Both xMin and xMax must be specified");
    }
}


#endif
