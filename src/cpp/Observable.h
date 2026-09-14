#ifndef OBSERVABLE_H
#define OBSERVABLE_H

#include "TH1.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TObject.h"

class Observable : public TObject {
   private:
    TGraphErrors* fGObs;

   public:
    // Empty Contructor
    Observable() : fGObs(nullptr) {};

    // Standard Contructor
    Observable(TGraphErrors* gObs);

    // Contructor from a histogram: one point per bin at the bin center
    Observable(TH1* hObs);

    // Destructor
    ~Observable();


    TGraphErrors* GetGraph() {return this->fGObs; };

    // Draw
    void Draw(const char* opt = "") const;

    // Fit
    void Fit(TF1 *fFit, const char* opt = "", double xMin=std::nan(""), double xMax=std::nan("")) const;

    ClassDef(Observable, 2)
};

ClassImp(Observable);

// Empty Constructor
Observable::~Observable() { delete fGObs; }

// Standard Constructor
Observable::Observable(TGraphErrors* gObs) { this->fGObs = gObs; }

// Constructor from a histogram
Observable::Observable(TH1* hObs) {
    int nBins = hObs->GetNbinsX();
    this->fGObs = new TGraphErrors(nBins);
    this->fGObs->SetName(hObs->GetName());
    for (int iBin = 0; iBin < nBins; iBin++) {
        this->fGObs->SetPoint(iBin, hObs->GetBinCenter(iBin + 1), hObs->GetBinContent(iBin + 1));
        this->fGObs->SetPointError(iBin, 0, hObs->GetBinError(iBin + 1));
    }
}

// Draw
void Observable::Draw(const char* opt) const { this->fGObs->Draw(opt); }

// Fit
void Observable::Fit(TF1 *fFit, const char* opt, double xMin, double xMax) const {
    if (xMin == xMin && xMax == xMax) {
        this->fGObs->Fit(fFit, opt, "", xMin, xMax);
    } else if (xMin != xMin && xMax != xMax) {
        this->fGObs->Fit(fFit, opt);
    } else {
        throw std::invalid_argument("Both xMin and xMax must be specified");
    }
}


#endif
