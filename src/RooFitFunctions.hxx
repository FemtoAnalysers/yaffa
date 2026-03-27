/*
* Usage: must be compiled first via 
*     root -l -b -q -e '.L RooFitFunctions.hxx++'
* and then you must load the .so in python
*     gSystem.Load('/path/to/yaffa/src/RooFitFunctions_hxx.so')
*     from ROOT import RooSourceAAA
*/

#include <Riostream.h>

#include "TMath.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooRealProxy.h"

#include "Functions.hxx"

class RooSourceAAA : public RooAbsPdf {
   public:
    RooSourceAAA() {};
    RooSourceAAA(const char* name, const char* title, RooAbsReal& _hyperRadius, RooAbsReal& _rho0)
        : RooAbsPdf(name, title),
          hyperRadius("hyperRadius", "hyperRadius", this, _hyperRadius),
          rho0("rho0", "rho0", this, _rho0) {};
    RooSourceAAA(const RooSourceAAA& other, const char* name = 0)
        : RooAbsPdf(other, name),
          hyperRadius("hyperRadius", this, other.hyperRadius),
          rho0("rho0", this, other.rho0) {};
    virtual TObject* clone(const char* newname) const { return new RooSourceAAA(*this, newname); }

   protected:
    RooRealProxy hyperRadius;
    RooRealProxy rho0;

    Double_t evaluate() const { return _SourceAAA(hyperRadius, rho0); };

   private:
    ClassDef(RooSourceAAA, 1)
};

ClassImp(RooSourceAAA)
