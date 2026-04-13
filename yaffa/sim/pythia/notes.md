Alternative implementation of he Blast Wave
```
/ BWM function definition with Newton integration
Double_t BWM_proper_fct(Double_t *x, Double_t *par)
{ // Parametrization and values taken from https://www.researchgate.net/profile/Smbat-Grigoryan/publication/354646091_A_three_component_model_for_hadron_p_rm_T-spectra_in_pp_and_Pb-Pb_collisions_at_the_LHC/links/62002b23b44cbe42272873dc/A-three-component-model-for-hadron-p-rm-T-spectra-in-pp-and-Pb-Pb-collisions-at-the-LHC.pdf?origin=publication_detail
    // Blast-Wave-Model
    // x[0] = intVar
    // par[0] = mass = 0.14 GeV
    // par[1] = betaS = 0.7937
    // par[2] = exponent = 2.261
    // par[3] = temperature = 0.158 GeV
    // par[4] = pT
    Int_t nSteps = 1000; // For now not settable in order to avoid false config.
    Float_t pT = x[0];
    Float_t mT = abs(TMath::Sqrt(TMath::Power(par[0], 2) + TMath::Power(pT, 2)));
    Float_t T = par[3];
    Float_t r = 0; // runs from 0 to 1;
    Float_t dr = 1. / (Float_t)nSteps;
    Float_t betaTn = 0;
    Float_t gammaTn = 0;
    Float_t sum = 0;
    for (unsigned int i = 0; i < nSteps + 1; ++i)
    {
        r = (0 + i * dr);
        betaTn = par[1] * TMath::Power(r, par[2]);
        gammaTn = 1. / sqrt((1 - TMath::Power(betaTn, 2)));
        if (false)
            if (i == 0 || i == nSteps)
            {
                printf("At iteration %i/%i: \n", i, nSteps);
                printf("r (with bounds 0 to 1): %f\n", r);
                printf("betaTn (with bounds 0 to %.4f): %f\n", betaTn, par[1]);
                printf("gammaTn (with bounds 0 to 1): %f\n", gammaTn);
            }
        sum += r * TMath::BesselI0(gammaTn * betaTn * pT / T) * TMath::BesselK1(gammaTn * mT / T);
    }
    // Take care of constant factors
    sum *= pT;
    sum *= mT;
    sum *= dr;
    return sum;
}
```