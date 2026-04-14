#ifndef SUPERFITTERMULTITRIAL_H
#define SUPERFITTERMULTITRIAL_H

#include <cmath>
#include <map>
#include <string>
#include <vector>
#include <stdio.h>
#include <algorithm>

#include "gsl/gsl_sf_dawson.h"

#include "Observable.h"
#include "Riostream.h"
#include "TF1.h"
#include "TFormula.h"
#include "TH1.h"
#include "TObject.h"
#include "TCanvas.h"

#include "Observable.h"
#include "SuperFitter.h"
#include "SuperFitterMultitrial.h"

using par = std::tuple<std::string, double, double, double>;

// Class for advanced multitrial fitting
class SuperFitterMultitrial : public SuperFitter {
   private:
    std::vector<Observable*> fCFs;                                                    // Observable to be fitted
    std::map<std::string, std::vector<TH1*>> fTemplates = {};
    std::map<std::string, std::vector<std::string>> fFunctions = {};
    std::map<std::string, std::vector<par>> fPars = {};
    std::vector<std::pair<std::string, std::string>> fDrawRecipes;

//     TF1* fFit;                                                           // Total fit function
//     std::vector<par> fPars;  // List of fit pars: (name, init, min, max)
//     double fFitRangeMin;                                                 // Fit range minimum
//     double fFitRangeMax;                                                 // Fit range maximum
//     double fDrawRangeMin;                                                // Draw range minimum
//     double fDrawRangeMax;                                                // Draw range maximum

   public:
    // // Empty Contructor
    SuperFitterMultitrial() : SuperFitter(), fCFs({}) {};



    // // Standard Contructor
    // SuperFitterMultitrial(Observable* oObservable, double xMin, double xMax)
    //     : TObject(), fCFs(oObservable) {};

//     // Destructor
//     ~SuperFitterMultitrial() { delete fCFs; };

    void AddCF(Observable* obs) {
        fCFs.push_back(obs);
    }

    void SetPars(std::string name, std::vector<par> pars) {
        if (fPars.contains(name))  {
            throw std::runtime_error("Parameters for '" + name + "' already defined");
            exit(1);
        }

        fPars.insert({name, pars});

    }
    void Add(std::string name, TH1 * hTemplate) {
        if (!fTemplates.contains(name)) fTemplates[name] = {};
        fTemplates[name].push_back(hTemplate);
    }

    void Add(std::string name, std::string func) {
        if (!fFunctions.contains(name)) fFunctions[name] = {};
        fFunctions[name].push_back(func);
    }

    void SetDrawRecipes(std::vector<std::pair<std::string, std::string>> recipes) {
        this->fDrawRecipes = recipes;
    }

    // Fit
    void FitMultitrials(std::string model, const char* opt = "") {
        int counter = 0;
        TCanvas *cFit = new TCanvas("cFitMultitrials", "", 600, 600);
        cFit->SaveAs("cFitMultitrials.pdf[");

        for (int iCF = 0; iCF < fCFs.size(); iCF++) {

            // Load templates
            // for (auto &[name, hTemplates] : fTemplates) {
            //     // find corresponding parameters
            //     if (!fPars.contains(name)) {
            //         throw std::runtime_error("Parameters for '" + name + "' not found");
            //     }

            //     std::vector<par> pars = fPars[name];
            //     printf("name: %s\n", name.data());

            //     // Template variations
            //     for (auto &hTemplate : hTemplates) {

            //     }
            // }

            SuperFitter* fitter = new SuperFitter(fCFs[iCF], 0, 0.5);
            
            // fitter->SetDrawRange(0, 0.5);
            for (auto &[name, pars] : fPars) {
                printf("name: %s\n", name.data()); 
                if (fTemplates.contains(name)) {
                    auto hTemplates = fTemplates[name];
                    fitter->Add(name, hTemplates[0], pars);
                } else if (fFunctions.contains(name)) {
                    auto fFuncs = fFunctions[name];
                    fitter->Add(name, fFuncs[0], pars);
                } else {
                    throw std::runtime_error("Term '" + name + "' is not implemented.");
                }
            }

            fitter->Fit(model, "MR+");
            cFit->DrawFrame(0, 0.98, 0.5, 1.5);
            fitter->Draw(fDrawRecipes);

            cFit->SaveAs("cFitMultitrials.pdf");

            // delete fitter;
        }
        cFit->SaveAs("cFitMultitrials.pdf]");


    //         Observable cf = fCFs[iCF];
    //         SuperFitter fitter(cf, 0, 0.5);



    //          #         if templFileName := term.get('file'):
    // #             templFile = TFile(templFileName)
    // #             template = utils.io.Load(templFile, term['path'])
    // #             if isinstance(template, TH1):
    // #                 hTemplate = utils.analysis.ChangeUnits(template, term.get('unit_mult', 1))
    // #                 hTemplate.SetDirectory(0)
    // #                 fitter.Add(term['name'], hTemplate, term['params'])

    // #             # elif isinstance(template, TF1):
    // #                 # Ccnvert to hist
    // #                 # hTemplate = TH1D(f"h{iTerm}", "", 500, 0, 2)
    // #                 # for iBin in range(500):
    // #                 #     bc = template.Eval(hTemplate.GetBinCenter(iBin + 1) * term.get('unit_mult', 1))
    // #                 #     print("bc", bc)
    // #                 #     hTemplate.SetBinContent(iBin + 1, bc)
    // #                 # hTemplate.SetDirectory(0)

    // #                 # c = TCanvas()
    // #                 # hTemplate.Draw()

    // #                 # hTemplate.Write()
    // #                 # # oFile.cd()
    // #                 # hTemplate.Write(f'hCF_{iTerm}')
    // #                 # c.SaveAs(f'test{iTerm}.png')


    // #                 # unitMult = term.get('unit_mult', 1)
    // #                 # SetOwnership(template, False)
    // #                 # template.SetName(f'f{iTerm}')

    // #                 # print('zezez', template)
    // #                 # fitter.Add(term['name'], template, term['params'], term.get('unit_mult', 1))
    // #                 # fitter.Add(term['name'], hTemplate, term['params'])
    // #             else:
    // #                 print('Type not implemented. Exit!')
    // #                 sys.exit()
    // #             templFile.Close()
    // #         else:
    // #             fitter.Add(term['name'], term['func'], term['params'])
    };

//     // Add fit component
//     void Add(std::string name, std::string func, std::vector<par> pars) {

//         if (func == "pol0") {
//             functions.push_back({name, Pol0, 1});
//         } else if (func == "pol1") {
//             functions.push_back({name, Pol1, 2});
//         } else if (func == "pol2") {
//             functions.push_back({name, Pol2, 3});
//         } else if (func == "pol3") {
//             functions.push_back({name, Pol3, 4});
//         } else if (func == "pol4") {
//             functions.push_back({name, Pol4, 5});
//         } else if (func == "pol5") {
//             functions.push_back({name, Pol5, 6});
//         } else if (func == "pol6") {
//             functions.push_back({name, Pol6, 7});
//         } else if (func == "pol7") {
//             functions.push_back({name, Pol7, 8});
//         } else if (func == "pol8") {
//             functions.push_back({name, Pol8, 9});
//         } else if (func == "pol9") {
//             functions.push_back({name, Pol9, 10});
//         } else if (func == "gaus") {
//             functions.push_back({name, Gaus, 3});
//         } else if (func == "lednicky") {
//             functions.push_back({name, Lednicky, 7});
//         } else {
//             throw std::runtime_error("Function " + func + " with name " + name + " is not implemented");
//         }

//         // Save fit settings
//         for (const auto& par : pars) {
//                   std::get<2>(par), std::get<3>(par));
//             this->fPars.push_back(par);
//         }
//     };

//     // Add template function
//     void Add(std::string name, TH1* hTemplate, std::vector<par> pars) {

//         functions.push_back(
//             {name, [hTemplate](double* x, double* p) { return p[0] * hTemplate->Interpolate(x[0]); }, 1});

//         // Save fit settings
//         for (const auto& par : pars) {
//                   std::get<2>(par), std::get<3>(par));
//             this->fPars.push_back(par);
//         }
//     };

//     // Add TF1 function
//     void Add(std::string name, TF1* fTemplate, std::vector<par> pars,
//              double unitMult) {  // todo: remove units mult here and put in .py
//         std::cout << "kekek" << fTemplate << std::endl;

//         functions.push_back(
//             {name,
//              [&, fTemplate, unitMult, this](double* x, double* p) { return p[0] * fTemplate->Eval(x[0] * unitMult); },
//              1});

//         // Save fit settings
//         for (const auto& par : pars) {
//                   std::get<2>(par), std::get<3>(par));
//             this->fPars.push_back(par);
//         }
//     };

//     // Draw
//     void Draw(std::vector<std::pair<std::string, std::string>> recipes) const {
//         std::cout << "Draw functions" << std::endl;
//         for (const auto& [name, _, __] : functions) {
//             std::cout << name << std::endl;
//         }


//         TLegend* leg = new TLegend(0.6, 0.6, 0.9, 0.9);

//         // Draw the fitted observable
//         this->fCFs->Draw("hist same pe");
//         leg->AddEntry(this->fCFs, "data", "pe");

//         // Draw the final fit function
//         this->fFit->Draw("same");

//         leg->AddEntry(this->fFit, "Total");
//         // Draw components based on the draw recipes
//         for (int iRecipe = 0; iRecipe < recipes.size(); iRecipe++) {
//             std::string legend = recipes[iRecipe].first;
//             std::string recipe = recipes[iRecipe].second;

//             // Tokenization of the recipe
//             auto tokens = Tokenize(recipe);

//             // Count the number of parameters
//             std::set<int> paraList = {};
//             std::vector<int> nParsDraw = {};
//             std::set<std::string> used_tokens = {};
//             for (const auto& token : tokens) {

//                 if (!IsFunction(token)) {
//                     continue;
//                 }

//                 int counter = 0;
//                 for (const auto& [name, _, __] : functions) {
//                     if (name == token) break;
//                     counter++;
//                 }

//                 int offset = 0;
//                 for (int iFunc = 0; iFunc < counter; iFunc++) {
//                     offset += std::get<2>(functions[iFunc]);
//                 }

//                 // Determine the number of parameters
//                 for (int iFunc = 0; iFunc < functions.size(); iFunc++) {
//                     auto name = std::get<0>(functions[iFunc]);
//                     if (name == token && used_tokens.find(token) == used_tokens.end()) {
//                         int nPars = std::get<2>(functions[iFunc]);
//                         nParsDraw.push_back(nPars);
//                         // Determine the position of the function in the list of functions
//                         for (int iPar = 0; iPar < nPars; iPar++) {
//                             paraList.insert(offset + iPar);
//                         }
//                         break;
//                     }
//                 }
//             }

//             for (const auto& e : paraList) {
//             }

//             // Convert to Reverse Polish Notation
//             auto rpn = toRPN(tokens);

//             for (const int& d : nParsDraw) {
//             }

//             // The following lambda evaluates the fit function
//             auto lambda = [this, rpn, nParsDraw](double* x, double* p) -> double {
//                 std::stack<double> stack;

//                 if (0.12 < x[0] && x[0] < .16)
//                 std::vector<std::pair<std::string, int>> nParameters = {};  // token, npars
//                 int idx = 0;
//                 for (const std::string& token : rpn) {
//                     if (stack.size() == 0) {
//                     } else {
//                         if (0.12 < x[0] && x[0] < .16)
//                     }

//                     if (isdigit(token[0]) || token[0] == '.') {
//                         // Push numbers
//                         stack.push(std::stod(token));
//                     } else if (IsFunction(token)) {

//                         // inly insert if not already present -> avoid duplicates
//                         if (std::find(nParameters.begin(), nParameters.end(), std::pair(token, 1)) ==
//                             nParameters.end()) {
//                             for (const auto& [name, _, npar] : functions) {
//                                 if (name == token) {
//                                     nParameters.push_back({token, npar});
//                                 }
//                             }
//                         }

//                         // Compute offset
//                         int offset = 0;
//                         for (const auto& [name, np] : nParameters) {
//                             if (name == token) break;
//                             offset += np;
//                         }


//                         // Determine the position of the function in the list of functions
//                         int counter = 0;
//                         for (const auto& [name, _, __] : functions) {
//                             if (name == token) break;
//                             counter++;
//                         }
//                         auto func = std::get<1>(functions[counter]);
//                         double value = func(x, p + offset);

//                         if (0.12 < x[0] && x[0] < .16)
//                         if (0.12 < x[0] && x[0] < .16)

//                         stack.push(value);
//                     } else if (IsOperator(token)) {
//                         // Apply operator
//                         if (stack.size() < 2) throw std::runtime_error("Insufficient arguments for operator");
//                         double b = stack.top();
//                         stack.pop();
//                         double a = stack.top();
//                         stack.pop();

//                         if (token == "+")
//                             stack.push(a + b);
//                         else if (token == "-")
//                             stack.push(a - b);
//                         else if (token == "*")
//                             stack.push(a * b);
//                         else if (token == "/")
//                             stack.push(a / b);
//                         else
//                             throw std::runtime_error("Unknown operator");
//                     } else {
//                         throw std::runtime_error("Unknown token: " + token);
//                     }
//                     if (0.12 < x[0] && x[0] < .16)
//                               join(" ", stack_to_vector(stack)).data());
//                 }

//                 if (0.12 < x[0] && x[0] < .16)
//                           join(" ", stack_to_vector(stack)).data());

//                 if (stack.size() != 1) throw std::runtime_error("Invalid RPN expression");
//                 return stack.top();
//             };


//             TF1* fTerm =
//                 new TF1(Form("fTerm%d", iRecipe), lambda, this->fDrawRangeMin, this->fDrawRangeMax, paraList.size());
//             int color = iRecipe + 2;
//             if (color >= 5) color++;
//             fTerm->SetLineColor(color);
//             fTerm->SetLineWidth(3);
//             fTerm->SetNpx(1000);

//             int counter = 0;
//             for (const int& par : paraList) {
//                 fTerm->FixParameter(counter++, this->fFit->GetParameter(par));
//             }
//             fTerm->Draw("same");
//             leg->AddEntry(fTerm, legend.data());
//         }
//         leg->Draw("same");
//     };

//     void SetDrawRange(double xMin, double xMax) {
//         this->fDrawRangeMin = xMin;
//         this->fDrawRangeMax = xMax;
//     }

//     TF1* GetFitFunction() { return this->fFit; }

    ClassDef(SuperFitterMultitrial, 1)
};

ClassImp(SuperFitterMultitrial);

#endif
