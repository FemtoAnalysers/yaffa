#ifndef SUPERFITTER_H
#define SUPERFITTER_H

#include <stdio.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>
#include "gsl/gsl_sf_dawson.h"

#include "Observable.h"
#include "Riostream.h"
#include "TF1.h"
#include "TFormula.h"
#include "TH1.h"
#include "TObject.h"

#if DO_DEBUG || 1
#define DEBUG(scopes, msg, ...)                          \
    do {                                                 \
        printf("[DEBUG] %s: ", __FUNCTION__);            \
        for (int i = 0; i < scopes; i++) printf("    "); \
        printf(msg, ##__VA_ARGS__);                      \
        printf("\n");                                    \
    } while (0)
#else
#define DEBUG(msg, ...)
#endif

// Definition of constants ---------------------------------------------------------------------------------------------
#define TINY std::numeric_limits<double>::min()
const double FmToNu(5.067731237e-3);
const double Pi(3.141592653589793);
const std::complex<double> i(0, 1);

// Definition of types -------------------------------------------------------------------------------------------------
namespace sf {
using parameter = std::tuple<std::string, double, double, double>;
}

// Definition of variables ---------------------------------------------------------------------------------------------

// List of TF1-compatible functions that can be used in the fit
std::vector<std::vector<std::tuple<std::string, std::function<double(double*, double*)>, int>>> functions = {};

// Utils ---------------------------------------------------------------------------------------------------------------

// Concatenate the elements of a std::vector via a separator. Equivalent of python's `" ".join(mylist)`
template <typename T>
std::string join(const std::string& separator, const std::vector<T>& list) {
    std::string output = "";
    for (const auto& element : list) {
        // Convert element to string
        std::stringstream ss;
        ss << element;
        output += ss.str() + separator;
    }

    if (!output.empty()) {
        output.resize(output.size() - separator.size());  // Remove last separator
    }

    return output;
}

// Processing of formulas ----------------------------------------------------------------------------------------------

// Convert a stack to a vector
template <typename T>
std::vector<T> stack_to_vector(const std::stack<T>& stack) {
    std::stack<T> tempStack = stack;  // Copy the original stack to preserve it
    std::vector<T> result;

    // Transfer elements from the stack to the vector (in reverse order)
    while (!tempStack.empty()) {
        result.push_back(tempStack.top());
        tempStack.pop();
    }

    // Reverse the vector to restore the original order of the stack
    std::reverse(result.begin(), result.end());

    return result;
}

// Tokenize a formula
std::vector<std::string> Tokenize(const std::string& formula) {
    DEBUG(0, "Start tokenization: formula = %s", formula.data());
    std::vector<std::string> tokens;
    std::string token;
    for (size_t i = 0; i < formula.size(); ++i) {
        char c = formula[i];
        DEBUG(1, "Analyzing char=%c", c);
        if (isspace(c)) continue;  // Ignore spaces

        if (isdigit(c) || c == '.') {
            // Handle numbers
            token += c;
        } else if (isalpha(c) || c == '_') {
            // Handle functions or variable names
            token += c;
        } else {
            // Handle operators or parentheses
            if (!token.empty()) {
                tokens.push_back(token);
                token.clear();
            }
            tokens.push_back(std::string(1, c));
        }
    }

    if (!token.empty()) {
        tokens.push_back(token);
    }
    DEBUG(0, "Tokenization is complete, reconstructed formula is '%s'", join(" ", tokens).data());

    return tokens;
}

// Get operator precedence
int GetPrecedence(const std::string& op) {
    if (op == "+" || op == "-") return 1;
    if (op == "*" || op == "/") return 2;
    return 0;
}

// Check if token is an operator
bool IsOperator(const std::string& token) { return token == "+" || token == "-" || token == "*" || token == "/"; }

// Check if token is a function. "raw" is a special token used for the total fit function
bool IsFunction(const std::string& token) {
    for (const auto& functionList : functions) {
        for (const auto& [name, _, __] : functionList) {
            if (token == name || token == "raw") return true;
        }
    }
    return false;
}

// Convert a vector of tokens into Reverse Polish Notation (RPN)
std::vector<std::string> toRPN(const std::vector<std::string>& tokens) {
    std::vector<std::string> output;
    std::stack<std::string> operators;

    for (const std::string& token : tokens) {
        if (isdigit(token[0]) || token[0] == '.' || IsFunction(token)) {
            // Numbers go directly to output
            output.push_back(token);
        } else if (token == "(") {
            operators.push(token);
        } else if (token == ")") {
            // Pop until matching "("
            while (!operators.empty() && operators.top() != "(") {
                output.push_back(operators.top());
                operators.pop();
            }
            if (!operators.empty()) operators.pop();  // Pop "("
            if (!operators.empty() && IsFunction(operators.top())) {
                output.push_back(operators.top());
                operators.pop();
            }
        } else if (IsOperator(token)) {
            while (!operators.empty() && GetPrecedence(operators.top()) >= GetPrecedence(token)) {
                output.push_back(operators.top());
                operators.pop();
            }
            operators.push(token);
        } else {
            // std::cout << "List of functions:" << std::endl;
            // for (const auto& [fn, _, __] : functions) {
            //     std::cout << fn << std::endl;
            // }
            throw std::runtime_error("Unrecognized token '" + token + "'");
        }
    }

    // Pop remaining operators
    while (!operators.empty()) {
        output.push_back(operators.top());
        operators.pop();
    }

    return output;
}

// Fit functions -------------------------------------------------------------------------------------------------------

// Normalized Gaussian
double Gaus(double* x, double* p) {
    double xx = x[0];
    double norm = p[0];
    double mean = p[1];
    double sigma = p[2];

    double normFactor = norm / (std::sqrt(2 * M_PI) * sigma);
    double exponent = -0.5 * std::pow((xx - mean) / sigma, 2);
    return normFactor * std::exp(exponent);
}

// Polynomial of degree 0
double Pol0(double* x, double* p) { return p[0]; }

// Polynomial of degree 1
double Pol1(double* x, double* p) { return Pol0(x, p) + p[1] * pow(x[0], 1); }

// Polynomial of degree 2
double Pol2(double* x, double* p) { return Pol1(x, p) + p[2] * pow(x[0], 2); }

// Polynomial of degree 3
double Pol3(double* x, double* p) { return Pol2(x, p) + p[3] * pow(x[0], 3); }

// Polynomial of degree 4
double Pol4(double* x, double* p) { return Pol3(x, p) + p[4] * pow(x[0], 4); }

// Polynomial of degree 5
double Pol5(double* x, double* p) { return Pol4(x, p) + p[5] * pow(x[0], 5); }

// Polynomial of degree 6
double Pol6(double* x, double* p) { return Pol5(x, p) + p[6] * pow(x[0], 6); }

// Polynomial of degree 7
double Pol7(double* x, double* p) { return Pol6(x, p) + p[7] * pow(x[0], 7); }

// Polynomial of degree 8
double Pol8(double* x, double* p) { return Pol7(x, p) + p[8] * pow(x[0], 8); }

// Polynomial of degree 9
double Pol9(double* x, double* p) { return Pol8(x, p) + p[9] * pow(x[0], 9); }

// Breit Wigner
double BreitWigner(double* x, double* par) {
    double kstar = x[0];

    double yield = par[0];
    double mean = par[1];
    double gamma = par[2];

    return yield * TMath::BreitWigner(kstar, mean, gamma);
}

// General Lednicky
double GeneralLednicky(double kstar, const double& GaussR, const complex<double>& a0, const double& effRange) {
    // printf("led\n");
    // Taken from
    // https://github.com/dimihayl/DLM/blob/c40f03eac38006f89eac8e5fa1533c9e48f2b455/CATS_Extentions/DLM_CkModels.cpp#L215
    if (GaussR != GaussR) {
        printf(
            "\033[1;33mWARNING:\033[0m GeneralLednicky got a bad value for the Radius (nan). Returning default value "
            "of 1.\n");
        exit(1);
        return 1;
    }

    kstar *= 1000;                  // change units to GeV/c
    kstar = std::max(kstar, 1.e-6);  // avoid problems with k* = 0

    const double Radius = GaussR * FmToNu;
    const complex<double> IsLen1 = 1. / (a0 * FmToNu + 1e-64);
    const double eRan1 = effRange * FmToNu;

    double F1 = gsl_sf_dawson(2. * kstar * Radius) / (2. * kstar * Radius);
    double F2 = (1. - exp(-4. * kstar * kstar * Radius * Radius)) / (2. * kstar * Radius);
    complex<double> ScattAmplSin = pow(IsLen1 + 0.5 * eRan1 * kstar * kstar - i * kstar, -1.);

    double CkValue = 0.;
    CkValue += 0.5 * pow(abs(ScattAmplSin) / Radius, 2) * (1. - (eRan1) / (2 * sqrt(Pi) * Radius)) +
               2 * real(ScattAmplSin) * F1 / (sqrt(Pi) * Radius) - imag(ScattAmplSin) * F2 / Radius;
    CkValue += 1;

    return CkValue;
}

// Lednicky
double Lednicky(double* x, double* par) {
    // Taken from
    // https://github.com/dimihayl/DLM/blob/c40f03eac38006f89eac8e5fa1533c9e48f2b455/CATS_Extentions/DLM_CkModels.cpp#L504C8-L504C53
    double kStar = x[0];
    double potPar0 = par[0];     // real part of the scattering length
    double potPar1 = par[1];     // imaginary part of the scattering length
    double potPar2 = par[2];     // effective range
    double sourcePar0 = par[3];  // radius of the first gaussian
    double sourcePar1 = par[4];  // radius of the second gaussian
    double sourcePar2 = par[5];  // relative weight of the two gaussians
    double sourcePar3 = par[6];  // normalization of the gaussians
    complex<double> ScatLen(potPar0, potPar1);

    double ll1 = GeneralLednicky(kStar, sourcePar0, ScatLen, potPar2);
    double ll2 = GeneralLednicky(kStar, sourcePar1, ScatLen, potPar2);
    return sourcePar3 * (sourcePar2 * ll1 + (1 - sourcePar2) * ll2) + 1. - sourcePar3;
}

// Class for advanced fitting ------------------------------------------------------------------------------------------
class SuperFitter : public TObject {
   private:
    std::vector<Observable*> fObs;                  // Observable to be fitted
    std::vector<TF1*> fFit;                         // Total fit function
    std::vector<std::vector<sf::parameter>> fPars;  // List of fit pars: (name, init, min, max)
    std::vector<TF1*> fTerms;
    std::vector<std::pair<double, double>> fFitRange;  // Fit range as the union of different intervals
    double fDrawRangeMin;                              // Draw range minimum
    double fDrawRangeMax;                              // Draw range maximum

   public:
    // Empty Contructor
    SuperFitter() : TObject(), fObs({}), fFit({}), fPars({}), fFitRange({}) {};

    // Destructor
    ~SuperFitter();

    bool IsInFitRange(double x);

    // SetModel
    void SetModel(int idx, std::string model);

    // Fit
    void Fit(const char* opt = "");

    // Add fit component
    void Add(int idx, std::string name, std::string func, std::vector<sf::parameter> pars);

    // Add template function
    void Add(int idx, std::string name, TH1* hTemplate, std::vector<sf::parameter> pars);

    // Add TF1 function
    void Add(int idx, std::string name, TF1* fTemplate, std::vector<sf::parameter> pars, double unitMult);

    // Draw
    void Draw(std::vector<std::pair<std::string, std::string>> recipes);

    // Set observable
    void AddObservable(Observable *obs) { this->fObs.push_back(obs); }

    // Set Fit range
    void SetFitRange(std::vector<std::pair<double, double>> fitRange) { this->fFitRange = fitRange; }

    // Set Draw range
    void SetDrawRange(double xMin, double xMax) {
        this->fDrawRangeMin = xMin;
        this->fDrawRangeMax = xMax;
    }

    TF1* GetFitFunction(int idx = 0) { return this->fFit[idx]; }
    TH1D* GetGenuineCF(int idx, std::string recipe);
    std::vector<TF1*> GetTerms() { return this->fTerms; }

    ClassDef(SuperFitter, 2)
};

// Destructor
SuperFitter::~SuperFitter() {
    fTerms.clear();
    functions.clear();
};

// Check if value is in fit range
bool SuperFitter::IsInFitRange(double x) {
    return true;
    for (const auto& [xMin, xMax] : this->fFitRange) {
        if (xMin < x && x < xMax) {
            return true;
        }
    }
    return false;
}

// Add fit component
void SuperFitter::Add(int idx, std::string name, std::string func, std::vector<sf::parameter> pars) {
    DEBUG(0, "Adding a new function '%s' to the fitter", name.data());

    if (idx > functions.size()) {
        throw std::invalid_argument("Index is larger than current length of the function list.");
    }

    if (idx > fPars.size()) {
        throw std::invalid_argument("Index is larger than current length of the parameter list.");
    }

    if (idx == functions.size()) {
        functions.push_back({});
    }

    if (idx == fPars.size()) {
        fPars.push_back({});
    }

    if (func == "pol0") {
        functions[idx].push_back({name, Pol0, 1});
    } else if (func == "pol1") {
        functions[idx].push_back({name, Pol1, 2});
    } else if (func == "pol2") {
        functions[idx].push_back({name, Pol2, 3});
    } else if (func == "pol3") {
        functions[idx].push_back({name, Pol3, 4});
    } else if (func == "pol4") {
        functions[idx].push_back({name, Pol4, 5});
    } else if (func == "pol5") {
        functions[idx].push_back({name, Pol5, 6});
    } else if (func == "pol6") {
        functions[idx].push_back({name, Pol6, 7});
    } else if (func == "pol7") {
        functions[idx].push_back({name, Pol7, 8});
    } else if (func == "pol8") {
        functions[idx].push_back({name, Pol8, 9});
    } else if (func == "pol9") {
        functions[idx].push_back({name, Pol9, 10});
    } else if (func == "gaus") {
        functions[idx].push_back({name, Gaus, 3});
    } else if (func == "breit_wigner") {
        functions[idx].push_back({name, BreitWigner, 3});
    } else if (func == "lednicky") {
        functions[idx].push_back({name, Lednicky, 7});
    } else {
        throw std::runtime_error("Function " + func + " with name " + name + " is not implemented");
    }

    // Save fit settings
    DEBUG(0, "Parameters are:");
    for (const auto& par : pars) {
        DEBUG(1, "name: %s   init: %.3f   min: %.3f   max: %.3f", std::get<0>(par).data(), std::get<1>(par),
              std::get<2>(par), std::get<3>(par));
        this->fPars[idx].push_back(par);
    }
};

// SetModel
void SuperFitter::SetModel(int idx, std::string model) {
    if (fFitRange.size() == 0) {
        throw std::runtime_error("Fit range is not specified!");
    }

    printf("Starting fit in Fit Range: ");
    for (const auto& [xMin, xMax] : this->fFitRange) {
        printf("[%.3f %.3f] U ", xMin, xMax);
    }
    printf("\n");

    // Tokenization of the model
    auto tokens = Tokenize(model);
    DEBUG(0, "Expression in infix: %s", join(" ", tokens).data());

    // Convert to Reverse Polish Notation
    auto rpn = toRPN(tokens);
    DEBUG(0, "Expression in RPN: %s", join(" ", rpn).data());

    // The following lambda evaluates the fit function
    auto lambda = [this, rpn, idx](double* x, double* p) -> double {
        std::stack<double> stack;

        if (false) DEBUG(0, "Compute fit function from RPN: '%s'", join(" ", rpn).data());
        for (const std::string& token : rpn) {
            if (false) DEBUG(1, "Start processing the token '%s'", token.data());
            if (stack.size() == 0) {
                if (false) DEBUG(1, "Stack is empty");
            } else {
                if (false) DEBUG(1, "Stack is: '%s'", join(" ", stack_to_vector(stack)).data());
            }
            if (isdigit(token[0]) || token[0] == '.') {
                if (false) DEBUG(2, "Token '%s' is a number", token.data());
                // Push numbers
                stack.push(std::stod(token));
            } else if (IsFunction(token)) {
                if (false) DEBUG(2, "Token '%s' is a function", token.data());

                // Determine the position of the function in the list of functions
                int counter = 0;
                for (const auto& [name, _, __] : functions[idx]) {
                    if (name == token) break;
                    counter++;
                }

                int offset = 0;
                for (int iFunc = 0; iFunc < counter; iFunc++) {
                    offset += std::get<2>(functions[idx][iFunc]);
                }
                if (false)
                    DEBUG(2, "Function '%s' was inserted in position: %d ==> Skipping %d parameters", token.data(),
                          counter, offset);

                auto func = std::get<1>(functions[idx][counter]);
                double value = func(x, p + offset);
                if (false) DEBUG(2, "Pushing %s(x=%.3f, p) = %.3f", token.data(), x[0], value);

                stack.push(value);
            } else if (IsOperator(token)) {
                if (false) DEBUG(2, "Token '%s' is an operator", token.data());
                // Apply operator
                if (stack.size() < 2) throw std::runtime_error("Insufficient arguments for operator");
                double b = stack.top();
                stack.pop();
                double a = stack.top();
                stack.pop();

                if (token == "+")
                    stack.push(a + b);
                else if (token == "-")
                    stack.push(a - b);
                else if (token == "*")
                    stack.push(a * b);
                else if (token == "/")
                    stack.push(a / b);
                else
                    throw std::runtime_error("Unknown operator");
            } else {
                if (false) DEBUG(2, "Token '%s' is unknown", token.data());
                throw std::runtime_error("Unknown token: " + token);
            }
            if (false) DEBUG(1, "Stack after processing the token: '%s'", join(" ", stack_to_vector(stack)).data());
        }

        if (false) DEBUG(0, "End of evaluation, return value is: '%s'", join(" ", stack_to_vector(stack)).data());

        if (stack.size() != 1) throw std::runtime_error("Invalid RPN expression");

        // Reject points outside of the fit range
        if (!IsInFitRange(x[0])) {
            TF1::RejectPoint();
        }

        return stack.top();
    };

    int nPars = 0;
    for (int iFunc = 0; iFunc < functions[idx].size(); iFunc++) {
        nPars += std::get<2>(functions[idx][iFunc]);
    }
    this->fFit.push_back(new TF1(Form("fFit_%d", idx), lambda, this->fDrawRangeMin, this->fDrawRangeMax, nPars));
    this->fFit[idx]->SetNpx(1000);

    for (int iPar = 0; iPar < this->fPars[idx].size(); iPar++) {
        auto par = this->fPars[idx][iPar];
        std::string name = std::get<0>(par);
        double centr = std::get<1>(par);
        double min = std::get<2>(par);
        double max = std::get<3>(par);

        printf("idx=%d\t par=%d\t name=%s\n", idx, iPar, name.data());
        this->fFit[idx]->SetParName(iPar, name.data());

        if (min > max) {
            this->fFit[idx]->FixParameter(iPar, centr);
        } else {
            if (!(min < centr && centr < max)) {
                printf("\033[33mWARNING: parameter '%s' is outside the allowed range\033[0m\n", name.data());
                centr = (min + max) / 2;
            }

            this->fFit[idx]->SetParameter(iPar, centr);
            this->fFit[idx]->SetParLimits(iPar, min, max);
        }
    }

    printf("Fit n. %d has %d parameters\n", idx, nPars);
    for (int iPar = 0; iPar < nPars; iPar++) {
        printf("\t%d\t%s\n", iPar, this->fFit[idx]->GetParName(iPar));
    }
    // for (int iFit = 0; iFit < fFit.size(); iFit++) {
    //     int np = this->fFit[iFit]->GetNpar();

    //     printf("Fit lalala %d has %d parameters:\n", iFit, np);

    //     nPars += np;
    // }  
};

// definition of shared parameter
// background function
std::vector<std::vector<int>> iPars = {
    {
        0,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
    },
    {
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        28,
        29,
        30,
        31,
        32,
        33,
        34,
        35,
        36,
        37,
        38,
    },
};


struct GlobalChi2 {
    GlobalChi2(std::vector<ROOT::Fit::Chi2Function *> chi2) : fChi2(chi2) {}

    double operator()(const double *par) const { 
        std::vector<int> sizes = {18, 21};

        double *pars[fChi2.size()];
        pars[0] = new double[sizes[0]];
        pars[1] = new double[sizes[1]];
        
        double chi2 = 0;
        for (size_t iChi2 = 0; iChi2 < fChi2.size(); iChi2++) {
            for (int i = 0; i < sizes[iChi2]; ++i) {
                pars[iChi2][i] = par[iPars[iChi2][i]];
            }
            chi2 += (*fChi2[iChi2])(pars[iChi2]);
        }

        return chi2;
    }

    const std::vector<ROOT::Fit::Chi2Function *> fChi2;
};


// Fit
void SuperFitter::Fit(const char* option) {

    printf("--- %.3f\n", fFit[0]->Eval(0.1));
    
    
    printf("Performing Combined fit\n");
    ROOT::Fit::DataOptions opt;

    ROOT::Math::WrappedMultiTF1 wf0(*(fFit[0]), 1);
    ROOT::Math::WrappedMultiTF1 wf1(*(fFit[1]), 1);

    printf("range: [%.3f, %.3f]\n", fFitRange[0].first, fFitRange[fFitRange.size() - 1].second);

    ROOT::Fit::DataRange range0;
    range0.SetRange(fFitRange[0].first, fFitRange[fFitRange.size() - 1].second);
    ROOT::Fit::BinData data0(opt, range0);
    fObs[0]->GetHistogram()->SetName("lala1");
    ROOT::Fit::FillData(data0, fObs[0]->GetHistogram());

    ROOT::Fit::DataRange range1;
    range1.SetRange(fFitRange[0].first, fFitRange[fFitRange.size() - 1].second);
    ROOT::Fit::BinData data1(opt, range1);
    ROOT::Fit::FillData(data1, fObs[1]->GetHistogram());

    ROOT::Fit::Chi2Function *chi2_0 = new ROOT::Fit::Chi2Function(data0, wf0);
    ROOT::Fit::Chi2Function *chi2_1 = new ROOT::Fit::Chi2Function(data1, wf1);

    std::vector<ROOT::Fit::Chi2Function *> chi2 = {chi2_0, chi2_1};
    // GlobalChi2 globalChi2(chi2_0, chi2_1);
    GlobalChi2 globalChi2(chi2);

    ROOT::Fit::Fitter fitter;

    std::vector<double> pars = {};

    // Count the number of fit parameters
    int nPars = 0;
    for (int iFit = 0; iFit < fFit.size(); iFit++) {
        int np = this->fFit[iFit]->GetNpar();
        for (int iPar = 0; iPar < np; iPar++) {
            printf("\t%d\t%s value=%.3f\n", nPars + iPar, this->fFit[iFit]->GetParName(iPar), this->fFit[iFit]->GetParameter(iPar));
            pars.push_back(this->fFit[iFit]->GetParameter(iPar));
        }
        nPars += np;
    }

    fitter.Config().SetParamsSettings(nPars, pars.data());

    // Set Parameter limits and fix
    nPars = 0;
    for (int iFit = 0; iFit < fFit.size(); iFit++) {
        int np = this->fFit[iFit]->GetNpar();

        printf("Fit %d has %d parameters:\n", iFit, np);
        for (int iPar = 0; iPar < np; iPar++) {
            // pars.push_back(this->fFit[iFit]->GetParameter(iPar));

            // Set Par Limits
            auto par = this->fPars[iFit][iPar];
            std::string name = std::get<0>(par);
            double centr = std::get<1>(par);
            double min = std::get<2>(par);
            double max = std::get<3>(par);

            fitter.Config().ParSettings(nPars + iPar).SetName(name.data());
            fitter.Config().ParSettings(nPars + iPar).SetValue(centr);

            if (min > max) {
                fitter.Config().ParSettings(nPars + iPar).Fix();
            } else {
                if (!(min < centr && centr < max)) {
                    printf("\033[33mWARNING: parameter '%s' is outside the allowed range\033[0m\n", name.data());
                    centr = (min + max) / 2;
                }

                fitter.Config().ParSettings(nPars + iPar).SetLimits(min, max);
            }
            printf("\t%d\t%s value=%.3f\n", nPars + iPar, this->fFit[iFit]->GetParName(iPar), this->fFit[iFit]->GetParameter(iPar));
        }

        nPars += np;
    }

    // create before the parameter settings in order to fix or set range on them
    // fix 5-th parameter
    // fitter.Config().ParSettings(4).Fix();
    // // set limits on the third and 4-th parameter
    // fitter.Config().ParSettings(2).SetLimits(-10, -1.E-4);
    // fitter.Config().ParSettings(3).SetLimits(0, 10000);
    // fitter.Config().ParSettings(3).SetStepSize(5);

    for (int iPar = 0; iPar < nPars; iPar++) {
        printf(" -> %s %.3f\n", fitter.Config().ParSettings(iPar).Name().data(), fitter.Config().ParSettings(iPar).Value());
    }
    fitter.Config().MinimizerOptions().SetPrintLevel(0);
    fitter.Config().SetMinimizer("Minuit2", "Migrad");

    // fit FCN function directly
    // (specify optionally data size and flag to indicate that is a chi2 fit)
    fitter.FitFCN(nPars, globalChi2, nullptr, data0.Size() + data1.Size(), true);
    ROOT::Fit::FitResult result = fitter.Result();
    result.Print(std::cout);
    
    // todo: change
    // this->fObs[0]->Fit(this->fFit[0], opt, fFitRange[0].first, fFitRange[fFitRange.size() - 1].second);
}

// Add TF1 function // todo: remove units mult here and put in .py
void SuperFitter::Add(int idx, std::string name, TF1* fTemplate, std::vector<sf::parameter> pars, double unitMult) {
    DEBUG(0, "Adding the function '%s' to the fitter", name.data());
    std::cout << "kekek" << fTemplate << std::endl;

    functions[idx].push_back(
        {name, [&, fTemplate, unitMult, this](double* x, double* p) { return p[0] * fTemplate->Eval(x[0] * unitMult); },
         1});

    // Save fit settings
    DEBUG(0, "Parameters are:");
    for (const auto& par : pars) {
        DEBUG(1, "name: %s   init: %.3f   min: %.3f   max: %.3f", std::get<0>(par).data(), std::get<1>(par),
              std::get<2>(par), std::get<3>(par));
        this->fPars[idx].push_back(par);
    }
};

// Add template function
void SuperFitter::Add(int idx, std::string name, TH1* hTemplate, std::vector<sf::parameter> pars) {
    DEBUG(0, "Adding the template '%s' to the fitter", name.data());

    functions[idx].push_back({name, [hTemplate](double* x, double* p) { return p[0] * hTemplate->Interpolate(x[0]); }, 1});

    // Save fit settings
    DEBUG(0, "Parameters are:");
    for (const auto& par : pars) {
        DEBUG(1, "name: %s   init: %.3f   min: %.3f   max: %.3f", std::get<0>(par).data(), std::get<1>(par),
              std::get<2>(par), std::get<3>(par));
        this->fPars[idx].push_back(par);
    }
};

// Draw
void SuperFitter::Draw(std::vector<std::pair<std::string, std::string>> recipes) {
    this->fTerms = {};

    for (const auto& [name, _, __] : functions[0]) {
        std::cout << name << std::endl;
    }

    DEBUG(0, "Start drawing");

    TLegend* leg = new TLegend(0.6, 0.6, 0.9, 0.9);

    // Draw the fitted observable
    // todo: change
    this->fObs[0]->Draw("hist same pe");
    // todo: change
    leg->AddEntry(this->fObs[0], "data", "pe");

    // Draw the final fit function
    this->fFit[0]->Draw("same");

    leg->AddEntry(this->fFit[0], "Total");
    // Draw components based on the draw recipes
    for (int iRecipe = 0; iRecipe < recipes.size(); iRecipe++) {
        std::string legend = recipes[iRecipe].first;
        std::string recipe = recipes[iRecipe].second;
        DEBUG(0, "Drawing the recipe '%s'", recipe.data());

        // Tokenization of the recipe
        auto tokens = Tokenize(recipe);
        DEBUG(0, "[DRAW] recipe in infix: %s", join(" ", tokens).data());

        // Count the number of parameters
        std::set<int> paraList = {};
        std::vector<int> nParsDraw = {};
        std::set<std::string> used_tokens = {};
        for (const auto& token : tokens) {
            DEBUG(1, "Processing token '%s'", token.data());

            if (!IsFunction(token)) {
                DEBUG(2, "Token '%s' is not a function --> skip!", token.data());
                continue;
            }

            int counter = 0;
            for (const auto& [name, _, __] : functions[0]) {
                if (name == token) break;
                counter++;
            }

            int offset = 0;
            for (int iFunc = 0; iFunc < counter; iFunc++) {
                offset += std::get<2>(functions[0][iFunc]);
            }

            // Determine the number of parameters
            for (int iFunc = 0; iFunc < functions[0].size(); iFunc++) {
                auto name = std::get<0>(functions[0][iFunc]);
                DEBUG(2, "Comparing with function '%s'", name.data());
                if (name == token && used_tokens.find(token) == used_tokens.end()) {
                    int nPars = std::get<2>(functions[0][iFunc]);
                    DEBUG(3, "It's a match! Number of parameters: %d", nPars);
                    nParsDraw.push_back(nPars);
                    // Determine the position of the function in the list of functions
                    for (int iPar = 0; iPar < nPars; iPar++) {
                        paraList.insert(offset + iPar);
                    }
                    break;
                }
            }
        }

        DEBUG(0, "ParaList:");
        for (const auto& e : paraList) {
            DEBUG(1, "%d value: %.3f", e, this->fFit[0]->GetParameter(e));
        }

        // Convert to Reverse Polish Notation
        auto rpn = toRPN(tokens);
        DEBUG(0, "[DRAW] Expression in RPN: %s", join(" ", rpn).data());

        for (const int& d : nParsDraw) {
            DEBUG(1, "par draw: %d", d);
        }

        // The following lambda evaluates the fit function
        auto lambda = [this, rpn, nParsDraw](double* x, double* p) -> double {
            std::stack<double> stack;

            if (0.12 < x[0] && x[0] < .16)
                DEBUG(0, "[DRAW] Compute fit function from RPN: '%s'", join(" ", rpn).data());
            std::vector<std::pair<std::string, int>> nParameters = {};  // token, npars
            int idx = 0;
            for (const std::string& token : rpn) {
                if (0.12 < x[0] && x[0] < .16) DEBUG(1, "[DRAW] Start processing the token '%s'", token.data());
                if (stack.size() == 0) {
                    if (0.12 < x[0] && x[0] < .16) DEBUG(1, "[DRAW] Stack is empty");
                } else {
                    if (0.12 < x[0] && x[0] < .16)
                        DEBUG(1, "[DRAW] Stack is: '%s'", join(" ", stack_to_vector(stack)).data());
                }

                if (isdigit(token[0]) || token[0] == '.') {
                    if (0.12 < x[0] && x[0] < .16) DEBUG(2, "[DRAW] Token '%s' is a number", token.data());
                    // Push numbers
                    stack.push(std::stod(token));
                } else if (IsFunction(token)) {
                    if (0.12 < x[0] && x[0] < .16) DEBUG(2, "Inserting token; %s   npar: %d\n", token.data(), 1);

                    // inly insert if not already present -> avoid duplicates
                    if (std::find(nParameters.begin(), nParameters.end(), std::pair(token, 1)) == nParameters.end()) {
                        for (const auto& [name, _, npar] : functions[0]) {
                            if (name == token) {
                                nParameters.push_back({token, npar});
                            }
                        }
                    }

                    // Compute offset
                    int offset = 0;
                    for (const auto& [name, np] : nParameters) {
                        if (name == token) break;
                        offset += np;
                    }

                    if (0.12 < x[0] && x[0] < .16) DEBUG(2, "[DRAW] Token '%s' is a function", token.data());

                    // Determine the position of the function in the list of functions
                    int counter = 0;
                    for (const auto& [name, _, __] : functions[0]) {
                        if (name == token) break;
                        counter++;
                    }
                    auto func = std::get<1>(functions[0][counter]);
                    double value = func(x, p + offset);

                    if (0.12 < x[0] && x[0] < .16)
                        DEBUG(2, "Counter: %d/%d    Offset: %d", counter, functions[0].size(), offset);
                    if (0.12 < x[0] && x[0] < .16)
                        DEBUG(2, "[DRAW] Pushing %s(x=%.3f, p) = %.3f", token.data(), x[0], value);

                    stack.push(value);
                } else if (IsOperator(token)) {
                    if (0.12 < x[0] && x[0] < .16) DEBUG(2, "[DRAW] Token '%s' is an operator", token.data());
                    // Apply operator
                    if (stack.size() < 2) throw std::runtime_error("Insufficient arguments for operator");
                    double b = stack.top();
                    stack.pop();
                    double a = stack.top();
                    stack.pop();

                    if (token == "+")
                        stack.push(a + b);
                    else if (token == "-")
                        stack.push(a - b);
                    else if (token == "*")
                        stack.push(a * b);
                    else if (token == "/")
                        stack.push(a / b);
                    else
                        throw std::runtime_error("Unknown operator");
                } else {
                    if (0.12 < x[0] && x[0] < .16) DEBUG(2, "[DRAW] Token '%s' is unknown", token.data());
                    throw std::runtime_error("Unknown token: " + token);
                }
                if (0.12 < x[0] && x[0] < .16)
                    DEBUG(1, "[DRAW] Stack after processing the token: '%s'", join(" ", stack_to_vector(stack)).data());
            }

            if (0.12 < x[0] && x[0] < .16)
                DEBUG(0, "[DRAW] End of evaluation, return value is: '%s'", join(" ", stack_to_vector(stack)).data());

            if (stack.size() != 1) throw std::runtime_error("Invalid RPN expression");
            return stack.top();
        };

        DEBUG(0, "Term '%s' needs %d parameters", recipe.data(), paraList.size());

        TF1* fTerm =
            new TF1(Form("fTerm%d", iRecipe), lambda, this->fDrawRangeMin, this->fDrawRangeMax, paraList.size());

        int colors[12] = {kBlue + 2,   kRed + 1,   kGreen + 3, kMagenta + 2, kCyan + 3, kOrange + 7,
                          kViolet + 3, kAzure + 4, kPink + 4,  kSpring - 7,  kTeal + 2, kGray + 2};

        fTerm->SetLineColor(colors[iRecipe]);
        fTerm->SetLineWidth(3);
        fTerm->SetNpx(1000);

        int counter = 0;
        for (const int& par : paraList) {
            DEBUG(2, "parameter val: %.3f", this->fFit[0]->GetParameter(par));
            fTerm->FixParameter(counter++, this->fFit[0]->GetParameter(par));
        }
        fTerm->Draw("same");
        fTerms.push_back(fTerm);
        leg->AddEntry(fTerm, legend.data());
    }
    leg->DrawClone("same");
};

// Get genuine correlation function
TH1D* SuperFitter::GetGenuineCF(int idx, std::string recipe) {
    // todo: change
    TH1D* hRawCF = (TH1D*)this->fObs[0]->GetHistogram();
    TH1D* hGenCF = (TH1D*)hRawCF->Clone("hGenCF");
    hGenCF->Reset();

    // Draw the fitted observable
    // todo: change
    this->fObs[0]->Draw("hist same pe");

    // Draw the final fit function
    this->fFit[0]->Draw("same");

    // Tokenization of the recipe
    auto tokens = Tokenize(recipe);

    // Count the number of parameters
    std::set<int> paraList = {};
    std::vector<int> nParsDraw = {};
    std::set<std::string> used_tokens = {};
    for (const auto& token : tokens) {
        int counter = 0;
        for (const auto& [name, _, __] : functions[idx]) {
            if (name == token) break;
            counter++;
        }

        int offset = 0;
        for (int iFunc = 0; iFunc < counter; iFunc++) {
            offset += std::get<2>(functions[idx][iFunc]);
        }

        // Determine the number of parameters
        for (int iFunc = 0; iFunc < functions[idx].size(); iFunc++) {
            auto name = std::get<0>(functions[idx][iFunc]);
            if (name == token && used_tokens.find(token) == used_tokens.end()) {
                int nPars = std::get<2>(functions[idx][iFunc]);
                nParsDraw.push_back(nPars);
                // Determine the position of the function in the list of functions
                for (int iPar = 0; iPar < nPars; iPar++) {
                    paraList.insert(offset + iPar);
                }
                break;
            }
        }
    }

    // Convert to Reverse Polish Notation
    auto rpn = toRPN(tokens);

    // The following lambda evaluates the fit function
    for (int iBin = 0; iBin < hGenCF->GetNbinsX(); iBin++) {
        double x = hGenCF->GetBinCenter(iBin + 1);

        std::stack<double> stack;

        std::vector<std::pair<std::string, int>> nParameters = {};  // token, npars
        int idx = 0;
        for (const std::string& token : rpn) {
            if (isdigit(token[0]) || token[0] == '.') {
                // Push numbers
                stack.push(std::stod(token));
            } else if (token == "raw") {
                stack.push(hRawCF->GetBinContent(iBin + 1));
            } else if (IsFunction(token)) {
                // inly insert if not already present -> avoid duplicates
                // if (std::find(nParameters.begin(), nParameters.end(), std::pair(token, 1)) ==
                //     nParameters.end()) {
                //     for (const auto& [name, _, npar] : functions) {
                //         if (name == token) {
                //             nParameters.push_back({token, npar});
                //         }
                //     }
                // }

                // // Compute offset
                // int offset = 0;
                // for (const auto& [name, np] : nParameters) {
                //     printf("offset name token %d %s %s\n", offset, name.data(), token.data());
                //     if (name == token) break;
                //     offset += np;
                // }

                // Determine the position of the function in the list of functions
                int counter = 0;
                for (const auto& [name, _, __] : functions[idx]) {
                    if (name == token) break;
                    counter++;
                }

                int offset = 0;
                for (const auto& [name, _, npar] : functions[idx]) {
                    if (name == token) break;
                    offset += npar;
                }

                auto func = std::get<1>(functions[idx][counter]);
                double* pars = fFit[0]->GetParameters();
                double value = func(&x, pars + offset);
                stack.push(value);
            } else if (IsOperator(token)) {
                // Apply operator
                if (stack.size() < 2) throw std::runtime_error("Insufficient arguments for operator");
                double b = stack.top();
                stack.pop();
                double a = stack.top();
                stack.pop();

                if (token == "+")
                    stack.push(a + b);
                else if (token == "-")
                    stack.push(a - b);
                else if (token == "*")
                    stack.push(a * b);
                else if (token == "/")
                    stack.push(a / b);
                else
                    throw std::runtime_error("Unknown operator");
            } else {
                throw std::runtime_error("Unknown token: " + token);
            }
        }

        if (stack.size() != 1) throw std::runtime_error("Invalid RPN expression");

        double cf = stack.top();
        double cfUnc = hRawCF->GetBinError(iBin + 1) * 2;  //! Use proper uncertainty

        if (std::isfinite(cf) && std::isfinite(cfUnc)) {
            hGenCF->SetBinContent(iBin + 1, stack.top());
            hGenCF->SetBinError(iBin + 1, cfUnc);
        }
        // return stack.top();
    }

    return hGenCF;
}

ClassImp(SuperFitter);

#endif
