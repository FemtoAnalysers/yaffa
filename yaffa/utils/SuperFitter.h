#ifndef SUPERFITTER_H
#define SUPERFITTER_H

#include <stdio.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "Observable.h"
#include "Riostream.h"
#include "TF1.h"
#include "TFormula.h"
#include "TH1.h"
#include "TObject.h"
#include "gsl/gsl_sf_dawson.h"

#define DEBUG(level, indent, msg, ...)                       \
    do {                                                     \
        if (level <= DEBUG_LEVEL) {                          \
            printf("[DEBUG] %s: ", __FUNCTION__);            \
            for (int i = 0; i < indent; i++) printf("    "); \
            printf(msg, ##__VA_ARGS__);                      \
            printf("\n");                                    \
        }                                                    \
    } while (0)

// Definition of constants ---------------------------------------------------------------------------------------------
#define TINY std::numeric_limits<double>::min()
const double FmToNu(5.067731237e-3);
const double Pi(3.141592653589793);
const std::complex<double> i(0, 1);
int colors[12] = {kBlue + 2,   kRed + 1,   kGreen + 3, kMagenta + 2, kCyan + 3, kOrange + 7,
                  kViolet + 3, kAzure + 4, kPink + 4,  kSpring - 7,  kTeal + 2, kGray + 2};
// Definition of types -------------------------------------------------------------------------------------------------
namespace sf {
using parameter = std::tuple<std::string, double, double, double>;
using func = std::function<double(double*, double*)>;
}

// Definition of variables ---------------------------------------------------------------------------------------------

// List of TF1-compatible functions that can be used in the fit
std::vector<std::vector<std::tuple<std::string, sf::func, int>>> functions = {};


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
    DEBUG(40, 0, "Start tokenization: formula = %s", formula.data());
    std::vector<std::string> tokens;
    std::string token;
    for (size_t i = 0; i < formula.size(); ++i) {
        char c = formula[i];
        DEBUG(40, 1, "Analyzing char=%c", c);
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
    DEBUG(40, 0, "Tokenization is complete, reconstructed formula is '%s'", join(" ", tokens).data());

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

    kstar *= 1000;                   // change units to GeV/c
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
    std::vector<Observable*> fObs;                     // Observable to be fitted
    std::vector<TF1*> fFit;                            // Total fit function
    std::vector<std::vector<sf::parameter>> fPars;     // List of fit pars: (name, init, min, max)
    std::vector<TF1*> fTerms;                          // Each function to be drawn
    std::vector<std::pair<double, double>> fFitRange;  // Fit range as the union of different intervals
    std::map<std::string, int> fParIndeces;            // Indeces of parameters for combined fit
    double fDrawRangeMin;                              // Draw range minimum
    double fDrawRangeMax;                              // Draw range maximum

   public:
    // Empty Contructor
    SuperFitter() : TObject(), fObs({}), fFit({}), fPars({}), fFitRange({}) {};

    // Destructor
    ~SuperFitter();

    bool IsInFitRange(double x);

    bool IsParameterPresent(std::string name);

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
    void Draw(int iFit, std::vector<std::pair<std::string, std::string>> recipes, std::string dataLabel="Data");

    // Set observable
    void AddObservable(Observable* obs) { this->fObs.push_back(obs); }

    // Set Fit range
    void SetFitRange(std::vector<std::pair<double, double>> fitRange) { this->fFitRange = fitRange; }

    // Set Draw range
    void SetDrawRange(double xMin, double xMax) {
        this->fDrawRangeMin = xMin;
        this->fDrawRangeMax = xMax;
    }

    int GetN();
    int GetNShared();
    int GetNIndependent();
    int GetNIndependent(int iFit);
    std::vector<double> GetInitialParameters();

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

// Checks if a parameter is already known to the fitter (in case of combined fit)
bool SuperFitter::IsParameterPresent(std::string name) {
    auto it = this->fParIndeces.find(name);

    if (it != this->fParIndeces.end()) {
        return true;
    }

    return false;
}

// Add fit component
void SuperFitter::Add(int idx, std::string name, std::string func, std::vector<sf::parameter> pars) {
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
    printf("Adding '%s' function with parameters:\n", name.data());
    for (const auto& par : pars) {
        auto [name, centr, min, max] = par;
        printf("    name: %s   init: %.3f   min: %.3f   max: %.3f\n", name.data(), centr, min, max);
        if (!IsParameterPresent(name)) {
            this->fPars[idx].push_back(par);
        }
    }
};

// Process operator token
void ProcessOperatorToken(std::stack<double> &stack, std::string token) {
    DEBUG(53, 2, "Token '%s' is an operator", token.data());
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
}

// Get index of function with a given name
int GetIndex(std::vector<std::tuple<std::string, sf::func, int>> funcs, std::string name) {
    int counter = 0;
    for (const auto& [fn, _, __] : funcs) {
        if (fn == name) break;
        counter++;
    }
    return counter;
}

// Compute how many parameters should be skipped
int ComputeOffset(std::vector<std::tuple<std::string, sf::func, int>> funcs, int counter) {
    int offset = 0;
    for (int iFunc = 0; iFunc < counter; iFunc++) {
        offset += std::get<2>(funcs[iFunc]);
    }
    return offset;
}

// SetModel
void SuperFitter::SetModel(int idx, std::string model) {
    // Tokenization of the model
    auto tokens = Tokenize(model);
    DEBUG(50, 0, "Expression in infix: %s", join(" ", tokens).data());

    // Convert to Reverse Polish Notation
    auto rpn = toRPN(tokens);
    DEBUG(50, 0, "Expression in RPN: %s", join(" ", rpn).data());

    // The following lambda evaluates the fit function
    auto lambda = [this, rpn, idx](double* x, double* p) -> double {
        std::stack<double> stack;

        DEBUG(51, 0, "Compute fit function from RPN: '%s'", join(" ", rpn).data());
        for (const std::string& token : rpn) {
            DEBUG(52, 1, "Start processing the token '%s'", token.data());
            if (stack.size() == 0) {
                DEBUG(53, 1, "Stack is empty");
            } else {
                DEBUG(53, 1, "Stack is: '%s'", join(" ", stack_to_vector(stack)).data());
            }
            if (isdigit(token[0]) || token[0] == '.') {
                DEBUG(53, 2, "Token '%s' is a number", token.data());
                // Push numbers
                stack.push(std::stod(token));
            } else if (IsFunction(token)) {
                int counter = GetIndex(functions[idx], token);
                int offset = ComputeOffset(functions[idx], counter);
                auto func = std::get<1>(functions[idx][counter]);

                DEBUG(53, 2, "Function '%s' at pos: %d ==> Skipping %d parameters", token.data(), counter, offset);
                stack.push(func(x, p + offset));
            } else if (IsOperator(token)) {
                ProcessOperatorToken(stack, token);
            } else {
                DEBUG(53, 2, "Token '%s' is unknown", token.data());
                throw std::runtime_error("Unknown token: " + token);
            }
            DEBUG(52, 1, "Stack after processing the token: '%s'", join(" ", stack_to_vector(stack)).data());
        }

        DEBUG(51, 0, "End of evaluation, return value is: '%s'", join(" ", stack_to_vector(stack)).data());

        if (stack.size() != 1) throw std::runtime_error("Invalid RPN expression");

        // Reject points outside of the fit range
        if (!IsInFitRange(x[0])) {
            TF1::RejectPoint();
        }

        return stack.top();
    };

    // Count how many parameter the function has
    int nPars = 0;
    for (int iFunc = 0; iFunc < functions[idx].size(); iFunc++) {
        nPars += std::get<2>(functions[idx][iFunc]);
    }

    this->fFit.push_back(new TF1(Form("fFit_%d", idx), lambda, this->fDrawRangeMin, this->fDrawRangeMax, nPars));
    this->fFit[idx]->SetNpx(1000);

    for (int iPar = 0; iPar < this->fPars[idx].size(); iPar++) {
        auto [name, _, __, ___] = this->fPars[idx][iPar];
        this->fFit[idx]->SetParName(iPar, name.data());
    }
};

// Global Chi2
struct GlobalChi2 {
    GlobalChi2(std::vector<ROOT::Fit::Chi2Function*> chi2, std::vector<std::vector<int>> parIndeces)
        : fChi2(chi2), fParIndeces(parIndeces) {}

    double operator()(const double* par) const {
        double chi2 = 0;
        for (size_t iChi2 = 0; iChi2 < fChi2.size(); iChi2++) {
            double* pars = new double[fParIndeces[iChi2].size()];

            for (int iPar = 0; iPar < fParIndeces[iChi2].size(); iPar++) {
                pars[iPar] = par[fParIndeces[iChi2][iPar]];
            }
            chi2 += (*fChi2[iChi2])(pars);
        }

        return chi2;
    }

    const std::vector<ROOT::Fit::Chi2Function*> fChi2;
    std::vector<std::vector<int>> fParIndeces;
};

// return the indeces of the fit parameters, taking into account the shared ones
std::vector<std::vector<int>> GetParameterIndeces(std::vector<std::vector<sf::parameter>> fPars) {
    std::map<std::string, int> parIndeces = {};
    std::vector<std::vector<int>> iPars = {};
    int idx = 0;
    for (size_t iFit = 0; iFit < fPars.size(); iFit++) {
        std::vector<int> pars = {};
        for (size_t iPar = 0; iPar < fPars[iFit].size(); iPar++) {
            auto [name, _, __, ___] = fPars[iFit][iPar];

            // Find if the parameter was already inserted in the list of parameters
            auto it = parIndeces.find(name);
            if (it != parIndeces.end()) {
                pars.push_back(it->second);
            } else {
                pars.push_back(idx);
                parIndeces[name] = idx;
                idx++;
            }
        }
        iPars.push_back(pars);
    }

    return iPars;
}

// return the total number of fit parameters
int SuperFitter::GetN() {
    int nPars = 0;
    for (size_t iParList = 0; iParList < fPars.size(); iParList++) {
        nPars += fPars[iParList].size();
    }
    return nPars;
}

// Returns the number of shared fit parameters
int SuperFitter::GetNShared() { return GetN() - GetNIndependent(); }

// Returns the number of independent fit parameters
int SuperFitter::GetNIndependent() {
    std::vector<std::string> seen = {};
    int nIndependent = 0;
    for (size_t iFit = 0; iFit < fFit.size(); iFit++) {
        for (size_t iPar = 0; iPar < fFit[iFit]->GetNpar(); iPar++) {
            std::string name = fFit[iFit]->GetParName(iPar);
            if (std::find(seen.begin(), seen.end(), name) == seen.end()) {
                nIndependent++;
                seen.push_back(name);
            }
        }
    }
    return nIndependent;
}

// Returns the number of independent fit parameters for a specific fit
int SuperFitter::GetNIndependent(int nFit) {
    std::vector<std::string> seen = {};
    int nIndependent = 0;
    for (size_t iFit = 0; iFit < nFit; iFit++) {
        for (size_t iPar = 0; iPar < fFit[iFit]->GetNpar(); iPar++) {
            std::string name = fFit[iFit]->GetParName(iPar);
            if (std::find(seen.begin(), seen.end(), name) == seen.end()) {
                nIndependent++;
                seen.push_back(name);
            }
        }
    }
    return nIndependent;
}

// Returns the values of the fit parameters at initialization
std::vector<double> SuperFitter::GetInitialParameters() {
    std::vector<double> pars = {};
    std::vector<std::string> seen = {};
    for (int iFit = 0; iFit < fFit.size(); iFit++) {
        for (int iPar = 0; iPar < this->fFit[iFit]->GetNpar(); iPar++) {
            auto par = this->fPars[iFit][iPar];
            std::string name = std::get<0>(par);
            if (std::find(seen.begin(), seen.end(), name) == seen.end()) {
                double centr = std::get<1>(par);
                pars.push_back(centr);
                seen.push_back(name);
            }
        }
    }
    return pars;
}
// Fit
void SuperFitter::Fit(const char* option) {
    if (fFitRange.size() == 0) {
        throw std::runtime_error("Fit range is not specified!");
    }

    printf("Starting fit in Fit Range: ");
    for (size_t iRange = 0; iRange < this->fFitRange.size(); iRange++) {
        const auto& [xMin, xMax] = this->fFitRange[iRange];
        if (iRange) {
            printf(" U ");
        }
        printf("[%.3f %.3f]", xMin, xMax);
    }
    printf("\n");

    // Count the number of parameters:
    int nPars = GetN();
    int nShared = GetNShared();
    printf("\nPerforming %zu fits simultaneously with %d parameters of which %d are shared\n", fFit.size(), nPars,
           nShared);

    ROOT::Fit::DataOptions opt;
    ROOT::Fit::DataRange range;
    range.SetRange(fFitRange[0].first, fFitRange[fFitRange.size() - 1].second);

    std::vector<ROOT::Fit::BinData> data = {};
    std::vector<ROOT::Math::WrappedMultiTF1> wf = {};
    std::vector<ROOT::Fit::Chi2Function*> chi2Func = {};

    // Prepare machinery for custom global chi2
    for (size_t iFit = 0; iFit < fFit.size(); iFit++) {
        data.push_back(ROOT::Fit::BinData(opt, range));
        wf.push_back(ROOT::Math::WrappedMultiTF1(*(fFit[iFit]), 1));
        ROOT::Fit::FillData(data[iFit], fObs[iFit]->GetHistogram());
        chi2Func.push_back(new ROOT::Fit::Chi2Function(data[iFit], wf[iFit]));
    }

    auto iPars = GetParameterIndeces(this->fPars);
    GlobalChi2 globalChi2(chi2Func, iPars);

    ROOT::Fit::Fitter fitter;

    std::vector<double> pars = GetInitialParameters();
    fitter.Config().SetParamsSettings(nPars - nShared, pars.data());

    // Count the number of fit parameters
    int idx = 0;
    std::vector<std::string> seen = {};
    for (int iFit = 0; iFit < fFit.size(); iFit++) {
        for (int iPar = 0; iPar < this->fFit[iFit]->GetNpar(); iPar++) {
            auto par = this->fPars[iFit][iPar];
            std::string name = std::get<0>(par);
            if (std::find(seen.begin(), seen.end(), name) != seen.end()) {
                continue;
            }

            double centr = std::get<1>(par);
            double min = std::get<2>(par);
            double max = std::get<3>(par);
            seen.push_back(name);

            fitter.Config().ParSettings(idx).SetName(name.data());

            // Set Par Limits
            if (min > max) {
                fitter.Config().ParSettings(idx).SetValue(centr);
                fitter.Config().ParSettings(idx).Fix();
            } else {
                if (!(min < centr && centr < max)) {
                    printf("\033[33mWARNING: parameter '%s' is outside the allowed range\033[0m\n", name.data());
                    centr = (min + max) / 2;
                }

                fitter.Config().ParSettings(idx).SetValue(centr);
                fitter.Config().ParSettings(idx).SetLimits(min, max);
            }
            idx++;
        }
    }

    fitter.Config().MinimizerOptions().SetPrintLevel(0);
    fitter.Config().SetMinimizer("Minuit2", "Migrad");
    fitter.FitFCN(nPars - nShared, globalChi2, nullptr, data[0].Size() + data[1].Size(), true);
    ROOT::Fit::FitResult result = fitter.Result();
    result.Print(std::cout);
}

// Add TF1 function // todo: remove units mult here and put in .py
void SuperFitter::Add(int idx, std::string name, TF1* fTemplate, std::vector<sf::parameter> pars, double unitMult) {
    auto lambda = [&, fTemplate, unitMult, this](double* x, double* p) {
        return p[0] * fTemplate->Eval(x[0] * unitMult);
    };
    functions[idx].push_back({name, lambda, 1});

    // Save fit settings
    printf("Adding '%s' function with parameters:\n", name.data());
    for (const auto& par : pars) {
        auto [name, centr, min, max] = par;
        printf("    name: %s   init: %.3f   min: %.3f   max: %.3f\n", name.data(), centr, min, max);
        if (!IsParameterPresent(name)) {
            this->fPars[idx].push_back(par);
        }
    }
};

// Add template function
void SuperFitter::Add(int idx, std::string name, TH1* hTemplate, std::vector<sf::parameter> pars) {
    auto lambda = [hTemplate](double* x, double* p) { return p[0] * hTemplate->Interpolate(x[0]); };
    functions[idx].push_back({name, lambda, 1});

    // Save fit settings
    printf("Adding '%s' template with parameters:\n", name.data());
    for (const auto& par : pars) {
        auto [name, centr, min, max] = par;
        printf("    name: %s   init: %.3f   min: %.3f   max: %.3f\n", name.data(), centr, min, max);
        if (!IsParameterPresent(name)) {
            this->fPars[idx].push_back(par);
        }
    }
};

// Draw
void SuperFitter::Draw(int iFit, std::vector<std::pair<std::string, std::string>> recipes, std::string dataLabel) {
    this->fTerms = {};

    printf("Start drawing\n");

    TLegend* leg = new TLegend(0.55, 0.9 - 0.05 * recipes.size(), 0.9, 0.9);

    // Draw the fitted observable
    this->fObs[iFit]->Draw("hist same pe");
    leg->AddEntry(this->fObs[iFit], dataLabel.data(), "pe");

    // Draw the final fit function
    this->fFit[iFit]->Draw("same");

    leg->AddEntry(this->fFit[iFit], "Total");
    // Draw components based on the draw recipes
    for (int iRecipe = 0; iRecipe < recipes.size(); iRecipe++) {
        std::string legend = recipes[iRecipe].first;
        std::string recipe = recipes[iRecipe].second;
        DEBUG(60, 0, "Drawing the recipe '%s'", recipe.data());

        // Tokenization of the recipe
        auto tokens = Tokenize(recipe);
        DEBUG(60, 0, "[DRAW] recipe in infix: %s", join(" ", tokens).data());

        // Count the number of parameters
        std::set<int> paraList = {};
        std::vector<int> nParsDraw = {};
        std::set<std::string> used_tokens = {};
        for (const auto& token : tokens) {
            DEBUG(61, 1, "Processing token '%s'", token.data());

            if (!IsFunction(token)) {
                DEBUG(62, 2, "Token '%s' is not a function --> skip!", token.data());
                continue;
            }

            int counter = GetIndex(functions[iFit], token);
            int offset = ComputeOffset(functions[iFit], counter);

            // Determine the number of parameters
            for (int iFunc = 0; iFunc < functions[iFit].size(); iFunc++) {
                auto name = std::get<0>(functions[iFit][iFunc]);
                DEBUG(62, 2, "Comparing with function '%s'", name.data());
                if (name == token && used_tokens.find(token) == used_tokens.end()) {
                    int nPars = std::get<2>(functions[iFit][iFunc]);
                    DEBUG(63, 3, "It's a match! Number of parameters: %d", nPars);
                    nParsDraw.push_back(nPars);
                    // Determine the position of the function in the list of functions
                    for (int iPar = 0; iPar < nPars; iPar++) {
                        paraList.insert(offset + iPar);
                    }
                    break;
                }
            }
        }

        DEBUG(60, 0, "ParaList:");
        for (const auto& e : paraList) {
            DEBUG(61, 1, "%d value: %.3f", e, this->fFit[iFit]->GetParameter(e));
        }

        // Convert to Reverse Polish Notation
        auto rpn = toRPN(tokens);
        DEBUG(60, 0, "[DRAW] Expression in RPN: %s", join(" ", rpn).data());

        for (const int& d : nParsDraw) {
            DEBUG(61, 1, "par draw: %d", d);
        }

        // The following lambda evaluates the fit function
        auto lambda = [this, rpn, nParsDraw, iFit](double* x, double* p) -> double {
            std::stack<double> stack;

            DEBUG(60, 0, "[DRAW] Compute fit function from RPN: '%s'", join(" ", rpn).data());
            std::vector<std::pair<std::string, int>> nParameters = {};  // token, npars
            int idx = 0;
            for (const std::string& token : rpn) {
                DEBUG(61, 1, "[DRAW] Start processing the token '%s'", token.data());
                if (stack.size() == 0) {
                    DEBUG(61, 1, "[DRAW] Stack is empty");
                } else {
                    DEBUG(61, 1, "[DRAW] Stack is: '%s'", join(" ", stack_to_vector(stack)).data());
                }

                if (isdigit(token[0]) || token[0] == '.') {
                    DEBUG(62, 2, "[DRAW] Token '%s' is a number", token.data());
                    // Push numbers
                    stack.push(std::stod(token));
                } else if (IsFunction(token)) {
                    DEBUG(62, 2, "[DRAW] Token '%s' is a function with %d parameters", token.data(), 1);

                    // inly insert if not already present -> avoid duplicates
                    if (std::find(nParameters.begin(), nParameters.end(), std::pair(token, 1)) == nParameters.end()) {
                        for (const auto& [name, _, npar] : functions[iFit]) {
                            if (name == token) {
                                nParameters.push_back({token, npar});
                            }
                        }
                    }

                    // Compute offset
                    int shift = 0;
                    for (const auto& [name, np] : nParameters) {
                        if (name == token) break;
                        shift += np;
                    }
                    
                    int counter = GetIndex(functions[iFit], token);

                    // Determine the position of the function in the list of functions
                    auto func = std::get<1>(functions[iFit][counter]);
                    double value = func(x, p + shift);

                    DEBUG(62, 2, "[DRAW] Counter: %d/%zu    Offset: %d", counter, functions[iFit].size(), shift);
                    DEBUG(62, 2, "[DRAW] Pushing %s(x=%.3f, p) = %.3f", token.data(), x[0], value);

                    stack.push(value);
                } else if (IsOperator(token)) {
                    ProcessOperatorToken(stack, token);
                } else {
                    DEBUG(62, 2, "[DRAW] Token '%s' is unknown", token.data());
                    throw std::runtime_error("Unknown token: " + token);
                }

                DEBUG(61, 1, "[DRAW] Stack after processing the token: '%s'", join(" ", stack_to_vector(stack)).data());
            }

            DEBUG(60, 0, "[DRAW] End of evaluation, return value is: '%s'", join(" ", stack_to_vector(stack)).data());

            if (stack.size() != 1) throw std::runtime_error("Invalid RPN expression");
            return stack.top();
        };

        DEBUG(60, 0, "Term '%s' needs %lu parameters", recipe.data(), paraList.size());

        TF1* fTerm =
            new TF1(Form("fTerm%d", iRecipe), lambda, this->fDrawRangeMin, this->fDrawRangeMax, paraList.size());

        fTerm->SetLineColor(colors[iRecipe]);
        fTerm->SetLineWidth(3);
        fTerm->SetNpx(1000);

        int counter = 0;
        for (const int& par : paraList) {
            DEBUG(62, 2, "parameter val: %.3f", this->fFit[iFit]->GetParameter(par));
            fTerm->FixParameter(counter++, this->fFit[iFit]->GetParameter(par));
        }
        fTerm->Draw("same");
        fTerms.push_back(fTerm);
        leg->AddEntry(fTerm, legend.data(), "l");
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
