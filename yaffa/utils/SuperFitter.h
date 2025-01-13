#ifndef SUPERFITTER_H
#define SUPERFITTER_H

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

#if DO_DEBUG
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

// List of TF1-compatible functions that can be used in the fit
std::vector<std::tuple<std::string, std::function<double(double*, double*)>, int>> functions = {};

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

// Check if token is a function
bool IsFunction(const std::string& token) {
    for (const auto& [name, _, __] : functions) {
        if (token == name) return true;
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
            std::cout << "List of functions:" << std::endl;
            for (const auto &[fn, _, __] : functions) {
                std::cout << fn << std::endl;
            }
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

// Fit functions ---------------------------------------------------------------

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

const double FmToNu(5.067731237e-3);
const double Pi(3.141592653589793);
const std::complex<double> i(0,1);

double GeneralLednicky(double Momentum, const double &GaussR,
                       const complex<double> &ScattLenSin, const double &EffRangeSin) {

    // Taken from https://github.com/dimihayl/DLM/blob/c40f03eac38006f89eac8e5fa1533c9e48f2b455/CATS_Extentions/DLM_CkModels.cpp#L215
    if (GaussR != GaussR)
    {
        printf("\033[1;33mWARNING:\033[0m GeneralLednicky got a bad value for the Radius (nan). Returning default value of 1.\n");
        return 1;
    }

    Momentum *= 1000; // change units to GeV/c

    const double Radius = GaussR * FmToNu;
    const complex<double> IsLen1 = 1. / (ScattLenSin * FmToNu + 1e-64);
    const double eRan1 = EffRangeSin * FmToNu;

    double F1 = gsl_sf_dawson(2. * Momentum * Radius) / (2. * Momentum * Radius);
    double F2 = (1. - exp(-4. * Momentum * Momentum * Radius * Radius)) / (2. * Momentum * Radius);
    complex<double> ScattAmplSin = pow(IsLen1 + 0.5 * eRan1 * Momentum * Momentum - i * Momentum, -1.);

    double CkValue = 0.;
    CkValue += 0.5 * pow(abs(ScattAmplSin) / Radius, 2) * (1. - (eRan1) / (2 * sqrt(Pi) * Radius)) +
               2 * real(ScattAmplSin) * F1 / (sqrt(Pi) * Radius) - imag(ScattAmplSin) * F2 / Radius;
    CkValue += 1;

    return CkValue;
}

double Lednicky(double *x, double *par)
{
    // Taken from https://github.com/dimihayl/DLM/blob/c40f03eac38006f89eac8e5fa1533c9e48f2b455/CATS_Extentions/DLM_CkModels.cpp#L504C8-L504C53
    double kStar = x[0];
    double potPar0 = par[0];     // real part of the scattering length
    double potPar1 = par[1];     // imaginary part of the scattering length
    double potPar2 = par[2];     // effective range
    double sourcePar0 = par[3];  // radius of the first gaussian
    double sourcePar1 = par[4];  // radius of the second gaussian
    double sourcePar2 = par[5];  // relative weight of the two gaussians
    double sourcePar3 = par[6];  // normalization of the gaussians
    complex<double> ScatLen(potPar0, potPar1);
    return sourcePar3 * (sourcePar2 * GeneralLednicky(kStar, sourcePar0, ScatLen, potPar2) 
           + (1 - sourcePar2) * GeneralLednicky(kStar, sourcePar1, ScatLen, potPar2)) + 1. - sourcePar3;
}

// Class for advanced fitting
class SuperFitter : public TObject {
   private:
    Observable* fObs;                                                    // Observable to be fitted
    TF1* fFit;                                                           // Total fit function
    std::vector<std::tuple<std::string, double, double, double>> fPars;  // List of fit pars: (name, init, min, max)
    double fMin;                                                         // Fit range minimum
    double fMax;                                                         // Fit range maximum

   public:
    // Empty Contructor
    SuperFitter() : TObject(), fObs(nullptr), fFit(nullptr), fPars({}), fMin(0), fMax(1) {};

    // Standard Contructor
    SuperFitter(Observable* oObservable, double xMin, double xMax)
        : TObject(), fObs(oObservable), fFit(nullptr), fPars({}), fMin(xMin), fMax(xMax) {
    };

    // Destructor
    ~SuperFitter() { delete fObs; };

    // Fit
    void Fit(std::string model, const char* opt = "") {
        // Tokenization of the model
        auto tokens = Tokenize(model);
        DEBUG(0, "Expression in infix: %s", join(" ", tokens).data());

        // Convert to Reverse Polish Notation
        auto rpn = toRPN(tokens);
        DEBUG(0, "Expression in RPN: %s", join(" ", rpn).data());

        // The following lambda evaluates the fit function
        auto lambda = [this, rpn](double* x, double* p) -> double {
            std::stack<double> stack;

if(false)
            DEBUG(0, "Compute fit function from RPN: '%s'", join(" ", rpn).data());
            for (const std::string& token : rpn) {
                if(false)
                DEBUG(1, "Start processing the token '%s'", token.data());
                if (stack.size() == 0) {
                    if(false)
                    DEBUG(1, "Stack is empty");
                } else {
                    if(false)
                    DEBUG(1, "Stack is: '%s'", join(" ", stack_to_vector(stack)).data());
                }
                if (isdigit(token[0]) || token[0] == '.') {
                    if(false)
                    DEBUG(2, "Token '%s' is a number", token.data());
                    // Push numbers
                    stack.push(std::stod(token));
                } else if (IsFunction(token)) {
                    if(false)
                    DEBUG(2, "Token '%s' is a function", token.data());

                    // Determine the position of the function in the list of functions
                    int counter = 0;
                    for (const auto& [name, _, __] : functions) {
                        if (name == token) break;
                        counter++;
                    }

                    int offset = 0;
                    for (int iFunc = 0; iFunc < counter; iFunc++) {
                        offset += std::get<2>(functions[iFunc]);
                    }
                    if(false)
                    DEBUG(2, "Function '%s' was inserted in position: %d ==> Skipping %d parameters", token.data(),
                          counter, offset);

                    auto func = std::get<1>(functions[counter]);
                    double value = func(x, p + offset);
                    if(false)
                    DEBUG(2, "Pushing %s(x=%.3f, p) = %.3f", token.data(), x[0], value);

                    stack.push(value);
                } else if (IsOperator(token)) {
                    if(false)
                    DEBUG(2, "Token '%s' is an operator", token.data());
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
                    if(false)
                    DEBUG(2, "Token '%s' is unknown", token.data());
                    throw std::runtime_error("Unknown token: " + token);
                }
                if(false)
                DEBUG(1, "Stack after processing the token: '%s'", join(" ", stack_to_vector(stack)).data());
            }

if(false)
            DEBUG(0, "End of evaluation, return value is: '%s'", join(" ", stack_to_vector(stack)).data());

            if (stack.size() != 1) throw std::runtime_error("Invalid RPN expression");
            return stack.top();
        };

        int nPars = 0;
        for (int iFunc = 0; iFunc < functions.size(); iFunc++) {
          nPars += std::get<2>(functions[iFunc]);
        }
        this->fFit = new TF1("fFit", lambda, this->fMin, this->fMax, nPars);
        this->fFit->SetNpx(1000);
        
        for (int iPar = 0; iPar < this->fPars.size(); iPar++) {
            auto par = this->fPars[iPar];
            std::string name = std::get<0>(par);
            double centr = std::get<1>(par);
            double min = std::get<2>(par);
            double max = std::get<3>(par);

            this->fFit->SetParName(iPar, name.data());

            if (min > max) {
                this->fFit->FixParameter(iPar, centr);
            } else {
                if (!(min < centr && centr < max)) {
                    printf("\033[33mWARNING: parameter '%s' is outside the allowed range\033[0m\n", name.data());
                    centr = (min + max) / 2;
                }

                this->fFit->SetParameter(iPar, centr);
                this->fFit->SetParLimits(iPar, min, max);
            }
        }
        this->fObs->Fit(this->fFit, opt);
    };

    // Add fit component
    void Add(std::string name, std::string func, std::vector<std::tuple<std::string, double, double, double>> pars) {
        DEBUG(0, "Adding a new function '%s' to the fitter", name.data());

        if (func == "pol0") {
            functions.push_back({name, Pol0, 1});
        } else if (func == "pol1") {
            functions.push_back({name, Pol1, 2});
        } else if (func == "pol2") {
            functions.push_back({name, Pol2, 3});
        } else if (func == "pol3") {
            functions.push_back({name, Pol3, 4});
        } else if (func == "pol4") {
            functions.push_back({name, Pol4, 5});
        } else if (func == "pol5") {
            functions.push_back({name, Pol5, 6});
        } else if (func == "pol6") {
            functions.push_back({name, Pol6, 7});
        } else if (func == "pol7") {
            functions.push_back({name, Pol7, 8});
        } else if (func == "pol8") {
            functions.push_back({name, Pol8, 9});
        } else if (func == "pol9") {
            functions.push_back({name, Pol9, 10});
        } else if (func == "gaus") {
            functions.push_back({name, Gaus, 3});
} else if (func == "lednicky") {
            functions.push_back({name, Lednicky, 7});
        } else {
            throw std::runtime_error("Function " + func + " with name " + name + " is not implemented");
        }

        // Save fit settings
        DEBUG(0, "Parameters are:");
        for (const auto& par : pars) {
            DEBUG(1, "name: %s   init: %.3f   min: %.3f   max: %.3f", std::get<0>(par).data(), std::get<1>(par),
                  std::get<2>(par), std::get<3>(par));
            this->fPars.push_back(par);
        }
    };

    // Add template function
    void Add(std::string name, TH1* hTemplate, std::vector<std::tuple<std::string, double, double, double>> pars) {
        DEBUG(0, "Adding the template '%s' to the fitter", name.data());

        functions.push_back({name, [hTemplate](double* x, double* p) { return p[0] * hTemplate->Interpolate(x[0]); }, 1});

        // Save fit settings
        DEBUG(0, "Parameters are:");
        for (const auto& par : pars) {
            DEBUG(1, "name: %s   init: %.3f   min: %.3f   max: %.3f", std::get<0>(par).data(), std::get<1>(par),
                  std::get<2>(par), std::get<3>(par));
            this->fPars.push_back(par);
        }
    };

    // Add TF1 function
    void Add(std::string name, TF1* fTemplate, std::vector<std::tuple<std::string, double, double, double>> pars, double unitMult) { // todo: remove units mult here and put in .py
        DEBUG(0, "Adding the function '%s' to the fitter", name.data());
        std::cout << "kekek" << fTemplate << std::endl;

        functions.push_back({name, [&, fTemplate, unitMult, this](double* x, double* p) { return p[0] * fTemplate->Eval(x[0] * unitMult); }, 1});


        // Save fit settings
        DEBUG(0, "Parameters are:");
        for (const auto& par : pars) {
            DEBUG(1, "name: %s   init: %.3f   min: %.3f   max: %.3f", std::get<0>(par).data(), std::get<1>(par),
                  std::get<2>(par), std::get<3>(par));
            this->fPars.push_back(par);
        }
    };

    // Draw
    void Draw(std::vector<std::pair<std::string, std::string>> recipes) const {
        std::cout << "Draw functions" << std::endl;
        for (const auto& [name, _, __] : functions) {
            std::cout << name << std::endl;
        }

        DEBUG(0, "Start drawing");

        TLegend *leg = new TLegend(0.6, 0.6, 0.9, 0.9);

        // Draw the fitted observable
        this->fObs->Draw("hist same pe");
        leg->AddEntry(this->fObs, "data", "pe");

        // Draw the final fit function
        this->fFit->Draw("same");

        leg->AddEntry(this->fFit, "Total");
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
            for (const auto &token : tokens) {
                DEBUG(1, "Processing token '%s'", token.data());

                if (!IsFunction(token)) {
                    DEBUG(2, "Token '%s' is not a function --> skip!", token.data());
                    continue;
                }

                int counter = 0;
                for (const auto& [name, _, __] : functions) {
                    if (name == token) break;
                    counter++;
                }

                int offset = 0;
                for (int iFunc = 0; iFunc < counter; iFunc++) {
                    offset += std::get<2>(functions[iFunc]);
                }

                // Determine the number of parameters
                for (int iFunc = 0; iFunc < functions.size(); iFunc++) {
                    auto name = std::get<0>(functions[iFunc]);
                    DEBUG(2, "Comparing with function '%s'", name.data());
                    if (name == token && used_tokens.find(token) == used_tokens.end()) {
                        int nPars = std::get<2>(functions[iFunc]);
                        DEBUG(3, "It's a match! Number of parameters: %d", nPars);
                        nParsDraw.push_back(nPars);
                        // Determine the position of the function in the list of functions
                        for (int iPar = 0; iPar< nPars; iPar++) {
                            paraList.insert(offset + iPar);
                        }
                        break;
                    }
                }
            }

            DEBUG(0, "ParaList:");
            for (const auto &e : paraList) {

                DEBUG(1, "%d value: %.3f", e, this->fFit->GetParameter(e));
            }

            // Convert to Reverse Polish Notation
            auto rpn = toRPN(tokens);
            DEBUG(0, "[DRAW] Expression in RPN: %s", join(" ", rpn).data());

            for (const int &d : nParsDraw) {
                DEBUG(1, "par draw: %d", d);
            }
            

            // The following lambda evaluates the fit function
            auto lambda = [this, rpn, nParsDraw](double* x, double* p) -> double {
                std::stack<double> stack;

if (0.12 < x[0] && x[0] < .16)
                DEBUG(0, "[DRAW] Compute fit function from RPN: '%s'", join(" ", rpn).data());
                std::vector<std::pair<std::string, int>> nParameters = {
                    // {"pol0", 1},
                    // {"ca", 1},
                    // {"nca", 1},
                    // // {"Sigma1385_direct", 1},
                    // {"Xi1530", 1},
                }; // token, npars
                int idx = 0;
                for (const std::string& token : rpn) {
                    if (0.12 < x[0] && x[0] < .16)
                    DEBUG(1, "[DRAW] Start processing the token '%s'", token.data());
                    if (stack.size() == 0) {
                        if (0.12 < x[0] && x[0] < .16)
                        DEBUG(1, "[DRAW] Stack is empty");
                    } else {
                        if (0.12 < x[0] && x[0] < .16)
                        DEBUG(1, "[DRAW] Stack is: '%s'", join(" ", stack_to_vector(stack)).data());
                    }

                    if (isdigit(token[0]) || token[0] == '.') {
                        if (0.12 < x[0] && x[0] < .16)
                        DEBUG(2, "[DRAW] Token '%s' is a number", token.data());
                        // Push numbers
                        stack.push(std::stod(token));
                    } else if (IsFunction(token)) {
                        if (0.12 < x[0] && x[0] < .16)
                        DEBUG(2, "Inserting token; %s   npar: %d\n", token.data(), 1);

                        // inly insert if not already present -> avoid duplicates
                        if (std::find(nParameters.begin(), nParameters.end(), std::pair(token, 1)) == nParameters.end()) {
                            for (const auto &[name, _, npar] : functions) {
                                if (name == token) {
                                    nParameters.push_back({token, npar});
                                }
                            }
                        }

                        // Compute offset
                        int offset = 0;
                        for (const auto &[name, np] : nParameters) {
                            if (name == token) break;
                            offset += np;
                        }
                        
                        if (0.12 < x[0] && x[0] < .16)
                        DEBUG(2, "[DRAW] Token '%s' is a function", token.data());

                        // Determine the position of the function in the list of functions
                        int counter = 0;
                        for (const auto& [name, _, __] : functions) {
                            if (name == token) break;
                            counter++;
                        }
                        // todo: why no offset?
                        auto func = std::get<1>(functions[counter]);
                        // int offset = std::accumulate(nParsDraw.begin(), nParsDraw.begin() + counter, 0);
                        // if (counter == 4) offset--;
                        // if (counter == 4) offset = 1;
                        double value = func(x, p + offset);
                        // value = func(x, p);

                        if (0.12 < x[0] && x[0] < .16)
                        DEBUG(2, "Counter: %d/%d    Offset: %d", counter, functions.size(), offset);
                        if (0.12 < x[0] && x[0] < .16)
                        DEBUG(2, "[DRAW] Pushing %s(x=%.3f, p) = %.3f", token.data(), x[0], value);

                        stack.push(value);
                    } else if (IsOperator(token)) {
                        if (0.12 < x[0] && x[0] < .16)
                        DEBUG(2, "[DRAW] Token '%s' is an operator", token.data());
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
                        if (0.12 < x[0] && x[0] < .16)
                        DEBUG(2, "[DRAW] Token '%s' is unknown", token.data());
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

            TF1 *fTerm = new TF1(Form("fTerm%d", iRecipe), lambda, this->fMin, this->fMax, paraList.size());
            int color = iRecipe + 2;
            if (color >= 5) color++;
            fTerm->SetLineColor(color);
            fTerm->SetNpx(1000);

            int counter = 0;
            for (const int &par : paraList) {
                DEBUG(2, "parameter val: %.3f",  this->fFit->GetParameter(par));
                fTerm->FixParameter(counter++, this->fFit->GetParameter(par));
            }
            fTerm->Draw("same");
            leg->AddEntry(fTerm, legend.data());
        }
        leg->Draw("same");
    };

    TF1* GetFitFunction() {
        return this->fFit;
    }

    ClassDef(SuperFitter, 1)
};

ClassImp(SuperFitter);

#endif
