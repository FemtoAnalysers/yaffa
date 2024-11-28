#ifndef SUPERFITTER_H
#define SUPERFITTER_H

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

#if 0
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
std::vector<std::pair<std::string, std::function<double(double*, double*)>>> functions = {};

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
    for (const auto& [name, _] : functions) {
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

// Class for advanced fitting
class SuperFitter : public TObject {
   private:
    Observable* fObs;
    TF1* fFit;
    std::vector<TF1*> fTerms;
    std::vector<int> fNPars;
    std::map<int, std::tuple<std::string, double, double, double>>
        fPars;  // List of fit parameters: (index, (name, init, min, max))
    double fMin;
    double fMax;

   public:
    // Empty Contructor
    SuperFitter() : TObject(), fObs(nullptr), fFit(nullptr), fPars({}), fNPars({}) {};

    // Standard Contructor
    SuperFitter(Observable* hObs, double xMin, double xMax) : TObject(), fPars({}) {
        this->fObs = hObs;
        this->fMin = xMin;
        this->fMax = xMax;
    };

    // Destructor
    ~SuperFitter() { delete fObs; };

    // Fit
    void Fit(std::string model, const char* opt = "") {
        // Tokenization of the model
        auto tokens = Tokenize(model);

        DEBUG(0, "Expression in infix: %s", join(" ", tokens).data());
        auto rpn = toRPN(tokens);
        DEBUG(0, "Expression in RPN: %s", join(" ", rpn).data());

        // The following lambda evaluates the fit function
        auto lambda = [this, rpn](double* x, double* p) -> double {
            std::stack<double> stack;

            DEBUG(0, "Compute fit function from RPN: '%s'", join(" ", rpn).data());
            for (const std::string& token : rpn) {
                DEBUG(1, "Start processing the token '%s'", token.data());
                if (stack.size() == 0) {
                    DEBUG(1, "Stack is empty");
                } else {
                    DEBUG(1, "Stack is: '%s'", join(" ", stack_to_vector(stack)).data());
                }
                if (isdigit(token[0]) || token[0] == '.') {
                    DEBUG(2, "Token '%s' is a number", token.data());
                    // Push numbers
                    stack.push(std::stod(token));
                } else if (IsFunction(token)) {
                    DEBUG(2, "Token '%s' is a function", token.data());

                    // Determine the position of the function in the list of functions
                    int counter = 0;
                    for (const auto& [name, _] : functions) {
                        if (name == token) break;
                        counter++;
                    }

                    int offset = std::accumulate(this->fNPars.begin(), this->fNPars.begin() + counter, 0);
                    DEBUG(2, "Function '%s' was inserted in position: %d ==> Skipping %d parameters", token.data(),
                          counter, offset);

                    auto func = functions[counter].second;
                    double value = func(x, p + offset);
                    DEBUG(2, "Pushing %s(x=%.3f, p) = %.3f", token.data(), x[0], value);

                    stack.push(value);
                } else if (IsOperator(token)) {
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
                    DEBUG(2, "Token '%s' is unknown", token.data());
                    throw std::runtime_error("Unknown token: " + token);
                }
                DEBUG(1, "Stack after processing the token: '%s'", join(" ", stack_to_vector(stack)).data());
            }

            DEBUG(0, "End of evaluation, return value is: '%s'", join(" ", stack_to_vector(stack)).data());

            if (stack.size() != 1) throw std::runtime_error("Invalid RPN expression");
            return stack.top();
        };

        int nPars = std::accumulate(this->fNPars.begin(), this->fNPars.end(), 0);
        this->fFit = new TF1("fFit", lambda, 0, 0.5, nPars);

        for (const auto& [iPar, par] : this->fPars) {
            this->fFit->SetParName(iPar, std::get<0>(par).data());

            if (std::get<2>(par) > std::get<3>(par)) {
                this->fFit->FixParameter(iPar, std::get<1>(par));
            } else {
                this->fFit->SetParameter(iPar, std::get<1>(par));
                this->fFit->SetParLimits(iPar, std::get<2>(par), std::get<3>(par));
            }
        }
        this->fObs->Fit(this->fFit, opt);
    };

    // Add fit component
    void Add(std::string name, std::vector<std::tuple<std::string, double, double, double>> pars) {
        DEBUG(0, "Adding a new function '%s' to the fitter", name.data());

        if (name == "pol1") {
            functions.push_back({name, Pol1});
            fNPars.push_back(2);
        } else if (name == "gaus") {
            functions.push_back({name, Gaus});
            fNPars.push_back(3);
        } else {
            throw std::runtime_error("Function " + name + " is not implemented");
        }

        // Save fit settings
        DEBUG(0, "Parameters are:");
        for (const auto& par : pars) {
            DEBUG(1, "name: %s   init: %.3f   min: %.3f   max: %.3f", std::get<0>(par).data(), std::get<1>(par),
                  std::get<2>(par), std::get<3>(par));
            this->fPars.insert({this->fPars.size(), par});
        }
    };

    // Add template function
    void Add(std::string name, TH1* hTemplate, std::vector<std::tuple<std::string, double, double, double>> pars) {
        DEBUG(0, "Adding the template '%s' to the fitter", name.data());

        functions.push_back({name, [hTemplate](double* x, double* p) { return p[0] * hTemplate->Interpolate(x[0]); }});
        fNPars.push_back(1);  // Only normalization

        // Save fit settings
        DEBUG(0, "Parameters are:");
        for (const auto& par : pars) {
            DEBUG(1, "name: %s   init: %.3f   min: %.3f   max: %.3f", std::get<0>(par).data(), std::get<1>(par),
                  std::get<2>(par), std::get<3>(par));
            this->fPars.insert({this->fPars.size(), par});
        }
    };

    // Draw
    void Draw(const char* opt = "") const { this->fObs->Draw(opt); };

    ClassDef(SuperFitter, 1)
};

ClassImp(SuperFitter);

#endif
