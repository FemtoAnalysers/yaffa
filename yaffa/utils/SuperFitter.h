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

#if true
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

// Lambda function to concatenate the elements of a std::vector via a separator. Equivalent of python's " ".join(mylist)
std::string join(std::string separator, std::vector<std::string> list) {
    std::string output = "";
    for (const auto &element : list) {
        output += element + separator;
    }

    if (output.size() > 0) output.resize(output.size() - separator.size()); // Remove last separator
    return output;
};

std::string join(std::string separator, std::vector<double> list) {
    std::string output = "";
    for (const auto &element : list) {
        output += std::to_string(element) + separator;
    }
    if (output.size() > 0) output.resize(output.size() - separator.size()); // Remove last separator
    return output;
};

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

double Gaus(double* x, double* p) {
    double xx = x[0];
    double norm = p[0];
    double mean = p[1];
    double sigma = p[2];

    double normFactor = norm / (std::sqrt(2 * M_PI) * sigma);
    double exponent = -0.5 * std::pow((xx - mean) / sigma, 2);
    return normFactor * std::exp(exponent);
}

double Poll0(double x, double p) { return p; }
double Pol0(double* x, double* p) { return p[0]; }
double Pol1(double* x, double* p) { return Pol0(x, p) + p[1] * pow(x[0], 1); }
double Pol2(double* x, double* p) { return Pol1(x, p) + p[2] * pow(x[0], 2); }
double Pol3(double* x, double* p) { return Pol2(x, p) + p[3] * pow(x[0], 3); }
double Pol4(double* x, double* p) { return Pol3(x, p) + p[4] * pow(x[0], 4); }
double Pol5(double* x, double* p) { return Pol4(x, p) + p[5] * pow(x[0], 5); }
double Pol6(double* x, double* p) { return Pol5(x, p) + p[6] * pow(x[0], 6); }
double Pol7(double* x, double* p) { return Pol6(x, p) + p[7] * pow(x[0], 7); }
double Pol8(double* x, double* p) { return Pol7(x, p) + p[8] * pow(x[0], 8); }
double Pol9(double* x, double* p) { return Pol8(x, p) + p[9] * pow(x[0], 9); }

class SuperFitter : public TObject {
   private:
    Observable* fObs;
    TF1* fFit;
    std::vector<TF1*> fTerms;
    std::map<int, std::tuple<std::string, double, double, double>> fPars;  // List of fit parameters: (index, (name, init, min, max))
    double fMin;
    double fMax;

   public:
    // Empty Contructor
    SuperFitter() : TObject(), fObs(nullptr), fFit(nullptr), fPars({}) {};

    // Standard Contructor
    SuperFitter(Observable* hObs, double xMin, double xMax);

    // Destructor
    ~SuperFitter();

    // Fit
    void Fit(std::string model, const char* opt = "");

    // Add fit component
    void Add(std::string name, std::vector<std::tuple<std::string, double, double, double>> pars);

    // Draw
    void Draw(const char* opt = "") const;

    ClassDef(SuperFitter, 1)
};

ClassImp(SuperFitter);

// Empty Constructor
SuperFitter::~SuperFitter() { delete fObs; }

// Standard Constructor
SuperFitter::SuperFitter(Observable* hObs, double xMin, double xMax) : TObject(), fPars({}) {
    this->fObs = hObs;
    this->fMin = xMin;
    this->fMax = xMax;
}

using Function = std::function<double(double*, double*)>;
std::map<std::string, Function> functions = {
    {"sin", [](double *x, double *p) { return sin(x[0]); }},
    {"pol1", Pol1},
    {"pol2", Pol2},
};

// Tokenizer function
std::vector<std::string> Tokenize(const std::string& str, const std::string& delimiters) {
    std::vector<std::string> tokens;
    size_t start = 0;
    size_t end;

    while ((end = str.find_first_of(delimiters, start)) != std::string::npos) {
        if (end != start) {  // Avoid empty tokens
            tokens.emplace_back(str.substr(start, end - start));
        }
        start = end + 1;
    }

    if (start < str.length()) {  // Add the last token
        tokens.emplace_back(str.substr(start));
    }

    return tokens;
}

// Tokenize the formula
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
        } else if (isalpha(c)) {
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

// Utility: Get operator precedence
int getPrecedence(const std::string& op) {
    if (op == "+" || op == "-") return 1;
    if (op == "*" || op == "/") return 2;
    return 0;
}

// Utility: Check if token is an operator
bool isOperator(const std::string& token) { return token == "+" || token == "-" || token == "*" || token == "/"; }

// Utility: Check if token is a function
bool isFunction(const std::string& token) { return functions.find(token) != functions.end(); }

// Convert to Reverse Polish Notation (RPN)
std::vector<std::string> toRPN(const std::vector<std::string>& tokens) {
    std::vector<std::string> output;
    std::stack<std::string> operators;

    for (const std::string& token : tokens) {
        if (isdigit(token[0]) || token[0] == '.') {
            // Numbers go directly to output
            output.push_back(token);
        } else if (isFunction(token)) {
            // Functions go to the operator stack
            operators.push(token);
        } else if (token == "(") {
            operators.push(token);
        } else if (token == ")") {
            // Pop until matching "("
            while (!operators.empty() && operators.top() != "(") {
                output.push_back(operators.top());
                operators.pop();
            }
            if (!operators.empty()) operators.pop();  // Pop "("
            if (!operators.empty() && isFunction(operators.top())) {
                output.push_back(operators.top());
                operators.pop();
            }
        } else if (isOperator(token)) {
            while (!operators.empty() && getPrecedence(operators.top()) >= getPrecedence(token)) {
                output.push_back(operators.top());
                operators.pop();
            }
            operators.push(token);
        }
    }

    // Pop remaining operators
    while (!operators.empty()) {
        output.push_back(operators.top());
        operators.pop();
    }

    return output;
}

// Fit
void SuperFitter::Fit(std::string model, const char* opt) {
    // Tokenization of the model
    auto tokens = Tokenize(model);

    auto rpn = toRPN(tokens);

    DEBUG(0, "Expression in RPN: %s", join(" ", rpn).data());
    auto lambda = [rpn](double* x, double* p) -> double {
        std::stack<double> stack;

        DEBUG(0, "Compute fit function");
        for (const std::string& token : rpn) {
            if (isdigit(token[0]) || token[0] == '.') {
                DEBUG(2, "Token '%s' is a digit", token.data());
                // Push numbers
                stack.push(std::stod(token));
            } else if (isFunction(token)) {
                DEBUG(2, "Token '%s' is a function", token.data());
                DEBUG(2, "Pushing %.3f into the stack", functions[token](x, p));

                // Call function
                stack.push(functions[token](x, p));

                DEBUG(2, "Stack is: '%s'", join(" ", stack_to_vector(stack)).data());

            } else if (isOperator(token)) {
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
        }

        DEBUG(0, "Stack is: '%s'", join(" ", stack_to_vector(stack)).data());

        if (stack.size() != 1) throw std::runtime_error("Invalid RPN expression");
        return stack.top();
    };

    this->fFit = new TF1("fFit", lambda, 0, 0.5, 4);
    this->fFit->SetNpx(15);

    for (const auto &[iPar, par] : this->fPars) {
        this->fFit->FixParameter(iPar, std::get<1>(par));
    }
    this->fObs->Fit(this->fFit, opt);
}

// Add fit component
void SuperFitter::Add(std::string name, std::vector<std::tuple<std::string, double, double, double>> pars) {
    DEBUG(0, "Adding a new function '%s' to the fitter", name.data());
    TF1* fFunc = nullptr;
    if (name == "pol0") {
        fFunc = new TF1("fFit", Pol0, this->fMin, this->fMax, 1);
    } else if (name == "pol1") {
        fFunc = new TF1("fFit", Pol1, this->fMin, this->fMax, 2);
    } else if (name == "pol2") {
        fFunc = new TF1("fFit", Pol2, this->fMin, this->fMax, 3);
    } else if (name == "pol3") {
        fFunc = new TF1("fFit", Pol3, this->fMin, this->fMax, 4);
    } else if (name == "pol4") {
        fFunc = new TF1("fFit", Pol4, this->fMin, this->fMax, 5);
    } else if (name == "pol5") {
        fFunc = new TF1("fFit", Pol5, this->fMin, this->fMax, 6);
    } else if (name == "pol6") {
        fFunc = new TF1("fFit", Pol6, this->fMin, this->fMax, 7);
    } else if (name == "pol7") {
        fFunc = new TF1("fFit", Pol7, this->fMin, this->fMax, 8);
    } else if (name == "pol8") {
        fFunc = new TF1("fFit", Pol8, this->fMin, this->fMax, 9);
    } else if (name == "pol9") {
        fFunc = new TF1("fFit", Pol9, this->fMin, this->fMax, 10);
    } else if (name == "gaus") {
        fFunc = new TF1("fFit", Gaus, this->fMin, this->fMax, 3);
    } else {
        throw std::runtime_error("Function " + name + " is not implemented");
    }

    this->fTerms.push_back(fFunc);

    // Save fit settings
    DEBUG(0, "Parameters are:");
    for (const auto &par : pars) {
        DEBUG(1, "name: %s   init: %.3f   min: %.3f   max: %.3f", std::get<0>(par).data(), std::get<1>(par), std::get<2>(par), std::get<3>(par));
        this->fPars.insert({this->fPars.size(), par});
    }
}

// Draw
void SuperFitter::Draw(const char* opt) const { this->fObs->Draw(opt); }

#endif
