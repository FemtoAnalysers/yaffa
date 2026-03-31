#ifndef UTILS_H
#define UTILS_H

#include <TH1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>


//**************Design of objects***********
void SetGraphInfo(TGraphErrors* graph, const int& color, const int& markerStyle, const double& markerSize, const double& xMin, const double& xMax, const char* name, const char* xTitle, const char* yTitle, const int& linestyle);
void SetGraphInfo(TGraph* graph, const int& color, const int& markerStyle, const double& markerSize, const double& xMin, const double& xMax, const char* name, const char* xTitle, const char* yTitle, const int& linestyle);

//**************Change object to graph***********
TGraphErrors* HistToGraph(TH1F* hist, double kmin, double kmax, int nBins);
TGraph* FitToGraph(TF1* fit, double kmin, double kmax, int nPoints);

#endif