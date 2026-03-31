#include "utils.h"

//**********Set the style of a Object***********
void SetGraphInfo(TGraph* graph, const int& color=kBlack, const int& markerStyle=kFullCircle, const double& markerSize=1.0, const double& xMin=0., const double& xMax=500., const char* name="name", const char* xTitle="x", const char* yTitle="y", const int& linestyle=1) {
    graph->SetMarkerColor(color);
    graph->SetLineColor(color);
    graph->SetLineStyle(linestyle);
    graph->SetLineWidth(3);
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerSize(markerSize);
    graph->SetTitle(name);
    graph->SetName(name);
    graph->GetXaxis()->SetRangeUser(xMin, xMax);
    graph->GetXaxis()->SetTitle(xTitle);
    graph->GetYaxis()->SetTitle(yTitle);
}

void SetGraphInfo(TGraphErrors* graph, const int& color=kBlack, const int& markerStyle=kFullCircle, const double& markerSize=1.0, const double& xMin=0., const double& xMax=500., const char* name="name", const char* xTitle="x", const char* yTitle="y", const int& linestyle=1) {
    graph->SetMarkerColor(color);
    graph->SetLineColor(color);
    graph->SetLineStyle(linestyle);
    graph->SetLineWidth(3);
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerSize(markerSize);
    graph->SetTitle(name);
    graph->SetName(name);
    graph->GetXaxis()->SetRangeUser(xMin, xMax);
    graph->GetXaxis()->SetTitle(xTitle);
    graph->GetYaxis()->SetTitle(yTitle);
}


// *******************change object to Graph *****************
TGraphErrors* HistToGraph(TH1F* hist, double kmin, double kmax, int nBins) {
    TGraphErrors* g = new TGraphErrors();
    for (int i = 0; i < nBins; i++) {
        double mon= hist->GetBinCenter(i+1);
        if(mon<kmin || mon>kmax) continue;
        g->SetPoint(i, mon, hist->GetBinContent(i+1));
        double xerror= (kmax-kmin)/(2*nBins);
        g->SetPointError(i, xerror, hist->GetBinError(i+1));
    }

    return g;
}

TGraph* FitToGraph(TF1* fit, double kmin, double kmax, int nPoints) {
    TGraph* g = new TGraph();

    for (int i = 0; i < nPoints; i++) {
        double mom = kmin + i*(kmax-kmin)/nPoints;
        g->SetPoint(i, mom, fit->Eval(mom));
    }

    return g;
}