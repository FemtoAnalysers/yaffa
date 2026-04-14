'''
A mix of utility functions
'''

from ROOT import TGraphErrors, TF1

def SliceVertically(hist, edges, name=None):
    '''
    Slice a TH2 vertically (ProjectionY) and return the list of slices
    '''

    slices = []
    lowEdges = edges[:-1]
    upEdges = edges[1:]
    if not name:
        name = hist.GetName()

    for lowEdge, upEdge in zip(lowEdges, upEdges):
        firstBin = hist.GetXaxis().FindBin(lowEdge * 1.0001)
        lastBin = hist.GetXaxis().FindBin(upEdge * 0.9999)
        slices.append(hist.ProjectionY(f'{name}{lowEdge:.0f}_{upEdge:.0f}', firstBin, lastBin))
    return slices

def FitHistList(hList:list, edges:list, func, nPar, range, name=None, normalize=False) -> TGraphErrors:
    gParameter = TGraphErrors()
    gParameter.SetName(name if name else 'gParameter')
    for iBin, hist in enumerate(hList):
        fFunc = TF1(f'f{iBin}', func, *range, nPar)
        fFunc.SetNpx(10000)
        for iPar in range(nPar):
            fFunc.SetParameter(iPar, 1)

        if normalize:
            hist.Scale(1./hist.GetEntries() / hist.GetBinWidth(1))
 
        hist.Fit(fFunc, 'MR+')
        hist.SetTitle(';#rho* (fm);Counts')
        hist.Write()
        fFunc.Write()

        # Fit results
        par = fFunc.GetParameter(0)
        parUnc = fFunc.GetParError(0)

        gParameter.SetPoint(iBin, (edges[iBin] + edges[iBin + 1]) / 2, par)
        gParameter.SetPointError(iBin, (edges[iBin] - edges[iBin + 1]) / 2, parUnc)
    # gParameter.SetTitle(';m_{T} (GeV);#rho_{0} (fm)')
    gParameter.Write()
    