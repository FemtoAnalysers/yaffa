'''
Module that contains various functions useful for the post-processing of the data.
'''

import numpy as np

from ROOT import TH1F, TH2D, TGraph  # pylint: disable=import-error

def GetSpread(objects):
    '''
    Given a list of TH1, returns a TH1 where the points are the averages of the bins and the width is the spread
    (standard deviation) of the bin entries.

    Args:
        objects (list[TH1]): list of histograms

    Returns:
        TH1: histogram containing the average and spread of the input histograms
    '''
    hSpread = objects[0].Clone('hSpread')
    hSpread.Reset()

    for iBin in range(hSpread.GetNbinsX()):
        yValues = np.array([obj.GetBinContent(iBin + 1) for obj in objects[1:]])
        spread = np.std(yValues)
        avg = np.average(yValues)
        hSpread.SetBinContent(iBin + 1, avg)
        hSpread.SetBinError(iBin + 1, spread)

    return hSpread


def ChangeUnits(hist, multiplier, name=None, title=''):
    '''
    Only for histogram with constant binwidth!
    '''
    nbins = hist.GetNbinsX()
    lowEdge = hist.GetBinLowEdge(1)
    upEdge = hist.GetBinLowEdge(nbins+1)
    if name is None:
        name = f'{hist.GetName()}_new'
    hNew = TH1F(name, title, nbins, lowEdge*multiplier, upEdge*multiplier)
    for i in range(0, nbins+2):
        hNew.SetBinContent(i, hist.GetBinContent(i))
        hNew.SetBinError(i, hist.GetBinError(i))
    return hNew


def ChangeUnits2D(hist, multiplier, name=None, title=''):
    '''
    Change the units of a 2D histogram. Useful to convert axes from/to MeV/GeV.
    Only for histogram with constant binwidth!

    Parameters
    ----------
    hist : TH2
        The histogram to be converted
    multiplier : float, (float, float)
        The multiplier to be applied to the axes. If only one is provided it is assumed that the same multiplier should
        be applied to the x and y axes. Provide a tuple (multX, multY) to have separate scaling factors.
    name : str, optional
        The name of the converted histogram. By default: `<old_name>_new`
    title : str, optional
        The title of the converted. By default: `''`

    Returns
    -------
    TH2D
        The histogram in the new units
    '''

    nbinsX = hist.GetNbinsX()
    nbinsY = hist.GetNbinsY()
    lowEdgeX = hist.GetXaxis().GetBinLowEdge(1)
    lowEdgeY = hist.GetXaxis().GetBinLowEdge(1)
    upEdgeX = hist.GetYaxis().GetBinLowEdge(nbinsX+1)
    upEdgeY = hist.GetYaxis().GetBinLowEdge(nbinsY+1)

    if name is None:
        name = f'{hist.GetName()}_new'

    if isinstance(multiplier, (int, float)):
        multX, multY = (multiplier, multiplier)
    else:
        multX, multY = multiplier

    hNew = TH2D(name, title, nbinsX, lowEdgeX * multX, upEdgeX * multX, nbinsX, lowEdgeY * multY, upEdgeY * multY)

    for iBinX in range(0, nbinsX+2):
        for iBinY in range(0, nbinsX+2):
            hNew.SetBinContent(iBinX, iBinY, hist.GetBinContent(iBinX, iBinY))
            hNew.SetBinError(iBinX, iBinY, hist.GetBinError(iBinX, iBinY))
    return hNew


def WeightedAverage(graph, weights):
    '''
    Compute the weighted average of a graph with an histogram (TH1)

    Parameters
    ----------
    graph : TGraph
        The graph to be reweighted
    weights : TH1
        The weights to be applied

    Returns
    -------
    TGraph
        The reweighted graph
    '''

    smeared = 0
    counts = weights.Integral(1, weights.GetNbinsX())
    for iBin in range(weights.GetNbinsX()):
        freq = weights.GetBinContent(iBin + 1) / counts
        y = graph.Eval(weights.GetBinCenter(iBin+1))
        smeared += freq * y
    return smeared


def SmearGraph(graph, matrix, name=None, title=''):
    '''
    Smear a graph with a smearing matrix that has:
     x axis: true variable
     y axis: reconstructed variable.

    Parameters
    ----------
    graph : TGraph
        The graph to be smeared
    matrix : TH2
        the smearing matrix
    name : str, optional
        The name of the smeared graph
    title : str, optional
        The title of the smeared graph


    Returns
    -------
    TGraph
        The smeared graph
    '''


    gSmeared = TGraph(1)
    if name:
        gSmeared.SetName(name)
    gSmeared.SetTitle(title)

    iPoint = 0
    for iBin in range(matrix.GetNbinsX()):
        hProj = matrix.ProjectionY(f'hProj_{iBin+1}', iBin+1, iBin+1)
        counts = hProj.Integral(1, hProj.GetNbinsX())

        if counts < 1:
            continue

        x = matrix.GetXaxis().GetBinCenter(iBin+1)
        if x > graph.GetPointX(graph.GetN() - 1):
            break

        ySmear = WeightedAverage(graph, hProj)
        gSmeared.SetPoint(iPoint, x, ySmear)
        iPoint += 1

    return gSmeared
