'''
Module that contains various functions useful for the post-processing of the data.
'''
import math

import numpy as np

from ROOT import TH1, TH1D, TH1F, TH1I, TH2D, TSpline3, TGraph, TH2, TGraphErrors, TF1  # pylint: disable=import-error

from yaffa import logger as log

def Convert(object, target_type):
    if type(object).__name__ == 'TH2D':
        if target_type == 'numpy_array':
            nx = object.GetNbinsX()
            ny = object.GetNbinsY()

            arr = np.empty((nx, ny))

            for ix in range(nx):
                for iy in range(ny):
                    arr[ix, iy] = object.GetBinContent(ix + 1, iy + 1)

            return arr

    log.critical(f'Conversion from type {type(object).__name__} to {target_type} is not implemented')

def ScaleGraph(graph, value, name=None):
    '''
    Scales a TGraphAsymmErrors by the specified value
    '''
    if name is None:
        name = f'{graph.GetName()}_scaled'

    gScaled = graph.Clone(name)

    scale = abs(value)

    for i in range(gScaled.GetN()):
        x = gScaled.GetPointX(i)
        y = gScaled.GetPointY(i)

        gScaled.SetPoint(i, x, y * value)

        gScaled.SetPointEYlow(i, gScaled.GetErrorYlow(i)  * scale)
        gScaled.SetPointEYhigh(i, gScaled.GetErrorYhigh(i) * scale)

    return gScaled

def SliceVertically(hist, edges=None, name=None):
    '''
    Slice a TH2 vertically (ProjectionY) and return the list of slices
    '''
    log.fatal("This function has an index bug")

    if edges == None:
        edges = [hist.GetXaxis().GetBinLowEdge(iBin + 1) for iBin in range(hist.GetNbinsX() + 1)]

    slices = []
    lowEdges = edges[:-1]
    upEdges = edges[1:]
    if not name:
        name = hist.GetName()

    for lowEdge, upEdge in zip(lowEdges, upEdges):
        firstBin = hist.GetXaxis().FindBin(lowEdge * 1.0001)
        lastBin = hist.GetXaxis().FindBin(upEdge * 0.9999)

        slices.append(hist.ProjectionY(f'{name}{lowEdge:.0f}_{upEdge:.0f}', firstBin, lastBin))

        # Exclude underflow and overflow
        if lastBin <= 0 or firstBin >= hist.GetNbinsX():
            slices[-1].Reset()

    return slices


def CopyHistInSubrange(hist, xMin, xMax):
    '''
    Crop a 1-dimentional constant-binning histogram in a specified subrange. Useful to ensure that histograms with
    different number of bins but same bin widths can still be compared and operations can be performed on them, like
    divisions.

    Parameters
    ----------
    hist : TH1
        The histogram to be cropped
    xMin : float
        lower limit. The histogram is copied from xMin * 1.0001 to avoid problems related to compiler precision which
        could result in selecting the wrong bin.
    cMax : float
        upper limit. The histogram is copied from xMax * 0.9999 to avoid problems related to compiler precision which
        could result in selecting the wrong bin.

    Returns
    -------
    TH1
        The cropped histogram
    '''
    # Check that the binning is consistent
    N = (xMax - xMin) / hist.GetBinWidth(1)
    print(f'N: {N:.3f} xmax: {xMax:.3f}  xmin: {xMin:.3f}  width: {hist.GetBinWidth(1):.3f}')
    if abs(round(N) - N) > 1.e-6:
        log.critical('Range is not a multiple of bin width')
    N = round(N)

    hNew = TH1D(f'{hist.GetName()}_new', hist.GetTitle(), N, xMin, xMax)
    for iBinNew in range(N):
        iBinOld = hist.FindBin(hNew.GetBinCenter(iBinNew + 1))

        if iBinOld < 1 or iBinOld > hist.GetNbinsX():
            hNew.SetBinContent(iBinNew + 1, 0)
            hNew.SetBinError(iBinNew + 1, 0)
        else:
            hNew.SetBinContent(iBinNew + 1, hist.GetBinContent(iBinOld))
            hNew.SetBinError(iBinNew + 1, hist.GetBinError(iBinOld))

    hNew.SetDirectory(0)
    return hNew

def GetSpread(objects):
    '''
    Given a list of TH1, returns a TH1 where the points are the averages of the bins and the width is the spread
    (standard deviation) of the bin entries.

    Args:
        objects (list[TH1]): list of histograms

    Returns:
        TH1: histogram containing the average and spread of the input histograms
    '''
    if len(objects) == 0:
        log.critical('No objects provided')

    if not objects or len(objects) == 0:
        raise ValueError('Empty list of objects')

    if all(type(obj) in (TF1, TSpline3) for obj in objects):  # pylint: disable=unidiomatic-typecheck
        xMin = objects[0].GetXmin()
        xMax = objects[0].GetXmax()
        gSpread = TGraphErrors(1)

        for iPoint, x in enumerate(np.linspace(xMin, xMax, num=1000)):
            yValues = np.array([obj.Eval(x) for obj in objects])
            spread = np.std(yValues)
            avg = np.average(yValues)
            gSpread.SetPoint(iPoint, x, avg)
            gSpread.SetPointError(iPoint, 0, spread)

        return gSpread

    if all(type(obj) in (TH1, TH1D, TH1F, TH1I) for obj in objects):
        hSpread = objects[0].Clone('hSpread')
        hSpread.Reset()

        for iBin in range(hSpread.GetNbinsX()):
            yValues = np.array([obj.GetBinContent(iBin + 1) for obj in objects[1:]])
            spread = np.std(yValues)
            avg = np.average(yValues)
            hSpread.SetBinContent(iBin + 1, avg)
            hSpread.SetBinError(iBin + 1, spread)

        return hSpread

    raise NotImplementedError('Spread only implemented for histograms and TF1')

def ChangeUnits(obj, multiplier, name=None, title=''):
    '''
    Only for objogram with constant binwidth!
    '''

    if isinstance(obj, TH1):
        nbins = obj.GetNbinsX()
        lowEdge = obj.GetBinLowEdge(1)
        upEdge = obj.GetBinLowEdge(nbins+1)
        if name is None:
            name = f'{obj.GetName()}_new'
        hNew = TH1F(name, title, nbins, lowEdge*multiplier, upEdge*multiplier)
        for i in range(0, nbins+2):
            hNew.SetBinContent(i, obj.GetBinContent(i))
            hNew.SetBinError(i, obj.GetBinError(i))
        return hNew

    if isinstance(obj, TGraphErrors):
        nPoints = obj.GetN()
        gNew = TGraphErrors(nPoints)
        gNew.SetName(name if name else f'{obj.GetName()}_stretch')

        for iPoint in range(nPoints):
            x = obj.GetPointX(iPoint)
            y = obj.GetPointY(iPoint)
            gNew.SetPoint(iPoint, x * multiplier, y)

            xUnc = obj.GetErrorX(iPoint)
            yUnc = obj.GetErrorY(iPoint)
            gNew.SetPointError(iPoint, xUnc, yUnc)

        return gNew

    raise NotImplementedError()

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


def IsBinningCompatible(*args): # pylint: disable=inconsistent-return-statements
    '''
    Checks that the bin width is compatible between:
     * x- and y-axes for a TH2
     * x-axes for two TH1

    Parameters
    ----------
    args : tuple
        If only one histogram is given it is assumed that it is a TH2, and this function checks if the bin width of the
        x- and y-axes are the same. If two histograms are given, one is expected to be a TH1 and the other can be a TH1
        or TH2. In this case it will be checked that the bin width of the TH1 matches the one of the other histogram.

    Returns
    -------
    bool
        True if the bin width is the same
    '''

    if len(args) == 1:
        hMatrix = args[0]
        if isinstance(hMatrix, TH2):
            bwx = hMatrix.GetXaxis().GetBinWidth(1)
            bwy = hMatrix.GetYaxis().GetBinWidth(1)
            return abs(bwx - bwy) < 1.e-6 * bwx
        log.critical('Binning check not implemented for class %s', type(hMatrix))
    elif len(args) == 2:
        h1 = args[0]
        h2 = args[0]

        if isinstance(h1, TH1) and isinstance(h2, TH1):
            bw1 = h1.GetBinWidth(1)
            bw2 = h2.GetBinWidth(1)
            return abs(bw1 - bw2) < 1.e-6 * bw1
        if isinstance(h1, TH2) and isinstance(h2, TH1):
            bw1 = hMatrix.GetXaxis().GetBinWidth(1)
            bw2 = h2.GetBinWidth(1)
            return IsBinningCompatible(h1) and abs(bw1 - bw2) > 1.e-6 * bw1
        if isinstance(h1, TH2) and isinstance(h2, TH1):
            return IsBinningCompatible(h2, h1)
        log.critical('Not implemented for types %s and %s', type(h1), type(h2))
    log.critical('Not implemented for more than two histograms')


def WeightedAverage(inObj, weights):
    '''
    Compute the weighted average of a graph with an histogram (TH1)

    Parameters
    ----------
    inObj : TGraph, TH1, TF1
        The graph to be reweighted
    weights : TH1
        The weights to be applied

    Returns
    -------
    float
        The weighted average
    '''

    avg = 0
    counts = weights.Integral(1, weights.GetNbinsX())
    for iBin in range(weights.GetNbinsX()):
        freq = weights.GetBinContent(iBin + 1) / counts
        y=0
        if isinstance(inObj, (TGraph, TF1)):
            y = inObj.Eval(weights.GetBinCenter(iBin+1))
        elif isinstance(inObj, TH1):
            if inObj.GetNbinsX() != weights.GetNbinsX():
                log.critical('Incompatible binning: %s and %s have %s and %s bins respectively', \
                          inObj, weights, inObj.GetNbinsX(), weights.GetNbinsX())
            iBin = inObj.FindBin(weights.GetBinCenter(iBin+1))
            y = inObj.GetBinContent(iBin)
        else:
            log.critical('Not implemented')
        avg += freq * y
    return avg


def SmearGraph(graph, matrix, name=None, title=''): # pylint: disable=inconsistent-return-statements
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

    if isinstance(graph, TGraph):
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

    if isinstance(graph, TH1):
        if not IsBinningCompatible(graph, matrix):
            bx = graph.GetNbinsX()
            xmin = graph.GetXaxis().GetXmin()
            xmax = graph.GetXaxis().GetXmax()
            bw = graph.GetBinWidth(1)
            log.error(f'Binning of {graph}: nbins=%s xmin=%s xmax=%s bw=%s', bx, xmin, xmax, bw)
            bx = matrix.GetNbinsX()
            by = matrix.GetNbinsY()
            xmin = matrix.GetXaxis().GetXmin()
            xmax = matrix.GetXaxis().GetXmax()
            ymin = matrix.GetYaxis().GetXmin()
            ymax = matrix.GetYaxis().GetXmax()
            bwx = matrix.GetXaxis().GetBinWidth(1)
            bwy = matrix.GetYaxis().GetBinWidth(1)
            log.error('Binning of %s: nbins=(%s, %s) xrange=(%s, %s) yrange=(%s, %s), bw=(%s, %s)', \
                      matrix, bx, by, xmin, xmax, ymin, ymax, bwx, bwy)
            log.critical('Binning is not compatible')

        hSmeared = graph.Clone(f'{graph}')
        hSmeared.Reset()
        if name:
            hSmeared.SetName(name)
        for iBin in range(matrix.GetNbinsX()):
            hProj = matrix.ProjectionY(f'hProj_{iBin+1}', iBin+1, iBin+1)
            hProj.Scale(1. / hProj.GetEntries())
            hSmeared += hProj * graph.GetBinContent(iBin + 1)
        for iBin in range(hSmeared.GetNbinsX()):
            bc = graph.GetBinContent(iBin+1)
            be = graph.GetBinError(iBin+1)
            hSmeared.SetBinError(iBin + 1, hSmeared.GetBinContent(iBin + 1) * be / bc)

        return hSmeared
    log.critical("Smearing for type %s is not implemented", type(graph))

def Divide(num, den, name=None): #pylint: disable=inconsistent-return-statements
    '''
    Divide two quantities.
    Implemented types:

    +---------------+-----+-------------+
    | den \\ num    | TH1 | TGraphErrors |
    +=========+=======+=================+
    | TH1          |  ✘  | ✔            |
    +---------+-------+-----------------+
    | TGraphErrors |  ✘  | ✘            |
    +---------+-------+-----------------+

    Parameters
    ----------
    num : ROOT object
        numerator.
    den : ROOT object
        denominator.

    Returns
    -------
    TH1
        the division between num and den.
    '''

    if not name:
        name = f'{num.GetName()}_ratio'

    if isinstance(den, TH1):
        if isinstance(num, TGraphErrors):
            nBins = den.GetNbinsX()
            ratio = den.Clone(name)
            ratio.Reset()

            if nBins != num.GetN():
                log.warning("You you are trying to divide two objects with different number of bins/points.")

                for iBin in range(nBins):
                    x = den.GetBinCenter(iBin + 1)
                    y = num.Eval(x)
                    ey = 0

                    binContent = den.GetBinContent(iBin + 1)
                    binError = den.GetBinError(iBin + 1)

                    if binContent > 0 and y > 0:
                        r = y / binContent
                        ratioUnc = r * math.sqrt((ey / y) ** 2 + (binError / binContent) ** 2)

                        ratio.SetBinContent(iBin + 1, r)
                        ratio.SetBinError(iBin + 1, ratioUnc)
                    else:
                        ratio.SetBinContent(iBin + 1, 0)
                        ratio.SetBinError(iBin + 1, 0)

                return ratio

            if den.FindBin(num.GetPointX(1)) != 1 or den.FindBin(num.GetN()) != den.GetNBins():
                log.critical("The binnings are not aligned.")

            for iBin in range(nBins):
                x = num.GetPointX(iBin)
                y = num.GetPointY(iBin)
                ey = num.GetErrorY(iBin)

                bin_idx = den.FindBin(x)
                binContent = den.GetBinContent(bin_idx)
                binError = den.GetBinError(bin_idx)

                if binContent > 0 and y > 0:
                    ratio = y / binContent
                    ratioUnc = ratio * math.sqrt((ey / y) ** 2 + (binError / binContent) ** 2)

                    ratio.SetBinContent(iBin + 1, ratio)
                    ratio.SetBinError(iBin + 1, ratioUnc)
                else:
                    ratio.SetBinContent(iBin + 1, 0)
                    ratio.SetBinError(iBin + 1, 0)

            return ratio
    elif isinstance(den, TF1):
        if isinstance(num, TH1):
            ratio = num.Clone(name)
            ratio.Reset()

            for iBin in range(num.GetNbinsX()):
                ratio.SetBinContent(iBin + 1, num.GetBinContent(iBin + 1) / den.Eval(num.GetBinCenter(iBin + 1)))
                ratio.SetBinError(iBin + 1, num.GetBinError(iBin + 1) / den.Eval(num.GetBinCenter(iBin + 1)))
            
            return ratio

    elif isinstance(den, TH2):
        if isinstance(num, TH2):
            ratio = num.Clone(name)
            ratio.Reset()
            ratio.GetZaxis().SetTitle('Ratio')
            for iBinX in range(num.GetNbinsX()):
                for iBinY in range(num.GetNbinsY()):
                    x = num.GetXaxis().GetBinCenter(iBinX + 1)
                    y = num.GetYaxis().GetBinCenter(iBinY + 1)
                    bc = num.GetBinContent(iBinX + 1, iBinY + 1)
                    den = den.Eval(x, y)
                    if den.Eval(x, y) > 0:
                        ratio = bc / den
                    else:
                        ratio = 0
                        
                    ratio.SetBinContent(iBinX + 1, iBinY + 1,  ratio)
            return ratio

    log.error("Division of %s by %s is not implemented aaaa.", num.ClassName(), den.ClassName())
    return None

def Bootstrap(obj):
    '''
    Returns a bootstrapped version of the input objects.
    Supperted types: TH1.

    Parameters
    ----------
    graph : TH1
        histogram to be bootstrapped

    returns:
        The bootstrapped histogram
    '''

    if type(obj) in (TH1, TH1D, TH1F, TH1I):
        hBS = obj.Clone(obj.GetName() + '_bs')
        for iBin in range(hBS.GetNbinsX()):
            bc = obj.GetBinContent(iBin + 1)
            bu = obj.GetBinError(iBin + 1)

            buNew = np.random.normal(loc=bc, scale=bu)

            hBS.SetBinContent(iBin + 1, buNew)

        return hBS

    print(type(obj))
    raise NotImplementedError('Bootstrap implemented only for TH1')
