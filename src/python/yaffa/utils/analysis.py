'''
Module that contains various functions useful for the post-processing of the data.
'''
import math

import numpy as np

from ROOT import TH1, TH1D, TH1F, TH1I, TH2D, TSpline3, TGraph, TH2, TGraphErrors, TGraphAsymmErrors, TGraphMultiErrors, TF1, TF2  # pylint: disable=import-error

from yaffa import logger as log

def EnforceMeV(obj):
    '''
    Detect if the axes of the object extend beyond 10. If so nothing happens, otherwise the object is assumed to be in 
    GeV and the axes is rescaled by a factor 1000 automatically.
    '''
    if isinstance(obj, TH1):
        xmax = obj.GetXaxis().GetXmax()
        if xmax < 10:
            log.warning('The object %s seems to be in GeV but MeV are required. Changing automatically to MeV. Provide ' \
            'a histogram in MeV to suppress this warning', obj.GetName())
        obj = ChangeUnits(obj, 1000, f'{obj.GetName()}_MeV')

        return obj

    log.critical('Unit check not implemented for objects of type %s', type(obj))


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
        if lastBin < 1 or firstBin > hist.GetNbinsX():
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
    Change the units of the x axis of a 1D histogram (constant bin width only) or of a graph. For graphs, the x errors
    are scaled too.

    Parameters
    ----------
    obj : TH1, TGraph, TGraphErrors, TGraphAsymmErrors, TGraphMultiErrors
        The object to be converted
    multiplier : float
        The multiplier to be applied to the x axis
    name : str, optional
        The name of the converted object. By default: `<old_name>_new` for histograms, `<old_name>_stretch` for graphs
    title : str, optional
        The title of the converted histogram. By default: `''`. Graphs keep their title

    Returns
    -------
    TH1F or graph of the same type as `obj`
        The object in the new units
    '''

    if isinstance(obj, TH1) and not isinstance(obj, TH2):
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

    if isinstance(obj, TGraph):
        gNew = obj.Clone(name if name else f'{obj.GetName()}_stretch')

        for iPoint in range(gNew.GetN()):
            gNew.SetPointX(iPoint, gNew.GetPointX(iPoint) * multiplier)
            if isinstance(gNew, TGraphErrors):
                gNew.SetPointError(iPoint, gNew.GetErrorX(iPoint) * multiplier, gNew.GetErrorY(iPoint))
            elif isinstance(gNew, TGraphAsymmErrors):
                gNew.SetPointError(iPoint, gNew.GetErrorXlow(iPoint) * multiplier, gNew.GetErrorXhigh(iPoint) * multiplier,
                                   gNew.GetErrorYlow(iPoint), gNew.GetErrorYhigh(iPoint))
            elif isinstance(gNew, TGraphMultiErrors):
                gNew.SetPointEX(iPoint, gNew.GetErrorXlow(iPoint) * multiplier, gNew.GetErrorXhigh(iPoint) * multiplier)

        return gNew

    raise NotImplementedError(f'Changing units is not implemented for {type(obj)}')

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
    lowEdgeY = hist.GetYaxis().GetBinLowEdge(1)
    upEdgeX = hist.GetXaxis().GetBinLowEdge(nbinsX+1)
    upEdgeY = hist.GetYaxis().GetBinLowEdge(nbinsY+1)

    if name is None:
        name = f'{hist.GetName()}_new'

    if isinstance(multiplier, (int, float)):
        multX, multY = (multiplier, multiplier)
    else:
        multX, multY = multiplier

    hNew = TH2D(name, title, nbinsX, lowEdgeX * multX, upEdgeX * multX, nbinsY, lowEdgeY * multY, upEdgeY * multY)

    for iBinX in range(0, nbinsX+2):
        for iBinY in range(0, nbinsY+2):
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
        for iBin in range(matrix.GetNbinsY()):
            hProj = matrix.ProjectionX(f'hProj_{iBin+1}', iBin+1, iBin+1)
            counts = hProj.Integral(1, hProj.GetNbinsX())

            if counts < 1:
                continue

            x = matrix.GetYaxis().GetBinCenter(iBin+1)
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

    +--------+-----+--------+-----+-----+-----+
    | num \\ den | TH1 | TGraph | TF1 | TH2 | TF2 |
    +========+=====+========+=====+=====+=====+
    | TH1    |  ✔  |   ✔    |  ✔  |  ✘  |  ✘  |
    +--------+-----+--------+-----+-----+-----+
    | TGraph |  ✔  |   ✔    |  ✔  |  ✘  |  ✘  |
    +--------+-----+--------+-----+-----+-----+
    | TH2    |  ✘  |   ✘    |  ✘  |  ✔  |  ✔  |
    +--------+-----+--------+-----+-----+-----+

    The denominator is evaluated at the x of the numerator and its uncertainty is propagated to
    the ratio. A TF1 denominator is exact, and so is a graph denominator that must be interpolated
    because x does not coincide with any of its points.

    Parameters
    ----------
    num : ROOT object
        numerator.
    den : ROOT object
        denominator.

    Returns
    -------
    TH1 or TGraphErrors
        the division between num and den. Graph numerators give a TGraphErrors, histograms a TH1.
    '''

    if not name:
        name = f'{num.GetName()}_ratio'

    def Evaluate(obj, x):
        '''
        Value and uncertainty of obj at x. Graphs carry an uncertainty only at their own points:
        when x falls between two of them the value is interpolated and taken as exact, as for TF1.
        '''
        if isinstance(obj, TH1):
            iBin = obj.FindBin(x)
            return obj.GetBinContent(iBin), obj.GetBinError(iBin)

        if isinstance(obj, TGraph):
            iPoint = next((i for i in range(obj.GetN()) if math.isclose(obj.GetPointX(i), x)), None)
            if iPoint is not None:
                return obj.GetPointY(iPoint), max(obj.GetErrorY(iPoint), 0)  # -1 without uncertainties

        return obj.Eval(x), 0

    if isinstance(num, TGraph) and isinstance(den, (TH1, TGraph, TF1)) and not isinstance(den, (TH2, TF2)):
        gRatio = TGraphErrors(1)
        gRatio.SetName(name)

        iRatioPoint = 0
        for iPoint in range(num.GetN()):
            x = num.GetPointX(iPoint)
            y = num.GetPointY(iPoint)
            yUnc = max(num.GetErrorY(iPoint), 0)  # graphs without uncertainties return -1
            d, dUnc = Evaluate(den, x)

            if d == 0:
                continue

            gRatio.SetPoint(iRatioPoint, x, y / d)
            gRatio.SetPointError(iRatioPoint, max(num.GetErrorX(iPoint), 0), math.hypot(yUnc, y * dUnc / d) / d)
            iRatioPoint += 1

        return gRatio

    if isinstance(num, TH1) and not isinstance(num, TH2) and isinstance(den, TGraph):
        hRatio = num.Clone(name)
        hRatio.Reset()

        for iBin in range(num.GetNbinsX()):
            y = num.GetBinContent(iBin + 1)
            d, dUnc = Evaluate(den, num.GetBinCenter(iBin + 1))
            if d == 0:
                continue

            hRatio.SetBinContent(iBin + 1, y / d)
            hRatio.SetBinError(iBin + 1, math.hypot(num.GetBinError(iBin + 1), y * dUnc / d) / d)

        return hRatio

    if isinstance(num, TH1) and isinstance(den, TH1) and not isinstance(num, TH2) and not isinstance(den, TH2):
        hRatio = num.Clone(name)
        if not hRatio.Divide(den):
            log.error('The division of %s by %s failed. Check that the binnings match', num.GetName(), den.GetName())
            return None

        return hRatio

    if isinstance(den, TF2):
        if isinstance(num, TH2):
            hRatio = num.Clone(name)
            hRatio.Reset()
            hRatio.GetZaxis().SetTitle('Ratio')

            for iBinX in range(num.GetNbinsX()):
                for iBinY in range(num.GetNbinsY()):
                    x = num.GetXaxis().GetBinCenter(iBinX + 1)
                    y = num.GetYaxis().GetBinCenter(iBinY + 1)
                    d = den.Eval(x, y)
                    if d == 0:
                        continue

                    hRatio.SetBinContent(iBinX + 1, iBinY + 1, num.GetBinContent(iBinX + 1, iBinY + 1) / d)

            return hRatio

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
            hRatio = num.Clone(name)
            hRatio.Reset()
            hRatio.GetZaxis().SetTitle('Ratio')
            for iBinX in range(num.GetNbinsX()):
                for iBinY in range(num.GetNbinsY()):
                    bcDen = den.GetBinContent(iBinX + 1, iBinY + 1)
                    if bcDen == 0:
                        continue

                    hRatio.SetBinContent(iBinX + 1, iBinY + 1, num.GetBinContent(iBinX + 1, iBinY + 1) / bcDen)
            return hRatio

    log.error("Division of %s by %s is not implemented.", num.ClassName(), den.ClassName())
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
