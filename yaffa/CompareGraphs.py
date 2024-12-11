'''
Compare histograms, graphs and functions.
'''

import argparse
import os
import yaml
from rich import print  # pylint: disable=redefined-builtin
import numpy as np
import numexpr

# pylint: disable=no-name-in-module
from ROOT import TFile, TCanvas, TLegend, TLine, TH1, TGraph, TGraphErrors, TGraphAsymmErrors, TH1D, TF1, gROOT

from yaffa import logger as log
from yaffa import utils

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


parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('cfg')
parser.add_argument('-b', action='store_true', default=False, help='Run in batch mode')
args = parser.parse_args()

gROOT.SetBatch(args.b)

# Load configuration file
with open(args.cfg, "r") as stream:
    try:
        cfg = yaml.safe_load(stream)
    except yaml.YAMLError:
        log.critical('Yaml file not loaded')

utils.style.SetStyle()

for plot in cfg:
    plot = plot["plot"]

    if 'root' in plot["opt"]["ext"]:
        oFile = TFile(f'{os.path.splitext(plot["output"])[0]}.root', 'recreate')

    panels = {'default': 1}
    if plot['ratio']['enable']:
        panels['ratio'] = len(panels)+1
    if plot['spread']['enable']:
        panels['spread'] = len(panels)+1
    if plot['relunc']['enable']:
        panels['relunc'] = len(panels)+1
    if plot.get('pulls').get('enable'):
        panels['pulls'] = len(panels)+1

    # Frame Coordinates
    fx1 = plot['opt']['rangex'][0]
    fy1 = plot['opt']['rangey'][0]
    fx2 = plot['opt']['rangex'][1]
    fy2 = plot['opt']['rangey'][1]

    # Load the objects to draw
    inObjs = []
    legends = []
    drawOpts = []
    for inputCfg in plot["input"]:
        inFile = TFile(inputCfg['file'])

        inObj = utils.io.Load(inFile, inputCfg['name'])

        if isinstance(inObj, TH1):
            inObj.SetDirectory(0)
            inObj.Rebin(inputCfg['rebin'])
            inObj = CopyHistInSubrange(inObj, fx1, fx2)

            if inputCfg['normalize']:
                inObj.Scale(1./inObj.Integral())
            if inputCfg['scale']:
                inObj.Scale(numexpr.evaluate(str(inputCfg['scale'])))
            if inputCfg['normalizecf']:
                inObj.Scale(inputCfg['normalizecf'])

        inObj.SetLineColor(utils.style.GetColor(inputCfg['color']))
        inObj.SetFillColorAlpha(utils.style.GetColor(inputCfg['color']), inputCfg['fillalpha'])
        inObj.SetMarkerColor(utils.style.GetColor(inputCfg['color']))
        inObj.SetLineWidth(inputCfg.get('thickness', 1))
        if isinstance(inObj, TH1):
            default = 'pe'
        elif isinstance(inObj, TGraph):
            default = 'p'
        elif isinstance(inObj, (TGraphErrors, TGraphAsymmErrors)):
            default = 'pe'
        elif isinstance(inObj, TF1):
            default = 'l'
        else:
            default = ''
        drawOpts.append(inputCfg.get('drawopt', default))
        inObj.SetMarkerStyle(inputCfg['markerstyle'])
        inObjs.append(inObj)
        legends.append(inputCfg['legend'])

    # Define the canvas
    nPanelsX, nPanelsY = utils.style.GetNPanels(len(panels))
    cPlot = TCanvas("cPlot", "cPlot", 600*nPanelsX, 600*nPanelsY)
    cPlot.Divide(nPanelsX, nPanelsY)
    pad = cPlot.cd(1)
    pad.SetLogx(plot["opt"]["logx"])
    pad.SetLogy(plot["opt"]["logy"])
    pad.DrawFrame(fx1, fy1, fx2, fy2, utils.style.SmartLabel(plot['opt']['title']))

    legx1 = plot['opt']['leg']['posx'][0]
    legy1 = plot['opt']['leg']['posy'][0]
    legx2 = plot['opt']['leg']['posx'][1]
    legy2 = plot['opt']['leg']['posy'][1]
    leg = TLegend(legx1, legy1, legx2, legy2)

    if 'root' in plot["opt"]["ext"]:
        oFile.cd()

    for iObj, (inObj, legend, drawOpt) in enumerate(zip(inObjs, legends, drawOpts)):
        if 'root' in plot["opt"]["ext"]:
            inObj.Write()

        inObj.Draw('same' + drawOpt)

        # Compute statistics for hist in the displayed range
        if isinstance(inObj, TH1):
            firstBin = inObj.FindBin(plot['opt']['rangex'][0]*1.0001)
            lastBin = inObj.FindBin(plot['opt']['rangex'][1]*0.9999)
            inObj.GetXaxis().SetRange(firstBin, lastBin)
            print(f'{legend}: mean = {inObj.GetMean()} sigma = {inObj.GetStdDev()}')
            if plot['opt']['leg']['yield']:
                legend += f';  Y={inObj.GetEntries():.1e}'
            if plot['opt']['leg']['relyield']:
                legend += f';  ({inObj.GetEntries() / inObjs[0].GetEntries() * 100:.1f}%)'
            if plot['opt']['leg']['mean']:
                legend += f';  #mu={inObj.GetMean():.3f}'
            if plot['opt']['leg']['sigma']:
                legend += f';  #sigma={inObj.GetStdDev():.3f}'
        if legend:
            leg.AddEntry(inObj, legend, 'lp')

    for line in plot['opt']['lines']:
        x1 = plot['opt']['rangex'][0] if(line['coordinates'][0] == 'min') else line['coordinates'][0]
        y1 = plot['opt']['rangey'][0] if(line['coordinates'][1] == 'min') else line['coordinates'][1]
        x2 = plot['opt']['rangex'][1] if(line['coordinates'][2] == 'max') else line['coordinates'][2]
        y2 = plot['opt']['rangey'][1] if(line['coordinates'][3] == 'max') else line['coordinates'][3]
        inputline = TLine(x1, y1, x2, y2)
        inputline.SetLineColor(utils.style.GetColor(line['color']))
        inputline.SetLineWidth(line['thickness'])
        inputline.Draw("same")
        leg.AddEntry(inputline, utils.style.SmartLabel(line['legendtag']),"l")

    leg.SetHeader(utils.style.SmartLabel(plot['opt']['leg']['header']), 'C')
    leg.Draw()

    # Compute ratio wrt the first obj
    while plot['ratio']['enable']: # use if-equivallent while scope to be able to control when to exit
        if len(inObjs) < 2:
            log.error("Not enough objects for making a ratio. Skipping ratio plot")
            break

        pad = cPlot.cd(panels['ratio'])
        pad.SetLogx(plot['ratio']['logx'])
        pad.SetLogy(plot['ratio']['logy'])
        y1 = plot['ratio']['rangey'][0]
        y2 = plot['ratio']['rangey'][1]
        frame = pad.DrawFrame(fx1, y1, fx2, y2, utils.style.SmartLabel(plot['opt']['title']))
        frame.GetYaxis().SetTitle('Ratio')
        hDen = inObjs[0].Clone()
        if isinstance(hDen, TH1):
            hDen.Rebin(plot['ratio']['rebin'])
            hDen.Sumw2()

        if isinstance(inObj, TH1):
            for inObj in inObjs[1:]:
                hRatio = inObj.Clone(f'{inObj.GetName()}_ratio')
                hRatio.Rebin(plot['ratio']['rebin'])
                hRatio.Divide(hDen)
                hRatio.Draw('same pe')
        elif isinstance(inObj, TGraphErrors):
            for inObj in inObjs[1:]:
                hRatio = utils.analysis.Divide(inObj, hDen)
                hRatio.SetName(f'{inObj.GetName()}_ratio')
                hRatio.Draw('same pe')
        else:
            log.error('Ratio for type %s is not implemented. Skipping this object', type(inObj))
            continue

        if 'root' in plot["opt"]["ext"]:
            hRatio.Write()

        line = TLine(plot['opt']['rangex'][0], 1, plot['opt']['rangex'][1], 1)
        line.SetLineColor(13)
        line.SetLineStyle(7)
        line.Draw('same pe')

        break

    if plot['spread']['enable']:
        # Calcualate the spread around the first oject
        if not isinstance(inObj, TH1):
            log.error('Spread for type %s is not implemented. Skipping this object', type(inObj))
            continue

        y1 = plot['spread']['rangey'][0]
        y2 = plot['spread']['rangey'][1]

        pad = cPlot.cd(panels['spread'])
        pad.SetLogx(plot['spread']['logx'])
        pad.SetLogy(plot['spread']['logy'])
        frame = pad.DrawFrame(fx1, y1, fx2, y2, utils.style.SmartLabel(plot['opt']['title']))
        title = 'relative spread #sigma/#mu'
        if plot['spread']['mode'] == 'percentage':
            title += ' (%)'
        frame.GetYaxis().SetTitle(title)

        hSpread = inObjs[0].Clone('hSpread')
        hSpread.Reset()

        for iBin in range(hSpread.GetXaxis().FindBin(x1 * 1.0001), hSpread.GetXaxis().FindBin(x2 * 0.9999)):
            yValues = np.array([obj.GetBinContent(iBin + 1) for obj in inObjs[1:]])
            spread = np.std(yValues) / inObjs[0].GetBinContent(iBin + 1)
            if plot['spread']['mode'] == 'relative':
                pass
            elif plot['spread']['mode'] == 'percentage':
                spread *= 100
            else:
                log.critical('Not implemented')
            hSpread.SetBinContent(iBin + 1, spread)

        hSpread.Draw('same')

    # Compute the relative uncertainties
    if plot['relunc']['enable']:
        pad = cPlot.cd(panels['relunc'])
        pad.SetGridx(plot['relunc']['gridx'])
        pad.SetGridy(plot['relunc']['gridy'])
        pad.SetLogx(plot["relunc"]["logx"])
        pad.SetLogy(plot["relunc"]["logy"])
        y1 = plot['relunc']['rangey'][0]
        y2 = plot['relunc']['rangey'][1]
        pad.SetLeftMargin(0.16)
        pad.SetRightMargin(0.1)
        frame = pad.DrawFrame(fx1, y1, fx2, y2, utils.style.SmartLabel(plot['opt']['title']))
        frame.GetYaxis().SetTitle('Relative uncertainty (%)')

        for iObj, inObj in enumerate(inObjs):
            if isinstance(inObj, TH1):
                nBins = inObj.GetNbinsX()
                hRelUnc = TH1D(f'hRelUnc_{iObj}', '', nBins, inObj.GetXaxis().GetXmin(), inObj.GetXaxis().GetXmax())

                for iBin in range(hRelUnc.GetNbinsX() + 1):
                    if inObj.GetBinContent(iBin) > 0:
                        hRelUnc.SetBinContent(iBin, 100 * inObj.GetBinError(iBin)/inObj.GetBinContent(iBin))
                        hRelUnc.SetBinError(iBin, 0)
            elif isinstance(inObj, TGraphErrors):
                hRelUnc = TGraphErrors(1)
                for iPoint in range(inObj.GetN()):
                    x = inObj.GetPointX(iPoint)
                    y = inObj.GetPointY(iPoint)
                    xUnc = inObj.GetErrorX(iPoint)
                    yUnc = inObj.GetErrorY(iPoint)

                    hRelUnc.SetPoint(iPoint, x,  100 * yUnc/y)
                    hRelUnc.SetPointError(iPoint, xUnc, 0)
            elif isinstance(inObj, TGraphAsymmErrors):
                hRelUnc = TGraphAsymmErrors(1)
                for iPoint in range(inObj.GetN()):
                    x = inObj.GetPointX(iPoint)
                    y = inObj.GetPointY(iPoint)
                    yUnc = inObj.GetErrorY(iPoint)  # Computes the average of the uncertainties
                    xUncUpper = inObj.GetErrorXhigh(iPoint)
                    xUncLower = inObj.GetErrorXlow(iPoint)

                    hRelUnc.SetPoint(iPoint, x,  100 * yUnc/y)
                    hRelUnc.SetPointError(iPoint, xUncLower, xUncUpper, 0, 0)
            else:
                log.error('Relative uncertainties for type %s are not implemented. Skipping this object', type(inObj))
                continue

            hRelUnc.SetLineColor(inObj.GetLineColor())
            hRelUnc.SetMarkerColor(inObj.GetMarkerColor())
            if isinstance(inObj, TH1):
                hRelUnc.DrawCopy('same hist')
            elif isinstance(inObj, TGraph):
                hRelUnc.Draw('same p')

            if 'root' in plot["opt"]["ext"]:
                hRelUnc.Write()

    # Compute the relative uncertainties
    if plot.get('pulls').get('enable'):
        pad = cPlot.cd(panels['pulls'])
        pad.SetGridx(plot['pulls']['gridx'])
        pad.SetGridy(plot['pulls']['gridy'])
        pad.SetLogx(plot["pulls"]["logx"])
        pad.SetLogy(plot["pulls"]["logy"])
        y1 = plot['pulls']['rangey'][0]
        y2 = plot['pulls']['rangey'][1]
        pad.SetLeftMargin(0.16)
        frame = pad.DrawFrame(fx1, y1, fx2, y2, utils.style.SmartLabel(plot['opt']['title']))
        frame.GetYaxis().SetTitle('Pulls')
        refObj = inObjs[0].Clone()
        if isinstance(inObjs[0], TF1):
            for iObj, inObj in enumerate(inObjs[1:]):
                if isinstance(inObj, TH1):
                    nBins = inObj.GetNbinsX()
                    hPulls = TH1D(f'hPulls_{iObj}', '', nBins, inObj.GetXaxis().GetXmin(), inObj.GetXaxis().GetXmax())

                    for iBin in range(hPulls.GetNbinsX() + 1):
                        if inObj.GetBinContent(iBin + 1) > 0:
                            delta = inObj.GetBinContent(iBin + 1) - refObj.Eval(inObj.GetBinCenter(iBin + 1))
                            pull = delta / inObj.GetBinError(iBin + 1)
                            hPulls.SetBinContent(iBin + 1, pull)
                            hPulls.SetBinError(iBin + 1, 0)
                elif isinstance(inObj, (TGraphAsymmErrors, TGraphErrors)) :
                    hPulls = TGraphErrors(1)
                    for iPoint in range(inObj.GetN()):
                        x = inObj.GetPointX(iPoint)
                        y = inObj.GetPointY(iPoint)
                        yUnc = inObj.GetErrorY(iPoint)

                        pull =  (y - refObj.Eval(x)) / yUnc

                        hPulls.SetPoint(iPoint, x,  pull)
                        hPulls.SetPointError(iPoint, 0, 0)
                else:
                    log.error('Pulls for type %s are not implemented. Skipping this object', type(inObj))
                    continue

                hPulls.SetLineColor(inObj.GetLineColor())
                hPulls.SetMarkerColor(inObj.GetMarkerColor())
                if isinstance(inObj, TH1):
                    hPulls.DrawCopy('same l')
                elif isinstance(inObj, TGraph):
                    hPulls.Draw('same l')
        else:
            log.error('Pulls not implemented for a reference of type %(type(refObj))')

    cPlot.Modified()
    cPlot.Update()

    # save canvas
    for ext in plot["opt"]["ext"]:
        if ext == 'root':
            oFile.Close()
            continue
        cPlot.SaveAs(f'{os.path.splitext(plot["output"])[0]}.{ext}')
