'''
Module that contains various functions useful for the post-processing of the data.
'''

from ROOT import TH1F  # pylint: disable=import-error

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
