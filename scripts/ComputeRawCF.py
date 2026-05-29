'''
Script to compute the raw 2B and 3B correlation functions.
'''

import argparse
from rich import print

from ROOT import TFile, TH1, TDirectory

from yaffa import logger as log

EVENT_TYPES = ['se', 'me']

# Which sub-systems should be analyzed for each analysis
SYSTEM_QUEUE = {
    'ppp' : ['ppp', 'apapap']
}

# How pairs/triplets should be combined with anti-pairs and anti-triplets
SUM_RECIPE = {
    'ppp' : {
        'ppp_apapap' : ['ppp', 'apapap']
    }
}

def SumHistograms(hists:list[TH1]) -> TH1:
    """Sum histograms

    Args:
        hists (list[TH1]): list of histograms to be summed

    Returns:
        TH1: the sum of the histograms
    """
    hist = hists[0].Clone()

    for h in hists[1:]:
        hist.Add(h)
    
    return hist

def SumDistr(distr:dict, system:str):
    """Sum the histograms according to the sum recipe

    Args:
        distr (dict): the input same- and mixed-event distributions
        system (str): the system to take into account e.g. ppp, ppL, etc
    """
    recipes = SUM_RECIPE.get(system)
    
    for key, recipe in recipes.items():
        distr[key] = {}

        for event in EVENT_TYPES:
            distr[key][event] = {}

            for hist in distr[recipe[0]][event]:
                to_sum = [dist[event][hist] for name, dist in distr.items() if name in recipe]
                distr[key][event][hist] = SumHistograms(to_sum)

def ComputeNormalization(hSE:TH1, hME:TH1, norm_region:list[float]=None) -> float:
    """Compute the normalization for the correlation function

    Args:
        hSE (TH1): same-event distribution
        hME (TH1): mixed-event distribution
        norm_region (list[float], optional): normalization region. If none, the normalization is computed via yields,
        i.e. the normalization region is [0, infinity]. Defaults to None.

    Returns:
        float: the normalization
    """    
    if norm_region is None:
        return hME.Integral() / hSE.Integral()
    
    firstBin = hSE.GetXaxis().FindBin(norm_region[0] * 1.0001)
    lastBin = hSE.GetXaxis().FindBin(norm_region[1] * 0.9999)

    return hME.Integral(firstBin, lastBin) / hSE.Integral(firstBin, lastBin)
    
def ProcessTriplets(distr:dict, dir:TDirectory, norm:list[float]=None) -> None:
    """Perform the full analysis of a 3B system

    Args:
        distr (dict): the same- and mixed-event distributions
        dir (TDirectory): the output file or directory where to save the results
        norm (list[float], optional): normalization region. If None, the CF is normalized via yields. Defaults to None.
    """    
    for system in distr:
        subdir = dir.mkdir(system)    
        subdir.cd()

        for event in EVENT_TYPES:
            subsubdir = dir.mkdir(f'{system}/{event}')
            subsubdir.cd()
            
            hQ3VsMult = distr[system][event]['q3vsmult']
            hQ3VsMult.Write()
 
            # mT and multiplicity integrated distributions
            distr[system][event]['q3'] = hQ3VsMult.ProjectionY('hQ3')
            distr[system][event]['q3'].Write()

        subsubdir = dir.mkdir(f'{system}/cf')
        subsubdir.cd()

        hSE = distr[system]['se']['q3']
        hME = distr[system]['me']['q3']
        
        hCF = hSE.Clone('hCF')
        hCF.Sumw2()
        hCF.Divide(hME)
        hCF.Scale(ComputeNormalization(hSE, hME, norm))
        hCF.Write()

    dir.cd()

def GetDistributions(file:str, system:str, bw:float=None) -> dict:
    """Load distributions from AnalysisResults.root

    Args:
        file (str): input file (AnalysisResults.root)
        system (str): the system to be studied
        bw (float, optional): bin width of the correlation function. The units must be the same as the SE and ME
         distributions in the input file. If None, the original binning is preserved. Defaults to None.

    Returns:
        dict: the distributions to be used for femto
    """    
    inFile = TFile.Open(file)

    distr = {}
    if queue := SYSTEM_QUEUE.get(system):
        for system in queue:
            distr[system] = {}

            for event in EVENT_TYPES:
                distr[system][event] = {}

                hTriplets = inFile.Get(f'femto-triplet-track-track-track_{system}/TrackTrackTrack/{event.upper()}/Analysis/hQ3VsMtVsMultVsCent')

                distr[system][event]['q3vsmult'] = hTriplets.Projection(0, 2)
                distr[system][event]['q3vsmult'].SetName('hQ3VsMult')
                distr[system][event]['q3vsmult'].SetDirectory(0)
                if bw is not None:
                    rebin = round(bw / distr[system][event]['q3vsmult'].GetYaxis().GetBinWidth(1))

                    if rebin < 1:
                        log.critical('Rebin value too small. Check your units')
                    if rebin >=  distr[system][event]['q3vsmult'].GetNbinsX():
                        log.critical('Very large rebin value. Check your units')

                    distr[system][event]['q3vsmult'].RebinY(rebin)

    else:
        log.error(f'System {system} is not implemented')

    return distr

def main(inFile:str, oFileName:str, system:str, norm:list[float]=None, bw:float=None) -> None:
    distr = GetDistributions(inFile, system, bw)

    SumDistr(distr, system)

    oFile = TFile(oFileName, 'recreate')
    oFile.cd()
    ProcessTriplets(distr, oFile, norm)
    oFile.Close()

    print(f'Output saved in {oFileName}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', default='AnalysisResults.root')
    parser.add_argument('ofile', nargs='?', default='RawCF.root')
    parser.add_argument('--system', choices=SYSTEM_QUEUE.keys())
    parser.add_argument('--norm', nargs=2, type=float, default=None)
    parser.add_argument('--bw', nargs='?', type=float, default=None, help='k* bin width of the CF')
    args = parser.parse_args()

    main(args.infile, args.ofile, args.system, args.norm, args.bw)
