'''
Script to compute the size of the source for mT-integrated and mT differential analyses
'''

import argparse

from ROOT import TFile, gROOT

from yaffa import utils

def PrintAvgMtInfo(hMt, mTBinLimits):
    '''Print average mT for integrated mT and for specific mT bins if available'''
    mTMins = mTBinLimits[:-1]
    mTMaxs = mTBinLimits[1:]

    print(f'\tmT (Integrated): {hMt.GetMean():.3f} GeV')
    for mTMin, mTMax in zip(mTMins, mTMaxs):
        hMt.GetXaxis().SetRangeUser(mTMin, mTMax)
        print(f'\tmT in [{mTMin:.3f}, {mTMax:.3f}]: {hMt.GetMean():.3f} GeV')


def PrintMtInfo(inFileName, combs, mTBinLimits):
    '''Compute the size of the source'''

    gROOT.SetBatch(True)

    suffix = 3001
    inFile = TFile(inFileName)

    for p1, p2 in combs:
        print(f'\nAverage mT for pair {p1}{p2} for k* < 200 MeV/c:')
        histName = f'HMResults{suffix}/HMResults{suffix}/Particle{p1}_Particle{p2}/MEmTDist_Particle{p1}_Particle{p2}'
        hMtVsKstar = utils.io.Load(inFile, histName)
        hMt = hMtVsKstar.ProjectionY("hMt", 1, hMtVsKstar.GetXaxis().FindBin(0.2 * 0.9999))
        PrintAvgMtInfo(hMt, mTBinLimits)


def main():
    '''main'''
    parser = argparse.ArgumentParser()
    parser.add_argument('inFile')
    parser.add_argument("--combs", type=float, nargs="+")
    parser.add_argument("--mT", type=float, nargs="+", default=[])
    args = parser.parse_args()

    combs = ['02', '13', '03', '12']

    PrintMtInfo(args.inFile, combs, args.mT)

if __name__ == '__main__':
    main()
