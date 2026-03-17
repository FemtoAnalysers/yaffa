from ROOT import TFile
import sys
sys.path.append('../../src')

from utils import SliceVertically

def main(cfg:dict):
    mTMins = cfg['mt_bins'][:-1]
    mTMaxs = cfg['mt_bins'][1:]

    # Open input file
    inFile = TFile(cfg['infile'])

    # Load hitograms from file
    hRhoVsMt = inFile.Get('hRhoVsMt')

    hMt = hRhoVsMt.ProjectionX('hMt')
    hRhos = SliceVertically(hRhoVsMt, cfg['mt_bins'], name='hRho')

    oFile = TFile(cfg['ofile'], 'recreate')
    for hRho in hRhos:
        hRho.Write()
    hMt.Write()

if __name__ == '__main__':
    import argparse
    import yaml
    
    parser = argparse.ArgumentParser()
    parser.add_argument('cfg')
    args = parser.parse_args()

    with open(args.cfg) as file:
        cfg = yaml.safe_load(file)
    
    main(cfg)
