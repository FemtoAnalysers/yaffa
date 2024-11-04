'''
Script to smear a graph
'''

import argparse
import yaml

from ROOT import TFile, gROOT, TH1
from yaffa.utils.io import Load
from yaffa.utils.analysis import SmearGraph
from yaffa import logger as log


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('cfgs', default='cfg_example_smear_graph.yml')
    parser.add_argument('-b', default=False, action='store_true')
    args = parser.parse_args()

    gROOT.SetBatch(args.b)

    # Load configuration file
    with open(args.cfgs, "r") as stream:
        try:
            cfgs = yaml.safe_load(stream)
        except yaml.YAMLError:
            log.critical('Yaml file not loaded')

    for cfg in cfgs:
        # Load graph to be smeared
        inFile = TFile(cfg['graph']['file'])
        inObj = Load(inFile, cfg['graph']['name'])
        inObj.SetDirectory(0)
        inFile.Close()

        if not isinstance(inObj, TH1):
            log.critical('Smearing for type %s is not implemented.', type(inObj))

        # Load smearing matrix
        inFile = TFile(cfg['smearing']['file'])
        hSmearingMatrix = Load(inFile, cfg['smearing']['name'])
        hSmearingMatrix.SetDirectory(0)
        inFile.Close()

        if not isinstance(inObj, TH1):
            log.critical('Smearing for type %s is not implemented.', type(inObj))

        # Perform the smearing
        oFile = TFile(cfg['output']['file'], 'recreate')
        hSmeared = SmearGraph(inObj, hSmearingMatrix, name=cfg['output']['name'])
        hSmeared.Write()
        inObj.Write()

        oFile.Close()
        print(f'Output saved in {cfg["output"]["file"]}')
