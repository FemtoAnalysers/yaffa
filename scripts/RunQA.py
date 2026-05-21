'''
Script to produce the QA plots.
'''

import os
import argparse

from ROOT import TFile, TCanvas, gPad, gROOT

from yaffa import utils
from yaffa import logger as log

utils.style.SetStyle()

def draw_objects(name, objects, drawopt='pe', normalize=False):
    c = TCanvas('c', '', 600, 600)

    for i, (leg, obj) in enumerate(objects.items()):
        obj.SetTitle(leg)
        obj.SetLineColor(i + 1)
        obj.SetLineWidth(2)
        
        obj.SetMarkerColor(i + 1)
        if normalize:
            obj.Scale(1./obj.Integral())

        obj.Draw(drawopt)
        
        if 'same' not in drawopt:
            drawopt += ' same'

    gPad.BuildLegend()

    c.SaveAs(f'qa/c{name}.pdf')

def do_triplet_qa(directory):
    se = directory.Get('SE/Analysis/hQ3VsMtVsMultVsCent')
    me = directory.Get('ME/Analysis/hQ3VsMtVsMultVsCent')

    if not se or not me:
        log.error('SE or ME THnSparse are not properly defined. Skipping triplet QA')
        return
    
    # Use a helper function to project THnSparse with name to avoid replacing existing histograms
    def proj(thn, axis, name):
        hist = thn.Projection(axis)
        hist.SetName(name)
        return hist

    draw_objects('Q3', {'SE': proj(se, 0, 'SE'), 'ME': proj(me, 0, 'ME')}, normalize=True)
    draw_objects('Mt', {'SE': proj(se, 1, 'SE'), 'ME': proj(me, 1, 'ME')}, normalize=True)
    draw_objects('Mult', {'SE': proj(se, 2, 'SE'), 'ME': proj(me, 2, 'ME')}, normalize=True)
    draw_objects('Cent', {'SE': proj(se, 3, 'SE'), 'ME': proj(me, 3, 'ME')}, normalize=True)

def process_combination(directory):
    for key in [k.GetName() for k in directory.GetListOfKeys()]:
        if key == 'TrackTrackTrack':
            do_triplet_qa(directory.Get(key))
            pass
        else:
            log.warning(f'QA not implemented for directory {key}')
            pass

def main(in_file : str):
    os.makedirs('qa', exist_ok=True)
    try:
        inFile = TFile.Open(in_file)
    except OSError:
        log.critical('Cannot open the input file')

    for key in [k.GetName() for k in inFile.GetListOfKeys()]:
        process_combination(inFile.Get(key))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('in_file', nargs='?', default='AnalysisResults.root', help='Input file (AnalysisResults.root)')
    parser.add_argument('-b', action='store_true', default=False, help='Set batch mode')
    args = parser.parse_args()

    gROOT.SetBatch(args.b)

    main(args.in_file)
