'''
Script to produce the QA plots.
'''

import os
import re
import argparse

from ROOT import TFile, TCanvas, TLatex, TLegend, SetOwnership, gPad, gROOT, kRed

from yaffa import utils
from yaffa import logger as log

utils.style.SetStyle()

# Mass (GeV/c^2) and charge of the particles that can be identified from the QA histograms. To analyze a new system,
# add here the corresponding (mass, charge) pair together with the label to be displayed in the plots
KNOWN_PARTICLES = {
    (0.938272, +1): 'p',
    (0.938272, -1): '#bar{p}',
}

# Maximum distance between the measured and the expected mass to consider a particle as identified
MASS_TOLERANCE = 0.001

# Size of the legend text. Set explicitly so that a long header can be shrunk to fit in one line
LEGEND_TEXT_SIZE = 0.04

# Multipage PDF files that are still open. Each of them collects the same quantity for the different systems
openPDFs = []

def save_canvas(canvas, name):
    '''Save the canvas as a new page of the multipage PDF file dedicated to the quantity called name.'''
    path = f'qa/c{name}.pdf'

    if path not in openPDFs:
        canvas.Print(f'{path}[')
        openPDFs.append(path)

    canvas.Print(path)

def close_pdfs():
    '''Close the multipage PDF files. Until this is done the files cannot be opened by a PDF viewer.'''
    canvas = TCanvas('cClose', '', 600, 600)

    for path in openPDFs:
        canvas.Print(f'{path}]')

    openPDFs.clear()

def draw_objects(name, objects, drawopt='pe', normalize=False, header=''):
    c = TCanvas('c', '', 600, 600)

    empty = []
    for i, (leg, obj) in enumerate(objects.items()):
        obj.SetTitle(leg if leg else '')
        obj.SetLineColor(i + 1)
        obj.SetLineWidth(2)

        obj.SetMarkerColor(i + 1)
        if normalize:
            integral = obj.Integral()
            if integral == 0:
                log.warning(f'The {leg if leg else name} histogram is empty. Skipping normalization')
                empty.append(f'EMPTY ({leg})' if leg else 'EMPTY')
            else:
                obj.Scale(1./integral)

        obj.Draw(drawopt)

        if 'same' not in drawopt:
            drawopt += ' same'

    draw_legend(objects, header)

    # Flag the empty histograms. Drawn after the legend so that they don't end up in it
    tl = TLatex()
    tl.SetNDC()
    tl.SetTextColor(kRed)
    tl.SetTextSize(0.06)
    for iEmpty, label in enumerate(empty):
        tl.DrawLatex(0.4, 0.5 - 0.07 * iEmpty, label)

    save_canvas(c, name)

def draw_legend(objects, header):
    '''Draw the legend of the objects in the current pad, with the analyzed system as header.'''
    # Make room at the top of the frame for the legend, keeping all the objects inside the drawing range. The range is
    # defined by the first object, the ones drawn on top of it with 'same' cannot enlarge it
    histograms = [obj for obj in objects.values() if obj.GetDimension() == 1]
    if histograms:
        first = list(objects.values())[0]
        first.SetMinimum(min(0, *[obj.GetMinimum() for obj in histograms]))
        first.SetMaximum(1.5 * max(obj.GetMaximum() for obj in histograms))

    # Update the pad, otherwise the axis ranges needed to compute the width of the header are not yet defined
    gPad.Modified()
    gPad.Update()

    xMin, xMax = gPad.GetLeftMargin(), 1 - gPad.GetRightMargin()
    if any(label for label in objects):
        legend = gPad.BuildLegend(xMin, 0.79, xMax, 0.91, header)
        legend.SetNColumns(max(len(objects), 1))
    else:  # None of the objects has a label: only the header is worth showing
        legend = TLegend(xMin, 0.79, xMax, 0.91, header)
        SetOwnership(legend, False)  # Let the pad own the legend, otherwise it is deleted before being painted
        legend.Draw()
    legend.SetFillStyle(0)

    # Set the text size explicitly, otherwise it is computed from the size of the box, and shrink it if the header
    # is too long to fit in one line
    text = TLatex(0, 0, header)
    text.SetTextSize(LEGEND_TEXT_SIZE)
    width = text.GetXsize() / (gPad.GetUxmax() - gPad.GetUxmin()) * (xMax - xMin)
    available = 0.90 * (xMax - xMin)
    legend.SetTextSize(LEGEND_TEXT_SIZE * min(1, available / width) if width > 0 else LEGEND_TEXT_SIZE)

    return legend

def get_particle(directory):
    '''Identify the particle analyzed in a combination directory from its charge and mass.'''
    tracks = sorted(
        [key.GetName() for key in directory.GetListOfKeys() if re.fullmatch(r'Track\d+', key.GetName())],
        key=lambda name: int(name.removeprefix('Track')))

    if not tracks:
        log.critical(f'No Track directory found in {directory.GetName()}. Cannot determine the analyzed system')

    hSign = directory.Get(f'{tracks[0]}/Analysis/hSign')
    hMass = directory.Get(f'{tracks[0]}/Analysis/hMass')

    if not hSign or not hMass:
        log.critical(f'hSign or hMass are missing in {directory.GetName()}/{tracks[0]}. '
                     'Cannot determine the analyzed system')

    charge = round(hSign.GetMean())
    mass = hMass.GetMean()

    for (knownMass, knownCharge), particle in KNOWN_PARTICLES.items():
        if charge == knownCharge and abs(mass - knownMass) < MASS_TOLERANCE:
            return particle

    log.critical(f'Unknown system: no known particle has charge {charge:+d} and mass {mass:.6f} GeV/c^2. '
                 'You may need to add your system to KNOWN_PARTICLES')
    return None

def do_triplet_qa(directory, header=''):
    se = directory.Get('SE/Analysis/hQ3VsMtVsMultVsCent')
    me = directory.Get('ME/Analysis/hQ3VsMtVsMultVsCent')

    if not se or not me:
        log.error('SE or ME THnSparse are not properly defined. Skipping triplet QA')
        return

    # Use a helper function to project THnSparse with name to avoid replacing existing histograms
    def proj(thn, axis, name):
        if not isinstance(axis, tuple):
            axis = (axis,)

        hist = thn.Projection(*axis)
        hist.SetName(name)
        return hist

    draw_objects('Q3', {'SE': proj(se, 0, 'SE'), 'ME': proj(me, 0, 'ME')}, normalize=True, header=header)
    draw_objects('Mt', {'SE': proj(se, 1, 'SE'), 'ME': proj(me, 1, 'ME')}, normalize=True, header=header)
    draw_objects('Mult', {'SE': proj(se, 2, 'SE'), 'ME': proj(me, 2, 'ME')}, normalize=True, header=header)
    draw_objects('Cent', {'SE': proj(se, 3, 'SE'), 'ME': proj(me, 3, 'ME')}, normalize=True, header=header)
    draw_objects('Q3VsMult', {None: proj(se, (0, 2), 'SE')}, drawopt='colz', normalize=True, header=header)

def process_combination(directory, particle):
    for key in [k.GetName() for k in directory.GetListOfKeys()]:
        if key == 'TrackTrackTrack':
            # For the time being the triplets are assumed to be made of three particles of the same type
            system = '-'.join([particle] * key.count('Track'))
            do_triplet_qa(directory.Get(key), f'{system} ({directory.GetName()})')
        else:
            log.warning(f'QA not implemented for directory {key}')

def main(in_file : str):
    os.makedirs('qa', exist_ok=True)
    try:
        inFile = TFile.Open(in_file)
    except OSError:
        log.critical('Cannot open the input file')

    directories = [inFile.Get(key.GetName()) for key in inFile.GetListOfKeys()]

    # Identify all the systems before drawing anything, so that an unknown one is reported before producing any plot
    particles = [get_particle(directory) for directory in directories]

    for directory, particle in zip(directories, particles):
        process_combination(directory, particle)

    close_pdfs()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('in_file', nargs='?', default='AnalysisResults.root', help='Input file (AnalysisResults.root)')
    parser.add_argument('-b', action='store_true', default=False, help='Set batch mode')
    args = parser.parse_args()

    gROOT.SetBatch(args.b)

    main(args.in_file)
