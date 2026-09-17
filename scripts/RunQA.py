'''
Script to produce the QA plots.
'''

import os
import re
import argparse

from ROOT import TFile, TCanvas, TDatabasePDG, TH2F, TLatex, TLegend, SetOwnership, gPad, gROOT, gStyle, kRed

from yaffa import utils
from yaffa import logger as log

utils.style.SetStyle()

# Print the bin contents drawn with the 'text' option as integers
gStyle.SetPaintTextFormat('.0f')

PARTICLES_LABELS = {
    2212: "p",
    -2212: "ap",
    321: "Kplus",
    -321: "Kminus",
}

PARTICLES_LATEX = {
    2212: "p",
    -2212: "#bar{p}",
    321: "K^{+}",
    -321: "K^{#minus}",
}


# Maximum distance between the measured and the expected mass to consider a particle as identified
MASS_TOLERANCE = 0.001

# Maximum Q3 (GeV/c) of the triplets counted in the summary plot of the yields
Q3_THRESHOLD = 0.6

# Size of the legend text. Set explicitly so that a long header can be shrunk to fit in one line
LEGEND_TEXT_SIZE = 0.04

# Multipage PDF files that are still open. Each of them collects the same quantity for the different systems
openPDFs = []

def save_canvas(canvas, name, subdir=''):
    '''Save the canvas as a new page of the multipage PDF file dedicated to the quantity called name.'''
    directory = os.path.join('qa', subdir)
    path = os.path.join(directory, f'c{name}.pdf')

    if path not in openPDFs:
        os.makedirs(directory, exist_ok=True)
        canvas.Print(f'{path}[')
        openPDFs.append(path)

    canvas.Print(path)

def close_pdfs():
    '''Close the multipage PDF files. Until this is done the files cannot be opened by a PDF viewer.'''
    canvas = TCanvas('cClose', '', 600, 600)

    for path in openPDFs:
        canvas.Print(f'{path}]')

    openPDFs.clear()

def draw_objects(name, objects, drawopt='pe', normalize=False, header='', title='', subdir=''):
    c = TCanvas('c', '', 600, 600)

    empty = []
    for i, (leg, obj) in enumerate(objects.items()):
        obj.SetTitle(leg if leg else title)
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

    save_canvas(c, name, subdir)

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

def get_particle_label(directory):
    '''Identify the particle analyzed in a combination directory from its charge and mass.'''

    hSign = directory.Get(f'Analysis/hSign')
    hMass = directory.Get(f'Analysis/hMass')

    if not hSign or not hMass:
        log.critical(f'hSign or hMass are missing in {directory.GetName()}/{tracks[0]}. '
                     'Cannot determine the analyzed system')

    charge = round(hSign.GetMean())
    mass = hMass.GetMean()

    database = TDatabasePDG.Instance()

    for pdg in [2212, 321]:
        if abs(database.GetParticle(pdg).Mass() - mass) < MASS_TOLERANCE:
            return PARTICLES_LABELS[charge * pdg]

    raise ValueError("Mass and charge combination not implemented")


def get_particle_latex(directory):
    '''Identify the particle analyzed in a combination directory from its charge and mass.'''

    hSign = directory.Get(f'Analysis/hSign')
    hMass = directory.Get(f'Analysis/hMass')

    if not hSign or not hMass:
        log.critical(f'hSign or hMass are missing in {directory.GetName()}/{tracks[0]}. '
                     'Cannot determine the analyzed system')

    charge = round(hSign.GetMean())
    mass = hMass.GetMean()

    database = TDatabasePDG.Instance()

    for pdg in [2212, 321]:
        if abs(database.GetParticle(pdg).Mass() - mass) < MASS_TOLERANCE:
            return PARTICLES_LATEX[charge * pdg]

    raise ValueError("Mass and charge combination not implemented")

def get_particles_latex(directory):
    '''Identify the three particles of the triplet analyzed in a combination directory.'''
    track_dirs = [directory.Get(key.GetName()) for key in directory.GetListOfKeys() if re.fullmatch(r'Track\d+', key.GetName())]

    if len(track_dirs) == 1:
        return [get_particle_latex(track_dirs[0])] * 3 # identical particles
    if len(track_dirs) == 2:
        return [get_particle_latex(track_dirs[0])] * 2 + [get_particle_latex(track_dirs[1])] # 2 identical, 1 different
    if len(track_dirs) == 3:
        return [get_particle_latex(track_dir) for track_dir in track_dirs]

    raise ValueError("Wrong number of particles")

def get_system_latex(directory):
    '''Name of the system analyzed in a combination directory, obtained by chaining its three particles.'''
    return ''.join(get_particles_latex(directory))

def do_track_qa(directory, header=''):
    '''Draw the QA of the tracks of one of the particle species of the combination.'''
    # Subdirectory of qa where the track plots are saved
    SUBDIR = 'triplet/track_' + get_particle_label(directory)

    for name in ['hPt', 'hEta', 'hPhi', 'hSign', 'hMass']:
        hist = directory.Get(f'Analysis/{name}')

        if not hist:
            log.error(f'Analysis/{name} is missing. Skipping it')
            continue

        draw_objects(name.removeprefix('h'), {None: hist}, header=header, subdir=SUBDIR)

def do_collision_qa(directory):
    '''Draw the event-level QA of the collisions.'''
    # Subdirectory of qa where the collision plots are saved
    SUBDIR = 'collision'

    for name in ['hPosZ', 'hMult', 'hCent', 'hMagField']:
        hist = directory.Get(f'Analysis/{name}')

        if not hist:
            log.error(f'Analysis/{name} is missing. Skipping it')
            continue

        draw_objects(name.removeprefix('h'), {None: hist}, subdir=SUBDIR)

def do_triplet_qa(directory, header=''):
    # Subdirectory of qa where the triplet plots are saved
    SUBDIR = 'triplet'

    se = directory.Get('SE/Analysis/hQ3VsMtVsMultVsCent')
    me = directory.Get('ME/Analysis/hQ3VsMtVsMultVsCent')

    if not se or not me:
        log.error('SE or ME THnSparse are not properly defined. Skipping triplet QA')
        return 0

    # Use a helper function to project THnSparse with name to avoid replacing existing histograms
    def proj(thn, axis, name):
        if not isinstance(axis, tuple):
            axis = (axis,)

        hist = thn.Projection(*axis)
        hist.SetName(name)
        return hist

    draw_objects('Q3', {'SE': proj(se, 0, 'SE'), 'ME': proj(me, 0, 'ME')}, normalize=True, header=header, subdir=SUBDIR)
    draw_objects('Mt', {'SE': proj(se, 1, 'SE'), 'ME': proj(me, 1, 'ME')}, normalize=True, header=header, subdir=SUBDIR)
    draw_objects('Mult', {'SE': proj(se, 2, 'SE'), 'ME': proj(me, 2, 'ME')}, normalize=True, header=header, subdir=SUBDIR)
    draw_objects('Cent', {'SE': proj(se, 3, 'SE'), 'ME': proj(me, 3, 'ME')}, normalize=True, header=header, subdir=SUBDIR)
    draw_objects('Q3VsMult', {None: proj(se, (0, 2), 'SE')}, drawopt='colz', normalize=True, header=header, subdir=SUBDIR)

    # Count the triplets in a dedicated projection, the one drawn above is normalized. The bin containing the threshold
    # is excluded so that only the triplets certainly below it are counted
    hQ3 = proj(se, 0, 'Q3Yield')
    return hQ3.Integral(1, hQ3.FindBin(Q3_THRESHOLD) - 1)

def draw_yields(yields):
    '''Draw the number of triplets below the Q3 threshold, one bin per system.'''
    # Subdirectory of qa where the triplet plots are saved
    SUBDIR = 'triplet'

    if not yields:
        log.warning('No triplet was analyzed. Skipping the plot of the yields')
        return

    # The first two particles of each system are on the y axis, the third one on the x axis
    pairs = sorted({pair for pair, _ in yields})
    singles = sorted({single for _, single in yields})

    hYields = TH2F('hYields', '', len(singles), 0, len(singles), len(pairs), 0, len(pairs))
    for iSingle, single in enumerate(singles):
        hYields.GetXaxis().SetBinLabel(iSingle + 1, single)
    for iPair, pair in enumerate(pairs):
        hYields.GetYaxis().SetBinLabel(iPair + 1, pair)

    for (pair, single), counts in yields.items():
        hYields.SetBinContent(singles.index(single) + 1, pairs.index(pair) + 1, counts)

    # Enlarge the text with the number of triplets, otherwise it is barely readable
    hYields.SetMarkerSize(1.5)

    gStyle.SetOptTitle(1)
    gStyle.SetPadTopMargin(0.1)

    title = f'Triplets with #it{{Q}}_{{3}} < {Q3_THRESHOLD * 1000:.0f} MeV/#it{{c}}'
    draw_objects('CountTriplets', {None: hYields}, drawopt='col text', title=title, subdir=SUBDIR)

    # Restore the style, so that the other canvases are not affected
    gStyle.SetOptTitle(0)
    gStyle.SetPadTopMargin(0.05)

def process_combination(directory, particle):
    '''Do the QA of a combination directory and return the triplets below the Q3 threshold of its system.'''
    yields = {}

    for key in [k.GetName() for k in directory.GetListOfKeys()]:
        if key == 'TrackTrackTrack':
            particles = get_particles_latex(directory)
            system = ''.join(particles)
            counts = do_triplet_qa(directory.Get(key), f'{system} ({directory.GetName()})')
            yields[(''.join(particles[:2]), particles[2])] = counts
        elif key == 'Collisions':
            do_collision_qa(directory.Get(key))
        elif re.fullmatch(r'Track\d+', key):
            track = directory.Get(key)
            do_track_qa(track, f'{get_particle_latex(track)} ({directory.GetName()})')
        else:
            log.warning(f'QA not implemented for directory {key}')

    return yields

def main(in_file : str):
    os.makedirs('qa', exist_ok=True)
    try:
        inFile = TFile.Open(in_file)
    except OSError:
        log.critical('Cannot open the input file')

    directories = [inFile.Get(key.GetName()) for key in inFile.GetListOfKeys()]

    # Identify all the systems before drawing anything, so that an unknown one is reported before producing any plot
    particles = [get_system_latex(directory) for directory in directories]

    yields = {}
    for directory, particle in zip(directories, particles):
        for system, counts in process_combination(directory, particle).items():
            yields[system] = yields.get(system, 0) + counts

    draw_yields(yields)

    close_pdfs()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('in_file', nargs='?', default='AnalysisResults.root', help='Input file (AnalysisResults.root)')
    parser.add_argument('-b', action='store_true', default=False, help='Set batch mode')
    args = parser.parse_args()

    gROOT.SetBatch(args.b)

    main(args.in_file)
