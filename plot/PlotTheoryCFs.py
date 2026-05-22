import numpy as np
from pathlib import Path


from ROOT import TCanvas, TGraph, TLegend

from yaffa import utils
from yaffa import logger as log

import os
from dotenv import load_dotenv
from pathlib import Path

env_path = Path(__file__).resolve().parent.parent / ".env"
print(f'Loading env from {env_path}')
if not load_dotenv(dotenv_path=env_path, verbose=True, override=True):
    print("Environment variables in .env not loaded")
YAFFA_PATH = os.getenv("YAFFA")
if not YAFFA_PATH:
    print("\033[33mWARNING: Path to yaffa is empty, something might break!\033[0m")

from ROOT import gInterpreter
gInterpreter.Declare(f'#include "{YAFFA_PATH}/src/cpp/RootFunctions.hxx"')
from ROOT import _SourceAAA


utils.style.SetStyle()

idx2Q3 = {
    1  :round( 53.05484208147264, 1),
    2  :round( 67.10965675120087, 1),
    3  :round( 91.89368207265395, 1),
    4  : round(106.10968416294529, 1),
    5  : round(140.36991799538106, 1),
    6  : round(150.06175444236283, 1),
    7  : round(167.77414187800218, 1),
    8  : round(183.7873641453079, 1),
    9  : round(198.51304177826705, 1),
    10 : round(212.21936832589057, 1),
    11 : round(225.09263166354424, 1),
    12 : round(237.26846685937852, 1),
    13 : round(248.84926743388257, 1),
    14 : round(259.9145829670971, 1),
    15 : round(280.7398359907621, 1),
    16 : round(300.12350888472565, 1),
    17 : round(335.54828375600437, 1),
    18 : round(375.15438610590707, 1),
    19 : round(410.96103963441595, 1),
    20 : round(474.53693371875704, 1),
    21 : round(581.1866758579381, 1),
    22 : round(671.0965675120087, 1),
    23 : round(750.3087722118141, 1),
    24 : round(821.9220792688319, 1),
}

def main():
    file = "../secrets/theory/wf/ppp.dat"
    path = Path(file)

    if not path.exists():
        log.error(f'File "{file}" does not exist.')
        return

    data = np.loadtxt(file)
    data = data[data[:, 0] > 1.2]
    
    c = TCanvas('c', '', 600, 600)
    c.DrawFrame(0, 0, 30, 2.5, ';#rho (fm);|#Psi|^{2}')
    leg = TLegend(0.5, 0.5, 0.85, 0.9)

    gWF = []
    hyp_rad = np.ascontiguousarray(data[:, 0], dtype='double')
    print(hyp_rad)
    source = [_SourceAAA(rho, 2.6) for rho in hyp_rad]    

    cf = []
    color = 1
    for iQ3, Q3 in idx2Q3.items():
        wf = np.ascontiguousarray(data[:, iQ3], dtype='double')
        gWF.append(TGraph(len(hyp_rad), hyp_rad, wf))  
        utils.style.SetObjectStyle(gWF[-1])
        gWF[-1].SetLineColor(color)

        cf.append(source @ wf / sum(source))

        if (iQ3 - 1) % 3 == 0:
            gWF[-1].DrawClone('l same')
            leg.AddEntry(gWF[-1], f'Q_{{3}} = {idx2Q3[iQ3]} GeV/#it{{c}}', 'l')
            color += 1
    leg.Draw('same')
    c.SaveAs('WaveFunctionsPPP.pdf')

    cCF = TCanvas('cCF', '', 600, 600)
    cCF.DrawFrame(0, 0, 799.9, 2.2, ';Q_{3} (GeV/#it{c});C_{ppp}')
    gCF = TGraph(len(idx2Q3), np.array(list(idx2Q3.values()), dtype='double'), np.array(cf, dtype='double'))
    gCF.Draw('lsame')
    leg = TLegend(0.2, 0.85, 0.5, 0.9)
    leg.AddEntry(gCF, 'E. Garrido et al., PLB 868 (2025) 139731')
    leg.Draw('same')

    cCF.SaveAs('CorrelationFunctionsPPP.pdf')

if __name__ == '__main__':
    main()
