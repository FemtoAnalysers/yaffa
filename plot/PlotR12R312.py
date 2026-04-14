from ROOT import TF2, TCanvas, TFile, gInterpreter, gROOT
gInterpreter.Declare('#include "../src/RootFunctions.hxx"')
from ROOT import SourceAAAJC

gROOT.SetBatch(True)

if __name__ == '__main__':
    fSourceAAAJC = TF2('fSourceAAAJC', SourceAAAJC, 0, 8, 0, 8, 1)
    fSourceAAAJC.SetParameter(0, 1.0)
    fSourceAAAJC.SetNpx(100)
    fSourceAAAJC.SetNpy(100)
    fSourceAAAJC.SetTitle(';r_{12} (fm);r_{3,12} (fm);(4#pi)^{2} r_{12}^{2} r_{3,12}^{2} S(r_{12}, r_{3,12})')
    fSourceAAAJC.GetZaxis().SetTitleOffset(1.7)

    cSource = TCanvas('cSource', '', 600, 600)
    cSource.SetTopMargin(0.02)
    cSource.SetRightMargin(0.2)
    fSourceAAAJC.Draw('colz')
    cSource.SaveAs('cR12R312.pdf')
    
    oFile = TFile('Source3B.root', 'recreate')
    fSourceAAAJC.Write()
    oFile.Close()
