from ROOT import TFile, TH1D

f = TFile('/home/daniel/an/LPi/systematics/fit/ancestors_0_2000.root')
oFile = TFile('/home/daniel/phsw/yaffa/yaffa/utils/ancestors_LPiminus.root', 'recreate')

for iAnc, anc in enumerate(["common_lpimin_gaussian_pol3", "noncommon_lpimin_pol3"]):
    func = f.Get(f'{anc}/fFit')

    # Ccnvert to hist
    hTemplate = TH1D(f'hCF_{iAnc}', "", 500, 0, 2)
    for iBin in range(500):
        bc = func.Eval(hTemplate.GetBinCenter(iBin + 1) * 1000)
        hTemplate.SetBinContent(iBin + 1, bc)
    hTemplate.SetDirectory(0)
    hTemplate.Draw()
    hTemplate.Write()
oFile.Close()
