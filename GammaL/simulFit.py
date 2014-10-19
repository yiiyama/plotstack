import sys
import array
import math
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

rootlogon = ROOT.gEnv.GetValue("Rint.Logon", "")
if rootlogon:
    ROOT.gROOT.Macro(rootlogon)

ROOT.gSystem.Load('libRooFit.so')
ROOT.gSystem.Load('/afs/cern.ch/user/y/yiiyama/src/Common/fitting/libCommonFitting.so')

output = ROOT.TFile.Open('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/simulFit.root', 'recreate')

templatePlot = 'DPhiLeptonMetMet4070'

leptons = ['el', 'mu']

sources = {
    'el': ROOT.TFile.Open('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/FloatingVGammaE.root'),
    'mu': ROOT.TFile.Open('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/FloatingVGammaM.root')
}

candSamples = {
    'el': 'PhotonAndElectron',
    'mu': 'PhotonAndMuon'
}

datasets = {
    'el': 'DataE',
    'mu': 'DataM'
}

fakeSamples = {
    'el': 'PhotonAndFakeElectron',
    'mu': 'PhotonAndFakeMuon'
}

targets = ROOT.TObjArray()
vgTemplates = ROOT.TObjArray()
vgShiftTemplates = dict([(ch, ROOT.TObjArray()) for ch in leptons])
qcdTemplates = dict([(ch, ROOT.TObjArray()) for ch in leptons])
scaleNorms = dict([(ch, ROOT.TObjArray()) for ch in leptons])

nbins = {}

for ch in leptons:
    source = sources[ch]
    target = source.Get('PreTemplateFit/target')
    targets.Add(target)

    nbins[ch] = target.GetNbinsX()

    emptyTemplate = target.Clone(ch + 'empty')
    emptyTemplate.Reset()

    

    nominal = None
    scaleUp = None
    scaleDown = None
    for sample in ['WGToLNuG_PtG-30-50', 'WGToLNuG_PtG-50-130', 'WGToLNuG_PtG-130', 'ZGToLLG_PtG-5-130', 'ZGToLLG_PtG-130']:
        basename = 'components/' + templatePlot + '/' + templatePlot + '_' + sample + '_' + candSamples[ch]

        if nominal is None:
            nominal = sources[ch].Get(basename)
        else:
            nominal.Add(sources[ch].Get(basename))
        if scaleUp is None:
            scaleUp = sources[ch].Get(basename + '_ScaleUp')
        else:
            scaleUp.Add(sources[ch].Get(basename + '_ScaleUp'))
        if scaleDown is None:
            scaleDown = sources[ch].Get(basename + '_ScaleDown')
        else:
            scaleDown.Add(sources[ch].Get(basename + '_ScaleDown'))

    vgTemplates.Add(nominal)
    scaleUp.Add(scaleDown, -1.)
    scaleUp.Scale(0.5)
    for iX in range(1, scaleUp.GetNbinsX() + 1):
        scaleUp.SetBinError(iX, 0.)

    for ch2 in leptons:
        scaleNorms[ch2].Add(emptyTemplate)
        if ch2 == ch:
            qcdTemplates[ch2].Add(sources[ch].Get('components/' + templatePlot + '/' + templatePlot + '_' + datasets[ch2] + '_' + fakeSamples[ch2]))
            vgShiftTemplates[ch2].Add(scaleUp)
        else:
            qcdTemplates[ch2].Add(emptyTemplate)
            vgShiftTemplates[ch2].Add(emptyTemplate)


corrScale = ROOT.TH1D('corrScale', 'Correction Scale', 1, 0., 1.)
corrScale.Fill(0.5)
corrScale.SetBinError(1, 0.)
corrScaleEmpty = ROOT.TH1D('corrScaleEmpty', 'Dummy', 1, 0., 1.)
corrScaleError = ROOT.TH1D('corrScaleError', 'Dummy chi2 denominator', 1, 0., 1.)
corrScaleError.SetBinError(1, 1.)
targets.Add(corrScaleError)
targets.Add(corrScaleError)
vgTemplates.Add(corrScaleEmpty)
vgTemplates.Add(corrScaleEmpty)
for ch in leptons:
    qcdTemplates[ch].Add(corrScaleEmpty)
    qcdTemplates[ch].Add(corrScaleEmpty)
    vgShiftTemplates[ch].Add(corrScaleEmpty)
    vgShiftTemplates[ch].Add(corrScaleEmpty)
    for ch2 in leptons:
        if ch2 == ch:
            scaleNorms[ch2].Add(corrScale)
        else:
            scaleNorms[ch2].Add(corrScaleEmpty)


vgScale = ROOT.RooRealVar('vg', 'vg', 1., 0., 5.)
effShifts = dict([(ch, ROOT.RooRealVar(ch + 'effShift', ch + 'effShift', 0., -100., 100.)) for ch in leptons])
effScales = dict([(ch, ROOT.RooFormulaVar(ch + 'eff', ch + 'eff', '@0 * @1', ROOT.RooArgList(vgScale, effShifts[ch]))) for ch in leptons])
qcdScales = dict([(ch, ROOT.RooRealVar(ch + 'qcd', ch + 'qcd', 0.1, 0., 10.)) for ch in leptons])

fitter = ROOT.TemplateChi2Fitter.singleton()
fitter.setTarget(targets)
fitter.addTemplate(vgTemplates, 'vg', vgScale)
for ch in leptons:
    fitter.addTemplate(vgShiftTemplates[ch], ch + 'vgShift', effScales[ch])
    fitter.addTemplate(qcdTemplates[ch], ch + 'qcd', qcdScales[ch])
    fitter.addTemplate(scaleNorms[ch], ch + 'eff', effShifts[ch])

directory = output.mkdir('PreTemplateFit')
fitter.plot(directory)

fitter.fit()

directory = output.mkdir('PostTemplateFit')
fitter.plot(directory)

central = vgScale.getVal()
fitErr = vgScale.getError()

qcdCentral = dict([(ch, qcdScales[ch].getVal()) for ch in leptons])
qcdErr = dict([(ch, qcdScales[ch].getError()) for ch in leptons])

print 'VGamma:', central, '+-', fitErr
print 'elQCD:', qcdCentral['el'], '+-', qcdErr['el']
print 'muQCD:', qcdCentral['mu'], '+-', qcdErr['mu']

#trees = dict([(ch, sources[ch].Get('TemplateFitError/toys')) for ch in leptons])
#
#ntotal = nbins['el'] + nbins['mu']
#targetC = array.array('d', [0.] * ntotal)
#targetE = array.array('d', [0.] * ntotal)
#vgC = array.array('d', [0.] * ntotal)
#vgE = array.array('d', [0.] * ntotal)
#qcdC = dict([(ch, array.array('d', [0.] * ntotal)) for ch in leptons])
#qcdE = dict([(ch, array.array('d', [0.] * ntotal)) for ch in leptons])
#
#for (ch, ind, offset) in [('el', 0, 0), ('mu', 1, nbins['el'])]:
#    for iX in range(nbins[ch]):
#        qcd = qcdTemplates[ch].At(ind)
#        qcdC[ch][iX + offset] = qcd.GetBinContent(iX + 1)
#        qcdE[ch][iX + offset] = qcd.GetBinError(iX + 1)
#
#vgScales = []
#
#iEntry = 0
#while trees['el'].GetEntry(iEntry) > 0 and trees['mu'].GetEntry(iEntry) > 0:
#    if iEntry % 100 == 0: print iEntry
#    iEntry += 1
#
#    for ch, offset in [('el', 0), ('mu', nbins['el'])]:
#        for iX in range(nbins[ch]):
#            targetC[iX + offset] = targetContents[ch][iX]
#            targetE[iX + offset] = targetErrors[ch][iX]
#            vgC[iX + offset] = vgContents[ch][iX]
#            vgE[iX + offset] = vgErrors[ch][iX]
#
#    fitter.setTarget(ntotal, targetC, targetE)
#    fitter.addTemplate(ntotal, vgC, vgE, 'vg', vgScale)
#    for ch in leptons:
#        fitter.addTemplate(ntotal, qcdC[ch], qcdE[ch], ch + 'qcd', qcdScales[ch])
#
#    if fitter.fit(-1) != 0: continue
#    
#    vgScales.append(vgScale.getVal())
#
#vgScales.sort()
#
#errHigh = vgScales[int(len(vgScales) * 0.84)] - central
#errLow = central - vgScales[int(len(vgScales) * 0.16)]
#scaleErr = max(errHigh, errLow)
#error = math.sqrt(scaleErr * scaleErr + fitErr * fitErr)
#
output.cd()
output.Write()
#
#vgGraph = ROOT.TGraphErrors(1)
#vgGraph.SetPoint(0, 0., central)
#vgGraph.SetPointError(0, 0., error)
#vgGraph.Write('VGamma')
#
#for ch in leptons:
#    graph = ROOT.TGraphErrors(1)
#    graph.SetPoint(0, 0., qcdCentral[ch])
#    graph.SetPointError(0, 0., qcdErr[ch])
#    graph.Write(ch + 'QCD')

#print 'Cross section:', central, '+-', error
