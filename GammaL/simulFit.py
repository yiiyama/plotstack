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

sources = {
    'el': ROOT.TFile.Open('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/FloatingVGammaE.root'),
    'mu': ROOT.TFile.Open('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/FloatingVGammaM.root')
}

leptons = ['el', 'mu']

targets = ROOT.TObjArray()
vgTemplates = ROOT.TObjArray()
qcdTemplates = dict([(ch, ROOT.TObjArray()) for ch in leptons])
vgQCDTemplates = dict([(ch, ROOT.TObjArray()) for ch in leptons])

nbins = {}

for ch in leptons:
    source = sources[ch]
    target = source.Get('PreTemplateFit/target')
    targets.Add(target)
    vgTemplates.Add(source.Get('PreTemplateFit/vg'))

    nbins[ch] = target.GetNbinsX()
    
qcdTemplates['el'].Add(sources['el'].Get('PreTemplateFit/qcd'))
temp = sources['mu'].Get('PreTemplateFit/qcd').Clone('elqcdEmpty')
temp.Reset()
qcdTemplates['el'].Add(temp)

temp = sources['el'].Get('PreTemplateFit/qcd').Clone('muqcdEmpty')
temp.Reset()
qcdTemplates['mu'].Add(temp)
qcdTemplates['mu'].Add(sources['mu'].Get('PreTemplateFit/qcd'))

vgQCDTemplates['el'].Add(sources['el'].Get('PreTemplateFit/vgQCD'))
temp = sources['mu'].Get('PreTemplateFit/vgQCD').Clone('elvgQCDEmpty')
temp.Reset()
vgQCDTemplates['el'].Add(temp)

temp = sources['el'].Get('PreTemplateFit/vgQCD').Clone('muvgQCDEmpty')
temp.Reset()
vgQCDTemplates['mu'].Add(temp)
vgQCDTemplates['mu'].Add(sources['mu'].Get('PreTemplateFit/vgQCD'))

vgScale = ROOT.RooRealVar('vg', 'vg', 1., -ROOT.RooNumber.infinity(), ROOT.RooNumber.infinity())
qcdScales = dict([(ch, ROOT.RooRealVar(ch + 'qcd', ch + 'qcd', 1., -ROOT.RooNumber.infinity(), ROOT.RooNumber.infinity())) for ch in leptons])
vgQCDScales = dict([(ch, ROOT.RooFormulaVar(ch + 'vgQCDScale', ch + 'vgQCDScale', '-1. * @0 * @1', ROOT.RooArgList(vgScale, qcdScales[ch]))) for ch in leptons])

fitter = ROOT.TemplateChi2Fitter.singleton()
fitter.setTarget(targets)
fitter.addTemplate(vgTemplates, 'vg', vgScale)
for ch in leptons:
    fitter.addTemplate(qcdTemplates[ch], ch + 'qcd', qcdScales[ch])
    fitter.addTemplate(vgQCDTemplates[ch], ch + 'vgQCD', vgQCDScales[ch])

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

trees = dict([(ch, sources[ch].Get('TemplateFitError/toys')) for ch in leptons])
targetContents = dict([(ch, array.array('d', [0.] * nbins[ch])) for ch in leptons])
targetErrors = dict([(ch, array.array('d', [0.] * nbins[ch])) for ch in leptons])
vgContents = dict([(ch, array.array('d', [0.] * nbins[ch])) for ch in leptons])
vgErrors = dict([(ch, array.array('d', [0.] * nbins[ch])) for ch in leptons])

for ch in leptons:
    trees[ch].SetBranchAddress('targetContents', targetContents[ch])
    trees[ch].SetBranchAddress('targetErrors', targetErrors[ch])
    trees[ch].SetBranchAddress('vgContents', vgContents[ch])
    trees[ch].SetBranchAddress('vgErrors', vgErrors[ch])

ntotal = nbins['el'] + nbins['mu']
targetC = array.array('d', [0.] * ntotal)
targetE = array.array('d', [0.] * ntotal)
vgC = array.array('d', [0.] * ntotal)
vgE = array.array('d', [0.] * ntotal)
qcdC = dict([(ch, array.array('d', [0.] * ntotal)) for ch in leptons])
qcdE = dict([(ch, array.array('d', [0.] * ntotal)) for ch in leptons])
vgQCDC = dict([(ch, array.array('d', [0.] * ntotal)) for ch in leptons])
vgQCDE = dict([(ch, array.array('d', [0.] * ntotal)) for ch in leptons])

for (ch, ind, offset) in [('el', 0, 0), ('mu', 1, nbins['el'])]:
    for iX in range(nbins[ch]):
        qcd = qcdTemplates[ch].At(ind)
        qcdC[ch][iX + offset] = qcd.GetBinContent(iX + 1)
        qcdE[ch][iX + offset] = qcd.GetBinError(iX + 1)
        vgQCD = vgQCDTemplates[ch].At(ind)
        vgQCDC[ch][iX + offset] = vgQCD.GetBinContent(iX + 1)
        vgQCDE[ch][iX + offset] = vgQCD.GetBinError(iX + 1)

vgScales = []

iEntry = 0
while trees['el'].GetEntry(iEntry) > 0 and trees['mu'].GetEntry(iEntry) > 0:
    iEntry += 1

    for ch, offset in [('el', 0), ('mu', nbins['el'])]:
        for iX in range(nbins[ch]):
            targetC[iX + offset] = targetContents[ch][iX]
            targetE[iX + offset] = targetErrors[ch][iX]
            vgC[iX + offset] = vgContents[ch][iX]
            vgE[iX + offset] = vgErrors[ch][iX]

    fitter.setTarget(ntotal, targetC, targetE)
    fitter.addTemplate(ntotal, vgC, vgE, 'vg', vgScale)
    for ch in leptons:
        fitter.addTemplate(ntotal, qcdC[ch], qcdE[ch], ch + 'qcd', qcdScales[ch])
        fitter.addTemplate(ntotal, vgQCDC[ch], vgQCDE[ch], ch + 'vgQCD', vgQCDScales[ch])

    if fitter.fit(-1) != 0: continue
    
    vgScales.append(vgScale.getVal())

vgScales.sort()

errHigh = vgScales[int(len(vgScales) * 0.84)] - central
errLow = central - vgScales[int(len(vgScales) * 0.16)]
scaleErr = max(errHigh, errLow)
error = math.sqrt(scaleErr * scaleErr + fitErr * fitErr)

output.cd()
output.Write()

vgGraph = ROOT.TGraphErrors(1)
vgGraph.SetPoint(0, 0., central)
vgGraph.SetPointError(0, 0., error)
vgGraph.Write('VGamma')

for ch in leptons:
    graph = ROOT.TGraphErrors(1)
    graph.SetPoint(0, 0., qcdCentral[ch])
    graph.SetPointError(0, 0., qcdErr[ch])
    graph.Write(ch + 'QCD')

print 'Cross section:', central, '+-', error
