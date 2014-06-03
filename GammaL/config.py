# Analysis-specific definitions
# Must define the following variables:
#  eventProcessor: event processor class
#  datasets: [Dataset]
#  weightCalc: {dataset name: {event class: weight calculator}}
#  stackConfigs: {stack name: StackConfig}
#  hdefs: [HDef]

import os
import ROOT

from dataset import Dataset
from histogram import HDef
from stack import Group, StackConfig
from floatfix import *

thisdir = os.path.dirname(os.path.abspath(__file__))

############## EVENT PROCESSOR ##############

ROOT.gROOT.LoadMacro(thisdir + '/GLSkimProcessor.cc+')

### EXPORTED ###
eventProcessor = ROOT.GLSkimProcessor


############## DATASET DEFINITIONS ##############

def realData(name, inputNames, L, filters):
    return Dataset(name, inputNames, True, L, 0., filters)

def mc(name, Leff, sigmaRelErr, filters):
    return Dataset(name, [name], False, Leff, sigmaRelErr, filters)

def addDataset(datasets, dataset):
    datasets[dataset.name] = dataset

data_E = ('PhotonAndElectron', 'ElePhotonAndElectron', 'FakePhotonAndElectron', 'PhotonAndFakeElectron')
data_Mu = ('PhotonAndMuon', 'ElePhotonAndMuon', 'FakePhotonAndMuon', 'PhotonAndFakeMuon')
All = data_E + data_Mu
MC = ('PhotonAndElectron', 'PhotonAndMuon')
ElePhoton = ('PhotonAndElectron', 'PhotonAndMuon', 'ElePhotonAndElectron', 'ElePhotonAndMuon')
FakePhoton = ('PhotonAndElectron', 'PhotonAndMuon', 'FakePhotonAndElectron', 'FakePhotonAndMuon')
EleFakePhoton = ('PhotonAndElectron', 'PhotonAndMuon', 'ElePhotonAndElectron', 'ElePhotonAndMuon', 'FakePhotonAndElectron', 'FakePhotonAndMuon')

### EXPORTED ###
datasets = {}
addDataset(datasets, realData('DataE', ['PhotonA', 'DoublePhotonB', 'DoublePhotonC', 'DoublePhotonD'], 876.225 + 4411.704 + 7054.732 + 7369.007, data_E))
addDataset(datasets, realData('DataM', ['MuEGA', 'MuEGB', 'MuEGC', 'MuEGD'], 876.225 + 4411.704 + 7054.732 + 7360.046, data_Mu))
addDataset(datasets, mc('WGToLNuG', 4802358. / 461.6, 0.15, ElePhoton))
addDataset(datasets, mc('WGToLNuG_PtG-30-50', 3000000. / 75.48, 0.15, ElePhoton)) # using parton-level cross section and total number of LHE events
addDataset(datasets, mc('WGToLNuG_PtG-50-130', 3000000. / 9.633, 0.15, ElePhoton)) # same (parton-level simple sum * matching eff / matched events = parton-level simple sum / all events)
addDataset(datasets, mc('WGToLNuG_PtG-130', 471458. / 0.2571, 0.15, ElePhoton)) # here PREP value is correctly set to parton-level xsec * matching eff
addDataset(datasets, mc('ZGToLLG', 6588161. / 132.6, 0.15, ElePhoton))
addDataset(datasets, mc('WWGJets', 304285. / 0.528, 0.15, ElePhoton))
addDataset(datasets, mc('TTGJets', 1791552. / 1.444, 0.15, ElePhoton))
addDataset(datasets, mc('WW', 10000431. / 56., 0.15, ElePhoton)) # NLO from twiki
addDataset(datasets, mc('TTJetsSemiLept', 24953451. / (245.8 * 0.324 * 0.676 * 2), 0.15, ElePhoton)) # NNLO from twiki * BR(W->lnu BR) * BR(W->had) * 2. Sudakov = 0.999
addDataset(datasets, mc('TTJetsFullLept', 12011428. / (245.8 * 0.324 * 0.324), 0.15, ElePhoton)) # NNLO from twiki * BR(W->lnu)^2. Sudakov = 0.999
addDataset(datasets, mc('WWNoGamma', 10000431. / 56., 0.15, EleFakePhoton)) # NLO from twiki
addDataset(datasets, mc('TTJetsSemiLeptNoGamma', 24953451. / (245.8 * 0.324 * 0.676 * 2), 0.15, FakePhoton)) # NNLO from twiki * BR(W->lnu BR) * BR(W->had) * 2. Sudakov = 0.999
addDataset(datasets, mc('TTJetsFullLeptNoGamma', 12011428. / (245.8 * 0.324 * 0.324), 0.15, EleFakePhoton)) # NNLO from twiki * BR(W->lnu)^2. Sudakov = 0.999
addDataset(datasets, mc('WJetsToLNu', 57709905. / 36703.2, 0.15, FakePhoton)) # NNLO from twiki. Sudakov = 0.997
addDataset(datasets, mc('WJetsToLNu_PtW-50To70', 48426609. / (36703.2 * 811.2 / 30400.), 0.15, FakePhoton)) # 811.2 / 30400: MG5 xsec ratio
addDataset(datasets, mc('WJetsToLNu_PtW-70To100', 22447541. / (36703.2 * 428.9 / 30400.), 0.15, FakePhoton)) # 428.9 / 30400: MG5 xsec ratio
addDataset(datasets, mc('WJetsToLNu_PtW-100', 12742382. / (36703.2 * 228.9 / 30400.), 0.15, FakePhoton)) # 228.9 / 30400: MG5 xsec ratio
addDataset(datasets, mc('DYJetsToLL', 30459503. / 3531.9, 0.15, EleFakePhoton)) # NNLO from twiki. Sudakov = 0.993
#addDataset(datasets, mc('TChiwg_300', 54908. / 0.146, 0.00661 / 0.0146, All))
#addDataset(datasets, mc('T5wg_500_425', 56217. / (4.525 * 0.5), 0.16, All)) # 0.5 from combinatorics (g/g->aW/aW)
addDataset(datasets, Dataset('T5wg_1000_425', ['T5wg/skim_1000_425'], False, 59387. / (0.0244 * 0.5), 0.26, All)) # 0.5 from combinatorics (g/g->aW/aW)


############## WEIGHT CALCULATORS ##############

### EXPORTED ###
weightCalc = {
    datasets['DataE'].PhotonAndElectron: ROOT.GLEventWeight(),
    datasets['DataE'].ElePhotonAndElectron: ROOT.ElePhotonFunctionalWeight(),
    datasets['DataE'].FakePhotonAndElectron: ROOT.JetPhotonHLTIsoWeight(),
    datasets['DataE'].PhotonAndFakeElectron: ROOT.GLEventWeight(),
    datasets['DataM'].PhotonAndMuon: ROOT.GLEventWeight(),
    datasets['DataM'].ElePhotonAndMuon: ROOT.ElePhotonFunctionalWeight(),
    datasets['DataM'].FakePhotonAndMuon: ROOT.JetPhotonWeight(),
    datasets['DataM'].PhotonAndFakeMuon: ROOT.GLEventWeight()
}

for name in ['WGToLNuG', 'WGToLNuG_PtG-30-50', 'WGToLNuG_PtG-50-130', 'WGToLNuG_PtG-130', 'ZGToLLG', 'WWGJets', 'TTGJets', 'WW', 'TTJetsSemiLept', 'TTJetsFullLept']:
    weightCalc[datasets[name].PhotonAndElectron] = ROOT.GLEventWeight()
    weightCalc[datasets[name].PhotonAndMuon] = ROOT.GLEventWeight()
    weightCalc[datasets[name].ElePhotonAndElectron] = ROOT.ElePhotonFunctionalWeight()
    weightCalc[datasets[name].ElePhotonAndMuon] = ROOT.ElePhotonFunctionalWeight()

weightCalc[datasets['T5wg_1000_425'].PhotonAndElectron] = ROOT.GLEventWeight()
weightCalc[datasets['T5wg_1000_425'].PhotonAndMuon] = ROOT.GLEventWeight()
weightCalc[datasets['T5wg_1000_425'].ElePhotonAndElectron] = ROOT.ElePhotonFunctionalWeight()
weightCalc[datasets['T5wg_1000_425'].ElePhotonAndMuon] = ROOT.ElePhotonFunctionalWeight()
weightCalc[datasets['T5wg_1000_425'].FakePhotonAndElectron] = ROOT.JetPhotonHLTIsoWeight()
weightCalc[datasets['T5wg_1000_425'].FakePhotonAndMuon] = ROOT.JetPhotonWeight()
weightCalc[datasets['T5wg_1000_425'].PhotonAndFakeElectron] = ROOT.GLEventWeight()
weightCalc[datasets['T5wg_1000_425'].PhotonAndFakeMuon] = ROOT.GLEventWeight()

weightCalcMC = {
    'PhotonAndElectron': ROOT.GLEventWeight(),
    'PhotonAndMuon': ROOT.GLEventWeight(),
    'PhotonAndFakeElectron': ROOT.GLEventWeight(),
    'PhotonAndFakeMuon': ROOT.GLEventWeight(),
    'ElePhotonAndElectron': ROOT.MCElePhotonFunctionalWeight(),
    'ElePhotonAndMuon': ROOT.MCElePhotonFunctionalWeight(),
    'FakePhotonAndElectron': ROOT.MCJetPhotonHLTIsoWeight(),
    'FakePhotonAndMuon': ROOT.MCJetPhotonWeight()
}

for s, c in weightCalcMC.items():
    weightCalc.update(dict((d.samples[s], c) for d in datasets.values() if not d.realData and s in d.samples and d.samples[s] not in weightCalc))

############## HISTOGRAMS ###############

metBinning1 = [0., 10., 20., 25., 30., 34., 38., 42., 46., 50., 54., 58., 62., 70., 80., 90.,100., 120.,200.,300.,400.]
metBinning2 = [0., 4., 8., 12., 16., 20., 24., 28., 32., 36., 40., 44., 48., 52., 56., 60., 65., 70., 75., 80., 90.,100., 110.,120.,160.,200.,300.,400.]
mtBinning1 = [0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 140., 200., 300., 400.]
mtBinning2 = [0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 140., 160., 200., 240., 280., 400.]
massBinning = [0.] + [x * 10. for x in range(1, 31)] + [300. + x * 20. for x in range(1, 11)] + [500. + x * 50. for x in range(1, 11)]
metSqrtHtBinning = [0., 4.] + [x + 5. for x in range(0, 8)] + [15., 40.]
dPhiBinning = (61, -3.15, 3.15)

MET = ('#slash{E}_{T}', 'GeV')
MT = ('M_{T}', 'GeV')
PT = ('P_{T}', 'GeV')
DR = '#DeltaR'
DPHI = '#Delta#phi'

### EXPORTED ###
hdefs = [
    HDef('MetHighMt', 'MET (M_{T} > 140 GeV)', metBinning1, xtitle = MET, overflowable = True, mask = ('-inf', 'inf'), logscale = True),
    HDef('MetHighPhotonPt', 'MET (P_{T}^{#gamma} > 180 GeV)', metBinning1, xtitle = MET, overflowable = True, mask = ('-inf', 'inf'), logscale = True),
    HDef('MetLowPhotonPt', 'MET (P_{T}^{#gamma} < 80 GeV)', metBinning2, xtitle = MET, overflowable = True, mask = (70., 'inf'), logscale = True),
    HDef('Met', 'MET', metBinning2, xtitle = MET, overflowable = True, mask = (70., 'inf'), logscale = True),
    HDef('MetLowPW', 'MET (P_{W} < 200 GeV)', metBinning2, xtitle = MET, overflowable = True, mask = (70., 'inf'), logscale = True),
    HDef('MetOverSqrtHtLowPhotonPt', 'MET / #sqrt{H_{T}} (P_{T}^{#gamma} < 80 GeV)', metSqrtHtBinning, xtitle = ('#slash{E}_{T}/#sqrt{H_{T}}', '#sqrt{GeV}')),
    HDef('MtHighMet', 'M_{T}^{e} (#slash{E}_{T} > 120 GeV)', mtBinning1, xtitle = MT, overflowable = True, mask = ('-inf', 'inf'), logscale = True),
    HDef('MtLowMetLowPhotonPt', 'M_{T}^{e} (#slash{E}_{T} < 70 GeV, P_{T}^{#gamma} < 80 GeV)', mtBinning2, xtitle = MT, overflowable = True, logscale = True),
    HDef('Mt', 'M_{T}^{e}', mtBinning2, xtitle = MT, overflowable = True, mask = (100., 'inf'), logscale = True),
    HDef('MetMt', 'MET v M_{T}', (50, 0., 400, 50, 0., 400.), xtitle = MET, ytitle = MT),
    HDef('Mass2', 'M_{l#gamma}', (60, 60., 120.), xtitle = ('M_{l#gamma}', 'GeV'), overflowable = True, logscale = True),
    HDef('Mass2Wide', 'M_{l#gamma}', massBinning, xtitle = ('M_{l#gamma}', 'GeV'), overflowable = True, logscale = True),
    HDef('Mass3', 'M_{ll#gamma}', (30, 60., 120.), xtitle = ('M_{ll#gamma}', 'GeV'), overflowable = True, logscale = True),
    HDef('Mass3Wide', 'M_{l#gamma}', massBinning, xtitle = ('M_{ll#gamma}', 'GeV'), overflowable = True, logscale = True),
    HDef('PhotonPt', 'Photon P_{T}', (100, 40., 440.), xtitle = PT, overflowable = True, mask = (80., 'inf'), logscale = True),
    HDef('PhotonPtHighMet', 'Photon P_{T} (#slash{E}_{T} > 120 GeV)', (100, 40., 400.), xtitle = PT, overflowable = True, mask = ('-inf', 'inf'), logscale = True),
    HDef('PhotonPtLowMet', 'Photon P_{T} (#slash{E}_{T} < 70 GeV)', (100, 40., 400.), xtitle = PT, overflowable = True, logscale = True),
    HDef('PhotonEta', 'Photon #eta', (60, -1.5, 1.5), xtitle = '#eta'),
    HDef('LeptonPt', 'Lepton P_{T}', (100, 25., 225.), xtitle = PT, mask = (120., 'inf'), logscale = True),
    HDef('LeptonPtHighMet', 'Lepton P_{T} (#slash{E}_{T} > 120 GeV)', (100, 25., 225.), xtitle = PT, mask = ('-inf', 'inf'), logscale = True),
    HDef('LeptonPtLowMet', 'Lepton P_{T} (#slash{E}_{T} < 70 GeV)', (100, 25., 225.), xtitle = PT, logscale = True),
    HDef('LeptonEta', 'Lepton #eta', (60, -2.5, 2.5), xtitle = '#eta'),
    HDef('DEtaPhotonLepton', '#Delta#eta(#gamma, l)', (40, -4., 4.), xtitle = '#Delta#eta'),
    HDef('DPhiPhotonLepton', '#Delta#phi(#gamma, l)', dPhiBinning, xtitle = DPHI),
    HDef('DPhiPhotonLeptonLowMet', '#Delta#phi(#gamma, l) (#slash{E}_{T} < 70 GeV)', dPhiBinning, xtitle = DPHI),
    HDef('DRPhotonLepton', '#DeltaR(#gamma, l)', (60, 0., 6.), xtitle = DR),
    HDef('DRPhotonLeptonLowMet', '#DeltaR(#gamma, l) (#slash{E}_{T} < 70 GeV)', (60, 0., 6.), xtitle = DR),
    HDef('DRPhotonLeptonLowPhotonPt', '#DeltaR(#gamma, l) (P_{T}^{#gamma} < 80 GeV)', (60, 0., 6.), xtitle = DR),
    HDef('DRPhotonJet', '#DeltaR(#gamma, j)', (40, 0.5, 5.), xtitle = DR),
    HDef('DRLeptonJet', '#DeltaR(l, j)', (40, 0.5, 5.), xtitle = DR),
    HDef('DPhiLeptonJetLowMet', '#Delta#phi(l, j) (#slash{E}_{T} < 70 GeV)', dPhiBinning, xtitle = DPHI),
    HDef('DPhiPhotonMet', '#Delta#phi(#gamma, #slash{E}_{T})', dPhiBinning, xtitle = DPHI),
    HDef('DPhiPhotonMetLowMet', '#Delta#phi(#gamma, #slash{E}_{T}) (#slash{E}_{T} < 70 GeV)', dPhiBinning, xtitle = DPHI),
    HDef('DPhiLeptonMet', '#Delta#phi(l, #slash{E}_{T})', dPhiBinning, xtitle = DPHI),
    HDef('DPhiLeptonMetLowMet', '#Delta#phi(l, #slash{E}_{T}) (#slash{E}_{T} < 70 GeV)', dPhiBinning, xtitle = DPHI),
    HDef('NPhoton', 'N_{#gamma}', (3, 0.5, 3.5), xtitle = 'N_{#gamma}', logscale = True),
    HDef('NLepton', 'N_{l}', (4, 0.5, 4.5), xtitle = 'N_{l}', logscale = True),
    HDef('NPhotonNLepton', 'N_{#gamma} vs N_{l}', (4, 0.5, 4.5, 3, 0.5, 3.5), xtitle = 'N_{l}', ytitle = 'N_{#gamma}', drawOption = 'BOX TEXT'),
    HDef('NJet', 'N_{jet} (P_{T}^{j} > 20 GeV, |#eta^{j}| < 3)', (11, -0.5, 10.5), xtitle = 'N_{jet}', mask = (4.5, 'inf'), logscale = True),
    HDef('PW', 'P_{W}', (50, 0., 400.), xtitle = ('P', 'GeV'), overflowable = True),
    HDef('MW', 'M_{e#slash{E}_{T}}', (50, 80., 280.), xtitle = ('M', 'GeV'), overflowable = True, mask = ('-inf', 'inf')),
#    HDef('NLW', '-ln(#Sigma)', (50, 20., 70.), xtitle = '-ln(#Sigma)', overflowable = True),
#    HDef('NLWHighMet', '-ln(#Sigma) (#slash{E}_{T} > 120 GeV)', (50, 20., 70.), xtitle = '-ln(#Sigma)', overflowable = True)
]

hdefIndices = dict([(hdefs[index].name, index) for index in range(len(hdefs))])

hdef = hdefs[hdefIndices['NPhoton']]
hdef.xlabels = map(str, range(1, 4))

hdef = hdefs[hdefIndices['NLepton']]
hdef.xlabels = map(str, range(1, 5))

hdef = hdefs[hdefIndices['NPhotonNLepton']]
hdef.xlabels = map(str, range(1, 4))
hdef.ylabels = map(str, range(1, 5))

hdef = hdefs[hdefIndices['NJet']]
hdef.xlabels = map(str, range(1, 10))

iH = 0
try:
    with open(thisdir + '/plots.h', 'r') as header:
        header.readline()
        while iH != len(hdefs):
            line = header.readline().strip()
            if line != 'k' + hdefs[iH].name + ',':
                break

            iH += 1
        
except IOError:
    pass

if iH != len(hdefs):
    with open(thisdir + '/plots.h', 'w') as header:
        header.write('enum Plot {\n')
        for hdef in hdefs:
            header.write('  k' + hdef.name + ',\n')
        header.write('  nPlots\n')
        header.write('};\n')


############## PLOT MAKERS #############

ROOT.gROOT.LoadMacro(thisdir + '/GLPlotMaker.cc+')

plotMakerE = ROOT.GLPlotMaker(0)
plotMakerM = ROOT.GLPlotMaker(1)
plotMakerETrue = ROOT.GLPlotMaker(0)
plotMakerETrue.matchTruePhoton()
plotMakerETrue.matchTrueLepton()
plotMakerMTrue = ROOT.GLPlotMaker(1)
plotMakerMTrue.matchTruePhoton()
plotMakerMTrue.matchTrueLepton()
plotMakerEZ = ROOT.GLPlotMaker(0)
plotMakerEZ.useZ()
plotMakerEE = ROOT.GLPlotMaker(0)
plotMakerEE.matchElePhoton()
plotMakerEE.useZ()
plotMakerME = ROOT.GLPlotMaker(1)
plotMakerME.matchElePhoton()
plotMakerEJ = ROOT.GLPlotMaker(0)
plotMakerEJ.matchFakePhoton()
plotMakerEJ.useZ()
plotMakerMJ = ROOT.GLPlotMaker(1)
plotMakerMJ.matchFakePhoton()

############## STACK CONFIGURATIONS ##############
### EXPORTED ###
stackConfigs = {}

######## Electron Channel Search ########

gObservedE = Group('Observed', 'Observed', ROOT.kBlack, Group.OBSERVED, [datasets['DataE'].PhotonAndElectron])
gEGFakeE = Group('EGFake', 'e#rightarrow#gamma fake', ROOT.kOrange, Group.BACKGROUND, [datasets['DataE'].ElePhotonAndElectron])
gJGFakeE = Group('JGFake', 'j#rightarrow#gamma fake', ROOT.kGreen, Group.BACKGROUND, [datasets['DataE'].FakePhotonAndElectron])
gJLFakeE = Group('JLFake', 'QCD', ROOT.kRed, Group.BACKGROUND, [datasets['DataE'].PhotonAndFakeElectron])
gVGammaE = Group('VGamma', 'V#gamma (MC)', ROOT.kMagenta, Group.BACKGROUND, [
    (datasets['WGToLNuG_PtG-30-50'].PhotonAndElectron, 1.),
    (datasets['WGToLNuG_PtG-30-50'].ElePhotonAndElectron, -1.),
    (datasets['WGToLNuG_PtG-50-130'].PhotonAndElectron, 1.),
    (datasets['WGToLNuG_PtG-50-130'].ElePhotonAndElectron, -1.),
    (datasets['WGToLNuG_PtG-130'].PhotonAndElectron, 1.),
    (datasets['WGToLNuG_PtG-130'].ElePhotonAndElectron, -1.),
    (datasets['ZGToLLG'].PhotonAndElectron, 1.),
    (datasets['ZGToLLG'].ElePhotonAndElectron, -1.)
])
gEWKE = Group('EWK', 't#bar{t}#gamma/W^{+}W^{-}#gamma (MC)', ROOT.kCyan, Group.BACKGROUND, [
    (datasets['WWGJets'].PhotonAndElectron, 1.),
    (datasets['WWGJets'].ElePhotonAndElectron, -1.),
    (datasets['WW'].PhotonAndElectron, 1.),
    (datasets['WW'].ElePhotonAndElectron, -1.),
    (datasets['TTGJets'].PhotonAndElectron, 1.),
    (datasets['TTGJets'].ElePhotonAndElectron, -1.),
    (datasets['TTJetsSemiLept'].PhotonAndElectron, 1.),
    (datasets['TTJetsSemiLept'].ElePhotonAndElectron, -1.),
    (datasets['TTJetsFullLept'].PhotonAndElectron, 1.),
    (datasets['TTJetsFullLept'].ElePhotonAndElectron, -1.)
])

gT5wg_1000_425E = Group('T5wg_1000_425', 'T5wg_1000_425', ROOT.kGreen + 3, Group.SIGNAL, [
    (datasets['T5wg_1000_425'].PhotonAndElectron, 1.),
    (datasets['T5wg_1000_425'].ElePhotonAndElectron, -1.),
    (datasets['T5wg_1000_425'].FakePhotonAndElectron, -1.),
    (datasets['T5wg_1000_425'].PhotonAndFakeElectron, -1.)
])

stackConfigs['FloatingVGammaE'] = StackConfig(plotMakerE)
stackConfigs['FloatingVGammaE'].groups = [
    gObservedE,
    gEGFakeE,
    gJLFakeE,
    gJGFakeE,
    gVGammaE,
    gEWKE,
    gT5wg_1000_425E
]
stackConfigs['FloatingVGammaE'].specialPlotMakers = {
    gVGammaE: plotMakerETrue,
    gEWKE: plotMakerETrue
}
stackConfigs['FloatingVGammaE'].floatFixer = Chi2FitFixer([gJLFakeE, gVGammaE], [hdefs[hdefIndices['DRPhotonLeptonLowMet']]])

######## Muon Channel Search ########

gObservedM = Group('Observed', 'Observed', ROOT.kBlack, Group.OBSERVED, [datasets['DataM'].PhotonAndMuon])
gEGFakeM = Group('EGFake', 'e#rightarrow#gamma fake', ROOT.kOrange, Group.BACKGROUND, [datasets['DataM'].ElePhotonAndMuon])
gJGFakeM = Group('JGFake', 'j#rightarrow#gamma fake', ROOT.kGreen, Group.BACKGROUND, [datasets['DataM'].FakePhotonAndMuon])
gJLFakeM = Group('JLFake', 'QCD', ROOT.kRed, Group.BACKGROUND, [datasets['DataM'].PhotonAndFakeMuon])
gVGammaM = Group('VGamma', 'V#gamma (MC)', ROOT.kMagenta, Group.BACKGROUND, [
    (datasets['WGToLNuG_PtG-30-50'].PhotonAndMuon, 1.),
#    (datasets['WGToLNuG_PtG-30-50'].ElePhotonAndMuon, -1.),
    (datasets['WGToLNuG_PtG-50-130'].PhotonAndMuon, 1.),
#    (datasets['WGToLNuG_PtG-50-130'].ElePhotonAndMuon, -1.),
    (datasets['WGToLNuG_PtG-130'].PhotonAndMuon, 1.),
#    (datasets['WGToLNuG_PtG-130'].ElePhotonAndMuon, -1.),
    (datasets['ZGToLLG'].PhotonAndMuon, 1.),
#    (datasets['ZGToLLG'].ElePhotonAndMuon, -1.)
])
gEWKM = Group('EWK', 't#bar{t}#gamma/W^{+}W^{-}#gamma (MC)', ROOT.kCyan, Group.BACKGROUND, [
    (datasets['WWGJets'].PhotonAndMuon, 1.),
#    (datasets['WWGJets'].ElePhotonAndMuon, -1.),
    (datasets['WW'].PhotonAndMuon, 1.),
#    (datasets['WW'].ElePhotonAndMuon, -1.),
    (datasets['TTGJets'].PhotonAndMuon, 1.),
#    (datasets['TTGJets'].ElePhotonAndMuon, -1.),
    (datasets['TTJetsSemiLept'].PhotonAndMuon, 1.),
#    (datasets['TTJetsSemiLept'].ElePhotonAndMuon, -1.),
    (datasets['TTJetsFullLept'].PhotonAndMuon, 1.),
#    (datasets['TTJetsFullLept'].ElePhotonAndMuon, -1.)
])
gT5wg_1000_425M = Group('T5wg_1000_425', 'T5wg_1000_425', ROOT.kGreen + 3, Group.SIGNAL, [
    (datasets['T5wg_1000_425'].PhotonAndMuon, 1.),
#    (datasets['T5wg_1000_425'].ElePhotonAndMuon, -1.),
#    (datasets['T5wg_1000_425'].FakePhotonAndMuon, -1.),
#    (datasets['T5wg_1000_425'].PhotonAndFakeMuon, -1.)
])


stackConfigs['FloatingVGammaM'] = StackConfig(plotMakerM)
stackConfigs['FloatingVGammaM'].groups = [
    gObservedM,
    gEGFakeM,
    gJLFakeM,
    gJGFakeM,
    gVGammaM,
    gEWKM,
    gT5wg_1000_425M
]
stackConfigs['FloatingVGammaM'].specialPlotMakers = {
#    gVGammaM: plotMakerMTrue,
#    gEWKM: plotMakerMTrue
}
stackConfigs['FloatingVGammaM'].floatFixer = Chi2FitFixer([gJLFakeM, gVGammaM], [hdefs[hdefIndices['DPhiLeptonMetLowMet']]])

######### E->gamma Fake Closure (Electron Channel) ########

gDYEGFakeE = Group('DYEGFake', 'Fake (MC truth)', ROOT.kBlack, Group.OBSERVED, [datasets['DYJetsToLL'].PhotonAndElectron])
gEWKEGFakeE = Group('EWKEGFake', 'Fake (MC truth)', ROOT.kBlack, Group.OBSERVED, [
    datasets['TTJetsFullLeptNoGamma'].PhotonAndElectron,
    datasets['WWNoGamma'].PhotonAndElectron
])
gDYEleProxyE = Group('DYEleProxy', 'DY Proxy (MC)', ROOT.kOrange, Group.BACKGROUND, [datasets['DYJetsToLL'].ElePhotonAndElectron])
gEWKEleProxyE = Group('EWKEleProxy', 'EWK Proxy (MC)', ROOT.kBlack, Group.OBSERVED, [
    datasets['TTJetsFullLeptNoGamma'].ElePhotonAndElectron,
    datasets['WWNoGamma'].ElePhotonAndElectron
])


stackConfigs['EGFakeClosureE'] = StackConfig(plotMakerEZ)
stackConfigs['EGFakeClosureE'].groups = [
    gDYEGFakeE,
    gEWKEGFakeE,
    gDYEleProxyE,
    gEWKEleProxyE
]
stackConfigs['EGFakeClosureE'].specialPlotMakers = {
    gDYEGFakeE: plotMakerEE,
    gEWKEGFakeE: plotMakerEE
}

######### E->gamma Fake Closure (Muon Channel) ########

gDYEGFakeM = Group('DYEGFake', 'Fake (MC truth)', ROOT.kBlack, Group.OBSERVED, [datasets['DYJetsToLL'].PhotonAndMuon])
gEWKEGFakeM = Group('EWKEGFake', 'Fake (MC truth)', ROOT.kBlack, Group.OBSERVED, [
    datasets['TTJetsFullLeptNoGamma'].PhotonAndMuon,
    datasets['WWNoGamma'].PhotonAndMuon
])
gDYEleProxyM = Group('DYEleProxy', 'DY Proxy (MC)', ROOT.kOrange, Group.BACKGROUND, [datasets['DYJetsToLL'].ElePhotonAndMuon])
gEWKEleProxyM = Group('EWKEleProxy', 'EWK Proxy (MC)', ROOT.kBlack, Group.OBSERVED, [
    datasets['TTJetsFullLeptNoGamma'].ElePhotonAndMuon,
    datasets['WWNoGamma'].ElePhotonAndMuon
])


stackConfigs['EGFakeClosureM'] = StackConfig(plotMakerM)
stackConfigs['EGFakeClosureM'].groups = [
    gDYEGFakeM,
    gEWKEGFakeM,
    gDYEleProxyM,
    gEWKEleProxyM
]
stackConfigs['EGFakeClosureM'].specialPlotMakers = {
    gDYEGFakeM: plotMakerME,
    gEWKEGFakeM: plotMakerME
}

######### J->gamma Fake Closure (Electron Channel) ########

gWJJGFakeE = Group('WJJGFake', 'Fake (MC truth)', ROOT.kBlack, Group.OBSERVED, [
    datasets['WJetsToLNu'].PhotonAndElectron,
    datasets['WJetsToLNu_PtW-50To70'].PhotonAndElectron,
    datasets['WJetsToLNu_PtW-70To100'].PhotonAndElectron,
    datasets['WJetsToLNu_PtW-100'].PhotonAndElectron
])
gDYJGFakeE = Group('DYJGFake', 'Fake (MC truth)', ROOT.kBlack, Group.OBSERVED, [datasets['DYJetsToLL'].PhotonAndElectron])
gEWKJGFakeE = Group('EWKJGFake', 'Fake (MC truth)', ROOT.kBlack, Group.OBSERVED, [
    datasets['TTJetsSemiLeptNoGamma'].PhotonAndElectron,
    datasets['TTJetsFullLeptNoGamma'].PhotonAndElectron,
    datasets['WWNoGamma'].PhotonAndElectron    
])
gWJJetProxyE = Group('WJJetProxy', 'WJets Proxy (MC)', ROOT.kGreen, Group.BACKGROUND, [
    datasets['WJetsToLNu'].FakePhotonAndElectron,
    datasets['WJetsToLNu_PtW-50To70'].FakePhotonAndElectron,
    datasets['WJetsToLNu_PtW-70To100'].FakePhotonAndElectron,
    datasets['WJetsToLNu_PtW-100'].FakePhotonAndElectron
])
gDYJetProxyE = Group('DYJetProxy', 'DY Proxy (MC)', ROOT.kSpring + 9, Group.BACKGROUND, [datasets['DYJetsToLL'].FakePhotonAndElectron])
gEWKJetProxyE = Group('EWKJetProxy', 'EWK Proxy (MC)', ROOT.kGreen + 8, Group.BACKGROUND, [
    datasets['TTJetsSemiLeptNoGamma'].FakePhotonAndElectron,
    datasets['TTJetsFullLeptNoGamma'].FakePhotonAndElectron,
    datasets['WWNoGamma'].FakePhotonAndElectron    
])


stackConfigs['JGFakeClosureE'] = StackConfig(plotMakerE)
stackConfigs['JGFakeClosureE'].groups = [
    gWJJGFakeE,
    gDYJGFakeE,
    gEWKJGFakeE,
    gWJJetProxyE,
    gDYJetProxyE,
    gEWKJetProxyE
]
stackConfigs['JGFakeClosureE'].specialPlotMakers = {
    gWJJGFakeE: plotMakerEJ,
    gDYJGFakeE: plotMakerEJ,
    gEWKJGFakeE: plotMakerEJ,
}

######### J->gamma Fake Closure (Muon Channel) ########

gWJJGFakeM = Group('WJJGFake', 'Fake (MC truth)', ROOT.kBlack, Group.OBSERVED, [
    datasets['WJetsToLNu'].PhotonAndMuon,
    datasets['WJetsToLNu_PtW-50To70'].PhotonAndMuon,
    datasets['WJetsToLNu_PtW-70To100'].PhotonAndMuon,
    datasets['WJetsToLNu_PtW-100'].PhotonAndMuon
])
gDYJGFakeM = Group('DYJGFake', 'Fake (MC truth)', ROOT.kBlack, Group.OBSERVED, [datasets['DYJetsToLL'].PhotonAndMuon])
gEWKJGFakeM = Group('EWKJGFake', 'Fake (MC truth)', ROOT.kBlack, Group.OBSERVED, [
    datasets['TTJetsSemiLeptNoGamma'].PhotonAndMuon,
    datasets['TTJetsFullLeptNoGamma'].PhotonAndMuon,
    datasets['WWNoGamma'].PhotonAndMuon    
])
gWJJetProxyM = Group('WJJetProxy', 'WJets Proxy (MC)', ROOT.kGreen, Group.BACKGROUND, [
    datasets['WJetsToLNu'].FakePhotonAndMuon,
    datasets['WJetsToLNu_PtW-50To70'].FakePhotonAndMuon,
    datasets['WJetsToLNu_PtW-70To100'].FakePhotonAndMuon,
    datasets['WJetsToLNu_PtW-100'].FakePhotonAndMuon
])
gDYJetProxyM = Group('DYJetProxy', 'DY Proxy (MC)', ROOT.kSpring + 9, Group.BACKGROUND, [datasets['DYJetsToLL'].FakePhotonAndMuon])
gEWKJetProxyM = Group('EWKJetProxy', 'EWK Proxy (MC)', ROOT.kGreen + 8, Group.BACKGROUND, [
    datasets['TTJetsSemiLeptNoGamma'].FakePhotonAndMuon,
    datasets['TTJetsFullLeptNoGamma'].FakePhotonAndMuon,
    datasets['WWNoGamma'].FakePhotonAndMuon    
])


stackConfigs['JGFakeClosureM'] = StackConfig(plotMakerM)
stackConfigs['JGFakeClosureM'].groups = [
    gWJJGFakeM,
    gDYJGFakeM,
    gEWKJGFakeM,
    gWJJetProxyM,
    gDYJetProxyM,
    gEWKJetProxyM
]
stackConfigs['JGFakeClosureM'].specialPlotMakers = {
    gWJJGFakeM: plotMakerMJ,
    gDYJGFakeM: plotMakerMJ,
    gEWKJGFakeM: plotMakerMJ,
}
