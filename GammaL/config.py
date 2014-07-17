# Analysis-specific definitions
# Must define the following variables:
#  eventProcessor: event processor class
#  datasets: [Dataset]
#  weightCalc: {dataset name: {event class: weight calculator}}
#  stackConfigs: {stack name: StackConfig}

import os
import ROOT

from dataset import Dataset
from histogram import HDef
from stack import Group, StackConfig

ROOT.gSystem.Load('libRooFit.so')
ROOT.gSystem.Load('/afs/cern.ch/user/y/yiiyama/src/Common/fitting/libCommonFitting.so')

thisdir = os.path.dirname(os.path.abspath(__file__))

############## EVENT PROCESSOR ##############

ROOT.gROOT.LoadMacro(thisdir + '/GLSkimProcessor.cc+')

### EXPORTED ###
eventProcessor = ROOT.GLSkimProcessor


############## DATASET DEFINITIONS ##############

def realdata(name, inputNames, L, filters):
    return Dataset(name, inputNames, Dataset.REALDATA, L, 0., filters)

def fullsim(name, Leff, sigmaRelErr, filters, inputNames = None):
    if not inputNames:
        return Dataset(name, [name], Dataset.FULLSIM, Leff, sigmaRelErr, filters)
    else:
        return Dataset(name, inputNames, Dataset.FULLSIM, Leff, sigmaRelErr, filters)

def fastsim(name, Leff, sigmaRelErr, filters, inputNames = None):
    if not inputNames:
        return Dataset(name, [name], Dataset.FASTSIM, Leff, sigmaRelErr, filters)
    else:
        return Dataset(name, inputNames, Dataset.FASTSIM, Leff, sigmaRelErr, filters)

def addDataset(datasets, dataset):
    datasets[dataset.name] = dataset

data_E = ('PhotonAndElectron', 'ElePhotonAndElectron', 'FakePhotonAndElectron', 'PhotonAndFakeElectron')
data_Mu = ('PhotonAndMuon', 'ElePhotonAndMuon', 'FakePhotonAndMuon', 'PhotonAndFakeMuon', 'PhotonAndDimuon', 'ElePhotonAndDimuon', 'FakePhotonAndDimuon')
All = data_E + data_Mu
MC = ('PhotonAndElectron', 'PhotonAndMuon')
FakePhoton = ('PhotonAndElectron', 'PhotonAndMuon', 'FakePhotonAndElectron', 'FakePhotonAndMuon', 'PhotonAndDimuon', 'FakePhotonAndDimuon')
EleFakePhoton = ('PhotonAndElectron', 'PhotonAndMuon', 'ElePhotonAndElectron', 'ElePhotonAndMuon', 'FakePhotonAndElectron', 'FakePhotonAndMuon', 'PhotonAndDimuon', 'ElePhotonAndDimuon', 'FakePhotonAndDimuon')

### EXPORTED ###
datasets = {}
addDataset(datasets, realdata('DataE', ['PhotonA', 'DoublePhotonB', 'DoublePhotonC', 'DoublePhotonD'], 876.225 + 4411.704 + 7054.732 + 7369.007, data_E)) # 19712
addDataset(datasets, realdata('DataEA', ['PhotonA'], 876.225, data_E)) # 19712
addDataset(datasets, realdata('DataM', ['MuEGA', 'MuEGB', 'MuEGC', 'MuEGD'], 876.225 + 4411.704 + 7054.732 + 7360.046, data_Mu))
addDataset(datasets, fullsim('WGToLNuG_PtG-30-50', 3000000. / 75.48, 0., All)) # using parton-level cross section and total number of LHE events
addDataset(datasets, fullsim('WGToLNuG_PtG-50-130', 3000000. / 9.633, 0., All)) # same (parton-level simple sum * matching eff / matched events = parton-level simple sum / all events)
addDataset(datasets, fullsim('WGToLNuG_PtG-130', 471458. / 0.2571, 0., All)) # here PREP value is correctly set to parton-level xsec * matching eff
addDataset(datasets, fullsim('ZGToLLG_PtG-5-130', 6588161. / 132.6, 0., All))
addDataset(datasets, fullsim('ZGToLLG_PtG-130', 497474. / 0.0478, 0., All))
addDataset(datasets, fullsim('WWGJets', 304285. / 0.528, 0.5, All))
addDataset(datasets, fullsim('TTGJets', 1791552. / (245.8 * 0.0107), 0.5, All))
addDataset(datasets, fullsim('WW', 10000431. / 56., 0.5, All)) # NLO from twiki
addDataset(datasets, fullsim('TTJetsSemiLept', 24953451. / (245.8 * 0.324 * 0.676 * 2), 0.5, All)) # NNLO from twiki * BR(W->lnu BR) * BR(W->had) * 2. Sudakov = 0.999
addDataset(datasets, fullsim('TTJetsFullLept', 12011428. / (245.8 * 0.324 * 0.324), 0.5, All)) # NNLO from twiki * BR(W->lnu)^2. Sudakov = 0.999
addDataset(datasets, fullsim('WJetsToLNu', 57709905. / 36703.2, 0.15, FakePhoton)) # NNLO from twiki. Sudakov = 0.997
addDataset(datasets, fullsim('WJetsToLNu_PtW-50To70', 48426609. / (36703.2 * 811.2 / 30400.), 0.15, FakePhoton)) # 811.2 / 30400: MG5 xsec ratio
addDataset(datasets, fullsim('WJetsToLNu_PtW-70To100', 22447541. / (36703.2 * 428.9 / 30400.), 0.15, FakePhoton)) # 428.9 / 30400: MG5 xsec ratio
addDataset(datasets, fullsim('WJetsToLNu_PtW-100', 12742382. / (36703.2 * 228.9 / 30400.), 0.15, FakePhoton)) # 228.9 / 30400: MG5 xsec ratio
addDataset(datasets, fullsim('DYJetsToLL', 30459503. / 3531.9, 0.15, EleFakePhoton)) # NNLO from twiki. Sudakov = 0.993
#addDataset(datasets, fastsim('TChiwg_300', 54908. / 0.146, 0.00661 / 0., All))
#addDataset(datasets, fastsim('T5wg_500_425', 56217. / (4.525 * 0.5), 0., All)) # 0.5 from combinatorics (g/g->aW/aW)
addDataset(datasets, fastsim('T5wg_1000_425', 59387. / (0.0244 * 0.5), 0., All, inputNames = ['T5wg/skim_1000_425'])) # 0.5 from combinatorics (g/g->aW/aW)

datasets['TTJetsSemiLeptMC'] = fullsim('TTJetsSemiLeptMC', 24953451. / (245.8 * 0.324 * 0.676 * 2), 0.5, ('FakePhotonAndElectron', 'FakePhotonAndMuon'), inputNames = ['TTJetsSemiLept'])
datasets['TTJetsFullLeptMC'] = fullsim('TTJetsFullLeptMC', 12011428. / (245.8 * 0.324 * 0.324), 0.5, ('ElePhotonAndElectron', 'ElePhotonAndMuon', 'FakePhotonAndElectron', 'FakePhotonAndMuon'), inputNames = ['TTJetsFullLept'])
datasets['WWMC'] = fullsim('WWMC', 10000431. / 56., 0.5, ('ElePhotonAndElectron', 'ElePhotonAndMuon', 'FakePhotonAndElectron', 'FakePhotonAndMuon'), inputNames = ['WW'])

############## WEIGHT CALCULATORS ##############

eventWeights = {}
eventWeights['default'] = ROOT.GLEventWeight()
eventWeights['egHLT'] = ROOT.ElePhotonFunctionalWeight()
eventWeights['eg'] = ROOT.ElePhotonFunctionalWeight(True)
eventWeights['jgHLT'] = ROOT.JetPhotonHLTIsoWeight()
eventWeights['jg'] = ROOT.JetPhotonWeight()
eventWeights['mcegHLT'] = ROOT.MCElePhotonFunctionalWeight()
eventWeights['mceg'] = ROOT.MCElePhotonFunctionalWeight(True)
eventWeights['mcjgHLT'] = ROOT.MCJetPhotonHLTIsoWeight()
eventWeights['mcjg'] = ROOT.MCJetPhotonWeight()

### EXPORTED ###
weightCalc = {
    datasets['DataE'].PhotonAndElectron: eventWeights['default'],
    datasets['DataE'].ElePhotonAndElectron: eventWeights['egHLT'],
    datasets['DataE'].FakePhotonAndElectron: eventWeights['jgHLT'],
    datasets['DataE'].PhotonAndFakeElectron: eventWeights['default'],
    datasets['DataM'].PhotonAndMuon: eventWeights['default'],
    datasets['DataM'].ElePhotonAndMuon: eventWeights['eg'],
    datasets['DataM'].FakePhotonAndMuon: eventWeights['jg'],
    datasets['DataM'].PhotonAndFakeMuon: eventWeights['default'],
    datasets['DataM'].PhotonAndDimuon: eventWeights['default'],
    datasets['DataM'].ElePhotonAndDimuon: eventWeights['eg'],
    datasets['DataM'].FakePhotonAndDimuon: eventWeights['jg']
}

for name in ['WGToLNuG_PtG-30-50', 'WGToLNuG_PtG-50-130', 'WGToLNuG_PtG-130', 'ZGToLLG_PtG-5-130', 'ZGToLLG_PtG-130', 'WWGJets', 'TTGJets', 'WW', 'TTJetsSemiLept', 'TTJetsFullLept', 'T5wg_1000_425']:
    weightCalc[datasets[name].ElePhotonAndElectron] = eventWeights['egHLT']
    weightCalc[datasets[name].ElePhotonAndMuon] = eventWeights['eg']
    weightCalc[datasets[name].ElePhotonAndDimuon] = eventWeights['eg']
    weightCalc[datasets[name].FakePhotonAndElectron] = eventWeights['jgHLT']
    weightCalc[datasets[name].FakePhotonAndMuon] = eventWeights['jg']

weightCalcMC = {
    'PhotonAndElectron': eventWeights['default'],
    'PhotonAndMuon': eventWeights['default'],
    'PhotonAndFakeElectron': eventWeights['default'],
    'PhotonAndFakeMuon': eventWeights['default'],
    'ElePhotonAndElectron': eventWeights['mcegHLT'],
    'ElePhotonAndMuon': eventWeights['mceg'],
    'FakePhotonAndElectron': eventWeights['mcjgHLT'],
    'FakePhotonAndMuon': eventWeights['mcjg'],
    'PhotonAndDimuon': eventWeights['default'],
    'ElePhotonAndDimuon': eventWeights['mceg'],
    'FakePhotonAndDimuon': eventWeights['mcjg']
}

for s, c in weightCalcMC.items():
    weightCalc.update(dict((d.samples[s], c) for d in datasets.values() if d.type != Dataset.REALDATA and s in d.samples and d.samples[s] not in weightCalc))

############## HISTOGRAMS ###############

metBinning1 = [0., 10., 20., 25., 30., 34., 38., 42., 46., 50., 54., 58., 62., 70., 80., 90.,100., 120.,200.,300.,400.]
metBinning2 = [0., 4., 8., 12., 16., 20., 24., 28., 32., 36., 40., 44., 48., 52., 56., 60., 65., 70., 75., 80., 90.,100., 110.,120.,160.,200.,300.,400.]
mtBinning1 = [0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 140., 200., 300., 400.]
mtBinning2 = [0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 140., 160., 200., 240., 280., 400.]
htBinning1 = [0., 10.] + [30. + x * 10. for x in range(12)] + [150. + x * 20. for x in range(15)] + [450. + x * 50. for x in range(7)] + [800. + x * 100. for x in range(5)]
htBinning2 = [0., 10.] + [30. + x * 40. for x in range(8)] + [350., 400., 450.] + [500. + x * 100. for x in range(8)]
massBinning = [0.] + [x * 10. for x in range(1, 31)] + [300. + x * 20. for x in range(1, 11)] + [500. + x * 50. for x in range(1, 11)]
photonPtBinning = [40. + 4. * x for x in range(30)] + [160. + 8. * x for x in range(10)] + [240., 260., 300., 350., 400., 500.]
leptonPtBinning = [25. + 5. * x for x in range(15)] + [100. + 10. * x for x in range(6)] + [160., 180., 220., 260., 300., 400.]
dPhiBinning = (61, -3.15, 3.15)
dEtaBinning = (40, -4., 4.)
dRBinning = (60, 0., 6.)

MET = ('#slash{E}_{T}', 'GeV')
MT = ('M_{T}', 'GeV')
HT = ('H_{T}', 'GeV')
PT = ('P_{T}', 'GeV')
DR = '#DeltaR'
DPHI = '#Delta#phi'
DETA = '#Delta#eta'

### EXPORTED ###
hdefs = [
    HDef('MetHighMtLowPhotonPt', 'MET (M_{T} > 100 GeV, 40 GeV < P_{T}^{#gamma} < 80 GeV)', metBinning1, xtitle = MET, overflowable = True, mask = ('-inf', 'inf'), logscale = True, vrange = (1.e-3, 1.e+2)),
    HDef('MetHighMtHighPhotonPt', 'MET (M_{T} > 100 GeV, P_{T}^{#gamma} > 80 GeV)', metBinning1, xtitle = MET, overflowable = True, mask = ('-inf', 'inf'), logscale = True, vrange = (1.e-3, 1.e+2)),
    HDef('MetLowMtLowPhotonPt', 'MET (M_{T} < 100 GeV, P_{T}^{#gamma} < 80 GeV)', metBinning1, xtitle = MET, overflowable = True, mask = (120., 'inf'), logscale = True),
    HDef('Met', 'MET', metBinning2, xtitle = MET, overflowable = True, mask = (70., 'inf'), logscale = True, vrange = (1.e-2, 5.e+3)),
    HDef('MtHighMetHighPhotonPt', 'M_{T} (#slash{E}_{T} > 120 GeV, P_{T}^{#gamma} > 80 GeV)', mtBinning1, xtitle = MT, overflowable = True, mask = ('-inf', 'inf'), logscale = True),
    HDef('MtLowMetLowPhotonPt', 'M_{T} (#slash{E}_{T} < 70 GeV, P_{T}^{#gamma} < 80 GeV)', mtBinning2, xtitle = MT, overflowable = True, logscale = True, vrange = (1.e-3, 1.e+3)),
    HDef('Mt', 'M_{T}', mtBinning2, xtitle = MT, overflowable = True, mask = (100., 'inf'), logscale = True),
    HDef('MetMt', 'MET v M_{T}', (50, 0., 400, 50, 0., 400.), xtitle = MET, ytitle = MT),
    HDef('Ht', 'H_{T}', htBinning1, xtitle = HT, overflowable = True, mask = (600., 'inf'), logscale = True),
    HDef('HtHighMtHighMetHighPhotonPt', 'H_{T} (M_{T} > 100 GeV, #slash{E}_{T} > 120 GeV, P_{T}^{#gamma} > 80 GeV)', htBinning2, xtitle = HT, overflowable = True, mask = ('-inf', 'inf'), logscale = True),
    HDef('HtLowMetLowPhotonPt', 'H_{T} (#slash{E}_{T} < 70 GeV, P_{T}^{#gamma} < 80 GeV)', htBinning2, xtitle = HT, overflowable = True, logscale = True),
    HDef('Mass2', 'M_{l#gamma}', (60, 60., 120.), xtitle = ('M_{l#gamma}', 'GeV')),
    HDef('Mass2Wide', 'M_{l#gamma}', massBinning, xtitle = ('M_{l#gamma}', 'GeV'), overflowable = True, logscale = True),
    HDef('Mass3', 'M_{ll#gamma}', (30, 60., 120.), xtitle = ('M_{ll#gamma}', 'GeV')),
    HDef('Mass3Wide', 'M_{l#gamma}', massBinning, xtitle = ('M_{ll#gamma}', 'GeV'), overflowable = True, logscale = True),
    HDef('Mll', 'M_{ll}', (60, 60., 120.), xtitle = ('M_{ll}', 'GeV'), overflowable = True),
    HDef('PhotonPt', 'Photon P_{T}', photonPtBinning, xtitle = PT, overflowable = True, mask = (80., 'inf'), logscale = True),
    HDef('PhotonPtZoomLowMet', 'Photon P_{T} (#slash{E}_{T} < 70 GeV)', (50, 40., 140.), xtitle = PT, overflowable = True, logscale = True),
    HDef('PhotonPtHighMet', 'Photon P_{T} (#slash{E}_{T} > 120 GeV)', photonPtBinning, xtitle = PT, overflowable = True, mask = ('-inf', 'inf'), logscale = True),
    HDef('PhotonPtLowMet', 'Photon P_{T} (#slash{E}_{T} < 70 GeV)', photonPtBinning, xtitle = PT, overflowable = True, logscale = True, vrange = (5.e-2, 1.e+4)),
    HDef('PhotonPtLowMetOnZ', 'Photon P_{T} (#slash{E}_{T} < 70 GeV, 81 GeV < M_{ll} < 101 GeV)', photonPtBinning, xtitle = PT, overflowable = True, logscale = True),
    HDef('PhotonEta', 'Photon #eta', (60, -1.5, 1.5), xtitle = '#eta'),
    HDef('LeptonPt', 'Lepton P_{T}', leptonPtBinning, xtitle = PT, overflowable = True, mask = (120., 'inf'), logscale = True),
    HDef('LeptonPtHighMet', 'Lepton P_{T} (#slash{E}_{T} > 120 GeV)', leptonPtBinning, xtitle = PT, overflowable = True, mask = ('-inf', 'inf'), logscale = True),
    HDef('LeptonPtLowMet', 'Lepton P_{T} (#slash{E}_{T} < 70 GeV)', leptonPtBinning, xtitle = PT, overflowable = True, logscale = True, vrange = (5.e-1, 3.e+3)),
    HDef('LeptonEta', 'Lepton #eta', (60, -2.5, 2.5), xtitle = '#eta'),
    HDef('DEtaPhotonLepton', '#Delta#eta(#gamma, l)', dEtaBinning, xtitle = DETA),
    HDef('DPhiPhotonLepton', '#Delta#phi(#gamma, l)', dPhiBinning, xtitle = DPHI),
    HDef('DPhiPhotonLeptonOnZ', '#Delta#phi(#gamma, l) (81 GeV < M_{ll} < 101 GeV)', dPhiBinning, xtitle = DPHI),
    HDef('DRPhotonLepton', '#DeltaR(#gamma, l)', dRBinning, xtitle = DR),
    HDef('DRPhotonLeptonLowMet1L', '#DeltaR(#gamma, l) (#slash{E}_{T} < 70 GeV, N_{l} = 1)', dRBinning, xtitle = DR),
    HDef('DEtaPhotonJet', '#Delta#eta(#gamma, j)', dEtaBinning, xtitle = DETA),
    HDef('DPhiPhotonJet', '#Delta#phi(#gamma, j)', dPhiBinning, xtitle = DPHI),
    HDef('DRPhotonJet', '#DeltaR(#gamma, j)', dRBinning, xtitle = DR),
    HDef('DEtaLeptonJet', '#Delta#eta(l, j)', dEtaBinning, xtitle = DETA),
    HDef('DPhiLeptonJet', '#Delta#phi(l, j)', dPhiBinning, xtitle = DPHI),
    HDef('DRLeptonJet', '#DeltaR(l, j)', dRBinning, xtitle = DR),
    HDef('DPhiPhotonMet', '#Delta#phi(#gamma, #slash{E}_{T})', dPhiBinning, xtitle = DPHI),
    HDef('DPhiLeptonMet', '#Delta#phi(l, #slash{E}_{T})', dPhiBinning, xtitle = DPHI),
    HDef('DPhiLeptonMetLowMet1L', '#Delta#phi(l, #slash{E}_{T}) (#slash{E}_{T} < 70 GeV, N_{l} = 1)', dPhiBinning, xtitle = DPHI),
    HDef('NPhoton', 'N_{#gamma}', (3, 0.5, 3.5), xtitle = 'N_{#gamma}', xlabels = [str(i) for i in range(1, 4)], logscale = True),
    HDef('NLepton', 'N_{l}', (4, 0.5, 4.5), xtitle = 'N_{l}', xlabels = [str(i) for i in range(1, 5)], logscale = True),
    HDef('NPhotonNLepton', 'N_{#gamma} vs N_{l}', (4, 0.5, 4.5, 3, 0.5, 3.5), xtitle = 'N_{l}', xlabels = [str(i) for i in range(1, 4)], ytitle = 'N_{#gamma}', ylabels = [str(i) for i in range(1, 5)], drawOption = 'BOX TEXT'),
    HDef('NJet', 'N_{jet} (P_{T}^{j} > 30 GeV, |#eta^{j}| < 3)', (11, -0.5, 10.5), xtitle = 'N_{jet}', xlabels = [str(i) for i in range(1, 10)], mask = (4.5, 'inf'), logscale = True),
    HDef('NJetLowMet', 'N_{jet} (P_{T}^{j} > 30 GeV, |#eta^{j}| < 3, #slash{E}_{T} < 70 GeV)', (11, -0.5, 10.5), xtitle = 'N_{jet}', xlabels = [str(i) for i in range(1, 10)], logscale = True),
    HDef('NVtx', 'N_{vtx}', (51, -0.5, 50.5), xtitle = 'N_{vtx}'),
    HDef('PW', 'P_{W}', (50, 0., 400.), xtitle = ('P', 'GeV'), overflowable = True),
    HDef('MW', 'M_{e#slash{E}_{T}}', (50, 80., 280.), xtitle = ('M', 'GeV'), overflowable = True, mask = ('-inf', 'inf')),
#    HDef('NLW', '-ln(#Sigma)', (50, 20., 70.), xtitle = '-ln(#Sigma)', overflowable = True),
#    HDef('NLWHighMet', '-ln(#Sigma) (#slash{E}_{T} > 120 GeV)', (50, 20., 70.), xtitle = '-ln(#Sigma)', overflowable = True)
]

hdefList = dict([(h.name, h) for h in hdefs])

############## PLOT MAKERS #############

ROOT.gROOT.LoadMacro(thisdir + '/GLPlotMaker.cc+')

############## STACK CONFIGURATIONS ##############
### EXPORTED ###
stackConfigs = {}

######## Search config with floating VGamma ########

from floatingVGamma import FloatingVGammaSearch

######## Electron Channel Search ########

gObservedE = Group('Observed', 'Observed', ROOT.kBlack, Group.OBSERVED, [datasets['DataE'].PhotonAndElectron])
gEGFakeE = Group('EGFake', 'e#rightarrow#gamma fake', ROOT.kOrange, Group.BACKGROUND, [datasets['DataE'].ElePhotonAndElectron])
gJGFakeE = Group('JGFake', 'j#rightarrow#gamma fake', ROOT.kGreen, Group.BACKGROUND, [datasets['DataE'].FakePhotonAndElectron])
gJLFakeE = Group('JLFake', 'QCD', ROOT.kRed, Group.BACKGROUND, [datasets['DataE'].PhotonAndFakeElectron])
gVGammaE = Group('VGamma', 'V#gamma (MC)', ROOT.kMagenta, Group.BACKGROUND, [
    datasets['WGToLNuG_PtG-30-50'].PhotonAndElectron,
    datasets['WGToLNuG_PtG-30-50'].ElePhotonAndElectron,
    datasets['WGToLNuG_PtG-30-50'].FakePhotonAndElectron,
    datasets['WGToLNuG_PtG-30-50'].PhotonAndFakeElectron,
    datasets['WGToLNuG_PtG-50-130'].PhotonAndElectron,
    datasets['WGToLNuG_PtG-50-130'].ElePhotonAndElectron,
    datasets['WGToLNuG_PtG-50-130'].FakePhotonAndElectron,
    datasets['WGToLNuG_PtG-50-130'].PhotonAndFakeElectron,
    datasets['WGToLNuG_PtG-130'].PhotonAndElectron,
    datasets['WGToLNuG_PtG-130'].ElePhotonAndElectron,
    datasets['WGToLNuG_PtG-130'].FakePhotonAndElectron,
    datasets['WGToLNuG_PtG-130'].PhotonAndFakeElectron,
    datasets['ZGToLLG_PtG-5-130'].PhotonAndElectron,
    datasets['ZGToLLG_PtG-5-130'].ElePhotonAndElectron,
    datasets['ZGToLLG_PtG-5-130'].FakePhotonAndElectron,
    datasets['ZGToLLG_PtG-5-130'].PhotonAndFakeElectron,
    datasets['ZGToLLG_PtG-130'].PhotonAndElectron,
    datasets['ZGToLLG_PtG-130'].ElePhotonAndElectron,
    datasets['ZGToLLG_PtG-130'].FakePhotonAndElectron,
    datasets['ZGToLLG_PtG-130'].PhotonAndFakeElectron
])
gEWKE = Group('EWK', 't#bar{t}#gamma/W^{+}W^{-}#gamma (MC)', ROOT.kCyan, Group.BACKGROUND, [
    datasets['WWGJets'].PhotonAndElectron,
    datasets['WWGJets'].ElePhotonAndElectron,
    datasets['WWGJets'].FakePhotonAndElectron,
    datasets['WWGJets'].PhotonAndFakeElectron,
    datasets['WW'].PhotonAndElectron,
    datasets['WW'].ElePhotonAndElectron,
    datasets['WW'].FakePhotonAndElectron,
    datasets['WW'].PhotonAndFakeElectron,
    datasets['TTGJets'].PhotonAndElectron,
    datasets['TTGJets'].ElePhotonAndElectron,
    datasets['TTGJets'].FakePhotonAndElectron,
    datasets['TTGJets'].PhotonAndFakeElectron,
    datasets['TTJetsSemiLept'].PhotonAndElectron,
    datasets['TTJetsSemiLept'].ElePhotonAndElectron,
    datasets['TTJetsSemiLept'].FakePhotonAndElectron,
    datasets['TTJetsSemiLept'].PhotonAndFakeElectron,
    datasets['TTJetsFullLept'].PhotonAndElectron,
    datasets['TTJetsFullLept'].ElePhotonAndElectron,
    datasets['TTJetsFullLept'].FakePhotonAndElectron,
    datasets['TTJetsFullLept'].PhotonAndFakeElectron
])
gT5wg_1000_425E = Group('T5wg_1000_425', 'T5wg_1000_425', ROOT.kGreen + 3, Group.SIGNAL, [
    datasets['T5wg_1000_425'].PhotonAndElectron,
    datasets['T5wg_1000_425'].ElePhotonAndElectron,
    datasets['T5wg_1000_425'].FakePhotonAndElectron,
    datasets['T5wg_1000_425'].PhotonAndFakeElectron
])

plotMakerE = ROOT.GLPlotMaker(0)

stackConfigs['FloatingVGammaE'] = FloatingVGammaSearch(plotMakerE, 'Electron', next(hdef for hdef in hdefs if hdef.name == 'DRPhotonLeptonLowMet1L'))
stackConfigs['FloatingVGammaE'].groups = [
    gObservedE,
    gEGFakeE,
    gJLFakeE,
    gJGFakeE,
    gVGammaE,
    gEWKE,
    gT5wg_1000_425E
]
stackConfigs['FloatingVGammaE'].hdefs = hdefs

######## Muon Channel Search ########

gObservedM = Group('Observed', 'Observed', ROOT.kBlack, Group.OBSERVED, [datasets['DataM'].PhotonAndMuon])
gEGFakeM = Group('EGFake', 'e#rightarrow#gamma fake', ROOT.kOrange, Group.BACKGROUND, [datasets['DataM'].ElePhotonAndMuon])
gJGFakeM = Group('JGFake', 'j#rightarrow#gamma fake', ROOT.kGreen, Group.BACKGROUND, [datasets['DataM'].FakePhotonAndMuon])
gJLFakeM = Group('JLFake', 'QCD', ROOT.kRed, Group.BACKGROUND, [datasets['DataM'].PhotonAndFakeMuon])
gVGammaM = Group('VGamma', 'V#gamma (MC)', ROOT.kMagenta, Group.BACKGROUND, [
    datasets['WGToLNuG_PtG-30-50'].PhotonAndMuon,
    datasets['WGToLNuG_PtG-30-50'].ElePhotonAndMuon,
    datasets['WGToLNuG_PtG-30-50'].FakePhotonAndMuon,
    datasets['WGToLNuG_PtG-30-50'].PhotonAndFakeMuon,
    datasets['WGToLNuG_PtG-50-130'].PhotonAndMuon,
    datasets['WGToLNuG_PtG-50-130'].ElePhotonAndMuon,
    datasets['WGToLNuG_PtG-50-130'].FakePhotonAndMuon,
    datasets['WGToLNuG_PtG-50-130'].PhotonAndFakeMuon,
    datasets['WGToLNuG_PtG-130'].PhotonAndMuon,
    datasets['WGToLNuG_PtG-130'].ElePhotonAndMuon,
    datasets['WGToLNuG_PtG-130'].FakePhotonAndMuon,
    datasets['WGToLNuG_PtG-130'].PhotonAndFakeMuon,
    datasets['ZGToLLG_PtG-5-130'].PhotonAndMuon,
    datasets['ZGToLLG_PtG-5-130'].ElePhotonAndMuon,
    datasets['ZGToLLG_PtG-5-130'].FakePhotonAndMuon,
    datasets['ZGToLLG_PtG-5-130'].PhotonAndFakeMuon,
    datasets['ZGToLLG_PtG-130'].PhotonAndMuon,
    datasets['ZGToLLG_PtG-130'].ElePhotonAndMuon,
    datasets['ZGToLLG_PtG-130'].FakePhotonAndMuon,
    datasets['ZGToLLG_PtG-130'].PhotonAndFakeMuon
])
gEWKM = Group('EWK', 't#bar{t}#gamma/W^{+}W^{-}#gamma (MC)', ROOT.kCyan, Group.BACKGROUND, [
    datasets['WWGJets'].PhotonAndMuon,
    datasets['WWGJets'].ElePhotonAndMuon,
    datasets['WWGJets'].FakePhotonAndMuon,
    datasets['WWGJets'].PhotonAndFakeMuon,
    datasets['WW'].PhotonAndMuon,
    datasets['WW'].ElePhotonAndMuon,
    datasets['WW'].FakePhotonAndMuon,
    datasets['WW'].PhotonAndFakeMuon,
    datasets['TTGJets'].PhotonAndMuon,
    datasets['TTGJets'].ElePhotonAndMuon,
    datasets['TTGJets'].FakePhotonAndMuon,
    datasets['TTGJets'].PhotonAndFakeMuon,
    datasets['TTJetsSemiLept'].PhotonAndMuon,
    datasets['TTJetsSemiLept'].ElePhotonAndMuon,
    datasets['TTJetsSemiLept'].FakePhotonAndMuon,
    datasets['TTJetsSemiLept'].PhotonAndFakeMuon,
    datasets['TTJetsFullLept'].PhotonAndMuon,
    datasets['TTJetsFullLept'].ElePhotonAndMuon,
    datasets['TTJetsFullLept'].FakePhotonAndMuon,
    datasets['TTJetsFullLept'].PhotonAndFakeMuon
])
gT5wg_1000_425M = Group('T5wg_1000_425', 'T5wg_1000_425', ROOT.kGreen + 3, Group.SIGNAL, [
    datasets['T5wg_1000_425'].PhotonAndMuon,
    datasets['T5wg_1000_425'].ElePhotonAndMuon,
    datasets['T5wg_1000_425'].FakePhotonAndMuon,
    datasets['T5wg_1000_425'].PhotonAndFakeMuon
])

plotMakerM = ROOT.GLPlotMaker(1)

stackConfigs['FloatingVGammaM'] = FloatingVGammaSearch(plotMakerM, 'Muon', next(hdef for hdef in hdefs if hdef.name == 'DPhiLeptonMetLowMet1L'))
stackConfigs['FloatingVGammaM'].groups = [
    gObservedM,
    gEGFakeM,
    gJLFakeM,
    gJGFakeM,
    gVGammaM,
    gEWKM,
    gT5wg_1000_425M
]
stackConfigs['FloatingVGammaM'].hdefs = hdefs

######### Post Simul-fit results ###########

from fixedScales import FixedScalesSearch

stackConfigs['FixedScalesE'] = FixedScalesSearch(stackConfigs['FloatingVGammaE'], (1.490, 0.103), 0.1697)
stackConfigs['FixedScalesM'] = FixedScalesSearch(stackConfigs['FloatingVGammaM'], (1.490, 0.103), 0.05418)

######### Dimuon ##########

plotMakerMM = ROOT.GLPlotMaker(1)
plotMakerMM.dileptonOnly()
plotMakerMMTrue = ROOT.GLPlotMaker(1)
plotMakerMMTrue.dileptonOnly()
plotMakerMMTrue.matchTruePhoton()
plotMakerMMTrue.matchTrueLepton()

class DimuonControl(StackConfig):
    def __init__(self, sourceName):
        StackConfig.__init__(self, plotMakerMM)
        self.sourceName = sourceName

    def scalePlots(self, outputDir):
        source = ROOT.TFile.Open(self.sourceName)
        gr = source.Get('TemplateFitError/VGamma')
        scale = gr.GetY()[0]
        errUp = gr.GetErrorYhigh(0)
        errDown = gr.GetErrorYlow(0)
        vgamma = next(g for g in self.groups if g.name == 'VGamma')
        for sample in vgamma.samples:
            if sample.eventClass == 'PhotonAndMuon':
                sign = 1.
            else:
                sign = -1.
                sample.counter.SetBinContent(3, 0.)
            
            for histogram in sample.histograms.values():
                histogram.hWeighted.Scale(sign * scale)
                histogram.hScaleUp.Scale(sign * (scale + errUp))
                histogram.hScaleDown.Scale(sign * (scale - errDown))

            count = sample.counter.GetBinContent(1)
            sample.counter.SetBinContent(1, count * scale)
            sample.counter.SetBinContent(2, count * (scale + max(errUp, errDown)))

        source.Close()


gVGammaMM = Group('VGamma', 'V#gamma (MC)', ROOT.kMagenta, Group.BACKGROUND, [
    datasets['ZGToLLG_PtG-5-130'].PhotonAndMuon,
    datasets['ZGToLLG_PtG-5-130'].ElePhotonAndMuon,
    datasets['ZGToLLG_PtG-5-130'].FakePhotonAndMuon,
    datasets['ZGToLLG_PtG-5-130'].PhotonAndFakeMuon,
    datasets['ZGToLLG_PtG-130'].PhotonAndMuon,
    datasets['ZGToLLG_PtG-130'].ElePhotonAndMuon,
    datasets['ZGToLLG_PtG-130'].FakePhotonAndMuon,
    datasets['ZGToLLG_PtG-130'].PhotonAndFakeMuon
])
gEWKMM = Group('EWK', 't#bar{t}#gamma/W^{+}W^{-}#gamma (MC)', ROOT.kCyan, Group.BACKGROUND, [
    datasets['WWGJets'].PhotonAndMuon,
    datasets['WW'].PhotonAndMuon,
    datasets['TTGJets'].PhotonAndMuon,
    datasets['TTJetsFullLept'].PhotonAndMuon
])

stackConfigs['DimuonControl'] = DimuonControl('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/FloatingVGammaM.root')
stackConfigs['DimuonControl'].groups = [
    gObservedM,
    gEGFakeM,
    gJGFakeM,
    gVGammaMM,
    gEWKMM
]
stackConfigs['DimuonControl'].hdefs = [
    HDef('Mass3LowMet', 'M_{#mu#mu#gamma} (#slash{E}_{T} < 70 GeV)', massBinning, xtitle = ('M_{#mu#mu#gamma}', 'GeV'), overflowable = True, logscale = True),
    HDef('MuonPtLowMet', 'P_{T}^{#mu} (#slash{E}_{T} < 70 GeV, N_{#mu} == 2)', metBinning1, xtitle = ('P_{T}', 'GeV'), overflowable = True, logscale = True)
]
stackConfigs['DimuonControl'].specialPlotMakers = {
    gVGammaM: plotMakerMMTrue,
    gEWKM:plotMakerMMTrue
}

######### E->gamma Fake Closure (Electron Channel) ########

gEGFakeE = Group('DYEGFake', 'Fake (MC truth)', ROOT.kBlack, Group.OBSERVED, [
    datasets['DYJetsToLL'].PhotonAndElectron,
    datasets['TTJetsFullLept'].PhotonAndElectron,
    datasets['WW'].PhotonAndElectron
])
gDYEleProxyE = Group('DYEleProxy', 'DY Proxy (MC)', ROOT.kOrange, Group.BACKGROUND, [datasets['DYJetsToLL'].ElePhotonAndElectron])
gEWKEleProxyE = Group('EWKEleProxy', 'EWK Proxy (MC)', ROOT.kOrange + 10, Group.BACKGROUND, [
    datasets['TTJetsFullLeptMC'].ElePhotonAndElectron,
    datasets['WWMC'].ElePhotonAndElectron
])

plotMakerEL = ROOT.GLPlotMaker(0)
plotMakerEL.noEffCorrection()
plotMakerEE = ROOT.GLPlotMaker(0)
plotMakerEE.matchElePhoton()
plotMakerEE.noEffCorrection()

stackConfigs['EGFakeClosureE'] = StackConfig(plotMakerEL)
stackConfigs['EGFakeClosureE'].groups = [
    gEGFakeE,
    gDYEleProxyE,
    gEWKEleProxyE
]
stackConfigs['EGFakeClosureE'].hdefs = [
    HDef('Met', 'MET', metBinning2, xtitle = MET, logscale = True),
    HDef('Mt', 'M_{T}', mtBinning2, xtitle = MT, logscale = True),
    HDef('PhotonPt', 'Photon P_{T}', [40. + x * 5. for x in range(8)] + [80., 90., 100., 120., 140., 200.], xtitle = PT, logscale = True),
    HDef('LeptonPt', 'Lepton P_{T}', [25. + x * 5. for x in range(11)] + [80., 90., 100., 140., 200.], xtitle = PT, logscale = True),
    hdefList['NVtx']
]
stackConfigs['EGFakeClosureE'].specialPlotMakers = {
    gEGFakeE: plotMakerEE
}

######### E->gamma Fake Closure (Muon Channel) ########

gEGFakeM = Group('DYEGFake', 'Fake (MC truth)', ROOT.kBlack, Group.OBSERVED, [
    datasets['DYJetsToLL'].PhotonAndMuon,
    datasets['TTJetsFullLept'].PhotonAndMuon,
    datasets['WW'].PhotonAndMuon
])
gDYEleProxyM = Group('DYEleProxy', 'DY Proxy (MC)', ROOT.kOrange, Group.BACKGROUND, [datasets['DYJetsToLL'].ElePhotonAndMuon])
gEWKEleProxyM = Group('EWKEleProxy', 'EWK Proxy (MC)', ROOT.kBlack, Group.OBSERVED, [
    datasets['TTJetsFullLeptMC'].ElePhotonAndMuon,
    datasets['WWMC'].ElePhotonAndMuon
])

plotMakerML = ROOT.GLPlotMaker(1)
plotMakerML.matchTrueLepton()
plotMakerML.noEffCorrection()
plotMakerME = ROOT.GLPlotMaker(1)
plotMakerME.matchTrueLepton()
plotMakerME.matchElePhoton()
plotMakerME.noEffCorrection()

stackConfigs['EGFakeClosureM'] = StackConfig(plotMakerM)
stackConfigs['EGFakeClosureM'].groups = [
    gEGFakeM,
    gDYEleProxyM,
    gEWKEleProxyM
]
stackConfigs['EGFakeClosureM'].specialPlotMakers = {
    gEGFakeM: plotMakerME
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
    datasets['TTJetsSemiLept'].PhotonAndElectron,
    datasets['TTJetsFullLept'].PhotonAndElectron,
    datasets['WW'].PhotonAndElectron    
])
gWJJetProxyE = Group('WJJetProxy', 'WJets Proxy (MC)', ROOT.kGreen, Group.BACKGROUND, [
    datasets['WJetsToLNu'].FakePhotonAndElectron,
    datasets['WJetsToLNu_PtW-50To70'].FakePhotonAndElectron,
    datasets['WJetsToLNu_PtW-70To100'].FakePhotonAndElectron,
    datasets['WJetsToLNu_PtW-100'].FakePhotonAndElectron
])
gDYJetProxyE = Group('DYJetProxy', 'DY Proxy (MC)', ROOT.kSpring + 9, Group.BACKGROUND, [datasets['DYJetsToLL'].FakePhotonAndElectron])
gEWKJetProxyE = Group('EWKJetProxy', 'EWK Proxy (MC)', ROOT.kGreen + 8, Group.BACKGROUND, [
    datasets['TTJetsSemiLept'].FakePhotonAndElectron,
    datasets['TTJetsFullLept'].FakePhotonAndElectron,
    datasets['WW'].FakePhotonAndElectron    
])

plotMakerEJ = ROOT.GLPlotMaker(0)
plotMakerEJ.matchTrueLepton()
plotMakerEJ.matchFakePhoton()

stackConfigs['JGFakeClosureE'] = StackConfig(plotMakerEL)
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
    datasets['TTJetsSemiLept'].PhotonAndMuon,
    datasets['TTJetsFullLept'].PhotonAndMuon,
    datasets['WW'].PhotonAndMuon    
])
gWJJetProxyM = Group('WJJetProxy', 'WJets Proxy (MC)', ROOT.kGreen, Group.BACKGROUND, [
    datasets['WJetsToLNu'].FakePhotonAndMuon,
    datasets['WJetsToLNu_PtW-50To70'].FakePhotonAndMuon,
    datasets['WJetsToLNu_PtW-70To100'].FakePhotonAndMuon,
    datasets['WJetsToLNu_PtW-100'].FakePhotonAndMuon
])
gDYJetProxyM = Group('DYJetProxy', 'DY Proxy (MC)', ROOT.kSpring + 9, Group.BACKGROUND, [datasets['DYJetsToLL'].FakePhotonAndMuon])
gEWKJetProxyM = Group('EWKJetProxy', 'EWK Proxy (MC)', ROOT.kGreen + 8, Group.BACKGROUND, [
    datasets['TTJetsSemiLept'].FakePhotonAndMuon,
    datasets['TTJetsFullLept'].FakePhotonAndMuon,
    datasets['WW'].FakePhotonAndMuon    
])

plotMakerMJ = ROOT.GLPlotMaker(1)
plotMakerMJ.matchTrueLepton()
plotMakerMJ.matchFakePhoton()

stackConfigs['JGFakeClosureM'] = StackConfig(plotMakerML)
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

######## Test ########

#from fixedVGamma import FixedVGammaSearch
#from floatingVGamma import FloatingWGammaSearch

gWGammaE = Group('WGamma', 'W#gamma (MC)', ROOT.kMagenta, Group.BACKGROUND, [
    datasets['WGToLNuG_PtG-30-50'].PhotonAndElectron,
    datasets['WGToLNuG_PtG-30-50'].ElePhotonAndElectron,
    datasets['WGToLNuG_PtG-30-50'].FakePhotonAndElectron,
    datasets['WGToLNuG_PtG-30-50'].PhotonAndFakeElectron,
    datasets['WGToLNuG_PtG-50-130'].PhotonAndElectron,
    datasets['WGToLNuG_PtG-50-130'].ElePhotonAndElectron,
    datasets['WGToLNuG_PtG-50-130'].FakePhotonAndElectron,
    datasets['WGToLNuG_PtG-50-130'].PhotonAndFakeElectron,
    datasets['WGToLNuG_PtG-130'].PhotonAndElectron,
    datasets['WGToLNuG_PtG-130'].ElePhotonAndElectron,
    datasets['WGToLNuG_PtG-130'].FakePhotonAndElectron,
    datasets['WGToLNuG_PtG-130'].PhotonAndFakeElectron,
])
gZGammaE = Group('ZGamma', 'Z#gamma (MC)', ROOT.kViolet + 10, Group.BACKGROUND, [
    datasets['ZGToLLG_PtG-5-130'].PhotonAndElectron,
    datasets['ZGToLLG_PtG-5-130'].ElePhotonAndElectron,
    datasets['ZGToLLG_PtG-5-130'].FakePhotonAndElectron,
    datasets['ZGToLLG_PtG-5-130'].PhotonAndFakeElectron,
    datasets['ZGToLLG_PtG-130'].PhotonAndElectron,
    datasets['ZGToLLG_PtG-130'].ElePhotonAndElectron,
    datasets['ZGToLLG_PtG-130'].FakePhotonAndElectron,
    datasets['ZGToLLG_PtG-130'].PhotonAndFakeElectron
])

stackConfigs['PlotVG'] = StackConfig(plotMakerE)
stackConfigs['PlotVG'].groups = [
    gWGammaE,
    gZGammaE
]
stackConfigs['PlotVG'].hdefs = [HDef('Met', 'MET', metBinning2, xtitle = MET, logscale = True)]

plotMakerETrue = ROOT.GLPlotMaker(0)
plotMakerETrue.matchTruePhoton()
plotMakerETrue.matchTrueLepton()

from floatingWGamma import FloatingWGammaSearch

stackConfigs['TestE'] = FloatingWGammaSearch(plotMakerE, 'Electron', next(hdef for hdef in hdefs if hdef.name == 'DRPhotonLeptonLowMet1L'))
stackConfigs['TestE'].groups = [
    gObservedE,
    gEGFakeE,
    gJLFakeE,
    gJGFakeE,
    gWGammaE,
    gZGammaE,
    gEWKE,
    gT5wg_1000_425E
]
stackConfigs['TestE'].hdefs = hdefs
stackConfigs['TestE'].specialPlotMakers = {
    gWGammaE: plotMakerETrue,
    gZGammaE: plotMakerETrue,
    gEWKE: plotMakerETrue
}

gWGammaM = Group('WGamma', 'W#gamma (MC)', ROOT.kPink + 10, Group.BACKGROUND, [
    datasets['WGToLNuG_PtG-30-50'].PhotonAndMuon,
    datasets['WGToLNuG_PtG-30-50'].ElePhotonAndMuon,
    datasets['WGToLNuG_PtG-30-50'].FakePhotonAndMuon,
    datasets['WGToLNuG_PtG-30-50'].PhotonAndFakeMuon,
    datasets['WGToLNuG_PtG-50-130'].PhotonAndMuon,
    datasets['WGToLNuG_PtG-50-130'].ElePhotonAndMuon,
    datasets['WGToLNuG_PtG-50-130'].FakePhotonAndMuon,
    datasets['WGToLNuG_PtG-50-130'].PhotonAndFakeMuon,
    datasets['WGToLNuG_PtG-130'].PhotonAndMuon,
    datasets['WGToLNuG_PtG-130'].ElePhotonAndMuon,
    datasets['WGToLNuG_PtG-130'].FakePhotonAndMuon,
    datasets['WGToLNuG_PtG-130'].PhotonAndFakeMuon,
])
gZGammaM = Group('ZGamma', 'Z#gamma (MC)', ROOT.kViolet + 10, Group.BACKGROUND, [
    datasets['ZGToLLG_PtG-5-130'].PhotonAndMuon,
    datasets['ZGToLLG_PtG-5-130'].ElePhotonAndMuon,
    datasets['ZGToLLG_PtG-5-130'].FakePhotonAndMuon,
    datasets['ZGToLLG_PtG-5-130'].PhotonAndFakeMuon,
    datasets['ZGToLLG_PtG-130'].PhotonAndMuon,
    datasets['ZGToLLG_PtG-130'].ElePhotonAndMuon,
    datasets['ZGToLLG_PtG-130'].FakePhotonAndMuon,
    datasets['ZGToLLG_PtG-130'].PhotonAndFakeMuon
])

plotMakerMTrue = ROOT.GLPlotMaker(1)
plotMakerMTrue.matchTruePhoton()
plotMakerMTrue.matchTrueLepton()

stackConfigs['TestM'] = FloatingWGammaSearch(plotMakerM, 'Muon', next(hdef for hdef in hdefs if hdef.name == 'DPhiLeptonMetLowMet1L'))
stackConfigs['TestM'].groups = [
    gObservedM,
    gEGFakeM,
    gJLFakeM,
    gJGFakeM,
    gWGammaM,
    gZGammaM,
    gEWKM,
    gT5wg_1000_425M
]
stackConfigs['TestM'].hdefs = hdefs
stackConfigs['TestM'].specialPlotMakers = {
    gWGammaM: plotMakerMTrue,
    gZGammaM: plotMakerMTrue,
    gEWKM: plotMakerMTrue
}
