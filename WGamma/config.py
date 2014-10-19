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

thisdir = os.path.dirname(os.path.abspath(__file__))

blindingPrescale = 20

############## EVENT PROCESSOR ##############

ROOT.gROOT.LoadMacro(thisdir + '/WGammaProcessor.cc+')

### EXPORTED ###
eventProcessor = ROOT.WGammaProcessor

############## DATASET DEFINITIONS ##############

def realdata(name, inputNames, L, filters):
    return Dataset(name, inputNames, Dataset.REALDATA, L / blindingPrescale, 0., filters)

def fullsim(name, Leff, sigmaRelErr, filters, inputNames = None):
    if not inputNames:
        return Dataset(name, [name], Dataset.FULLSIM, Leff / blindingPrescale, sigmaRelErr, filters)
    else:
        return Dataset(name, inputNames, Dataset.FULLSIM, Leff / blindingPrescale, sigmaRelErr, filters)

def fastsim(name, Leff, sigmaRelErr, filters, inputNames = None):
    if not inputNames:
        return Dataset(name, [name], Dataset.FASTSIM, Leff / blindingPrescale, sigmaRelErr, filters)
    else:
        return Dataset(name, inputNames, Dataset.FASTSIM, Leff / blindingPrescale, sigmaRelErr, filters)

def addDataset(datasets, dataset):
    datasets[dataset.name] = dataset

data_E = ('SoftPhotonAndElectron', 'SoftElePhotonAndElectron', 'SoftFakePhotonAndElectron', 'SoftPhotonAndFakeElectron')
data_Mu = ('SoftPhotonAndMuon', 'SoftElePhotonAndMuon', 'SoftFakePhotonAndMuon', 'SoftPhotonAndFakeMuon')
All = data_E + data_Mu
MC = ('SoftPhotonAndElectron', 'SoftPhotonAndMuon')
FakePhoton = ('SoftPhotonAndElectron', 'SoftPhotonAndMuon', 'SoftFakePhotonAndElectron', 'SoftFakePhotonAndMuon')
EleFakePhoton = ('SoftPhotonAndElectron', 'SoftPhotonAndMuon', 'SoftElePhotonAndElectron', 'SoftElePhotonAndMuon', 'SoftFakePhotonAndElectron', 'SoftFakePhotonAndMuon')
OnlyFakes = ('SoftElePhotonAndElectron', 'SoftElePhotonAndMuon', 'SoftFakePhotonAndElectron', 'SoftFakePhotonAndMuon', 'SoftPhotonAndFakeElectron', 'SoftPhotonAndFakeMuon')

### EXPORTED ###

datasets = {}
addDataset(datasets, realdata('DataE', ['SingleElectronA', 'SingleElectronB', 'SingleElectronC', 'SingleElectronD'], 876.225 + 4411.704 + 7054.732 + 7369.007, data_E)) # 19712
addDataset(datasets, realdata('DataM', ['SingleMuA', 'SingleMuB', 'SingleMuC', 'SingleMuD'], 876.225 + 4411.704 + 7054.732 + 7360.046, data_Mu))
addDataset(datasets, fullsim('WGToLNuG_PtG-30-50', 3000000. / 75.48, 0., MC)) # using parton-level cross section and total number of LHE events
addDataset(datasets, fullsim('WGToLNuG_PtG-50-130', 3000000. / 9.633, 0., MC)) # same (parton-level simple sum * matching eff / matched events = parton-level simple sum / all events)
addDataset(datasets, fullsim('WGToLNuG_PtG-130', 471458. / 0.2571, 0., MC)) # here PREP value is correctly set to parton-level xsec * matching eff
addDataset(datasets, fullsim('ZGToLLG_PtG-5-130', 6588161. / 132.6, 0., MC))
addDataset(datasets, fullsim('ZGToLLG_PtG-130', 497474. / 0.0478, 0., MC))
addDataset(datasets, fullsim('WWGJets', 304285. / 0.528, 0.5, MC))
addDataset(datasets, fullsim('TTGJets', 1719954. / (1.444 * 2.), 0.5, MC)) # factor 0.453765 = 453765 events in phase space of TOP-13-011 / 1M events studied
addDataset(datasets, fullsim('WWJetsTo2L2Nu', 10000431. / 56., 0.5, All)) # NLO from twiki
addDataset(datasets, fullsim('WZJetsTo2L2Q', 10000283. / 33.21 / 0.67, 0.5, All)) # NLO from twiki
addDataset(datasets, fullsim('WZJetsTo3LNu', 2133868. / 33.21 / 0.33, 0.5, All)) # NLO from twiki
addDataset(datasets, fullsim('TTJetsSemiLept', 24953451. / (245.8 * 0.324 * 0.676 * 2), 0.5, MC)) # NNLO from twiki * BR(W->lnu BR) * BR(W->had) * 2. Sudakov = 0.999
addDataset(datasets, fullsim('TTJetsFullLept', 12011428. / (245.8 * 0.324 * 0.324), 0.5, MC)) # NNLO from twiki * BR(W->lnu)^2. Sudakov = 0.999
addDataset(datasets, fullsim('WJetsToLNu', 57709905. / 36703.2, 0.15, FakePhoton)) # NNLO from twiki. Sudakov = 0.997
addDataset(datasets, fullsim('WJetsToLNu_PtW-50To70', 48426609. / (36703.2 * 811.2 / 30400.), 0.15, FakePhoton)) # 811.2 / 30400: MG5 xsec ratio
addDataset(datasets, fullsim('WJetsToLNu_PtW-70To100', 22447541. / (36703.2 * 428.9 / 30400.), 0.15, FakePhoton)) # 428.9 / 30400: MG5 xsec ratio
addDataset(datasets, fullsim('WJetsToLNu_PtW-100', 12742382. / (36703.2 * 228.9 / 30400.), 0.15, FakePhoton)) # 228.9 / 30400: MG5 xsec ratio
addDataset(datasets, fullsim('DYJetsToLL', 30459503. / 3531.9, 0.15, EleFakePhoton)) # NNLO from twiki. Sudakov = 0.993

datasets['DataE'].prescale = blindingPrescale
datasets['DataM'].prescale = blindingPrescale

############## WEIGHT CALCULATORS ##############

eventWeights = {}
eventWeights['default'] = ROOT.WGEventWeight()
eventWeights['egHLT'] = ROOT.ElePhotonFunctionalWeight()
eventWeights['eg'] = ROOT.ElePhotonFunctionalWeight(True)
eventWeights['jgHLT'] = ROOT.JetPhotonHLTIsoWeight()
eventWeights['jg'] = ROOT.JetPhotonWeight()
eventWeights['mcegHLT'] = ROOT.MCElePhotonFunctionalWeight()
eventWeights['mceg'] = ROOT.MCElePhotonFunctionalWeight(True)
eventWeights['mcjgHLT'] = ROOT.MCJetPhotonWeight(0)
eventWeights['mcjg'] = ROOT.MCJetPhotonWeight(1)

### EXPORTED ###
weightCalc = {
    datasets['DataE'].SoftPhotonAndElectron: eventWeights['default'],
    datasets['DataE'].SoftElePhotonAndElectron: eventWeights['egHLT'],
    datasets['DataE'].SoftFakePhotonAndElectron: eventWeights['jgHLT'],
    datasets['DataE'].SoftPhotonAndFakeElectron: eventWeights['default'],
    datasets['DataM'].SoftPhotonAndMuon: eventWeights['default'],
    datasets['DataM'].SoftElePhotonAndMuon: eventWeights['eg'],
    datasets['DataM'].SoftFakePhotonAndMuon: eventWeights['jg'],
    datasets['DataM'].SoftPhotonAndFakeMuon: eventWeights['default'],
}

weightCalcMC = {
    'SoftPhotonAndElectron': eventWeights['default'],
    'SoftPhotonAndMuon': eventWeights['default'],
    'SoftPhotonAndFakeElectron': eventWeights['default'],
    'SoftPhotonAndFakeMuon': eventWeights['default'],
    'SoftElePhotonAndElectron': eventWeights['mcegHLT'],
    'SoftElePhotonAndMuon': eventWeights['mceg'],
    'SoftFakePhotonAndElectron': eventWeights['mcjgHLT'],
    'SoftFakePhotonAndMuon': eventWeights['mcjg'],
}

for s, c in weightCalcMC.items():
    weightCalc.update(dict((d.samples[s], c) for d in datasets.values() if d.dataType != Dataset.REALDATA and s in d.samples and d.samples[s] not in weightCalc))

############## HISTOGRAMS ###############

metBinning1 = [0., 4., 8., 12., 16., 20., 24., 28., 32., 36., 40., 44., 48., 52., 56., 60., 65., 70., 75., 80., 90.,100., 110.,120.,130.,140.,150.,160.,180.,200.,250.,300.,400.]
metBinning2 = [0., 10., 20., 25., 30., 34., 38., 42., 46., 50., 54., 58., 62., 70., 80., 90.,100.,110.,120.,140.,160.,200.,250.,300.,400.]
metBinning3 = [0., 10., 20., 30., 40., 50., 60., 70., 80., 90.,100., 120.,200.,300.,400.]
mtBinning1 = [0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 140., 200., 300., 400.]
mtBinning2 = [0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 140., 160., 200., 240., 280., 400.]
htBinning1 = [0.] + [30. + x * 10. for x in range(12)] + [150. + x * 20. for x in range(15)] + [450. + x * 50. for x in range(7)] + [800. + x * 100. for x in range(5)]
htBinning2 = [0.] + [30. + x * 80. for x in range(4)] + [350., 450., 600., 1000.]
massBinning = [0.] + [x * 10. for x in range(1, 31)] + [300. + x * 20. for x in range(1, 11)] + [500. + x * 50. for x in range(1, 11)]
photonPtBinning = [40. + 4. * x for x in range(30)] + [160. + 8. * x for x in range(10)] + [240., 260., 300., 350., 400., 500.]
photonPtBinning2 = [40. + 10. * x for x in range(10)] + [140. + 20. * x for x in range(3)] + [200., 300.]
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
    HDef('MetHighMtHighHtHighPhotonPt', 'MET (M_{T} > 100 GeV, H_{T} > 400 GeV, P_{T}^{#gamma} > 80 GeV)', metBinning3, xtitle = MET, overflowable = True, mask = ('-inf', 'inf'), logscale = True, vrange = (1.e-3, 1.e+2)),
    HDef('MetHighMtMidHtHighPhotonPt', 'MET (M_{T} > 100 GeV, 100 GeV < H_{T} < 400 GeV, P_{T}^{#gamma} > 80 GeV)', metBinning2, xtitle = MET, overflowable = True, mask = ('-inf', 'inf'), logscale = True, vrange = (1.e-3, 1.e+2)),
    HDef('MetHighMtLowHtHighPhotonPt', 'MET (M_{T} > 100 GeV, H_{T} < 100 GeV, P_{T}^{#gamma} > 80 GeV)', metBinning2, xtitle = MET, overflowable = True, mask = ('-inf', 'inf'), logscale = True, vrange = (1.e-3, 1.e+2)),
    HDef('MetHighMtHighHtLowPhotonPt', 'MET (M_{T} > 100 GeV, H_{T} > 400 GeV, P_{T}^{#gamma} < 80 GeV)', metBinning3, xtitle = MET, overflowable = True, mask = ('-inf', 'inf'), logscale = True, vrange = (1.e-3, 1.e+2)),
    HDef('MetHighMtMidHtLowPhotonPt', 'MET (M_{T} > 100 GeV, 100 GeV < H_{T} < 400 GeV, P_{T}^{#gamma} < 80 GeV)', metBinning2, xtitle = MET, overflowable = True, mask = ('-inf', 'inf'), logscale = True, vrange = (1.e-3, 1.e+2)),
    HDef('MetHighMtLowHtLowPhotonPt', 'MET (M_{T} > 100 GeV, H_{T} < 100 GeV, P_{T}^{#gamma} < 80 GeV)', metBinning2, xtitle = MET, overflowable = True, mask = ('-inf', 'inf'), logscale = True, vrange = (1.e-3, 1.e+2)),
    HDef('MetDiLepHighHtHighPhotonPt', 'MET (N_{l} #geq 2, H_{T} > 400 GeV, P_{T}^{#gamma} > 80 GeV)', metBinning3, xtitle = MET, overflowable = True, mask = ('-inf', 'inf'), logscale = True, vrange = (1.e-3, 1.e+2)),
    HDef('MetDiLepMidHtHighPhotonPt', 'MET (N_{l} #geq 2, 100 GeV < H_{T} < 400 GeV, P_{T}^{#gamma} > 80 GeV)', metBinning2, xtitle = MET, overflowable = True, mask = ('-inf', 'inf'), logscale = True, vrange = (1.e-3, 1.e+2)),
    HDef('MetDiLepLowHtHighPhotonPt', 'MET (N_{l} #geq 2, H_{T} < 100 GeV, P_{T}^{#gamma} > 80 GeV)', metBinning2, xtitle = MET, overflowable = True, mask = ('-inf', 'inf'), logscale = True, vrange = (1.e-3, 1.e+2)),
    HDef('Met', 'MET', metBinning1, xtitle = MET, overflowable = True, mask = (70., 'inf'), logscale = True, vrange = (1.e-2, 5.e+3)),
    HDef('MetHighMt', 'MET (M_{T} > 100 GeV)', metBinning2, xtitle = MET, overflowable = True, mask = (70., 'inf'), logscale = True, vrange = (1.e-2, 5.e+3)),
    HDef('Mt', 'M_{T}', mtBinning2, xtitle = MT, overflowable = True, mask = (100., 'inf'), logscale = True),
    HDef('MetMt', 'MET v M_{T}', (50, 0., 400, 50, 0., 400.), xtitle = MET, ytitle = MT),
    HDef('Ht', 'H_{T}', htBinning1, xtitle = HT, overflowable = True, mask = (600., 'inf'), logscale = True),
    HDef('HtHighMetHighMt', 'H_{T} (#slash{E}_{T} > 120 GeV, M_{T} > 100 GeV)', htBinning2, xtitle = HT, overflowable = True, mask = (600., 'inf'), logscale = True),
    HDef('Mass2', 'M_{l#gamma}', (60, 60., 120.), xtitle = ('M_{l#gamma}', 'GeV')),
    HDef('Mass2Wide', 'M_{l#gamma}', massBinning, xtitle = ('M_{l#gamma}', 'GeV'), overflowable = True, logscale = True),
    HDef('Mass3', 'M_{ll#gamma}', (30, 60., 120.), xtitle = ('M_{ll#gamma}', 'GeV')),
    HDef('Mass3Wide', 'M_{ll#gamma}', massBinning, xtitle = ('M_{ll#gamma}', 'GeV'), overflowable = True, logscale = True),
    HDef('Mll', 'M_{ll}', (60, 60., 120.), xtitle = ('M_{ll}', 'GeV'), overflowable = True),
    HDef('PhotonPt', 'Photon P_{T}', photonPtBinning, xtitle = PT, overflowable = True, mask = (80., 'inf'), logscale = True),
    HDef('PhotonPtZoomLowMet', 'Photon P_{T} (#slash{E}_{T} < 70 GeV)', (50, 40., 140.), xtitle = PT, overflowable = True, logscale = True),
    HDef('PhotonPtHighMet', 'Photon P_{T} (#slash{E}_{T} > 120 GeV)', photonPtBinning, xtitle = PT, overflowable = True, mask = ('-inf', 'inf'), logscale = True),
    HDef('PhotonPtHighMetHighMt', 'Photon P_{T} (#slash{E}_{T} > 120 GeV, M_{T} > 100 GeV)', photonPtBinning2, xtitle = PT, overflowable = True, mask = ('-inf', 'inf'), logscale = True),
    HDef('PhotonPtLowMet', 'Photon P_{T} (#slash{E}_{T} < 70 GeV)', photonPtBinning, xtitle = PT, overflowable = True, logscale = True, vrange = (5.e-2, 1.e+4)),
    HDef('PhotonPtLowMetOnZ', 'Photon P_{T} (#slash{E}_{T} < 70 GeV, 81 GeV < M_{ll} < 101 GeV)', photonPtBinning, xtitle = PT, overflowable = True, logscale = True),
    HDef('PhotonEta', 'Photon #eta', (60, -1.5, 1.5), xtitle = '#eta'),
    HDef('LeptonPt', 'Lepton P_{T}', leptonPtBinning, xtitle = PT, overflowable = True, mask = (120., 'inf'), logscale = True),
    HDef('LeptonPtHighMet', 'Lepton P_{T} (#slash{E}_{T} > 120 GeV)', leptonPtBinning, xtitle = PT, overflowable = True, mask = ('-inf', 'inf'), logscale = True),
    HDef('LeptonPtLowMet', 'Lepton P_{T} (#slash{E}_{T} < 70 GeV)', leptonPtBinning, xtitle = PT, overflowable = True, logscale = True, vrange = (5.e-2, 3.e+3)),
    HDef('LeptonEta', 'Lepton #eta', (60, -2.5, 2.5), xtitle = '#eta'),
    HDef('DEtaPhotonLepton', '#Delta#eta(#gamma, l)', dEtaBinning, xtitle = DETA),
    HDef('DPhiPhotonLepton', '#Delta#phi(#gamma, l)', dPhiBinning, xtitle = DPHI),
    HDef('DPhiPhotonLeptonMet4070', '#Delta#phi(#gamma, l) (40 GeV < #slash{E}_{T} < 70 GeV)', dPhiBinning, xtitle = DPHI),
    HDef('DPhiPhotonLeptonOnZ', '#Delta#phi(#gamma, l) (81 GeV < M_{ll} < 101 GeV)', dPhiBinning, xtitle = DPHI),
    HDef('DRPhotonLepton', '#DeltaR(#gamma, l)', dRBinning, xtitle = DR),
    HDef('DRPhotonLeptonMet4070', '#DeltaR(#gamma, l) (40 GeV < #slash{E}_{T} < 70 GeV)', dRBinning, xtitle = DR),
    HDef('DRPhotonTrailLepton', '#DeltaR(#gamma, l_{2})', dRBinning, xtitle = DR),
    HDef('DEtaPhotonJet', '#Delta#eta(#gamma, j)', dEtaBinning, xtitle = DETA),
    HDef('DPhiPhotonJet', '#Delta#phi(#gamma, j)', dPhiBinning, xtitle = DPHI),
    HDef('DRPhotonJet', '#DeltaR(#gamma, j)', dRBinning, xtitle = DR),
    HDef('DEtaLeptonJet', '#Delta#eta(l, j)', dEtaBinning, xtitle = DETA),
    HDef('DPhiLeptonJet', '#Delta#phi(l, j)', dPhiBinning, xtitle = DPHI),
    HDef('DRLeptonJet', '#DeltaR(l, j)', dRBinning, xtitle = DR),
    HDef('DPhiPhotonMet', '#Delta#phi(#gamma, #slash{E}_{T})', dPhiBinning, xtitle = DPHI),
    HDef('DPhiPhotonMetMet4070', '#Delta#phi(#gamma, #slash{E}_{T}) (40 GeV < #slash{E}_{T} < 70 GeV)', dPhiBinning, xtitle = DPHI),
    HDef('DPhiLeptonMet', '#Delta#phi(l, #slash{E}_{T})', dPhiBinning, xtitle = DPHI),
    HDef('DPhiLeptonMetMet4070', '#Delta#phi(l, #slash{E}_{T}) (40 GeV < #slash{E}_{T} < 70 GeV)', dPhiBinning, xtitle = DPHI),
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

hdefsClosure = [
    HDef('Met', 'MET', metBinning1, xtitle = MET, logscale = True),
    HDef('Mt', 'M_{T}', mtBinning2, xtitle = MT, logscale = True),
    HDef('PhotonPt', 'Photon P_{T}', [40. + x * 5. for x in range(8)] + [80., 90., 100., 120., 140., 200.], xtitle = PT, logscale = True),
    HDef('LeptonPt', 'Lepton P_{T}', [25. + x * 5. for x in range(11)] + [80., 90., 100., 140., 200.], xtitle = PT, logscale = True),
    hdefList['NVtx']
]

############## PLOT MAKERS #############

#ROOT.gROOT.LoadMacro(thisdir + '/WGPlotMaker.cc+')

############## STACK CONFIGURATIONS ##############
### EXPORTED ###
stackConfigs = {}

######## Electron Channel Search ########

gObservedE = Group('Observed', 'Observed', ROOT.kBlack, Group.OBSERVED, [datasets['DataE'].SoftPhotonAndElectron])
gEGFakeE = Group('EGFake', 'e#rightarrow#gamma fake', ROOT.kOrange, Group.BACKGROUND, [datasets['DataE'].SoftElePhotonAndElectron])
gJGFakeE = Group('JGFake', 'j#rightarrow#gamma fake', ROOT.kGreen, Group.BACKGROUND, [datasets['DataE'].SoftFakePhotonAndElectron])
gJLFakeE = Group('JLFake', 'QCD', ROOT.kRed, Group.BACKGROUND, [datasets['DataE'].SoftPhotonAndFakeElectron])
gWGammaE = Group('WGamma', 'W#gamma (MC)', ROOT.kMagenta, Group.BACKGROUND, [
    datasets['WGToLNuG_PtG-30-50'].SoftPhotonAndElectron,
    datasets['WGToLNuG_PtG-50-130'].SoftPhotonAndElectron,
    datasets['WGToLNuG_PtG-130'].SoftPhotonAndElectron
])
gZGammaE = Group('ZGamma', 'Z#gamma (MC)', ROOT.kBlue - 10, Group.BACKGROUND, [
    datasets['ZGToLLG_PtG-5-130'].SoftPhotonAndElectron,
    datasets['ZGToLLG_PtG-130'].SoftPhotonAndElectron
])
gMultiBosonE = Group('MultiBoson', 'WW#gamma/WZ#gamma (MC)', ROOT.kCyan, Group.BACKGROUND, [
    datasets['WWGJets'].SoftPhotonAndElectron,
    datasets['WWJetsTo2L2Nu'].SoftPhotonAndElectron,
    datasets['WZJetsTo2L2Q'].SoftPhotonAndElectron,
    datasets['WZJetsTo3LNu'].SoftPhotonAndElectron
])
gTTE = Group('TT', 't#bar{t}#gamma (MC)', ROOT.kRed - 7, Group.BACKGROUND, [
    datasets['TTGJets'].SoftPhotonAndElectron,
    datasets['TTJetsSemiLept'].SoftPhotonAndElectron,
    datasets['TTJetsFullLept'].SoftPhotonAndElectron
])

######## Muon Channel Search ########

gObservedM = Group('Observed', 'Observed', ROOT.kBlack, Group.OBSERVED, [datasets['DataM'].SoftPhotonAndMuon])
gEGFakeM = Group('EGFake', 'e#rightarrow#gamma fake', ROOT.kOrange, Group.BACKGROUND, [datasets['DataM'].SoftElePhotonAndMuon])
gJGFakeM = Group('JGFake', 'j#rightarrow#gamma fake', ROOT.kGreen, Group.BACKGROUND, [datasets['DataM'].SoftFakePhotonAndMuon])
gJLFakeM = Group('JLFake', 'QCD', ROOT.kRed, Group.BACKGROUND, [datasets['DataM'].SoftPhotonAndFakeMuon])
gWGammaM = Group('WGamma', 'W#gamma (MC)', ROOT.kMagenta, Group.BACKGROUND, [
    datasets['WGToLNuG_PtG-30-50'].SoftPhotonAndMuon,
    datasets['WGToLNuG_PtG-50-130'].SoftPhotonAndMuon,
    datasets['WGToLNuG_PtG-130'].SoftPhotonAndMuon
])
gZGammaM = Group('ZGamma', 'Z#gamma (MC)', ROOT.kBlue - 10, Group.BACKGROUND, [
    datasets['ZGToLLG_PtG-5-130'].SoftPhotonAndMuon,
    datasets['ZGToLLG_PtG-130'].SoftPhotonAndMuon
])
gMultiBosonM = Group('MultiBoson', 'WW#gamma/WZ#gamma (MC)', ROOT.kCyan, Group.BACKGROUND, [
    datasets['WWGJets'].SoftPhotonAndMuon,
    datasets['WWJetsTo2L2Nu'].SoftPhotonAndMuon,
    datasets['WZJetsTo2L2Q'].SoftPhotonAndMuon,
    datasets['WZJetsTo3LNu'].SoftPhotonAndMuon
])
gTTM = Group('TT', 't#bar{t}#gamma (MC)', ROOT.kRed - 7, Group.BACKGROUND, [
    datasets['TTGJets'].SoftPhotonAndMuon,
    datasets['TTJetsSemiLept'].SoftPhotonAndMuon,
    datasets['TTJetsFullLept'].SoftPhotonAndMuon
])
