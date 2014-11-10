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

ROOT.gSystem.Load('/afs/cern.ch/user/y/yiiyama/src/Common/fitting/libCommonFitting.so')

thisdir = os.path.dirname(os.path.abspath(__file__))

############## EVENT PROCESSOR ##############

ROOT.gROOT.LoadMacro(thisdir + '/GLSkimProcessor.cc+')

### EXPORTED ###
eventProcessor = ROOT.GLSkimProcessor


############## DATASET DEFINITIONS ##############

def realdata(name, inputNames, L, filters):
    return Dataset(name, inputNames, Dataset.REALDATA, 1, 1. / L, 0., filters)

def fullsim(name, nEvents, sigma, sigmaRelErr, filters, inputNames = None):
    if not inputNames:
        dataset = Dataset(name, [name], Dataset.FULLSIM, nEvents, sigma, sigmaRelErr, filters)
    else:
        dataset = Dataset(name, inputNames, Dataset.FULLSIM, nEvents, sigma, sigmaRelErr, filters)

    dataset.entryList = 'hardPhotonList'
    return dataset

def fastsim(name, nEvents, sigma, sigmaRelErr, filters, inputNames = None):
    if not inputNames:
        return Dataset(name, [name], Dataset.FASTSIM, nEvents, sigma, sigmaRelErr, filters)
    else:
        return Dataset(name, inputNames, Dataset.FASTSIM, nEvents, sigma, sigmaRelErr, filters)

def addDataset(datasets, dataset):
    datasets[dataset.name] = dataset

data_E = ('PhotonAndElectron', 'ElePhotonAndElectron', 'FakePhotonAndElectron', 'PhotonAndFakeElectron', 'ElePhotonAndFakeElectron', 'FakePhotonAndFakeElectron')
data_Mu = ('PhotonAndMuon', 'ElePhotonAndMuon', 'FakePhotonAndMuon', 'PhotonAndFakeMuon', 'PhotonAndDimuon', 'ElePhotonAndDimuon', 'FakePhotonAndDimuon', 'ElePhotonAndFakeMuon', 'FakePhotonAndFakeMuon')
All = data_E + data_Mu
MC = ('PhotonAndElectron', 'PhotonAndMuon', 'PhotonAndDimuon')
FakePhoton = ('PhotonAndElectron', 'PhotonAndMuon', 'FakePhotonAndElectron', 'FakePhotonAndMuon', 'PhotonAndDimuon', 'FakePhotonAndDimuon')
FakeLepton = ('PhotonAndElectron', 'PhotonAndMuon', 'PhotonAndFakeElectron', 'PhotonAndFakeMuon')
EleFakePhoton = ('PhotonAndElectron', 'PhotonAndMuon', 'ElePhotonAndElectron', 'ElePhotonAndMuon', 'FakePhotonAndElectron', 'FakePhotonAndMuon', 'PhotonAndDimuon', 'ElePhotonAndDimuon', 'FakePhotonAndDimuon')
OnlyFakes = ('ElePhotonAndElectron', 'ElePhotonAndMuon', 'FakePhotonAndElectron', 'FakePhotonAndMuon', 'PhotonAndFakeElectron', 'PhotonAndFakeMuon', 'FakePhotonAndDimuon')

### EXPORTED ###
datasets = {}
addDataset(datasets, realdata('DataE', ['PhotonA', 'DoublePhotonB', 'DoublePhotonC', 'DoublePhotonD'], 876.225 + 4411.704 + 7054.732 + 7369.007, data_E)) # 19712
addDataset(datasets, realdata('DataEA', ['PhotonA'], 876.225, data_E)) # 19712
addDataset(datasets, realdata('DataM', ['MuEGA', 'MuEGB', 'MuEGC', 'MuEGD'], 876.225 + 4411.704 + 7054.732 + 7360.046, data_Mu))
addDataset(datasets, realdata('DataMA', ['MuEGA'], 876.225, data_Mu))
addDataset(datasets, fullsim('WGToLNuG_PtG-30-50', 869591, 75.48 * 869591. / 3000000., 0., MC)) # using parton-level cross section and total number of LHE events
addDataset(datasets, fullsim('WGToLNuG_PtG-50-130', 1135698, 9.633 * 1135698. / 3000000., 0., MC)) # same (parton-level simple sum * matching eff / matched events = parton-level simple sum / all events)
addDataset(datasets, fullsim('WGToLNuG_PtG-130', 471458, 0.2571, 0., MC)) # here PREP value is correctly set to parton-level xsec * matching eff
addDataset(datasets, fullsim('ZGToLLG_PtG-5-130', 6588161, 132.6, 0., MC))
addDataset(datasets, fullsim('ZGToLLG_PtG-130', 497474, 0.0478, 0., MC))
addDataset(datasets, fullsim('WWGJets', 304285, 0.528, 0.5, MC))
addDataset(datasets, fullsim('TTGJets', 1719954. 1.444 * 2., 0.5, MC)) # factor 0.453765 = 453765 events in phase space of TOP-13-011, 1M events studied
addDataset(datasets, fullsim('ttA', 1074860, 228.4 * 0.0107, 0.5, MC))
addDataset(datasets, fullsim('WW', 10000431, 56., 0.5, All)) # NLO from twiki
addDataset(datasets, fullsim('WZ', 10000283, 33.21, 0.5, All)) # NLO from twiki
addDataset(datasets, fullsim('TTJetsSemiLept', 24953451, 245.8 * 0.324 * 0.676 * 2, 0.5, MC)) # NNLO from twiki * BR(W->lnu BR) * BR(W->had) * 2. Sudakov = 0.999
addDataset(datasets, fullsim('TTJetsSemiLeptFake', 24953451, 245.8 * 0.324 * 0.676 * 2, 0.5, OnlyFakes, inputNames = ['TTJetsSemiLept']))
addDataset(datasets, fullsim('TTJetsFullLept', 12011428, 245.8 * 0.324 * 0.324, 0.5, MC)) # NNLO from twiki * BR(W->lnu)^2. Sudakov = 0.999
addDataset(datasets, fullsim('TTJetsFullLeptFake', 12011428, 245.8 * 0.324 * 0.324, 0.5, OnlyFakes, inputNames = ['TTJetsFullLept']))
addDataset(datasets, fullsim('WJetsToLNu', 57709905, 36703.2, 0.15, FakePhoton)) # NNLO from twiki. Sudakov = 0.997
addDataset(datasets, fullsim('WJetsToLNu_PtW-50To70', 48426609, 36703.2 * 811.2 / 30400., 0.15, FakePhoton)) # 811.2 / 30400: MG5 xsec ratio
addDataset(datasets, fullsim('WJetsToLNu_PtW-70To100', 22447541, 36703.2 * 428.9 / 30400., 0.15, FakePhoton)) # 428.9 / 30400: MG5 xsec ratio
addDataset(datasets, fullsim('WJetsToLNu_PtW-100', 12742382, 36703.2 * 228.9 / 30400., 0.15, FakePhoton)) # 228.9 / 30400: MG5 xsec ratio
addDataset(datasets, fullsim('DYJetsToLL', 30459503, 3531.9, 0.15, EleFakePhoton)) # NNLO from twiki. Sudakov = 0.993
addDataset(datasets, fullsim('GJets_HT-40To100', 19857930, 20930.0, 0.15, FakeLepton))
addDataset(datasets, fullsim('GJets_HT-100To200', 9612703, 5212.0, 0.15, FakeLepton))
addDataset(datasets, fullsim('GJets_HT-200To400', 10494617, 960.5, 0.15, FakeLepton))
addDataset(datasets, fullsim('GJets_HT-400ToInf', 9539562, 107.5, 0.15, FakeLepton))
#addDataset(datasets, fastsim('T5wg_500_425', 56217, 4.525 * 0.5, 0., All)) # 0.5 from combinatorics (g/g->aW/aW)
addDataset(datasets, fastsim('T5wg_1000_425', 59387, 0.0244 * 0.5, 0., All, inputNames = ['T5wg/skim_1000_425'])) # 0.5 from combinatorics (g/g->aW/aW)
addDataset(datasets, fastsim('TChiwg_300', 55240, 0.146, 0.045, All, inputNames = ['TChiwg/skim_300']))
addDataset(datasets, fastsim('Spectra_gW_M3_1015_M2_455_gg', 250000, 0.0176, 0.00093, All, inputNames = ['Spectra_gW/skim_M3_1015_M2_455_gg']))
addDataset(datasets, fastsim('Spectra_gW_M3_1515_M2_305_nc', 125000, 0.0913 + 0.0372, 0.00062, All, inputNames = ['Spectra_gW/skim_M3_1515_M2_305_ncp', 'Spectra_gW/skim_M3_1515_M2_305_ncm']))

datasets['TTJetsFullLeptFake'].prescale = 100
datasets['TTJetsSemiLeptFake'].prescale = 100

############## WEIGHT CALCULATORS ##############

eventWeights = {}
eventWeights['default'] = ROOT.GLEventWeight()
eventWeights['egHLT'] = ROOT.ElePhotonFunctionalWeight()
eventWeights['eg'] = ROOT.ElePhotonFunctionalWeight(True)
eventWeights['jgHLT'] = ROOT.JetPhotonHLTIsoWeight()
eventWeights['jg'] = ROOT.JetPhotonWeight()
eventWeights['mcegHLT'] = ROOT.MCElePhotonFunctionalWeight()
eventWeights['mceg'] = ROOT.MCElePhotonFunctionalWeight(True)
eventWeights['mcjgHLT'] = ROOT.MCJetPhotonWeight(0)
eventWeights['mcjg'] = ROOT.MCJetPhotonWeight(1)

#eventWeights['geniso'] = ROOT.PhotonGenIsoWeight()

### EXPORTED ###
weightCalc = {
    datasets['DataE'].PhotonAndElectron: eventWeights['default'],
    datasets['DataE'].ElePhotonAndElectron: eventWeights['egHLT'],
    datasets['DataE'].FakePhotonAndElectron: eventWeights['jgHLT'],
    datasets['DataE'].PhotonAndFakeElectron: eventWeights['default'],
    datasets['DataE'].ElePhotonAndFakeElectron: eventWeights['egHLT'],
    datasets['DataE'].FakePhotonAndFakeElectron: eventWeights['jgHLT'],
    datasets['DataM'].PhotonAndMuon: eventWeights['default'],
    datasets['DataM'].ElePhotonAndMuon: eventWeights['eg'],
    datasets['DataM'].FakePhotonAndMuon: eventWeights['jg'],
    datasets['DataM'].PhotonAndFakeMuon: eventWeights['default'],
    datasets['DataM'].ElePhotonAndFakeMuon: eventWeights['eg'],
    datasets['DataM'].FakePhotonAndFakeMuon: eventWeights['jg'],
    datasets['DataM'].PhotonAndDimuon: eventWeights['default'],
    datasets['DataM'].ElePhotonAndDimuon: eventWeights['eg'],
    datasets['DataM'].FakePhotonAndDimuon: eventWeights['jg']
}

for name in ['T5wg_1000_425', 'TChiwg_300']:
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

MET = ('E_{T}^{miss}', 'GeV')
MT = ('M_{T}', 'GeV')
HT = ('H_{T}', 'GeV')
PHOTONPT = ('Photon P_{T}', 'GeV')
LEPTONPT = ('Lepton P_{T}', 'GeV')

### EXPORTED ###
hdefs = [
    HDef('MetHighMtHighHtHighPhotonPt', metBinning3, cond = ['M_{T} > 100 GeV', 'H_{T} > 400 GeV', 'P_{T}^{#gamma} > 80 GeV'], xtitle = MET, overflowable = True, mask = ('-inf', 'inf'), logscale = True, vrange = (1.e-3, 1.e+2)),
    HDef('MetHighMtMidHtHighPhotonPt', metBinning2, cond = ['M_{T} > 100 GeV', '100 GeV < H_{T} < 400 GeV', 'P_{T}^{#gamma} > 80 GeV'], xtitle = MET, overflowable = True, mask = ('-inf', 'inf'), logscale = True, vrange = (1.e-3, 1.e+2)),
    HDef('MetHighMtLowHtHighPhotonPt', metBinning2, cond = ['M_{T} > 100 GeV', 'H_{T} < 100 GeV', 'P_{T}^{#gamma} > 80 GeV'], xtitle = MET, overflowable = True, mask = ('-inf', 'inf'), logscale = True, vrange = (1.e-3, 1.e+2)),
    HDef('MetHighMtHighHtLowPhotonPt', metBinning3, cond = ['M_{T} > 100 GeV', 'H_{T} > 400 GeV', 'P_{T}^{#gamma} < 80 GeV'], xtitle = MET, overflowable = True, mask = ('-inf', 'inf'), logscale = True, vrange = (1.e-3, 1.e+2)),
    HDef('MetHighMtMidHtLowPhotonPt', metBinning2, cond = ['M_{T} > 100 GeV', '100 GeV < H_{T} < 400 GeV', 'P_{T}^{#gamma} < 80 GeV'], xtitle = MET, overflowable = True, mask = ('-inf', 'inf'), logscale = True, vrange = (1.e-3, 1.e+2)),
    HDef('MetHighMtLowHtLowPhotonPt', metBinning2, cond = ['M_{T} > 100 GeV', 'H_{T} < 100 GeV', 'P_{T}^{#gamma} < 80 GeV'], xtitle = MET, overflowable = True, mask = ('-inf', 'inf'), logscale = True, vrange = (1.e-3, 1.e+2)),
    HDef('MetDiLepHighHtHighPhotonPt', metBinning3, cond = ['N_{l} #geq 2', 'H_{T} > 400 GeV', 'P_{T}^{#gamma} > 80 GeV'], xtitle = MET, overflowable = True, mask = ('-inf', 'inf'), logscale = True, vrange = (1.e-3, 1.e+2)),
    HDef('MetDiLepMidHtHighPhotonPt', metBinning2, cond = ['N_{l} #geq 2', '100 GeV < H_{T} < 400 GeV', 'P_{T}^{#gamma} > 80 GeV'], xtitle = MET, overflowable = True, mask = ('-inf', 'inf'), logscale = True, vrange = (1.e-3, 1.e+2)),
    HDef('MetDiLepLowHtHighPhotonPt', metBinning2, cond = ['N_{l} #geq 2', 'H_{T} < 100 GeV', 'P_{T}^{#gamma} > 80 GeV'], xtitle = MET, overflowable = True, mask = ('-inf', 'inf'), logscale = True, vrange = (1.e-3, 1.e+2)),
    HDef('Met', metBinning1, xtitle = MET, overflowable = True, mask = (70., 'inf'), logscale = True, vrange = (1.e-2, 5.e+3)),
    HDef('MetHighMt', metBinning2, cond = 'M_{T} > 100 GeV', xtitle = MET, overflowable = True, mask = (70., 'inf'), logscale = True, vrange = (1.e-2, 5.e+3)),
    HDef('Mt', mtBinning2, xtitle = MT, overflowable = True, mask = (100., 'inf'), logscale = True),
    HDef('MetMt', (50, 0., 400, 50, 0., 400.), xtitle = MET, ytitle = MT),
    HDef('Ht', htBinning1, cond = ['P_{T}^{j} > 30 GeV, |#eta^{j}| < 3'], xtitle = HT, overflowable = True, mask = (600., 'inf'), logscale = True, vrange = (1.e-2, 2.e+3)),
    HDef('HtHighMetHighMt', htBinning2, cond = ['P_{T}^{j} > 30 GeV, |#eta^{j}| < 3', 'E_{T}^{miss} > 120 GeV', 'M_{T} > 100 GeV'], xtitle = HT, overflowable = True, mask = (600., 'inf'), logscale = True),
    HDef('Mass2', (60, 60., 120.), xtitle = ('M_{l#gamma}', 'GeV')),
    HDef('Mass2Wide', massBinning, xtitle = ('M_{l#gamma}', 'GeV'), overflowable = True, logscale = True),
    HDef('Mass3', (30, 60., 120.), xtitle = ('M_{ll#gamma}', 'GeV')),
    HDef('Mass3Wide', massBinning, xtitle = ('M_{ll#gamma}', 'GeV'), overflowable = True, logscale = True),
    HDef('Mll', (60, 60., 120.), xtitle = ('M_{ll}', 'GeV'), overflowable = True),
    HDef('PhotonPt', photonPtBinning, xtitle = PHOTONPT, overflowable = True, mask = (80., 'inf'), logscale = True),
    HDef('PhotonPtZoomLowMet', (50, 40., 140.), cond = 'E_{T}^{miss} < 70 GeV', xtitle = PHOTONPT, overflowable = True, logscale = True),
    HDef('PhotonPtHighMet', photonPtBinning, cond = 'E_{T}^{miss} > 120 GeV', xtitle = PHOTONPT, overflowable = True, mask = ('-inf', 'inf'), logscale = True),
    HDef('PhotonPtHighMetHighMt', photonPtBinning2, cond = ['E_{T}^{miss} > 120 GeV', 'M_{T} > 100 GeV'], xtitle = PHOTONPT, overflowable = True, mask = ('-inf', 'inf'), logscale = True),
    HDef('PhotonPtLowMet', photonPtBinning, cond = 'E_{T}^{miss} < 70 GeV', xtitle = PHOTONPT, overflowable = True, logscale = True, vrange = (5.e-2, 1.e+4)),
    HDef('PhotonPtLowMetOnZ', photonPtBinning, cond = ['E_{T}^{miss} < 70 GeV', '81 GeV < M_{ll} < 101 GeV'], xtitle = PHOTONPT, overflowable = True, logscale = True),
    HDef('PhotonEta', (60, -1.5, 1.5), xtitle = 'Photon #eta'),
    HDef('LeptonPt', leptonPtBinning, xtitle = LEPTONPT, overflowable = True, mask = (120., 'inf'), logscale = True),
    HDef('LeptonPtHighMet', leptonPtBinning, cond = 'E_{T}^{miss} > 120 GeV', xtitle = LEPTONPT, overflowable = True, mask = ('-inf', 'inf'), logscale = True),
    HDef('LeptonPtLowMet', leptonPtBinning, cond = 'E_{T}^{miss} < 70 GeV', xtitle = LEPTONPT, overflowable = True, logscale = True, vrange = (5.e-2, 3.e+3)),
    HDef('LeptonEta', (60, -2.5, 2.5), xtitle = 'Lepton #eta'),
    HDef('DEtaPhotonLepton', dEtaBinning, xtitle = '#Delta#eta(#gamma, l)'),
    HDef('DPhiPhotonLepton', dPhiBinning, xtitle = '#Delta#phi(#gamma, l)'),
    HDef('DPhiPhotonLeptonMet4070', dPhiBinning, xtitle = '#Delta#phi(#gamma, l)', cond = '40 GeV < E_{T}^{miss} < 70 GeV'),
    HDef('DPhiPhotonLeptonOnZ', dPhiBinning, xtitle = '#Delta#phi(#gamma, l)', cond = '81 GeV < M_{ll} < 101 GeV'),
    HDef('DRPhotonLepton', dRBinning, xtitle = '#DeltaR(#gamma, l)'),
    HDef('DRPhotonLeptonMet4070', dRBinning, xtitle = '#DeltaR(#gamma, l)', cond = '40 GeV < E_{T}^{miss} < 70 GeV'),
    HDef('DRPhotonTrailLepton', dRBinning, xtitle = '#DeltaR(#gamma, l_{2})'),
    HDef('DEtaPhotonJet', dEtaBinning, xtitle = '#Delta#eta(#gamma, j)'),
    HDef('DPhiPhotonJet', dPhiBinning, xtitle = '#Delta#phi(#gamma, j)'),
    HDef('DRPhotonJet', dRBinning, xtitle = '#DeltaR(#gamma, j)'),
    HDef('DEtaLeptonJet', dEtaBinning, xtitle = '#Delta#eta(l, j)'),
    HDef('DPhiLeptonJet', dPhiBinning, xtitle = '#Delta#phi(l, j)'),
    HDef('DRLeptonJet', dRBinning, xtitle = '#DeltaR(l, j)'),
    HDef('DPhiPhotonMet', dPhiBinning, xtitle = '#Delta#phi(#gamma, E_{T}^{miss})'),
    HDef('DPhiPhotonMetMet4070', dPhiBinning, xtitle = '#Delta#phi(#gamma, E_{T}^{miss})', cond = '40 GeV < E_{T}^{miss} < 70 GeV'),
    HDef('DPhiLeptonMet', dPhiBinning, xtitle = '#Delta#phi(l, E_{T}^{miss})'),
    HDef('DPhiLeptonMetMet4070', dPhiBinning, xtitle = '#Delta#phi(l, E_{T}^{miss})', cond = '40 GeV < E_{T}^{miss} < 70 GeV'),
    HDef('NPhoton', (3, 0.5, 3.5), xtitle = ('N_{#gamma}', 'NoUnit'), xlabels = [str(i) for i in range(1, 4)], logscale = True),
    HDef('NLepton', (4, 0.5, 4.5), xtitle = ('N_{l}', 'NoUnit'), xlabels = [str(i) for i in range(1, 5)], logscale = True),
    HDef('NPhotonNLepton', (4, 0.5, 4.5, 3, 0.5, 3.5), xtitle = ('N_{l}', 'NoUnit'), xlabels = [str(i) for i in range(1, 4)], ytitle = 'N_{#gamma}', ylabels = [str(i) for i in range(1, 5)], drawOption = 'BOX TEXT'),
    HDef('NJet', (11, -0.5, 10.5), cond = ['P_{T}^{j} > 30 GeV, |#eta^{j}| < 3'], xtitle = ('N_{jet}', 'NoUnit'), xlabels = [str(i) for i in range(1, 10)], mask = (4.5, 'inf'), logscale = True, vrange = (1.e-2, 1.e+6)),
    HDef('NJetLowMet', (11, -0.5, 10.5), cond = ['P_{T}^{j} > 30 GeV, |#eta^{j}| < 3', 'E_{T}^{miss} < 70 GeV'], xtitle = ('N_{jet}', 'NoUnit'), xlabels = [str(i) for i in range(1, 10)], logscale = True),
    HDef('NVtx', (51, -0.5, 50.5), xtitle = ('N_{vtx}', 'NoUnit')),
    HDef('PW', (50, 0., 400.), xtitle = ('P_{W}', 'GeV'), overflowable = True),
    HDef('MW', (50, 80., 280.), xtitle = ('M_{l#slash{E}_T}', 'GeV'), overflowable = True, mask = ('-inf', 'inf')),
#    HDef('NLW', '-ln(#Sigma)', (50, 20., 70.), xtitle = '-ln(#Sigma)', overflowable = True),
#    HDef('NLWHighMet', '-ln(#Sigma) (#slash{E}_{T} > 120 GeV)', (50, 20., 70.), xtitle = '-ln(#Sigma)', overflowable = True)
]

hdefList = dict([(h.name, h) for h in hdefs])

hdefsClosure = [
    HDef('Met', metBinning1, xtitle = MET, logscale = True),
    HDef('Mt', mtBinning2, xtitle = MT, logscale = True),
    HDef('PhotonPt', [40. + x * 5. for x in range(8)] + [80., 90., 100., 120., 140., 200.], xtitle = PHOTONPT, logscale = True),
    HDef('LeptonPt', [25. + x * 5. for x in range(11)] + [80., 90., 100., 140., 200.], xtitle = LEPTONPT, logscale = True),
    hdefList['NVtx']
]

############## PLOT MAKERS #############

ROOT.gROOT.LoadMacro(thisdir + '/GLPlotMaker.cc+')

############## STACK CONFIGURATIONS ##############
### EXPORTED ###
stackConfigs = {}
searchPlots = [
    hdefList['MetHighMtLowHtHighPhotonPt'],
    hdefList['MetHighMtMidHtHighPhotonPt'],
    hdefList['MetHighMtHighHtHighPhotonPt'],
    hdefList['MetHighMtLowHtLowPhotonPt'],
    hdefList['MetHighMtMidHtLowPhotonPt'],
    hdefList['MetHighMtHighHtLowPhotonPt']
]

hdefsE = map(lambda h: h.clone(), hdefs)
for h in hdefsE: h.conditions = ['e + #gamma selection'] + h.conditions
hdefsM = map(lambda h: h.clone(), hdefs)
for h in hdefsM: h.conditions = ['#mu + #gamma selection'] + h.conditions
searchPlotsE = map(lambda h: h.clone(), searchPlots)
for h in searchPlotsE: h.conditions = ['e + #gamma selection'] + h.conditions
searchPlotsM = map(lambda h: h.clone(), searchPlots)
for h in searchPlotsM: h.conditions = ['#mu + #gamma selection'] + h.conditions

######## Search config with floating VGamma ########

from floatingVGamma import FloatingVGammaSearch

######## Systematic evaluation using lepton Pt reweighting #########

leptonPtSource = ROOT.TFile.Open('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/dimuonPt.root')
ptlWeight = leptonPtSource.Get('ratio')

######## Histogram Colors #########

cEWK = ROOT.TColor.GetColor(255, 168, 0)
cEGFake = ROOT.TColor.GetColor(255, 101, 102)
cJLFake = ROOT.TColor.GetColor(108, 156, 254)
cJGFake = ROOT.TColor.GetColor(107, 255, 153)
cVGamma = ROOT.TColor.GetColor(220, 56, 252)

######## Electron Channel Search ########

gObservedE = Group('Observed', 'Observed', ROOT.kBlack, Group.OBSERVED, [datasets['DataE'].PhotonAndElectron])
gEGFakeE = Group('EGFake', 'e#rightarrow#gamma fake', cEGFake, Group.BACKGROUND, [datasets['DataE'].ElePhotonAndElectron])
gJGFakeE = Group('JGFake', 'j#rightarrow#gamma fake', cJGFake, Group.BACKGROUND, [datasets['DataE'].FakePhotonAndElectron])
gJLFakeE = Group('JLFake', 'QCD', cJLFake, Group.BACKGROUND, [datasets['DataE'].PhotonAndFakeElectron])
gVGammaE = Group('VGamma', 'V#gamma (MC)', cVGamma, Group.BACKGROUND, [
    datasets['WGToLNuG_PtG-30-50'].PhotonAndElectron,
    datasets['WGToLNuG_PtG-50-130'].PhotonAndElectron,
    datasets['WGToLNuG_PtG-130'].PhotonAndElectron,
    datasets['ZGToLLG_PtG-5-130'].PhotonAndElectron,
    datasets['ZGToLLG_PtG-130'].PhotonAndElectron
])
gEWKE = Group('EWK', 't#bar{t}#gamma/W^{+}W^{-}#gamma (MC)', cEWK, Group.BACKGROUND, [
    datasets['WWGJets'].PhotonAndElectron,
    datasets['WW'].PhotonAndElectron,
    datasets['WZ'].PhotonAndElectron,
#    datasets['ttA'].PhotonAndElectron,
    datasets['TTGJets'].PhotonAndElectron,
    datasets['TTJetsSemiLept'].PhotonAndElectron,
    datasets['TTJetsFullLept'].PhotonAndElectron
])
gT5wg_1000_425E = Group('T5wg_1000_425', 'T5wg_1000_425', ROOT.kGreen + 3, Group.SIGNAL, [
    datasets['T5wg_1000_425'].PhotonAndElectron,
    datasets['T5wg_1000_425'].ElePhotonAndElectron,
    datasets['T5wg_1000_425'].FakePhotonAndElectron,
    datasets['T5wg_1000_425'].PhotonAndFakeElectron
])
gTChiwg_300E = Group('TChiwg_300', 'TChiwg_300', ROOT.kBlue, Group.SIGNAL, [
    datasets['TChiwg_300'].PhotonAndElectron,
    datasets['TChiwg_300'].ElePhotonAndElectron,
    datasets['TChiwg_300'].FakePhotonAndElectron,
    datasets['TChiwg_300'].PhotonAndFakeElectron
])

plotMakerE = ROOT.GLPlotMaker(0)

stackConfigs['FloatingVGammaE'] = FloatingVGammaSearch(plotMakerE, 'Electron', hdefList['DPhiLeptonMetMet4070'])
stackConfigs['FloatingVGammaE'].groups = [
    gObservedE,
    gEWKE,
    gEGFakeE,
    gJLFakeE,
    gJGFakeE,
    gVGammaE,
    gT5wg_1000_425E,
    gTChiwg_300E
]
stackConfigs['FloatingVGammaE'].hdefs = hdefsE

######## Muon Channel Search ########

gObservedM = Group('Observed', 'Observed', ROOT.kBlack, Group.OBSERVED, [datasets['DataM'].PhotonAndMuon])
gEGFakeM = Group('EGFake', 'e#rightarrow#gamma fake', cEGFake, Group.BACKGROUND, [datasets['DataM'].ElePhotonAndMuon])
gJGFakeM = Group('JGFake', 'j#rightarrow#gamma fake', cJGFake, Group.BACKGROUND, [datasets['DataM'].FakePhotonAndMuon])
gJLFakeM = Group('JLFake', 'QCD', cJLFake, Group.BACKGROUND, [datasets['DataM'].PhotonAndFakeMuon])
gVGammaM = Group('VGamma', 'V#gamma (MC)', cVGamma, Group.BACKGROUND, [
    datasets['WGToLNuG_PtG-30-50'].PhotonAndMuon,
    datasets['WGToLNuG_PtG-50-130'].PhotonAndMuon,
    datasets['WGToLNuG_PtG-130'].PhotonAndMuon,
    datasets['ZGToLLG_PtG-5-130'].PhotonAndMuon,
    datasets['ZGToLLG_PtG-130'].PhotonAndMuon
])
gEWKM = Group('EWK', 't#bar{t}#gamma/W^{+}W^{-}#gamma (MC)', cEWK, Group.BACKGROUND, [
    datasets['WWGJets'].PhotonAndMuon,
    datasets['WW'].PhotonAndMuon,
    datasets['WZ'].PhotonAndMuon,
#    datasets['ttA'].PhotonAndMuon,
    datasets['TTGJets'].PhotonAndMuon,
    datasets['TTJetsSemiLept'].PhotonAndMuon,
    datasets['TTJetsFullLept'].PhotonAndMuon
])
gT5wg_1000_425M = Group('T5wg_1000_425', 'T5wg_1000_425', ROOT.kGreen + 3, Group.SIGNAL, [
    datasets['T5wg_1000_425'].PhotonAndMuon,
    datasets['T5wg_1000_425'].ElePhotonAndMuon,
    datasets['T5wg_1000_425'].FakePhotonAndMuon,
    datasets['T5wg_1000_425'].PhotonAndFakeMuon
])
gTChiwg_300M = Group('TChiwg_300', 'TChiwg_300', ROOT.kBlue, Group.SIGNAL, [
    datasets['TChiwg_300'].PhotonAndMuon,
    datasets['TChiwg_300'].ElePhotonAndMuon,
    datasets['TChiwg_300'].FakePhotonAndMuon,
    datasets['TChiwg_300'].PhotonAndFakeMuon
])

plotMakerM = ROOT.GLPlotMaker(1)

stackConfigs['FloatingVGammaM'] = FloatingVGammaSearch(plotMakerM, 'Muon', hdefList['DPhiLeptonMetMet4070'])
stackConfigs['FloatingVGammaM'].groups = [
    gObservedM,
    gEWKM,
    gEGFakeM,
    gJLFakeM,
    gJGFakeM,
    gVGammaM,
    gT5wg_1000_425M,
    gTChiwg_300M
]
stackConfigs['FloatingVGammaM'].hdefs = hdefsM

######### Post Simul-fit results ###########

from fixedScales import FixedScalesSearch

stackConfigs['FixedScalesE'] = FixedScalesSearch(stackConfigs['FloatingVGammaE'], (1.60221674546, -7.97187e-02, 0.11049638913), 0.234907178758)
stackConfigs['FixedScalesM'] = FixedScalesSearch(stackConfigs['FloatingVGammaM'], (1.60221674546, -1.87756e+00, 0.11049638913), 0.00404211957313)

######### Dimuon ##########

plotMakerMM = ROOT.GLPlotMaker(1)
plotMakerMM.dileptonOnly()

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
        source.Close()

        zgamma = next(g for g in self.groups if g.name == 'ZGamma')
        for sample in zgamma.samples:
            for histogram in sample.histograms.values():
                histogram.hWeighted.Scale(scale)
                histogram.hScaleUp.Scale(scale + errUp)
                histogram.hScaleDown.Scale(scale - errDown)

            count = sample.counter.GetBinContent(1)
            sample.counter.SetBinContent(1, count * scale)
            sample.counter.SetBinContent(2, count * (scale + max(errUp, errDown)))


gObservedMM = Group('Observed', 'Observed', ROOT.kBlack, Group.OBSERVED, [datasets['DataM'].PhotonAndDimuon])
gEGFakeMM = Group('EGFake', 'e#rightarrow#gamma fake', ROOT.kOrange, Group.BACKGROUND, [datasets['DataM'].ElePhotonAndDimuon])
gJGFakeMM = Group('JGFake', 'j#rightarrow#gamma fake', ROOT.kGreen, Group.BACKGROUND, [datasets['DataM'].FakePhotonAndDimuon])
gZGammaMM = Group('ZGamma', 'Z#gamma (MC)', ROOT.kMagenta, Group.BACKGROUND, [
    datasets['ZGToLLG_PtG-5-130'].PhotonAndDimuon,
    datasets['ZGToLLG_PtG-130'].PhotonAndDimuon
])
gWWMM = Group('WW', 'W^{+}W^{-}#gamma (MC)', ROOT.kSpring, Group.BACKGROUND, [
    datasets['WWGJets'].PhotonAndDimuon,
    datasets['WW'].PhotonAndDimuon
])
gTTMM = Group('TT', 't#bar{t}#gamma (MC)', ROOT.kCyan, Group.BACKGROUND, [
    datasets['TTGJets'].PhotonAndDimuon,
    datasets['TTJetsFullLept'].PhotonAndDimuon
])

stackConfigs['DimuonControl'] = DimuonControl('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/FloatingVGammaM.root')
stackConfigs['DimuonControl'].groups = [
    gObservedMM,
    gEGFakeMM,
    gJGFakeMM,
    gZGammaMM,
    gWWMM,
    gTTMM
]
stackConfigs['DimuonControl'].hdefs = [
    HDef('MuonPtLowMet', [4. * x for x in range(101)], cond = ['E_{T}^{miss} < 70 GeV', 'N_{#mu} == 2'], xtitle = ('Muon P_{T}', 'GeV'), overflowable = True, logscale = True),
    hdefList['PhotonPtLowMet'],
    hdefList['LeptonPtLowMet'],
    hdefList['Mass2'],
    hdefList['Mass3'],
    hdefList['Mll']
]

######### E->gamma Fake Closure (Electron Channel) ########

gEGFakeTruthE = Group('EGFake', 'Fake (MC truth)', ROOT.kBlack, Group.OBSERVED, [
    datasets['DYJetsToLL'].PhotonAndElectron,
    datasets['TTJetsFullLept'].PhotonAndElectron,
    datasets['WW'].PhotonAndElectron
])
gDYEleProxyE = Group('DYEleProxy', 'DY Proxy (MC)', ROOT.kOrange, Group.BACKGROUND, [datasets['DYJetsToLL'].ElePhotonAndElectron])
gEWKEleProxyE = Group('EWKEleProxy', 'EWK Proxy (MC)', ROOT.kOrange + 10, Group.BACKGROUND, [
    datasets['TTJetsFullLeptFake'].ElePhotonAndElectron,
    datasets['WW'].ElePhotonAndElectron
])

plotMakerEL = ROOT.GLPlotMaker(0)
plotMakerEL.noEffCorrection()
plotMakerEE = ROOT.GLPlotMaker(0)
plotMakerEE.matchElePhoton()
plotMakerEE.noEffCorrection()

stackConfigs['EGFakeClosureE'] = StackConfig(plotMakerEL)
stackConfigs['EGFakeClosureE'].groups = [
    gEGFakeTruthE,
    gDYEleProxyE,
    gEWKEleProxyE
]
stackConfigs['EGFakeClosureE'].hdefs = hdefsClosure
stackConfigs['EGFakeClosureE'].specialPlotMakers = {
    gEGFakeTruthE: plotMakerEE
}

######### E->gamma Fake Closure (Muon Channel) ########

gEGFakeTruthM = Group('EGFake', 'Fake (MC truth)', ROOT.kBlack, Group.OBSERVED, [
    datasets['DYJetsToLL'].PhotonAndMuon,
    datasets['TTJetsFullLept'].PhotonAndMuon,
    datasets['WW'].PhotonAndMuon
])
gDYEleProxyM = Group('DYEleProxy', 'DY Proxy (MC)', ROOT.kOrange, Group.BACKGROUND, [datasets['DYJetsToLL'].ElePhotonAndMuon])
gEWKEleProxyM = Group('EWKEleProxy', 'EWK Proxy (MC)', ROOT.kOrange + 10, Group.BACKGROUND, [
    datasets['TTJetsFullLeptFake'].ElePhotonAndMuon,
    datasets['WW'].ElePhotonAndMuon
])

plotMakerML = ROOT.GLPlotMaker(1)
plotMakerML.noEffCorrection()
plotMakerME = ROOT.GLPlotMaker(1)
plotMakerME.matchElePhoton()
plotMakerME.noEffCorrection()

stackConfigs['EGFakeClosureM'] = StackConfig(plotMakerM)
stackConfigs['EGFakeClosureM'].groups = [
    gEGFakeTruthM,
    gDYEleProxyM,
    gEWKEleProxyM
]
stackConfigs['EGFakeClosureM'].hdefs = hdefsClosure
stackConfigs['EGFakeClosureM'].specialPlotMakers = {
    gEGFakeTruthM: plotMakerME
}

######### J->gamma Fake Closure (Electron Channel) ########

gJGFakeTruthE = Group('JGFake', 'Fake (MC truth)', ROOT.kBlack, Group.OBSERVED, [
#    datasets['WJetsToLNu'].PhotonAndElectron,
    datasets['WJetsToLNu_PtW-50To70'].PhotonAndElectron,
    datasets['WJetsToLNu_PtW-70To100'].PhotonAndElectron,
    datasets['WJetsToLNu_PtW-100'].PhotonAndElectron,
    datasets['DYJetsToLL'].PhotonAndElectron,
    datasets['TTJetsSemiLept'].PhotonAndElectron,
    datasets['TTJetsFullLept'].PhotonAndElectron,
    datasets['WW'].PhotonAndElectron
])
gWJJetProxyE = Group('WJJetProxy', 'WJets Proxy (MC)', ROOT.kGreen, Group.BACKGROUND, [
#    datasets['WJetsToLNu'].FakePhotonAndElectron,
    datasets['WJetsToLNu_PtW-50To70'].FakePhotonAndElectron,
    datasets['WJetsToLNu_PtW-70To100'].FakePhotonAndElectron,
    datasets['WJetsToLNu_PtW-100'].FakePhotonAndElectron
])
gDYJetProxyE = Group('DYJetProxy', 'DY Proxy (MC)', ROOT.kSpring + 9, Group.BACKGROUND, [datasets['DYJetsToLL'].FakePhotonAndElectron])
gEWKJetProxyE = Group('EWKJetProxy', 'EWK Proxy (MC)', ROOT.kGreen + 8, Group.BACKGROUND, [
    datasets['TTJetsSemiLeptFake'].FakePhotonAndElectron,
    datasets['TTJetsFullLeptFake'].FakePhotonAndElectron,
    datasets['WW'].FakePhotonAndElectron    
])

plotMakerEJ = ROOT.GLPlotMaker(0)
plotMakerEJ.noEffCorrection()
plotMakerEJ.matchFakePhoton()

stackConfigs['JGFakeClosureE'] = StackConfig(plotMakerEL)
stackConfigs['JGFakeClosureE'].groups = [
    gJGFakeTruthE,
    gWJJetProxyE,
    gDYJetProxyE,
    gEWKJetProxyE
]
stackConfigs['JGFakeClosureE'].hdefs = hdefsClosure
stackConfigs['JGFakeClosureE'].specialPlotMakers = {
    gJGFakeTruthE: plotMakerEJ
}

######### J->gamma Fake Closure (Muon Channel) ########

gJGFakeTruthM = Group('JGFake', 'Fake (MC truth)', ROOT.kBlack, Group.OBSERVED, [
#    datasets['WJetsToLNu'].PhotonAndMuon,
    datasets['WJetsToLNu_PtW-50To70'].PhotonAndMuon,
    datasets['WJetsToLNu_PtW-70To100'].PhotonAndMuon,
    datasets['WJetsToLNu_PtW-100'].PhotonAndMuon,
    datasets['DYJetsToLL'].PhotonAndMuon,
    datasets['TTJetsSemiLept'].PhotonAndMuon,
    datasets['TTJetsFullLept'].PhotonAndMuon,
    datasets['WW'].PhotonAndMuon
])
gWJJetProxyM = Group('WJJetProxy', 'WJets Proxy (MC)', ROOT.kGreen, Group.BACKGROUND, [
#    datasets['WJetsToLNu'].FakePhotonAndMuon,
    datasets['WJetsToLNu_PtW-50To70'].FakePhotonAndMuon,
    datasets['WJetsToLNu_PtW-70To100'].FakePhotonAndMuon,
    datasets['WJetsToLNu_PtW-100'].FakePhotonAndMuon
])
gDYJetProxyM = Group('DYJetProxy', 'DY Proxy (MC)', ROOT.kSpring + 9, Group.BACKGROUND, [datasets['DYJetsToLL'].FakePhotonAndMuon])
gEWKJetProxyM = Group('EWKJetProxy', 'EWK Proxy (MC)', ROOT.kGreen + 8, Group.BACKGROUND, [
    datasets['TTJetsSemiLeptFake'].FakePhotonAndMuon,
    datasets['TTJetsFullLeptFake'].FakePhotonAndMuon,
    datasets['WW'].FakePhotonAndMuon
])

plotMakerMJ = ROOT.GLPlotMaker(1)
plotMakerMJ.noEffCorrection()
plotMakerMJ.matchFakePhoton()

stackConfigs['JGFakeClosureM'] = StackConfig(plotMakerML)
stackConfigs['JGFakeClosureM'].groups = [
    gJGFakeTruthM,
    gWJJetProxyM,
    gDYJetProxyM,
    gEWKJetProxyM
]
stackConfigs['JGFakeClosureM'].hdefs = hdefsClosure
stackConfigs['JGFakeClosureM'].specialPlotMakers = {
    gJGFakeTruthM: plotMakerMJ
}

######## Test ########

stackConfigs['TestE'] = FloatingVGammaSearch(plotMakerE, 'Electron', hdefList['DPhiLeptonMetMet4070'])
stackConfigs['TestE'].groups = [
    gObservedE,
    gEWKE,
    gEGFakeE,
    gJLFakeE,
    gJGFakeE,
    gVGammaE,
    gT5wg_1000_425E,
    gTChiwg_300E
]
stackConfigs['TestE'].hdefs = hdefsE

stackConfigs['TestM'] = FloatingVGammaSearch(plotMakerM, 'Muon', hdefList['DPhiLeptonMetMet4070'])
stackConfigs['TestM'].groups = [
    gObservedM,
    gEWKM,
    gEGFakeM,
    gJLFakeM,
    gJGFakeM,
    gVGammaM,
    gT5wg_1000_425M,
    gTChiwg_300M
]
stackConfigs['TestM'].hdefs = hdefsM
