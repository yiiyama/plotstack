import sys
import os
import ROOT

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import rootconfig
from locations import sourceDir
from dataset import Dataset
from histogram import HDef
from stack import Group, StackConfig
from weightEvents import weightEvents
from makeStack import fillPlots, formGroups
from config import hdefList

thisdir = os.path.dirname(os.path.abspath(__file__))

class SignalStack(StackConfig):

    def __init__(self, plotMaker, lepton):
        StackConfig.__init__(self, plotMaker)

        self.jlClass = 'PhotonAndFake' + lepton
        self.candClass = 'PhotonAnd' + lepton

        source = ROOT.TFile('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/simulFit.root')
        if lepton == 'Electron':
            self.jlScale = source.Get('elQCD').GetY()[0]
        elif lepton == 'Muon':
            self.jlScale = source.Get('muQCD').GetY()[0]

        self.hdefs = [
            hdefList['MetHighMtLowPhotonPt'],
            hdefList['MetHighMtHighPhotonPt']
        ]

    def scalePlots(self, outputDir):
        for group in self.groups: # should only be one
            for sample in group.samples:
                for histogram in sample.histograms.values():
                    if sample.eventClass == self.jlClass:
                        histogram.hWeighted.Scale(-self.jlScale)
                        histogram.hScaleUp.Scale(-self.jlScale)
                        histogram.hScaleDown.Scale(-self.jlScale)
                    elif sample.eventClass != self.candClass:
                        histogram.hWeighted.Scale(-1.)
                        histogram.hScaleUp.Scale(-1.)
                        histogram.hScaleDown.Scale(-1.)


def countSignal(model, mass):
    point = model + '_' + mass
    
    with open('/afs/cern.ch/user/y/yiiyama/src/GammaL/xsec/' + model + '.xsecs') as xsecsource:
        for line in xsecsource:
            p, c, u, n = line.strip().split()
            if p == point:
                xsec = float(c)
                nEvents = int(n)
                break
        else:
            raise RuntimeError('No point ' + point + ' found')
    
    eventWeights = {}
    eventWeights['default'] = ROOT.GLEventWeight()
    eventWeights['eg'] = ROOT.ElePhotonFunctionalWeight()
    eventWeights['jgHLT'] = ROOT.JetPhotonHLTIsoWeight()
    eventWeights['jg'] = ROOT.JetPhotonWeight()
    
    eventweights = {
        'PhotonAndElectron': eventWeights['default'],
        'ElePhotonAndElectron': eventWeights['eg'],
        'FakePhotonAndElectron': eventWeights['jgHLT'],
        'PhotonAndFakeElectron': eventWeights['default'],
        'PhotonAndMuon': eventWeights['default'],
        'ElePhotonAndMuon': eventWeights['eg'],
        'FakePhotonAndMuon': eventWeights['jg'],
        'PhotonAndFakeMuon': eventWeights['default']
    }
    
    dataset = Dataset(point, [model + '/skim_' + mass], Dataset.FASTSIM, nEvents / xsec, 0., eventweights.keys())
    
    weightCalc = dict([(dataset.samples[s], w) for s, w in eventweights.items()])
    
    weightEvents(dataset, ROOT.GLSkimProcessor, weightCalc, sourceDir, os.environ['TMPDIR'] + '/countSignal/trees')
    
    outputFile = ROOT.TFile(os.environ['TMPDIR'] + '/countSignal/plots/' + point + '.root', 'recreate')
    
    for iL, lepton, Lnorm in [(0, 'Electron', 876.225 + 4411.704 + 7054.732 + 7369.007), (1, 'Muon', 876.225 + 4411.704 + 7054.732 + 7360.046)]:
        plotMaker = ROOT.GLPlotMaker(iL)
    
        config = SignalStack(plotMaker, lepton)
        config.groups = [
            Group(point, point, ROOT.kBlack, Group.SIGNAL, [
                dataset.samples['PhotonAnd' + lepton],
                dataset.samples['ElePhotonAnd' + lepton],
                dataset.samples['FakePhotonAnd' + lepton],
                dataset.samples['PhotonAndFake' + lepton]
            ])
        ]

        directory = outputFile.mkdir(lepton)
    
        fillPlots(config, os.environ['TMPDIR'] + '/countSignal/trees', directory, integratedLumi = Lnorm)
        formGroups(config, directory)
    
    outputFile.Close()

if __name__ == '__main__':

    try:
        os.makedirs(os.environ['TMPDIR'] + '/countSignal/trees')
        os.makedirs(os.environ['TMPDIR'] + '/countSignal/plots')
    except OSError:
        pass

    for mglu in range(400, 1550, 50):
        if mglu >= 800 and mglu < 1000: continue
        if mglu >= 1350: continue
        for mchi in range(25, mglu, 50):
            pointName = '{mglu}_{mchi}'.format(mglu = mglu, mchi = mchi)
            print pointName
            countSignal('T5wg', pointName)
