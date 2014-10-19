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
from config import searchPlots

thisdir = os.path.dirname(os.path.abspath(__file__))

class SignalStack(StackConfig):

    def __init__(self, plotMaker, lepton):
        StackConfig.__init__(self, plotMaker)

        self.jlClass = 'PhotonAndFake' + lepton
        self.candClass = 'PhotonAnd' + lepton

        if lepton == 'Electron':
            source = ROOT.TFile('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/FloatingVGammaE.root')
            self.jlScale = source.Get('TemplateFitError/QCD').GetY()[0]
        elif lepton == 'Muon':
            source = ROOT.TFile('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/FloatingVGammaM.root')
            self.jlScale = source.Get('TemplateFitError/QCD').GetY()[0]

        self.hdefs = searchPlots

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


def getDataset(model, mass):
    point = model + '_' + mass

    with open('/afs/cern.ch/user/y/yiiyama/output/GammaL/limits/xsecs/' + model + '.xsecs') as xsecsource:
        for line in xsecsource:
            p, c, u, n = line.strip().split()
            if p == point:
                xsec = float(c)
                nEvents = int(n)
                break
        else:
            raise RuntimeError('No point ' + point + ' found')

    if nEvents == 0:
        return None

    eventClasses = [
        'PhotonAndElectron',
        'ElePhotonAndElectron',
        'FakePhotonAndElectron',
        'PhotonAndFakeElectron',
        'PhotonAndMuon',
        'ElePhotonAndMuon',
        'FakePhotonAndMuon',
        'PhotonAndFakeMuon'
    ]

    return Dataset(point, [model + '/skim_' + mass], Dataset.FASTSIM, nEvents / xsec, 0., eventClasses)

def countSignal(model, mass, treeDir):
    point = model + '_' + mass

    dataset = getDataset(model, mass)
    if not dataset: return
    
    eventWeights = {}
    eventWeights['default'] = ROOT.GLEventWeight()
    eventWeights['egHLT'] = ROOT.ElePhotonFunctionalWeight()
    eventWeights['eg'] = ROOT.ElePhotonFunctionalWeight(True)
    eventWeights['jgHLT'] = ROOT.JetPhotonHLTIsoWeight()
    eventWeights['jg'] = ROOT.JetPhotonWeight()
    
    weightCalcList = {
        'PhotonAndElectron': eventWeights['default'],
        'ElePhotonAndElectron': eventWeights['egHLT'],
        'FakePhotonAndElectron': eventWeights['jgHLT'],
        'PhotonAndFakeElectron': eventWeights['default'],
        'PhotonAndMuon': eventWeights['default'],
        'ElePhotonAndMuon': eventWeights['eg'],
        'FakePhotonAndMuon': eventWeights['jg'],
        'PhotonAndFakeMuon': eventWeights['default']
    }
    
    weightCalc = dict([(dataset.samples[s], w) for s, w in weightCalcList.items()])
        
    weightEvents(dataset, ROOT.GLSkimProcessor, weightCalc, sourceDir, treeDir)


if __name__ == '__main__':

    from optparse import OptionParser

    parser = OptionParser(usage = 'Usage: countSignal.py [options]')

    parser.add_option('-p', '--point', dest = 'point', default = '', help = 'Process single point')
    parser.add_option('-m', '--model', dest = 'model', default = '', help = 'Model to process')

    options, args = parser.parse_args()

    treeDir = os.environ['TMPDIR'] + '/countSignal'
    try:
        os.makedirs(treeDir)
    except OSError:
        pass

    models = []
    if options.model:
        if options.point:
            countSignal(options.model, options.point, treeDir)

        else:
            models.append(options.model)

    else:
        if options.point:
            raise RuntimeError('Point given without model')

        models = ['T5wg', 'TChiwg', 'Spectra_gW']

    if 'T5wg' in models:
        try:
            os.makedirs(treeDir + '/T5wg')
        except OSError:
            pass

        for mglu in range(400, 1550, 50):
            for mchi in range(25, mglu, 50):
                pointName = '{mglu}_{mchi}'.format(mglu = mglu, mchi = mchi)
                print 'T5wg_' + pointName
                countSignal('T5wg', pointName, treeDir + '/T5wg')

    if 'TChiwg' in models:
        try:
            os.makedirs(treeDir + '/TChiwg')
        except OSError:
            pass

        for mchi in range(100, 810, 10):
            pointName = str(mchi)
            print 'TChiwg_' + pointName
            countSignal('TChiwg', pointName, treeDir + '/TChiwg')
        
    if 'TChiwgSuppl' in models:
        try:
            os.makedirs(treeDir + '/TChiwg')
        except OSError:
            pass

        for mchi in range(125, 1625, 50):
            pointName = str(mchi)
            print 'TChiwg_' + pointName
            countSignal('TChiwg', pointName, treeDir + '/TChiwg')

    if 'Spectra_gW' in models:
        try:
            os.makedirs(treeDir + '/Spectra_gW')
        except OSError:
            pass

        for m3 in range(715, 1565, 50):
            for m2 in range(205, m3, 50):
                for proc in ['gg', 'ncp', 'ncm']:
                    pointName = 'M3_{m3}_M2_{m2}_{proc}'.format(m3 = m3, m2 = m2, proc = proc)
                    print 'Spectra_gW_' + pointName
                    countSignal('Spectra_gW', pointName, treeDir + '/Spectra_gW')

    print 'Signal trees and plots are created in', treeDir, '- do not forget to copy them back!'
