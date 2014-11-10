import sys
import os
import subprocess
import time
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
    parser.add_option('-c', '--copy', dest = 'copy', action = 'store_true', help = 'Copy files to ncmu40 at the end of execution')

    options, args = parser.parse_args()

    treeDir = os.environ['TMPDIR'] + '/countSignal'
    try:
        os.makedirs(treeDir)
    except OSError:
        pass

    if options.model:
        model = [options.model]

    else:
        models = ['T5wg', 'TChiwg', 'Spectra_gW']

    pointList = {}

    if options.point:
        if not options.model:
            raise RuntimeError('Point given without model')

        pointList[options.model] = [options.point]

    else:
        if 'T5wg' in models:
            pointList['T5wg'] = ['{mglu}_{mchi}'.format(mglu = mglu, mchi = mchi) for mglu in range(400, 1550, 50) for mchi in range(25, mglu, 50)]
    
        if 'TChiwg' in models:
            pointList['TChiwg'] = [str(mchi) for mchi in range(100, 810, 10)]
            
        if 'TChiwgSuppl' in models:
            pointList['TChiwgSuppl'] = [str(mchi) for mchi in range(125, 1625, 50)]
    
        if 'Spectra_gW' in models:
            pointList['Spectra_gW'] = ['M3_{m3}_M2_{m2}_{proc}'.format(m3 = m3, m2 = m2, proc = proc) for m3 in range(715, 1565, 50) for m2 in range(205, m3, 50) for proc in ['gg', 'ncp', 'ncm']]

    for model, points in pointList.items():
        try:
            os.makedirs(treeDir + '/' + model)
        except OSError:
            pass

        for pointName in points:
            print model + '_' + pointName
            countSignal(model, pointName, treeDir + '/' + model)

    if options.copy:
        for directory in os.listdir(treeDir):
            proc = subprocess.Popen(['tar', '-czf', treeDir + '/' + directory + '.tar.gz', treeDir + '/' + directory], stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
            proc.communicate()

            proc = subprocess.Popen(['scp', treeDir + '/' + directory + '.tar.gz', 'ncmu40:/store/countSignal/'], stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
            proc.communicate()

            proc = subprocess.Popen(['ssh', 'ncmu40', '"tar -xzf /store/countSignal/' + directory + '.tar.gz -C /store/countSignal/ && rm /store/countSignal/' + directory + '.tar.gz"'], stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
            proc.communicate()

    else:
        print 'Signal trees and plots are created in', treeDir, '- do not forget to copy them back!'
