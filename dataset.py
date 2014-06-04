import ROOT

from histogram import HistogramContainer

class Dataset(object):

    class Sample(HistogramContainer):
        def __init__(self, dataset, eventClass):
            HistogramContainer.__init__(self, dataset.name + '_' + eventClass)
            self.dataset = dataset
            self.tree = None
    
        def loadTree(self, inputDir):
            self._source = ROOT.TFile.Open(inputDir + '/' + self.name + '.root')
            self.tree = self._source.Get('eventList')

        def releaseTree(self):
            self._source.Close()
            del self._source
            self.tree = None

    REALDATA = 0
    FULLSIM = 1
    FASTSIM = 2

    def __init__(self, name, inputNames, type, Leff, sigmaRelErr, eventClasses):
        self.name = name
        self.inputNames = inputNames
        self.type = type
        self.Leff = Leff
        self.sigmaRelErr = sigmaRelErr

        self.samples = {}
        for c in eventClasses:
            sample = Dataset.Sample(self, c)
            self.samples[c] = sample
            setattr(self, c, sample)
        
