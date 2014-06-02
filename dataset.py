import ROOT

from histogram import Histogram

class Dataset(object):

    class Sample(object):
        def __init__(self, dataset, eventClass):
            self.name = dataset.name + '_' + eventClass
            self.dataset = dataset
            self.histograms = []
            self.tree = None
    
        def bookHistograms(self, hdefs, outputFile):
            self.histograms = []
            for hdef in hdefs:
                outputFile.cd(hdef.name)
                self.histograms.append(Histogram(hdef, suffix = self.name))

        def loadHistograms(self, hdefs, sourceFile):
            self.histograms = []
            for hdef in hdefs:
                directory = sourceFile.GetDirectory(hdef.name)
                self.histograms.append(Histogram(hdef, suffix = self.name, directory = directory))

        def getHistogram(self, name):
            return next(h for h in self.histograms if h.name == name + '_' + self.name)

        def loadTree(self, inputDir):
            self._source = ROOT.TFile.Open(inputDir + '/' + self.name + '.root')
            self.tree = self._source.Get('eventList')

        def releaseTree(self):
            self._source.Close()
            del self._source
            self.tree = None
    
        def postFill(self, applyMask):
            for h in self.histograms:
                h.postFill(applyMask)
        
        def setHistogramErrors(self):
            for h in self.histograms:
                h.setError()
    
        def scaleHistograms(self, scale):
            for h in self.histograms:
                h.scale(scale)


    def __init__(self, name, inputNames, realData, Leff, sigmaRelErr, eventClasses):
        self.name = name
        self.inputNames = inputNames
        self.realData = realData
        self.Leff = Leff
        self.sigmaRelErr = sigmaRelErr

        self.samples = {}
        for c in eventClasses:
            sample = Dataset.Sample(self, c)
            self.samples[c] = sample
            setattr(self, c, sample)
        
