import ROOT

from histogram import HDef, Histogram

class Dataset(object):

    class Sample(object):

        COUNTER = HDef('Counter', 'Counter', (3, 0., 3.), xlabels = ['Rate', 'ScaleUncert', 'Entries'])
        
        def __init__(self, dataset, eventClass):
            self.name = dataset.name + '_' + eventClass
            self.dataset = dataset
            self.eventClass = eventClass
            self.tree = None
            self.counter = None
            self.histograms = {}
    
        def loadTree(self, inputDir):
            self._source = ROOT.TFile.Open(inputDir + '/' + self.name + '.root')
            self.tree = self._source.Get('eventList')

        def releaseTree(self):
            self._source.Close()
            del self._source
            self.tree = None

        def bookHistograms(self, hdefs, outputDir):
            directory = outputDir.GetDirectory('counters')
            if not directory:
                directory = outputDir.mkdir('counters')
            directory.cd()

            self.counter = Dataset.Sample.COUNTER.generate(suffix = self.name)
            self.counter.Sumw2()
            
            self.histograms = {}
            for hdef in hdefs:
                directory = outputDir.GetDirectory(hdef.name)
                if not directory:
                    directory = outputDir.mkdir(hdef.name)
                directory.cd()
                   
                self.histograms[hdef.name] = Histogram(hdef, suffix = self.name)

        def loadHistograms(self, hdefs, sourceDir):
            self.counter = sourceDir.Get('counters/Counter_' + self.name)
            
            self.histograms = {}
            for hdef in hdefs:
                directory = sourceDir.GetDirectory(hdef.name)
                self.histograms[hdef.name]= Histogram(hdef, suffix = self.name, directory = directory)


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
        
