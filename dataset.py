import ROOT

from histogram import HDef, Histogram

class Dataset(object):

    class Sample(object):

        COUNTER = HDef('Counter', (3, 0., 3.), xlabels = ['Rate', 'ScaleUncert', 'Entries'])
        
        def __init__(self, dataset, eventClass):
            self.name = dataset.name + '_' + eventClass
            self.dataset = dataset
            self.eventClass = eventClass
            self.tree = None
            self.counter = None
            self.histograms = {}
    
        def loadTree(self, inputDir):
            fileName = self.name
            if self.dataset.prescale != 1: fileName += '_ps%d' % self.dataset.prescale
            fileName += '.root'
            self._source = ROOT.TFile.Open(inputDir + '/' + fileName)
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

        def postFillHistograms(self, applyMask = False):
            for iBin in [1, 2]:
                self.counter.SetBinContent(iBin, self.counter.GetBinContent(iBin) * self.dataset.prescale)

            for h in self.histograms.values():
                h.scale(self.dataset.prescale)
                h.postFill(applyMask)

    REALDATA = ROOT.Dataset.kRealData
    FULLSIM = ROOT.Dataset.kFullSim
    FASTSIM = ROOT.Dataset.kFastSim

    def __init__(self, name, inputNames, dataType, Leff, sigmaRelErr, eventClasses, prescale = 1):
        self.name = name
        self.inputNames = inputNames
        self.dataType = dataType
        self.Leff = Leff
        self.sigmaRelErr = sigmaRelErr
        self.prescale = prescale
        self.entryList = ''

        self.samples = {}
        for c in eventClasses:
            sample = Dataset.Sample(self, c)
            self.samples[c] = sample
            setattr(self, c, sample)
        
    def cppObject(self):
        return ROOT.Dataset(self.name, self.dataType, self.Leff, self.sigmaRelErr, self.prescale)
        
