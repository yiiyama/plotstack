import ROOT

from stack import Group

class FloatFixer(object):
    def __init__(self, floating):
        if len(floating) == 0:
            raise RuntimeError('No floating component')
        
        self.plotIndex = -1
        self.floating = floating
        self.scales = dict([(g, 0) for g in floating])
        
    def setHistogram(self, group, histograms):
        pass

    def calculate(self):
        pass

class MatrixFixer(FloatFixer):
    def __init__(self, floating):
        FloatFixer.__init__(self, floating)
        self.calculators = dict([(g, None) for g in self.floating])
        self.targets = dict([(g, 0.) for g in self.floating])
        self.integrals = dict([(gc, dict([(gt, 0.) for gt in self.floating])) for gc in self.floating])

    def setHistogram(self, group, histograms):
        for g in self.floating:
            calc = self.calculators[g]
            norm = calc.getNorm(histograms[calc.hIndex].hWeighted)
        
            if group.category == Group.OBSERVED:
                self.targets[g] += norm
            elif group in self.floating:
                self.integrals[g][group] += norm
            else:
                self.targets[g] -= norm

    def calculate(self):
        mIntegral = ROOT.TMatrixD(len(self.floating), len(self.floating))
        # sum {mIntegral[i][j] * scale[j]} = target[i]
        for i in len(self.floating):
            for j in len(self.floating):
                mIntegral[i][j] = self.integrals[self.floating[i]][self.floating[j]]

        vTarget = ROOT.TVectorD(len(self.floating))
        for i in len(self.floating):
            vTarget[i] = self.targets[self.floating[i]]

        mIntegral.Invert()
        vTarget *= mIntegral
        for j in len(self.floating):
            self.scales[self.floating[j]] = vTarget[j]


class Chi2FitFixer(FloatFixer):

    ROOT.gSystem.Load('libRooFit.so')
    ROOT.gSystem.Load('../../Common/fitting/libCommonFitting.so')
    
    def __init__(self, floating, index):
        FloatFixer.__init__(self, floating)
        #index can be an array/tuple of indices too
        self.plotIndex = index
        self.templates = dict([(g, {}) for g in self.floating])
        self.targets = {}

    def setHistogram(self, group, histograms):
        plots = []
        try:
            for index in self.plotIndex:
                plots.append((histograms[index].hdef, histograms[index].hWeighted))
        except TypeError:
            plots.append((histograms[self.plotIndex].hdef, histograms[self.plotIndex].hWeighted))
        
        if group in self.floating:
            templates = self.templates[group]
            for hdef, plot in plots:
                try:
                    template = templates[hdef.name]
                except KeyError:
                    template = hdef.generate('Template_' + group.name)
                    template.SetDirectory(0)
                    templates[hdef.name] = template
                    
                template.Add(plot)
                
        else:
            for hdef, plot in plots:
                try:
                    target = self.targets[hdef.name]
                except KeyError:
                    target = hdef.generate('Target')
                    target.SetDirectory(0)
                    self.targets[hdef.name] = target
                
                if group.category == Group.OBSERVED:
                    factor = 1.
                else:
                    factor = -1.

                for name, plot in plots:
                    target.Add(plot, factor)
            
    def calculate(self):
        fitter = ROOT.TemplateChi2Fitter.singleton()

        targs = ROOT.TObjArray()
        for name in sorted(self.targets.keys()):
            targs.Add(self.targets[name])

        fitter.setTarget(targs)
        
        for g in self.floating:
            temps = ROOT.TObjArray()
            for name in sorted(self.templates[g].keys()):
                temps.Add(self.templates[g][name])
                
            fitter.addTemplate(temps, g.name)

        fitter.fit()

        for g in self.floating:
            self.scales[g] = fitter.getScale(g.name)
