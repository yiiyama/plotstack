import ROOT

from stack import Group

class FloatFixer(object):
    def __init__(self, floating):
        if len(floating) == 0:
            raise RuntimeError('No floating component')
        
        self.floating = floating
        self.scales = dict([(g, 0) for g in floating])
        
    def setHistogram(self, group):
        pass

    def calculate(self):
        pass

class MatrixFixer(FloatFixer):
    def __init__(self, floating):
        FloatFixer.__init__(self, floating)
        self.calculators = dict([(g, None) for g in self.floating])
        self.targets = dict([(g, 0.) for g in self.floating])
        self.integrals = dict([(gc, dict([(gt, 0.) for gt in self.floating])) for gc in self.floating])

    def setHistogram(self, group):
        for g in self.floating:
            calc = self.calculators[g]
            norm = calc.getNorm(group.getHistogram(calc.hName).hWeighted)
        
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
    
    def __init__(self, floating, hdefs):
        FloatFixer.__init__(self, floating)
        self.hdefs = hdefs
        self.templates = dict([(g, {}) for g in self.floating])
        self.targets = {}

    def setHistogram(self, group):
        plots = [(hdef, group.getHistogram(hdef.name).hWeighted) for hdef in self.hdefs]
        
        if group in self.floating:
            for hdef, plot in plots:
                template = hdef.generate('Template_' + group.name)
                template.SetDirectory(0)
                template.Add(plot)
                self.templates[group][hdef.name] = template
                
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
            
    def calculate(self, outputFile):
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

        fitter.plot(outputFile.mkdir('PreTemplateFit'))

        fitter.initializeScales()
        for g in self.floating:
            self.scales[g] = fitter.getScale(g.name)

        nTrial = 0

        while True:
            fitter.fit()

            nTrial += 1

            if nTrial == 10: break
    
            for g in self.floating:
                if fitter.getScale(g.name) < 0.:
                    self.scales[g] *= 0.1
                    fitter.setScale(g.name, self.scales[g])
                    break
            else:
                break

        for g in self.floating:
            self.scales[g] = fitter.getScale(g.name)

        fitter.plot(outputFile.mkdir('PostTemplateFit'))


