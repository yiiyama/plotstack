import ROOT

from stack import Group, StackConfig

class TemplateFit(StackConfig):
    def __init__(self, plotMaker, targetName, floatGroups, hdef):
        StackConfig.__init__(self, plotMaker)

        self.targetName = targetName
        self.floatGroups = floatGroups
        self.tempDef = hdef

    def scalePlots(self, outputDir):
        groups = dict([(g.name, g) for g in self.groups])

        scales = dict([(groupName, ROOT.RooRealVar(groupName, groupName, 0.1, 0., 10.)) for groupName in self.floatGroups])

        fitter = ROOT.TemplateChi2Fitter.singleton()

        for sample in groups[self.targetName].samples:
            hist = self.tempDef.generate()
            hist.Add(sample.histograms[self.tempDef.name].hWeighted)
            fitter.setTarget(hist)
            hist.Delete()

        for groupName in self.floatGroups:
            hist = self.tempDef.generate()
            for sample in groups[groupName].samples:
                hist.Add(sample.histograms[self.tempDef.name].hWeighted)

            fitter.addTemplate(hist, groupName, scales[groupName])

        if fitter.fit() != 0:
            raise RuntimeError('Template fit did not converge')

        for groupName in self.floatGroups:
            scale = scales[groupName].getVal()
            for sample in groups[groupName].samples:
                for histogram in sample.histograms.values():
                    histogram.hWeighted.Scale(scale)
                    histogram.hScaleUp.Scale(scale)
                    histogram.hScaleDown.Scale(scale)

                sample.counter.SetBinContent(1, sample.counter.GetBinContent(1) * scale)
