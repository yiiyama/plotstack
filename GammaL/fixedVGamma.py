import sys
import math
import random
import array
import ROOT

from stack import Group, StackConfig
from histogram import Histogram

class FixedVGammaSearch(StackConfig):

    def __init__(self, plotMaker, lepton, hdef):
        StackConfig.__init__(self, plotMaker)

        self.lepton = lepton
        self.tempDef = hdef

    def scalePlots(self, outputDir):
        groups = dict([(g.name, g) for g in self.groups])
        floatGroups = ['JGFake', 'JLFake']
        fixedGroups = ['EGFake', 'VGamma']

        candClass = 'PhotonAnd' + self.lepton
        jlClass = 'PhotonAndFake' + self.lepton

        for sample in groups['VGamma'].samples:
            for histogram in sample.histograms.values():
                histogram.hWeighted.Scale(1.415)
                histogram.hScaleUp.Scale(1.459)
                histogram.hScaleDown.Scale(1.371)

            cont = sample.counter.GetBinContent(1)
            sample.counter.SetBinContent(1, cont * 1.415)
            sample.counter.SetBinContent(2, cont * 0.044)

        ### CREATE TEMPLATES FROM SAMPLE HISTOGRAMS IN GROUPS
        
        observed = self.tempDef.generate(suffix = 'observed')
        floatHistograms = dict([(gname, self.tempDef.generate(suffix = gname)) for gname in floatGroups])
        fixedHistograms = dict([(gname, Histogram(self.tempDef, suffix = gname)) for gname in fixedGroups])

        for sample in groups['Observed'].samples:
            observed.Add(sample.histograms[self.tempDef.name].hWeighted)
        for sample in groups['EWK'].samples:
            observed.Add(sample.histograms[self.tempDef.name].hWeighted, -1.)
        for gname in fixedGroups:
            for sample in groups[gname].samples:
                fixedHistograms[gname].add(sample.histograms[self.tempDef.name])
        for gname in floatGroups:
            for sample in groups[gname].samples:
                floatHistograms[gname].Add(sample.histograms[self.tempDef.name].hWeighted)

        ### FIT FOR CENTRAL VALUES

        fitter = ROOT.TemplateChi2Fitter.singleton()
        
        scales = dict([(gname, ROOT.RooRealVar(gname, gname, 1., -ROOT.RooNumber.infinity(), ROOT.RooNumber.infinity())) for gname in floatGroups])

        target = observed.Clone('target')
        for fixed in fixedHistograms.values():
            target.Add(fixed.hWeighted, -1.)

        fitter.setTarget(target)
        for gname, histo in floatHistograms.items():
            fitter.addTemplate(histo, gname, scales[gname])

        target.Delete()

        directory = outputDir.GetDirectory('PreTemplateFit')
        if not directory:
            directory = outputDir.mkdir('PreTemplateFit')
        directory.cd()

        fitter.plot(directory)
        directory.Write()

        if fitter.fit() != 0:
            raise RuntimeError('Template fit did not converge')

        directory = outputDir.GetDirectory('PostTemplateFit')
        if not directory:
            directory = outputDir.mkdir('PostTemplateFit')
        directory.cd()

        fitter.plot(directory)
        directory.Write()

        centrals = dict([(gname, scales[gname].getVal()) for gname in floatGroups])
        errors = dict([(gname, scales[gname].getError()) for gname in floatGroups])

        ### ERROR ESTIMATION

        directory = outputDir.GetDirectory('TemplateFitError')
        if not directory:
            directory = outputDir.mkdir('TemplateFitError')
        directory.cd()

        fixedRelErrors = {}
        for gname in fixedGroups:
            fixedRelErrors[gname] = self.tempDef.generate(suffix = gname + 'Errors')
            for iX in range(1, self.tempDef.nx + 1):
                cent = fixedHistograms[gname].hWeighted.GetBinContent(iX)
                if cent == 0.: continue
                high = fixedHistograms[gname].hScaleUp.GetBinContent(iX)
                low = fixedHistograms[gname].hScaleDown.GetBinContent(iX)
                fixedRelErrors[gname].SetBinContent(iX, max(high - cent, cent - low) / cent)

            fixedRelErrors[gname].Write()

        tree = ROOT.TTree('toys', 'Toys')
        sigma = {}
        fl = {}
        for gname in fixedGroups:
            sigma[gname] = array.array('d', [0.])
            tree.Branch('sigma' + gname, sigma[gname], 'sigma' + gname + '/D')
        for gname in floatGroups:
            fl[gname] = array.array('d', [0.])
            tree.Branch(gname, fl[gname], gname + '/D')
            
        usedScales = dict([(gname, []) for gname in floatGroups])

        print 'Error evaluation with 1000 toys:'
        for iToy in range(1000):
            if iToy % 100 == 0:
                sys.stdout.write('\r' + str(iToy))
                sys.stdout.flush()
                
            target = observed.Clone('target')
            for gname, histo in fixedHistograms.items():
                s = random.gauss(0., 1.)
                sigma[gname][0] = s
                relErr = fixedRelErrors[gname]
                subt = histo.hWeighted.Clone('subt')
                for iX in range(1, self.tempDef.nx + 1):
                    scale = 1. + relErr.GetBinContent(iX) * s
                    subt.SetBinContent(iX, histo.hWeighted.GetBinContent(iX) * scale)
                    subt.SetBinError(iX, histo.hWeighted.GetBinError(iX) * scale)
                    
                target.Add(subt, -1.)
                subt.Delete()
            
            fitter.changeTarget(target)

            target.Delete()

            if fitter.fit(-1) != 0: continue

            for gname in floatGroups:
                usedScales[gname].append(scales[gname].getVal())
                fl[gname] = scales[gname].getVal()

            tree.Fill()

        sys.stdout.write('\n')

        directory.cd()
        tree.Write()

        errHigh = {}
        errLow = {}

        for gname in floatGroups:
            scaleList = usedScales[gname]
            central = centrals[gname]
            error = errors[gname]

            scaleList.sort()
            err = scaleList[int(len(scaleList) * 0.84)] - central
            errHigh[gname] = math.sqrt(err * err + error * error)
            err = central - scaleList[int(len(scaleList) * 0.16)]
            errLow[gname] = math.sqrt(err * err + error * error)

            graph = ROOT.TGraphAsymmErrors(1)
            graph.SetPoint(0, 0., central)
            graph.SetPointEYhigh(0, errHigh[gname])
            graph.SetPointEYlow(0, errLow[gname])
            directory.cd()
            graph.Write(gname)

        ### SET SCALES AND ERRORS

        for gname in floatGroups:
            for sample in groups[gname].samples:
                for histogram in sample.histograms.values():
                    histogram.hWeighted.Scale(centrals[gname])
                    histogram.hScaleUp.Scale(centrals[gname] + errHigh[gname])
                    histogram.hScaleDown.Scale(centrals[gname] - errLow[gname])

                cont = sample.counter.GetBinContent(1)
                sample.counter.SetBinContent(1, cont * centrals[gname])
                sample.counter.SetBinContent(2, cont * max(errHigh[gname], errLow[gname]))
