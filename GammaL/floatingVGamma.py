import sys
import math
import random
import array
import ROOT

from stack import Group, StackConfig
from histogram import Histogram

class FloatingVGammaSearch(StackConfig):

    def __init__(self, plotMaker, lepton, hdef):
        StackConfig.__init__(self, plotMaker)

        self.lepton = lepton
        self.tempDef = hdef

    def scalePlots(self, outputDir):
        groups = dict([(g.name, g) for g in self.groups])
        fakeGroups = ['EGFake', 'JGFake']

        candClass = 'PhotonAnd' + self.lepton
        jlClass = 'PhotonAndFake' + self.lepton
        fakeClasses = ['ElePhotonAnd' + self.lepton, 'FakePhotonAnd' + self.lepton]

        ### CREATE TEMPLATES FROM SAMPLE HISTOGRAMS IN GROUPS
        
        observed = self.tempDef.generate(suffix = 'observed')
        fakeHistograms = dict([(gname, Histogram(self.tempDef, suffix = gname)) for gname in fakeGroups])
        vgHistogram = Histogram(self.tempDef, suffix = 'vg')
        qcdHistogram = Histogram(self.tempDef, suffix = 'qcd')

        for sample in groups['Observed'].samples:
            observed.Add(sample.histograms[self.tempDef.name].hWeighted)
        for sample in groups['EWK'].samples:
            observed.Add(sample.histograms[self.tempDef.name].hWeighted, -1.)
        for gname in fakeGroups:
            for sample in groups[gname].samples:
                fakeHistograms[gname].add(sample.histograms[self.tempDef.name])
        for sample in groups['VGamma'].getSamples(candClass):
            vgHistogram.add(sample.histograms[self.tempDef.name])
        for sample in groups['JLFake'].samples:
            qcdHistogram.add(sample.histograms[self.tempDef.name])

        ### FIT FOR CENTRAL VALUES

        fitter = ROOT.TemplateChi2Fitter.singleton()
        
        vgScale = ROOT.RooRealVar('vg', 'vg', 1., -ROOT.RooNumber.infinity(), ROOT.RooNumber.infinity())
        qcdScale = ROOT.RooRealVar('qcd', 'qcd', 1., -ROOT.RooNumber.infinity(), ROOT.RooNumber.infinity())

        target = observed.Clone('target')
        for histo in fakeHistograms.values():
            target.Add(histo.hWeighted, -1.)

        fitter.setTarget(target)
        fitter.addTemplate(vgHistogram.hWeighted, 'vg', vgScale)
        fitter.addTemplate(qcdHistogram.hWeighted, 'qcd', qcdScale)

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

        vgCentral = vgScale.getVal()
        vgCentralErr = vgScale.getError()
        qcdCentral = qcdScale.getVal()
        qcdCentralErr = qcdScale.getError()

        ### ERROR ESTIMATION

        directory = outputDir.GetDirectory('TemplateFitError')
        if not directory:
            directory = outputDir.mkdir('TemplateFitError')
        directory.cd()

        tree = ROOT.TTree('toys', 'Toys')
        vg = array.array('d', [0.])
        qcd = array.array('d', [0.])
        targetContents = array.array('d', [0.] * self.tempDef.nx)
        targetErrors = array.array('d', [0.] * self.tempDef.nx)
        vgContents = array.array('d', [0.] * self.tempDef.nx)
        vgErrors = array.array('d', [0.] * self.tempDef.nx)

        tree.Branch('vg', vg, 'vg/D')
        tree.Branch('qcd', qcd, 'qcd/D')
        tree.Branch('targetContents', targetContents, 'content[%d]/D' % self.tempDef.nx)
        tree.Branch('targetErrors', targetErrors, 'error[%d]/D' % self.tempDef.nx)
        tree.Branch('vgContents', vgContents, 'content[%d]/D' % self.tempDef.nx)
        tree.Branch('vgErrors', vgErrors, 'error[%d]/D' % self.tempDef.nx)
            
        vgScales = []
        qcdScales = []

        def modifyTemplate(histogram, sigma, name = 'template'):
            result = histogram.hWeighted.Clone(name)
            for iX in range(1, result.GetNbinsX() + 1):
                err = max(abs(histogram.hScaleUp.GetBinContent(iX) - result.GetBinContent(iX)), abs(result.GetBinContent(iX) - histogram.hScaleDown.GetBinContent(iX)))
                result.SetBinContent(iX, result.GetBinContent(iX) + err * sigma)

            return result

        print 'Error evaluation with 1000 toys:'
        for iToy in range(1000):
            if iToy % 100 == 0:
                sys.stdout.write('\r' + str(iToy))
                sys.stdout.flush()

            egScaleVar = random.gauss(0., 1.)
            jgScaleVar = random.gauss(0., 1.)
            effScaleVar = random.gauss(0., 1.)
                
            target = observed.Clone('target')
            egTemplate = modifyTemplate(fakeHistograms['EGFake'], egScaleVar)
            jgTemplate = modifyTemplate(fakeHistograms['JGFake'], jgScaleVar)
            target.Add(egTemplate, -1.)
            target.Add(jgTemplate, -1.)

            vgTemplate = modifyTemplate(vgHistogram, effScaleVar, 'vg')

            for iX in range(self.tempDef.nx):
                targetContents[iX] = target.GetBinContent(iX + 1)
                targetErrors[iX] = target.GetBinError(iX + 1)
                vgContents[iX] = vgTemplate.GetBinContent(iX + 1)
                vgErrors[iX] = vgTemplate.GetBinError(iX + 1)

            fitter.setTarget(target)
            fitter.addTemplate(vgTemplate, 'vg', vgScale)
            fitter.addTemplate(qcdHistogram.hWeighted, 'qcd', qcdScale)

            target.Delete()
            egTemplate.Delete()
            jgTemplate.Delete()
            vgTemplate.Delete()

            if fitter.fit(-1) != 0: continue

            vgScales.append(vgScale.getVal())
            qcdScales.append(qcdScale.getVal())

            vg[0] = vgScale.getVal()
            qcd[0] = qcdScale.getVal()

            tree.Fill()

        sys.stdout.write('\n')

        directory.cd()
        tree.Write()

        for scales, central, centralErr, name in [(vgScales, vgCentral, vgCentralErr, 'VGamma'), (qcdScales, qcdCentral, qcdCentralErr, 'QCD')]:
            scales.sort()
            err = scales[int(len(scales) * 0.84)] - central
            errHigh = math.sqrt(err * err + centralErr * centralErr)
            err = central - scales[int(len(scales) * 0.16)]
            errLow = math.sqrt(err * err + centralErr * centralErr)

            if name == 'VGamma':
                vgErrHigh = errHigh
                vgErrLow = errLow
        
            graph = ROOT.TGraphAsymmErrors(1)
            graph.SetPoint(0, 0., central)
            graph.SetPointEYhigh(0, errHigh)
            graph.SetPointEYlow(0, errLow)
            directory.cd()
            graph.Write(name)

        ### SET SCALES AND ERRORS

        for sample in groups['VGamma'].getSamples(candClass):
            count = sample.counter.GetBinContent(1)
#            scaleError = max(vgErrHigh, vgErrLow)
#            scaleHigh = vgCentral + vgErrHigh
#            scaleLow = vgCentral - vgErrLow
            scaleError = vgCentralErr
            scaleHigh = vgCentral + vgCentralErr
            scaleLow = vgCentral - vgCentralErr

            for histogram in sample.histograms.values():
                histogram.hWeighted.Scale(vgCentral)
                histogram.hScaleUp.Scale(scaleHigh)
                histogram.hScaleDown.Scale(scaleLow)

            sample.counter.SetBinContent(1, count * vgCentral)
            sample.counter.SetBinContent(2, count * scaleError)

        for sample in groups['JLFake'].samples:
            for histogram in sample.histograms.values():
                histogram.hWeighted.Scale(qcdCentral)
                histogram.hScaleUp.Scale(qcdCentral)
                histogram.hScaleDown.Scale(qcdCentral)

            cont = sample.counter.GetBinContent(1)
            sample.counter.SetBinContent(1, cont * qcdCentral)

        for gname in fakeGroups:
            for sample in groups[gname].samples:
                for histogram in sample.histograms.values():
                    name = histogram.hScaleUp.GetName()
                    histogram.hScaleUp.Delete()
                    histogram.hScaleUp = histogram.hWeighted.Clone(name)
                    name = histogram.hScaleDown.GetName()
                    histogram.hScaleDown.Delete()
                    histogram.hScaleDown = histogram.hWeighted.Clone(name)
                    
                sample.counter.SetBinContent(2, 0.)

        for group in self.groups:
            if group.name != 'EWK' and group.category != Group.SIGNAL: continue

            for sample in group.samples:
                count = sample.counter.GetBinContent(1)
                countError = sample.counter.GetBinContent(2)

                if sample.eventClass == jlClass:
                    for histogram in sample.histograms.values():
                        histogram.hWeighted.Scale(-qcdCentral)
                        histogram.hScaleUp.Scale(-qcdCentral)
                        histogram.hScaleDown.Scale(-qcdCentral)

                        sample.counter.SetBinContent(1, -count * qcdCentral)
                        sample.counter.SetBinContent(2, -countError * qcdCentral)

                elif sample.eventClass != candClass:
                    for histogram in sample.histograms.values():
                        histogram.hWeighted.Scale(-1.)
                        histogram.hScaleUp.Scale(-1.)
                        histogram.hScaleDown.Scale(-1.)

                        sample.counter.SetBinContent(1, -count)
                        sample.counter.SetBinContent(2, -countError)
