#!/usr/bin/env python

import os
import math

import ROOT

import plotflags
from stack import Group, Stack

class MakerWrapper(object):
    def __init__(self, core, tree):
        self._core = core
        self._core.eventList = tree

    def run(self, tree, histograms, lumi):
        for iH in range(len(histograms)):
            self._core.setHistogram(iH, histograms[iH]._core)

        self._core.input = tree
        if lumi > 0.:
            self._core.Lnorm = lumi
        else:
            self._core.Lnorm = 1.

        w = ROOT.Double()
        e2 = ROOT.Double()
        self._core.run(w, e2)

        return w, e2
    

def fillPlots(config, hdefs, eventListDir, outputPath):
    outputFile = ROOT.TFile.Open(outputPath, 'recreate')
    
    for hdef in hdefs:
        outputFile.mkdir(hdef.name)

    try:
        obsSamples = next(s for g, s in config.sampleSets if g.category == Group.OBSERVED and s[0].dataset.realData)
        integratedLumi = max([s.dataset.Leff for s in obsSamples])
    except StopIteration:
        integratedLumi = -1.

    groupWeights = {}
    
    for group, samples in config.sampleSets:
        groupWeight = 0.
        groupWeightErr2 = 0.

        applyMask = group.category == Group.OBSERVED and plotflags.HIDESENSITIVE

        for sample in samples:
            print sample.name
            
            sample.bookHistograms(hdefs, outputFile)

            sample.loadTree(eventListDir)

            try:
                plotMaker = MakerWrapper(config.specialPlotMakers[sample], sample.tree)
            except KeyError:
                try:
                    plotMaker = MakerWrapper(config.specialPlotMakers[group], sample.tree)
                except KeyError:
                    plotMaker = MakerWrapper(config.plotMaker, sample.tree)

            weight, weightErr2 = plotMaker.run(sample.tree, sample.histograms, integratedLumi)

            sample.releaseTree()

            groupWeight += weight
            groupWeightErr2 += weightErr2

            sample.postFill(applyMask)

            if config.floatFixer:
                if group not in config.floatFixer.floating:
                    sample.setHistogramErrors()
                    
                config.floatFixer.setHistogram(group, sample.histograms)

        groupWeights[group] = (groupWeight, groupWeightErr2)


    if config.floatFixer:
        config.floatFixer.calculate()

        for group, scale in config.floatFixer.scales.items():
            for sample in next(s for g, s in config.sampleSets if g == group):
                sample.scaleHistograms(scale)
    
            weight, weightErr2 = groupWeights[group]
            groupWeights[group] = (weight * scale, weightErr2 * scale)

    outputFile.cd()
    outputFile.Write()

    ROOT.TObjString('L = ' + str(integratedLumi)).Write()

    # release the histograms from the samples before closing the file
    for group, samples in config.sampleSets:
        for sample in samples:
            sample.histograms = []
    
    outputFile.Close()

    sumObs = 0.
    err2Obs = 0.
    sumBkg = 0.
    err2Bkg = 0.

    print 'Sum of weights'

    for group, samples in config.sampleSets:
        weight, weightErr2 = groupWeights[group]
        if group.category == Group.OBSERVED:
            spacer = ' +'
        elif group.category == Group.BACKGROUND:
            spacer = ' -'
        else:
            spacer = '  '
            
        print spacer + group.name + ':', '%.3e' % weight, '+-', '%.3e' % math.sqrt(weightErr2)

        if group.category == Group.OBSERVED:
            sumObs += weight
            err2Obs += weightErr2
        if group.category == Group.BACKGROUND:
            sumBkg += weight
            err2Bkg += weightErr2

    print "Obs total:", '%.3e' % sumObs, "+-", '%.3e' % math.sqrt(err2Obs)
    print "Bkg total:", '%.3e' % sumBkg, "+-", '%.3e' % math.sqrt(err2Bkg)


def makeStack(config, hdefs, inputPath, plotsDir):
    source = ROOT.TFile.Open(inputPath)

    keys = source.GetListOfKeys()
    for key in keys:
        if key.GetName().startswith('L = '):
            integratedLumi = float(key.GetName().split()[2])

    paves = []
    cmsPave = ROOT.TPaveText()
    cmsPave.SetY2NDC(0.04)
    cmsPave.SetTextFont(62)
    cmsPave.SetTextSize(0.03)
    cmsPave.SetTextAlign(13)
    cmsPave.SetBorderSize(0)
    cmsPave.SetFillStyle(0)
    paves.append(cmsPave)
    if integratedLumi > 0.: 
        cmsPave.AddText('CMS Preliminary 2014')
        lumiPave = ROOT.TPaveText()
        lumiPave.SetY2NDC(0.06)
        lumiPave.SetTextFont(62)
        lumiPave.SetTextSize(0.03)
        lumiPave.SetTextAlign(13)
        lumiPave.SetBorderSize(0)
        lumiPave.SetFillStyle(0)
        lumiPave.AddText('#int#it{L}dt = %.1f pb^{-1}' % integratedLumi)
        paves.append(lumiPave)
        arbitraryUnit = False
    else:
        cmsPave.AddText('CMS Simulation 2014')
        arbitraryUnit = True

    for group, samples in config.sampleSets:
        for sample in samples:
            sample.loadHistograms(hdefs, source)

    for iH in range(len(hdefs)):
        hdef = hdefs[iH]

        stack = Stack(hdef)

        for group, samples in config.sampleSets:
            stack.addGroup(group)
            for sample in samples:
                stack.add(group, sample.getHistogram(hdef.name).hWeighted)

        stack.draw(plotsDir, texts = paves, arbitraryUnit = arbitraryUnit, maskObserved = plotflags.HIDESENSITIVE)
    

if __name__ == '__main__':

    import sys
    import os
    from optparse import OptionParser

    import rootconfig
    import locations
    exec('from ' + locations.config + ' import hdefs, stackConfigs')

    parser = OptionParser(usage = 'Usage: makeStack.py [options] stackType')

    parser.add_option('-s', '--only-stack', action = 'store_true', dest = 'onlyStack', help = 'Only make stack from existing plots file.')

    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.print_usage()
        sys.exit(1)

    stackType = args[0]
    config = stackConfigs[stackType]

    rootFileName = locations.plotsOutputDir + '/' + stackType + '.root'

    if not options.onlyStack:
        fillPlots(config, hdefs, locations.eventListDir, rootFileName)

    try:
        os.makedirs(locations.plotsDir + '/' + stackType)
    except OSError:
        pass

    makeStack(config, hdefs, rootFileName, locations.plotsDir + '/' + stackType)

