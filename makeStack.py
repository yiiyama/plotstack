#!/usr/bin/env python

import os
import math

import ROOT
import rootconfig

import plotflags
from dataset import Dataset
from stack import Group, Stack

class MakerWrapper(object):

    def __init__(self, core):
        self._core = core

    def run(self, sample, hdefs, lumi):
        self._core.setCounter(sample.counter)

        iH = 0
        for hdef in hdefs:
            self._core.setHistogram(hdef.name, sample.histograms[hdef.name]._core)
            iH += 1

        self._core.eventList = sample.tree

        if lumi > 0.:
            self._core.Lnorm = lumi
        else:
            self._core.Lnorm = 1.

        self._core.run()


def fillPlots(config, eventListDir, outputFile, integratedLumi = -1.):

    if integratedLumi < 0.:
        try:
            observed = next(g for g in config.groups if g.category == Group.OBSERVED)
            integratedLumi = next(s.dataset.Leff for s in observed.samples if s.dataset.dataType == Dataset.REALDATA)
        except:
            pass

    outputDir = outputFile.GetDirectory('components')
    if not outputDir:
        outputDir = outputFile.mkdir('components')

    for group in config.groups:
        applyMask = group.category == Group.OBSERVED and plotflags.HIDESENSITIVE

        for sample in group.samples:
            print sample.name
            
            sample.bookHistograms(config.hdefs, outputDir)

            sample.loadTree(eventListDir)

            try:
                plotMaker = MakerWrapper(config.specialPlotMakers[sample])
            except KeyError:
                try:
                    plotMaker = MakerWrapper(config.specialPlotMakers[group])
                except KeyError:
                    plotMaker = MakerWrapper(config.plotMaker)

            plotMaker.run(sample, config.hdefs, integratedLumi)

            sample.releaseTree()

            sample.postFillHistograms(applyMask)

    outputDir.cd()
    outputDir.Write()
    ROOT.TObjString('L = ' + str(integratedLumi)).Write()


def formGroups(config, outputFile, printResult = False):

    componentsDir = outputFile.GetDirectory('components')
    for group in config.groups:
        for sample in group.samples:
            sample.loadHistograms(config.hdefs, componentsDir)
    
    config.scalePlots(outputFile)

    groupsDir = outputFile.GetDirectory('groups')
    if not groupsDir:
        groupsDir = outputFile.mkdir('groups')

    for group in config.groups:
        group.bookHistograms(config.hdefs, groupsDir)

    groupsDir.cd()
    groupsDir.Write()

    if not printResult: return

    sumObs = 0.
    err2Obs = 0.
    sumBkg = 0.
    err2Bkg = 0.

    print 'Sum of weights'

    for group in config.groups:
        if group.category == Group.OBSERVED:
            spacer = ' +'
        elif group.category == Group.BACKGROUND:
            spacer = ' -'
        else:
            spacer = '  '

        weight = group.counter.GetBinContent(1)
        statErr = group.counter.GetBinError(1)
        if weight > 0.:
            scaleErr = group.counter.GetBinContent(2) / weight
        else:
            scaelErr = 0.
        entries = int(group.counter.GetBinContent(3))
        if group.counter.GetBinContent(2) == 0.:
            print spacer + group.name + ':', '%.3e +- %.3e (%d entries)' % (weight, statErr, entries)
        else:
            print spacer + group.name + ':', '(%.3e +- %.3e) * (1 +- %.2f) (%d entries)' % (weight, statErr, scaleErr, entries)

        if group.category == Group.OBSERVED:
            sumObs += weight
            err2Obs += statErr * statErr
        if group.category == Group.BACKGROUND:
            sumBkg += weight
            err2Bkg += statErr * statErr + scaleErr * scaleErr * weight * weight

    print "Obs total:", '%.3e' % sumObs, "+-", '%.3e' % math.sqrt(err2Obs)
    print "Bkg total:", '%.3e' % sumBkg, "+-", '%.3e' % math.sqrt(err2Bkg)


def makeStack(config, source, plotsDir):
    ROOT.gErrorIgnoreLevel = 2000
    
    keys = source.GetDirectory('components').GetListOfKeys()
    for key in keys:
        if key.GetName().startswith('L = '):
            integratedLumi = float(key.GetName().split()[2])

    paves = []
    cmsPave = ROOT.TPaveText()
    cmsPave.SetY2NDC(0.045)
    cmsPave.SetTextFont(62)
    cmsPave.SetTextSize(0.03)
    cmsPave.SetTextAlign(12)
    cmsPave.SetBorderSize(0)
    cmsPave.SetFillStyle(0)
    paves.append(cmsPave)
    if integratedLumi > 0.: 
        cmsPave.AddText('CMS Preliminary 2014')
        lumiPave = ROOT.TPaveText()
        lumiPave.SetY2NDC(0.045)
        lumiPave.SetTextFont(62)
        lumiPave.SetTextSize(0.03)
        lumiPave.SetTextAlign(12)
        lumiPave.SetBorderSize(0)
        lumiPave.SetFillStyle(0)
        lumiPave.AddText('L = %.1f fb^{-1}' % (integratedLumi / 1000.))
        paves.append(lumiPave)
        arbitraryUnit = False
    else:
        cmsPave.AddText('CMS Simulation 2014')
        arbitraryUnit = True

    for group in config.groups:
        group.loadHistograms(config.hdefs, source.GetDirectory('groups'))

    for hdef in config.hdefs:
        stack = Stack(hdef)

        for group in config.groups:
            stack.addGroup(group)

        print stack.name
        stack.draw(plotsDir, texts = paves, arbitraryUnit = arbitraryUnit, maskObserved = plotflags.HIDESENSITIVE, drawEmpty = plotflags.DRAWEMPTY)
    

if __name__ == '__main__':

    import sys
    import os
    from optparse import OptionParser

    import locations
    exec('from ' + locations.analysis + '.config import stackConfigs')

    parser = OptionParser(usage = 'Usage: makeStack.py [options] stackType')

    parser.add_option('-i', '--input', dest = 'histoSource', default = '', help = 'Use pre-filled histogram file. Implies -g')
    parser.add_option('-g', '--form-group', action = 'store_true', dest = 'formGroup', help = 'Form group plots from existing component histograms.')
    parser.add_option('-s', '--only-stack', action = 'store_true', dest = 'onlyStack', help = 'Only make stack from existing group histograms.')
    parser.add_option('-f', '--set-flags', dest = 'setFlags', default = '', help = 'Set plot flags.')
    parser.add_option('-L', '--integrated-lumi', dest = 'integratedLumi', default = -1., help = 'Integrated luminosity to normalize to in 1/pb.')
    parser.add_option('-l', '--list-config', action = 'store_true', dest = 'listConfig', help = 'Print list of configurations and exit.')
    parser.add_option('-d', '--print-datasets', action = 'store_true', dest = 'printDatasets', help = 'Print list of datasets for the configuration and exit.')

    options, args = parser.parse_args()

    if options.listConfig:
        for conf in sorted(stackConfigs.keys()):
            print conf
        sys.exit(0)

    if len(args) != 1:
        parser.print_usage()
        sys.exit(1)

    if options.setFlags:
        for stmt in options.setFlags.split():
            exec('plotflags.' + stmt)
        
    stackName = args[0]
    config = stackConfigs[stackName]

    if options.printDatasets:
        dsList = []
        for group in config.groups:
            for sample in group.samples:
                dsList.append(sample.dataset)

        dsList = list(set(dsList))

        output = ''
        for dataset in dsList:
            output += dataset.name + ' '

        print output
        sys.exit(0)

    if options.histoSource:
        def copyDir(sourceDir, targetDir):
            keys = sourceDir.GetListOfKeys()
            for key in keys:
                obj = key.ReadObj()
                if obj.InheritsFrom(ROOT.TDirectory.Class()):
                    targ = targetDir.mkdir(obj.GetName())
                    copyDir(obj, targ)
                elif obj.InheritsFrom(ROOT.TTree.Class()):
                    targetDir.cd()
                    tree = obj.CloneTree(-1, 'fast')
                    tree.Write()
                else:
                    try:
                        obj.SetDirectory(targetDir)
                    except AttributeError:
                        pass
                    targetDir.cd()
                    obj.Write()
                    

        source = ROOT.TFile.Open(locations.plotsOutputDir + '/' + options.histoSource + '.root')
        if not source:
            raise IOError('Histogram input not found')

        outputFile = ROOT.TFile.Open(locations.plotsOutputDir + '/' + stackName + '.root', 'recreate')
        copyDir(source.GetDirectory('components'), outputFile.mkdir('components'))

        outputFile.Close()

        options.formGroup = True
        
    try:
        os.makedirs(locations.plotsDir + '/' + stackName)
    except OSError:
        pass

    if not options.onlyStack:
        if not options.formGroup:
            outputFile = ROOT.TFile.Open(locations.plotsOutputDir + '/' + stackName + '.root', 'recreate')
            fillPlots(config, locations.eventListDir, outputFile, float(options.integratedLumi))
            outputFile.Close()

        outputFile = ROOT.TFile.Open(locations.plotsOutputDir + '/' + stackName + '.root', 'update')
        formGroups(config, outputFile, printResult = True)
        outputFile.Close()

    source = ROOT.TFile.Open(locations.plotsOutputDir + '/' + stackName + '.root')
    makeStack(config, source, locations.plotsDir + '/' + stackName)
    source.Close()
