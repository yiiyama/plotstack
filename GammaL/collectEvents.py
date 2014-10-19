import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('event', default = '', mult = VarParsing.multiplicity.singleton, mytype = VarParsing.varType.string)
options.parseArguments()

process = cms.Process('DATA')

#files = []
#events = []
#
#with open('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/eventList.dat') as eventList:
#    for line in eventList:
#        if line.startswith('#'): continue
#        if not line.strip(): continue
#        
#        run, lumi, event, fname = line.strip().split()
#        events.append(':'.join((run, lumi, event)))
#        files.append('root://xrootd.ba.infn.it/' + fname)
#
#files = list(set(files))

process.source = cms.Source('PoolSource',
                            fileNames = cms.untracked.vstring(options.inputFiles)
)

if options.event:
    eventsToProcess = cms.untracked.VEventRange(options.event)

process.out = cms.OutputModule('PoolOutputModule',
                               fileName = cms.untracked.string(options.outputFile)
)

process.ep = cms.EndPath(process.out)
