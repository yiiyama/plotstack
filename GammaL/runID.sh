#!/bin/bash

SAMPLE=$1

CMSSW_BASE=/afs/cern.ch/user/y/yiiyama/cmssw/SLC6Ntuplizer5314
cd $CMSSW_BASE
eval `scram runtime -sh`

cd /afs/cern.ch/user/y/yiiyama/src/GammaL/plotstack/GammaL
root -l -b -q 'idEfficiency.cc+("rooth://ncmu40/rootd/data3/simpleTree/'$SAMPLE'/simpleTree*.root", "/afs/cern.ch/user/y/yiiyama/output/GammaL/main/idEfficiency/'$SAMPLE'/input/counts.root")'
$HOME/public/dcmutools/reducer.py -o ${SAMPLE}.root HEfficiencyAdder /afs/cern.ch/user/y/yiiyama/output/GammaL/main/idEfficiency/$SAMPLE
mv /afs/cern.ch/user/y/yiiyama/output/GammaL/main/idEfficiency/$SAMPLE/${SAMPLE}.root /afs/cern.ch/user/y/yiiyama/output/GammaL/main/idEfficiency/
