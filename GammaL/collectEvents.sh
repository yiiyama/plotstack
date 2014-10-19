#!/bin/bash

CWD=$HOME/src/GammaL/plotstack/GammaL

if [ $# -eq 0 ]; then
    while read line; do
        [[ $line =~ ^[0-9]+ ]] || continue
        echo $line
        read -ra arr <<< "$line"
        bsub -q 8nh "$CWD/collectEvents.sh ${arr[3]} ${arr[0]} ${arr[1]} ${arr[2]}"
    done < $HOME/output/GammaL/main/eventList.dat
else
    export CMSSW_BASE=$HOME/cmssw/SLC6Ntuplizer5314
    export X509_USER_PROXY=$HOME/tmp/x509up_u51268

    RUN=$2
    LUMI=$3
    EVENT=$4

    cd $CMSSW_BASE
    eval `scram runtime -sh`

    cd $TMPDIR
    cmsRun $HOME/src/GammaL/plotstack/GammaL/collectEvents.py inputFiles="root://xrootd.ba.infn.it/$1" outputFile=$TMPDIR/${RUN}_${LUMI}_${EVENT}.root event=${RUN}:${LUMI}:${EVENT}
    scp $TMPDIR/${RUN}_${LUMI}_${EVENT}.root ncmu40:/store/signalEvents/
fi
