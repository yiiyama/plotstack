#!/usr/bin/env python

def weightEvents(dataset, processorClass, weightCalc, inputDir, outputDir, eventClasses = []):
    processor = processorClass(dataset.name, dataset.type, dataset.Leff, dataset.sigmaRelErr * dataset.sigmaRelErr, outputDir)

    for name in dataset.inputNames:
        processor.addInput(inputDir + '/' + name + '(:?_[0-9]+|).root')

    for eventClass, sample in dataset.samples.items():
        if len(eventClasses) and eventClass not in eventClasses: continue
        processor.setOutput(eventClass, weightCalc[sample])

    processor.book()
        
    processor.process()

    processor.write()


if __name__ == '__main__':

    import sys
    import os
    import time
    import subprocess
    import re
    from optparse import OptionParser

    import rootconfig
    import locations
    exec('from ' + locations.analysis + '.config import eventProcessor, datasets, weightCalc')

    parser = OptionParser(usage = 'Usage:\n weightEvents.py [options] datasets\n weightEvetns.py -s STACK')

    parser.add_option('-s', '--stack', dest = 'stack', help = 'Process all datasets used in STACK', default = '', metavar = 'STACK')
    parser.add_option('-e', '--event-class', dest = 'eventClasses', help = 'Comma separated list of event classes to process', default = '', metavar = 'CLASSES')

    options, args = parser.parse_args()

    try:
        if options.stack:
            exec('from ' + locations.analysis + '.config import stackConfigs')
            dsList = []
            try:
                stackConfig = stackConfigs[options.stack]
            except KeyError:
                print 'No stack', options.stack, 'defined'
                raise

            for group in stackConfig.groups:
                for sample in group.samples:
                    dsList.append(sample.dataset)
    
            dsList = list(set(dsList))
        else:
            dsList = []
            for name in args:
                if '*' in name:
                    for dsName, ds in datasets.items():
                        if re.match(name.replace('*', '.*'), dsName):
                            dsList.append(ds)
                else:
                    dsList.append(datasets[name])

        if len(dsList) == 0:
            raise Exception()

    except:
        parser.print_usage()
        sys.exit(1)

    try:
        outputDir = os.environ['TMPDIR'] + '/WE' + str(int(time.time()))
        os.mkdir(outputDir)
    except KeyError:
        outputDir = locations.eventListDir

    if options.eventClasses.strip():
        eventClasses = options.eventClasses.strip().split(',')
    else:
        eventClasses = []

    for ds in dsList:
        weightEvents(ds, eventProcessor, weightCalc, locations.sourceDir, outputDir, eventClasses = eventClasses)

    if outputDir != locations.eventListDir:
        targetDir = locations.eventListDir
        if targetDir.startswith('rooth://'):
            targetDir = targetDir.replace('rooth://', '')
            host = targetDir[:targetDir.find('/')]
            remotePath = targetDir[targetDir.find('/') + 1:]
            targetDir = host + ':' + remotePath

        for file in os.listdir(outputDir):
            proc = subprocess.Popen(['scp', outputDir + '/' + file, targetDir + '/'])
            while proc.poll() is None: time.sleep(1)

