#!/usr/bin/env python

def weightEvents(dataset, processorClass, weightCalc, inputDir, outputDir):
    processor = processorClass(dataset.name, dataset.type, dataset.Leff, dataset.sigmaRelErr * dataset.sigmaRelErr, outputDir)

    for name in dataset.inputNames:
        processor.addInput(inputDir + '/' + name + '(:?_[0-9]+|).root')

    for eventClass, sample in dataset.samples.items():
        processor.setFilter(eventClass, weightCalc[sample])

    processor.book()
        
    processor.process()

    processor.write()


if __name__ == '__main__':

    import sys
    import os
    import time
    import subprocess
    from optparse import OptionParser

    import rootconfig
    import locations
    exec('from ' + locations.analysis + '.config import eventProcessor, datasets, weightCalc')

    parser = OptionParser(usage = 'Usage:\n weightEvents.py [options] datasets\n weightEvetns.py -s STACK')

    parser.add_option('-s', '--stack', dest = 'stack', help = 'Process all datasets used in STACK', default = '', metavar = 'STACK')

    options, args = parser.parse_args()

    try:
        if options.stack:
            exec('from ' + locations.analysis + '.config import stackConfigs')
            dsList = []
            for group in stackConfigs[options.stack].groups:
                for sample, factor in group.content:
                    dsList.append(sample.dataset)
    
            dsList = list(set(dsList))
        else:
            dsList = []
            for name in args:
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

    for ds in dsList:
        weightEvents(ds, eventProcessor, weightCalc, locations.sourceDir, outputDir)

    if outputDir != locations.eventListDir:
        targetDir = locations.eventListDir
        if targetDir.startswith('rooth://'):
            targetDir = targetDir.replace('rooth://', '')
            targetDir = targetDir[:targetDir.find('/')] + ':./' + targetDir[targetDir.find('/') + 1:]

        for file in os.listdir(outputDir):
            proc = subprocess.Popen(['scp', outputDir + '/' + file, targetDir + '/'])
            while proc.poll() is None: time.sleep(1)

