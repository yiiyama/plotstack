#!/usr/bin/env python

def weightEvents(dataset, processorClass, weightCalc, inputDir, outputDir):
    processor = processorClass(dataset.name, dataset.Leff, dataset.sigmaRelErr * dataset.sigmaRelErr, outputDir)

    for name in dataset.inputNames:
        processor.addInput(inputDir + '/' + name + '(:?_[0-9]+|).root')

    for c in dataset.samples.keys():
        processor.setFilter(c, weightCalc[c])

    processor.book()
        
    processor.process()

    processor.write()


if __name__ == '__main__':

    import sys
    import os
    import time
    import subprocess

    import rootconfig
    import locations
    exec('from ' + locations.config + ' import eventProcessor, datasets, weightCalc')

    try:
        datasetNames = sys.argv[1:]
    except:
        print 'Usage: weightEvents.py dataset'
        sys.exit(1)

    try:
        outputDir = os.environ['TMPDIR'] + '/WE' + str(int(time.time()))
        os.mkdir(outputDir)
    except KeyError:
        outputDir = locations.eventListDir

    for name in datasetNames:
        dataset = datasets[name]
        weightEvents(dataset, eventProcessor, weightCalc[name], locations.sourceDir, outputDir)

    if outputDir != locations.eventListDir:
        targetDir = locations.eventListDir
        if targetDir.startswith('rooth://'):
            targetDir = targetDir.replace('rooth://', '')
            targetDir = targetDir[:targetDir.find('/')] + ':./' + targetDir[targetDir.find('/') + 1:]

        for file in os.listdir(outputDir):
            proc = subprocess.Popen(['scp', outputDir + '/' + file, targetDir + '/'])
            while proc.poll() is None: time.sleep(1)

