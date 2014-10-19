import sys
import ROOT

sys.path.append('/afs/cern.ch/user/y/yiiyama/work/ext/bin')
import das_client

outputFile = open('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/eventList.dat', 'w')

for lepton in ['E', 'M']:
    if lepton == 'E':
        lep = 'el'
        sample = 'PhotonAndElectron'
        addCut = ' && (mass2 < 81. || mass2 > 101.)'
        datasets = ['/Photon/Run2012A-22Jan2013-v1/AOD', '/DoublePhoton/Run2012B-22Jan2013-v1/AOD', '/DoublePhoton/Run2012C-22Jan2013-v2/AOD', '/DoublePhoton/Run2012D-22Jan2013-v1/AOD']
    else:
        lep = 'mu'
        sample = 'PhotonAndMuon'
        addCut = ''
        datasets = ['/MuEG/Run2012A-22Jan2013-v1/AOD', '/MuEG/Run2012B-22Jan2013-v1/AOD', '/MuEG/Run2012C-22Jan2013-v1/AOD', '/MuEG/Run2012D-22Jan2013-v1/AOD']

    source = ROOT.TFile.Open('rooth://ncmu40//store/glweighted/Data' + lepton + '_' + sample + '.root')
    tree = source.Get('eventList')

    tree.SetEstimate(tree.GetEntries() + 1)

    for pt in [('photon.pt[0] <= 80.', 'LowPt'), ('photon.pt[0] > 80.', 'HighPt')]:
        for ht in [('ht <= 100.', 'LowHt'), ('ht > 100. && ht <= 400.', 'MidHt'), ('ht > 400.', 'HighHt')]:
            for met in [('met > 120. && met <= 200.', '120'), ('met > 200. && met <= 300.', '200'), ('met > 300.', '300')]:

                channelName = lep + pt[1] + ht[1] + met[1]
                print channelName
                outputFile.write('# ' + channelName + '\n')

                cut = 'mt > 100. && ' + pt[0] + ' && ' + ht[0] + ' && ' + met[0] + addCut
                nEntries = tree.Draw('run:lumi:event', cut, 'goff')
                run = tree.GetV1()
                lumi = tree.GetV2()
                event = tree.GetV3()

                for iEntry in range(nEntries):
                    for dataset in datasets:
                        query = 'file dataset=%s run=%d lumi=%d' % (dataset, int(run[iEntry]), int(lumi[iEntry]))
                        dasData = das_client.get_data('https://cmsweb.cern.ch', query, 0, 10, 0, 300, '', '')
                        files = dasData['data'][0]['file']
                        if len(files):
                            fileName = files[0]['name']
                            break

                    outputFile.write('%d %d %d %s\n' % (int(run[iEntry]), int(lumi[iEntry]), int(event[iEntry]), fileName))

                outputFile.write('\n')

outputFile.close()
