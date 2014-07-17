import sys
import pickle
import array
import math

import ROOT

photonPtBins = array.array('d', [40., 50., 8000.])
photonEtaBins = array.array('d', [0., 0.8, 1.4442])
electronPtBins = array.array('d', [20., 30., 40., 50., 8000.])
electronEtaBins = array.array('d', [0., 0.8, 1.442, 1.556, 2., 2.5])
muonPtBins = array.array('d', [25., 30., 35., 40., 50., 60., 90., 140., 8000.])
muonEtaBins = array.array('d', [0., 0.9, 1.2, 2.1, 2.4])

idSFOutputFile = ROOT.TFile('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/idEfficiency/scalefactors.root', 'recreate')

photon = ROOT.TH2D('photon', 'Photon ID scale factor', len(photonPtBins) - 1, photonPtBins, len(photonEtaBins) - 1, photonEtaBins)
electron = ROOT.TH2D('electron', 'Electron ID scale factor', len(electronPtBins) - 1, electronPtBins, len(electronEtaBins) - 1, electronEtaBins)
muon = ROOT.TH2D('muon', 'Muon ID scale factor', len(muonPtBins) - 1, muonPtBins, len(muonEtaBins) - 1, muonEtaBins)
muonData = ROOT.TH2D('muonData', 'Muon ID efficiency data', len(muonPtBins) - 1, muonPtBins, len(muonEtaBins) - 1, muonEtaBins)
muonMC = ROOT.TH2D('muonMC', 'Muon ID efficiency MC', len(muonPtBins) - 1, muonPtBins, len(muonEtaBins) - 1, muonEtaBins)

photonIDBaseSource = ROOT.TFile('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/idEfficiency/Photon_ID_CSEV_SF_Jan22rereco_Full2012_S10_MC_V01.root')
photonIDBase = photonIDBaseSource.Get('PhotonIDSF_LooseWP_Jan22rereco_Full2012_S10_MC_V01')
photonIDVetoSource = ROOT.TFile('/afs/cern.ch/user/y/yiiyama/output/GammaL/photonID/scaleFactorZmumug_pt.root')
photonIDVeto = photonIDVetoSource.Get('scale_nochnogsf_all')

for iP in range(photonIDVeto.GetN()):
    iX = photon.GetXaxis().FindFixBin(photonIDVeto.GetX()[iP])
    if iX == 0: continue
    if iX > photon.GetNbinsX(): break

    for iY in range(1, photon.GetNbinsY() + 1):
        bin = photon.GetBin(iX, iY)

        photon.SetBinContent(bin, photonIDVeto.GetY()[iP])
        photon.SetBinError(bin, photonIDVeto.GetErrorY(iP))

for iX in range(1, photon.GetNbinsX() + 1):
    sX = photonIDBase.GetXaxis().FindFixBin(photon.GetXaxis().GetBinLowEdge(iX))
    
    for iY in range(1, photon.GetNbinsY() + 1):
        sY = photonIDBase.GetYaxis().FindFixBin(photon.GetYaxis().GetBinCenter(iY))
        
        bin = photon.GetBin(iX, iY)
        sbin = photonIDBase.GetBin(sX, sY)

        sf = photon.GetBinContent(bin) * photonIDBase.GetBinContent(sbin)
        relErrVeto = photon.GetBinError(bin) / photon.GetBinContent(bin)
        relErrBase = photonIDBase.GetBinError(sbin) / photonIDBase.GetBinContent(sbin)

        photon.SetBinContent(bin, sf)
        photon.SetBinError(bin, sf * math.sqrt(relErrVeto * relErrVeto + relErrBase * relErrBase))


photonIDBaseSource.Close()
photonIDVetoSource.Close()

#[20., 30., 40., 50., 8000.] x [0., 0.8, 1.442, 1.556, 2., 2.5]
electronSFs = [
    [(0.986, 0.002), (1.002, 0.001), (1.005, 0.001), (1.004, 0.001)],
    [(0.959, 0.003), (0.980, 0.001), (0.988, 0.001), (0.988, 0.002)],
    [(0.967, 0.007), (0.950, 0.007), (0.958, 0.005), (0.966, 0.009)],
    [(0.941, 0.005), (0.967, 0.003), (0.992, 0.002), (1.000, 0.003)],
    [(1.020, 0.003), (1.021, 0.003), (1.019, 0.002), (1.022, 0.004)]
]
electronSFSystematics = [
    [1.40, 0.28, 0.14, 0.41],
    [1.40, 0.28, 0.14, 0.41],
    [5.70, 2.40, 0.28, 0.43],
    [2.20, 0.59, 0.30, 0.53],
    [2.20, 0.59, 0.30, 0.53]
]

for iEta in range(len(electronSFs)):
    row = electronSFs[iEta]
    for iPt in range(len(row)):
        bin = electron.GetBin(iPt + 1, iEta + 1)
        sf, err = row[iPt]
        syst = electronSFSystematics[iEta][iPt] * 0.01
        electron.SetBinContent(bin, sf)
        electron.SetBinError(bin, math.sqrt(err * err + syst * syst))

with open('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/idEfficiency/MuonEfficiencies_Run2012ReReco_53X.pkl') as source:
    muonID = pickle.load(source)['Tight']

with open('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/idEfficiency/MuonEfficiencies_ISO_Run_2012ReReco_53X.pkl') as source:
    muonISO = pickle.load(source)['combRelIsoPF04dBeta<012_Tight']

for iPt in range(len(muonPtBins) - 1):
    if muonPtBins[iPt] < 140.:
        ptbinName = str(int(muonPtBins[iPt])) + '_' + str(int(muonPtBins[iPt + 1]))
    else:
        ptbinName = '140_300'

    for iEta in range(len(muonEtaBins) - 1):
        if muonEtaBins[iEta] < 0.5:
            etabinName = 'ptabseta<0.9'
        else:
            etabinName = 'ptabseta%.1f-%.1f' % (muonEtaBins[iEta], muonEtaBins[iEta + 1])

        idset = tuple([muonID[etabinName][ptbinName]['data/mc'][k] for k in ['efficiency_ratio', 'err_low', 'err_hi']])
        isoset = tuple([muonISO[etabinName][ptbinName]['data/mc'][k] for k in ['efficiency_ratio', 'err_low', 'err_hi']])

        idRelErr = max(idset[1:]) / idset[0]
        isoRelErr = max(isoset[1:]) / isoset[0]

        muon.SetBinContent(iPt + 1, iEta + 1, idset[0] * isoset[0])
        muon.SetBinError(iPt + 1, iEta + 1, idset[0] * isoset[0] * math.sqrt(math.pow(idRelErr, 2.) + math.pow(isoRelErr, 2.)))

        idset = tuple([muonID[etabinName][ptbinName]['data'][k] for k in ['efficiency', 'err_low', 'err_hi']])
        isoset = tuple([muonISO[etabinName][ptbinName]['data'][k] for k in ['efficiency', 'err_low', 'err_hi']])

        idRelErr = max(idset[1:]) / idset[0]
        isoRelErr = max(isoset[1:]) / isoset[0]

        muonData.SetBinContent(iPt + 1, iEta + 1, idset[0] * isoset[0])
        muonData.SetBinError(iPt + 1, iEta + 1, idset[0] * isoset[0] * math.sqrt(math.pow(idRelErr, 2.) + math.pow(isoRelErr, 2.)))

        idset = tuple([muonID[etabinName][ptbinName]['mc'][k] for k in ['efficiency', 'err_low', 'err_hi']])
        isoset = tuple([muonISO[etabinName][ptbinName]['mc'][k] for k in ['efficiency', 'err_low', 'err_hi']])

        idRelErr = max(idset[1:]) / idset[0]
        isoRelErr = max(isoset[1:]) / isoset[0]

        muonMC.SetBinContent(iPt + 1, iEta + 1, idset[0] * isoset[0])
        muonMC.SetBinError(iPt + 1, iEta + 1, idset[0] * isoset[0] * math.sqrt(math.pow(idRelErr, 2.) + math.pow(isoRelErr, 2.)))

idSFOutputFile.cd()
idSFOutputFile.Write()
idSFOutputFile.Close()

samples = ['WGToLNuG', 'ZGToLLG', 'TTGJets']

for sample in samples:
    effOutputFile = ROOT.TFile.Open('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/idEfficiency/' + sample + '.root', 'recreate')
    
    source = ROOT.TFile.Open('rooth://ncmu40//store/idEfficiency/' + sample + '.root')
    tree = source.Get('idTree')

    photon_eff = ROOT.TProfile2D('photon_eff', 'Photon ID efficiency', len(photonPtBins) - 1, photonPtBins, len(photonEtaBins) - 1, photonEtaBins)
    electron_eff = ROOT.TProfile2D('electron_eff', 'Electron ID efficiency', len(electronPtBins) - 1, electronPtBins, len(electronEtaBins) - 1, electronEtaBins)
    muon_eff = ROOT.TProfile2D('muon_eff', 'Muon ID efficiency', len(muonPtBins) - 1, muonPtBins, len(muonEtaBins) - 1, muonEtaBins)

    tree.Draw('recoId>0:TMath::Abs(eta):pt>>photon_eff', 'weight * (pdgId == 22 && genIso < 5.)', 'prof goff')
    tree.Draw('recoId>1:TMath::Abs(eta):pt>>electron_eff', 'weight * (TMath::Abs(pdgId) == 11)', 'prof goff')
    tree.Draw('recoId>1:TMath::Abs(eta):pt>>muon_eff', 'weight * (TMath::Abs(pdgId) == 13)', 'prof goff')

    photon_eff.SetDirectory(effOutputFile)
    electron_eff.SetDirectory(effOutputFile)
    muon_eff.SetDirectory(effOutputFile)

    source.Close()

    effOutputFile.cd()
    effOutputFile.Write()
    effOutputFile.Close()


hltSFOutputFile = ROOT.TFile('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/hltEfficiency/scalefactors.root', 'recreate')

photonDataSource = ROOT.TFile('/afs/cern.ch/user/y/yiiyama/output/GammaL/triggers/photon/data.root')
#photonDataL1 = photonDataSource.Get('SingleEG22_etaNVtx_eff')
#photonDataHLT = photonDataSource.Get('Ph36IdIso_etaNVtx_eff')
photonDataL1 = photonDataSource.Get('SingleEG22_ptEta_eff')
photonDataHLT = photonDataSource.Get('Ph36IdIso_ptEta_eff')

photonMCSource = ROOT.TFile('/afs/cern.ch/user/y/yiiyama/output/GammaL/triggers/photon/dye.root')
#photonMCL1 = photonMCSource.Get('SingleEG22_etaNVtx_eff')
#photonMCHLT = photonMCSource.Get('Ph36IdIso_etaNVtx_eff')
photonMCL1 = photonMCSource.Get('SingleEG22_ptEta_eff')
photonMCHLT = photonMCSource.Get('Ph36IdIso_ptEta_eff')

hltSFOutputFile.cd()
photon_e = photonDataHLT.Clone('photon_e')

for iX in range(1, photon_e.GetNbinsX() + 1):
    for iY in range(1, photon_e.GetNbinsY() + 1):
        bin = photon_e.GetBin(iX, iY)

        sf = photonDataL1.GetBinContent(bin) * photonDataHLT.GetBinContent(bin) / photonMCL1.GetBinContent(bin) / photonMCHLT.GetBinContent(bin)
        relErr2 = 0.
        relErr = photonDataL1.GetBinError(bin) / photonDataL1.GetBinContent(bin)
        relErr2 += relErr * relErr
        relErr = photonDataHLT.GetBinError(bin) / photonDataHLT.GetBinContent(bin)
        relErr2 += relErr * relErr
        relErr = photonMCL1.GetBinError(bin) / photonMCL1.GetBinContent(bin)
        relErr2 += relErr * relErr
        relErr = photonMCHLT.GetBinError(bin) / photonMCHLT.GetBinContent(bin)
        relErr2 += relErr * relErr

        photon_e.SetBinContent(bin, sf)
        photon_e.SetBinError(bin, sf * math.sqrt(relErr2))

photonDataSource.Close()
photonMCSource.Close()

photonDataSource = ROOT.TFile('/afs/cern.ch/user/y/yiiyama/output/GammaL/triggers/mueg_ph/data.root')
photonDataHLT = photonDataSource.Get('Ph22CaloIdL_inclusive_eff')

photonMCSource = ROOT.TFile('/afs/cern.ch/user/y/yiiyama/output/GammaL/triggers/mueg_ph/zg.root')
photonMCHLT = photonMCSource.Get('Ph22CaloIdL_inclusive_eff')

photon_mu = ROOT.TGraphErrors(1)
photon_mu.SetName('photon_mu')

sf = photonDataHLT.GetY()[0] / photonMCHLT.GetY()[0]
relErr2 = 0.
relErr = photonDataHLT.GetErrorY(0) / photonDataHLT.GetY()[0]
relErr2 += relErr * relErr
relErr = photonMCHLT.GetErrorY(0) / photonMCHLT.GetY()[0]
relErr2 += relErr * relErr

photon_mu.SetPoint(0, 0., sf)
photon_mu.SetPointError(0, 0., sf * math.sqrt(relErr2))

photonDataSource.Close()
photonMCSource.Close()

hltSFOutputFile.cd()
photon_mu.Write()

electronDataSource = ROOT.TFile('/afs/cern.ch/user/y/yiiyama/output/GammaL/triggers/electron/data.root')
#electronDataHLT = electronDataSource.Get('Ph22IdIso_etaNVtx_eff')
electronDataHLT = electronDataSource.Get('Ph22IdIso_ptEta_eff')

electronMCSource = ROOT.TFile('/afs/cern.ch/user/y/yiiyama/output/GammaL/triggers/electron/dy.root')
#electronMCHLT = electronMCSource.Get('Ph22IdIso_etaNVtx_eff')
electronMCHLT = electronMCSource.Get('Ph22IdIso_ptEta_eff')

hltSFOutputFile.cd()
electron = electronDataHLT.Clone('electron')

for iX in range(1, electron.GetNbinsX() + 1):
    for iY in range(1, electron.GetNbinsY() + 1):
        bin = electron.GetBin(iX, iY)

        sf = electronDataHLT.GetBinContent(bin) / electronMCHLT.GetBinContent(bin)
        relErr2 = 0.
        try:
            relErr = electronDataHLT.GetBinError(bin) / electronDataHLT.GetBinContent(bin)
        except ZeroDivisionError:
            relErr = 0.
        relErr2 += relErr * relErr
        relErr = electronMCHLT.GetBinError(bin) / electronMCHLT.GetBinContent(bin)
        relErr2 += relErr * relErr

        electron.SetBinContent(bin, sf)
        electron.SetBinError(bin, sf * math.sqrt(relErr2))

electronDataSource.Close()
electronMCSource.Close()

muonDataSource = ROOT.TFile('/afs/cern.ch/user/y/yiiyama/output/GammaL/triggers/mueg_mu/data.root')
muonDataHLT = muonDataSource.Get('Mu22_ptEta_eff')

muonMCSource = ROOT.TFile('/afs/cern.ch/user/y/yiiyama/output/GammaL/triggers/mueg_mu/zg.root')
muonMCHLT = muonMCSource.Get('Mu22_ptEta_eff')

hltSFOutputFile.cd()
muon = muonDataHLT.Clone('muon')

for iX in range(1, muon.GetNbinsX() + 1):
    for iY in range(1, muon.GetNbinsY() + 1):
        bin = muon.GetBin(iX, iY)

        sf = muonDataHLT.GetBinContent(bin) / muonMCHLT.GetBinContent(bin)
        relErr2 = 0.
        relErr = muonDataHLT.GetBinError(bin) / muonDataHLT.GetBinContent(bin)
        relErr2 += relErr * relErr
        relErr = muonMCHLT.GetBinError(bin) / muonMCHLT.GetBinContent(bin)
        relErr2 += relErr * relErr

        muon.SetBinContent(bin, sf)
        muon.SetBinError(bin, sf * math.sqrt(relErr2))

muonDataSource.Close()
muonMCSource.Close()

muegDataSource = ROOT.TFile('/afs/cern.ch/user/y/yiiyama/output/GammaL/triggers/mueg_cross/data.root')
muegDataHLT = muegDataSource.Get('Mu3p5EG12_inclusive_eff')

muegMCSource = ROOT.TFile('/afs/cern.ch/user/y/yiiyama/output/GammaL/triggers/mueg_cross/zg.root')
muegMCHLT = muegMCSource.Get('Mu3p5EG12_inclusive_eff')

mueg = ROOT.TGraphErrors(1)
mueg.SetName('mueg')

sf = muegDataHLT.GetY()[0] / muegMCHLT.GetY()[0]
relErr2 = 0.
relErr = muegDataHLT.GetErrorY(0) / muegDataHLT.GetY()[0]
relErr2 += relErr * relErr
relErr = muegMCHLT.GetErrorY(0) / muegMCHLT.GetY()[0]
relErr2 += relErr * relErr

mueg.SetPoint(0, 0., sf)
mueg.SetPointError(0, 0., sf * math.sqrt(relErr2))

muegDataSource.Close()
muegMCSource.Close()

hltSFOutputFile.cd()
mueg.Write()
hltSFOutputFile.Write()
hltSFOutputFile.Close()


samples = [('WGToLNuG', 'wg'), ('ZGToLLG', 'zg'), ('TTGJets', 'ttg')]
#plots = [('photon_e', 'photon', ('SingleEG22', 'Ph36IdIso'), 'etaNVtx'), ('photon_mu', 'mueg_ph', 'Ph22CaloIdL', 'etaNVtx'), ('electron', 'electron', 'Ph22IdIso', 'etaNVtx'), ('muon', 'mueg_mu', 'Mu22', 'ptEta')]
plots = [('photon_e', 'photon', ('SingleEG22', 'Ph36IdIso'), 'ptEta'), ('photon_mu', 'mueg_ph', 'Ph22CaloIdL', 'inclusive'), ('electron', 'electron', 'Ph22IdIso', 'ptEta'), ('muon', 'mueg_mu', 'Mu22', 'ptEta')]

for fullname, shortname in samples:
    effOutputFile = ROOT.TFile('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/hltEfficiency/' + fullname + '.root', 'recreate')
    for obj, measName, filterName, var in plots:
        source = ROOT.TFile('/afs/cern.ch/user/y/yiiyama/output/GammaL/triggers/' + measName + '/' + shortname + '.root')
        eff = None
        if type(filterName) is tuple:
            if var == 'inclusive':
                raise RuntimeError('Cannot multiply graphs')
                
            for filt in filterName:
                if not eff:
                    eff = source.Get(filt + '_' + var + '_eff')
                else:
                    eff.Multiply(source.Get(filt + '_' + var + '_eff'))
        else:
            eff = source.Get(filterName + '_' + var + '_eff')

        effOutputFile.cd()
        eff.Write(obj + '_eff')

        source.Close()

    effOutputFile.Close()
