import ROOT

ROOT.gSystem.Load('libRooFit.so')

ROOT.gStyle.SetOptStat(0)

canvas = ROOT.TCanvas('c1', 'c1')
canvas.SetLogy()

pt = ROOT.RooRealVar('pt', 'P_{T}', 0., 1000., 'GeV')
wgt = ROOT.RooRealVar('wgt', 'wgt', 0., 1000.)
c1bkg = ROOT.RooRealVar('c1bkg', 'c1bkg', -0.1, -10., 0.)
c2bkg = ROOT.RooRealVar('c2bkg', 'c2bkg', -0.01, -10., 0.)
fbkg = ROOT.RooRealVar('fbkg', 'fbkg', 0.5, 0., 1.)
c1obs = ROOT.RooRealVar('c1obs', 'c1obs', -0.1, -10., 0.)
c2obs = ROOT.RooRealVar('c2obs', 'c2obs', -0.01, -10., 0.)
fobs = ROOT.RooRealVar('fobs', 'fobs', 0.5, 0., 1.)
bkgPdf = ROOT.RooGenericPdf('bkgPdf', 'bkgPdf', '@3 * TMath::Exp(@1 * @0) + (1. - @3) * TMath::Exp(@2 * @0)', ROOT.RooArgList(pt, c1bkg, c2bkg, fbkg))
obsPdf = ROOT.RooGenericPdf('obsPdf', 'obsPdf', '@3 * TMath::Exp(@1 * @0) + (1. - @3) * TMath::Exp(@2 * @0)', ROOT.RooArgList(pt, c1obs, c2obs, fobs))

pt.setRange('fitRange', 50., 400.)
pt.setRange('plotRange', 0., 400.)
pt.setRange('normRange', 50., 70.)
pt.setBins(80, 'plotRange')

ptset = ROOT.RooArgSet(pt)

bkgFit = ROOT.TF1('bkgFit', '[3] * ([2] * TMath::Exp([0] * x) + (1. - [2]) * TMath::Exp([1] * x))', 50., 400.)
obsFit = ROOT.TF1('obsFit', '[3] * ([2] * TMath::Exp([0] * x) + (1. - [2]) * TMath::Exp([1] * x))', 50., 400.)

outputFile = ROOT.TFile.Open('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/dimuonPt.root', 'recreate')

# FILL BACKGROUND DATASET

bkgDataset = ROOT.RooDataSet('bkgDataset', 'bkgDataset', ROOT.RooArgSet(pt, wgt), 'wgt')

bkgTree = ROOT.TChain('eventList')
bkgTree.Add('rooth://ncmu40//store/glweighted/DataM_ElePhotonAndDimuon.root')
bkgTree.Add('rooth://ncmu40//store/glweighted/DataM_FakePhotonAndDimuon.root')
bkgTree.Add('rooth://ncmu40//store/glweighted/WWGJets_PhotonAndDimuon.root')
bkgTree.Add('rooth://ncmu40//store/glweighted/WW_PhotonAndDimuon.root')
bkgTree.Add('rooth://ncmu40//store/glweighted/TTGJets_PhotonAndDimuon.root')
bkgTree.Add('rooth://ncmu40//store/glweighted/TTJetsFullLept_PhotonAndDimuon.root')

bkgTree.SetEstimate(bkgTree.GetEntries() + 1)
nEntries = bkgTree.Draw('muon.pt:eventSigma*puWeight*effScale*19712.', 'met < 70. && muon.pt > 50.', 'goff')

ptArr = bkgTree.GetV1()
wgtArr = bkgTree.GetV2()
for iE in range(nEntries):
    if ptArr[iE] > 1000.: continue
    pt.setVal(ptArr[iE])
    bkgDataset.add(ptset, wgtArr[iE])

zgTree = ROOT.TChain('eventList')
zgTree.Add('rooth://ncmu40//store/glweighted/ZGToLLG_PtG-5-130_PhotonAndDimuon.root')
zgTree.Add('rooth://ncmu40//store/glweighted/ZGToLLG_PtG-130_PhotonAndDimuon.root')

scale = ROOT.TFile.Open('/afs/cern.ch/user/y/yiiyama/output/GammaL/main/FloatingVGammaM.root').Get('TemplateFitError/VGamma').GetY()[0]

nEntries = zgTree.Draw('muon.pt:eventSigma*puWeight*effScale*19712.*' + str(scale), 'met < 70. && muon.pt > 50.', 'goff')

ptArr = zgTree.GetV1()
wgtArr = zgTree.GetV2()
for iE in range(nEntries):
    if ptArr[iE] > 1000.: continue
    pt.setVal(ptArr[iE])
    bkgDataset.add(ptset, wgtArr[iE])

# FIT BACKGROUND DISTRIBUTION

bkgPdf.fitTo(bkgDataset, ROOT.RooFit.Range('fitRange'), ROOT.RooFit.SumW2Error(True))

bkgFit.SetParameters(c1bkg.getVal(), c2bkg.getVal(), fbkg.getVal(), 1.)
bkgFit.SetParameter(3, bkgDataset.sumEntries() / bkgFit.Integral(50., 400.))

outputFile.cd()
bkgFit.Write('bkgFit')

# PLOT BACKGROUND DISTRIBUTION

bkgFrame = pt.frame(ROOT.RooFit.Range('plotRange'), ROOT.RooFit.Bins(80), ROOT.RooFit.Title('P_{T}^{#mu} (#slash{E}_{T} < 70 GeV)'))

bkgDataset.plotOn(bkgFrame, ROOT.RooFit.MarkerStyle(24), ROOT.RooFit.MarkerColor(ROOT.kRed - 8), ROOT.RooFit.LineColor(ROOT.kRed - 8), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
bkgPdf.plotOn(bkgFrame, ROOT.RooFit.LineColor(ROOT.kRed))

bkgFrame.Draw()

canvas.Print('/afs/cern.ch/user/y/yiiyama/www/plots/dimuonPt/bkg.pdf')

# FILL OBSERVED DATASET

obsDataset = ROOT.RooDataSet('obsDataset', 'obsDataset', ROOT.RooArgSet(pt))

obsTree = ROOT.TChain('eventList')
obsTree.Add('rooth://ncmu40//store/glweighted/DataM_PhotonAndDimuon.root')

obsTree.SetEstimate(obsTree.GetEntries() + 1)
nEntries = obsTree.Draw('muon.pt', 'met < 70. && muon.pt > 50.', 'goff')

ptArr = obsTree.GetV1()
for iE in range(nEntries):
    if ptArr[iE] > 1000.: continue
    pt.setVal(ptArr[iE])
    obsDataset.add(ptset)

# FIT OBSERVED

obsPdf.fitTo(obsDataset, ROOT.RooFit.Range('fitRange'))

obsFit.SetParameters(c1obs.getVal(), c2obs.getVal(), fobs.getVal(), 1.)
obsFit.SetParameter(3, obsDataset.sumEntries() / obsFit.Integral(50., 400.))

outputFile.cd()
obsFit.Write('obsFit')

# PLOT OBSERVED

obsFrame = pt.frame(ROOT.RooFit.Range('plotRange'), ROOT.RooFit.Bins(80), ROOT.RooFit.Title('P_{T}^{#mu} (#slash{E}_{T} < 70 GeV)'))

obsDataset.plotOn(obsFrame, ROOT.RooFit.MarkerStyle(24), ROOT.RooFit.MarkerColor(ROOT.kBlue - 8), ROOT.RooFit.LineColor(ROOT.kBlue - 8))
obsPdf.plotOn(obsFrame, ROOT.RooFit.LineColor(ROOT.kBlue))

obsFrame.Draw()

canvas.Print('/afs/cern.ch/user/y/yiiyama/www/plots/dimuonPt/obs.pdf')

canvas = ROOT.TCanvas('c2', 'c2', 800, 400)
canvas.Divide(2, 1)
canvas.cd(1).SetLogy()

overFrame = pt.frame(ROOT.RooFit.Range('plotRange'), ROOT.RooFit.Bins(80), ROOT.RooFit.Title('P_{T}^{#mu} (#slash{E}_{T} < 70 GeV)'))

bkgDataset.plotOn(overFrame, ROOT.RooFit.MarkerStyle(24), ROOT.RooFit.MarkerColor(ROOT.kRed - 8), ROOT.RooFit.LineColor(ROOT.kRed - 8), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
obsDataset.plotOn(overFrame, ROOT.RooFit.MarkerStyle(24), ROOT.RooFit.MarkerColor(ROOT.kBlue - 8), ROOT.RooFit.LineColor(ROOT.kBlue - 8), ROOT.RooFit.DrawOption('SAME EP'))

overFrame.Draw()

bkgFit.SetParameter(3, bkgFit.GetParameter(3) * 5.) # 400 / 80
obsFit.SetParameter(3, obsFit.GetParameter(3) * 5.) # 400 / 80

bkgFit.SetLineColor(ROOT.kRed)
obsFit.SetLineColor(ROOT.kBlue)

bkgFit.Draw('SAME')
obsFit.Draw('SAME')

legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.SetFillStyle(0)
legend.SetBorderSize(0)

legend.AddEntry(bkgFit, 'Background', 'L')
legend.AddEntry(obsFit, 'Data', 'L')

legend.Draw()

canvas.cd(2).SetGrid()
ratio = ROOT.TF1('ratio', '([3] * ([2] * TMath::Exp([0] * x) + (1. - [2]) * TMath::Exp([1] * x))) / ([7] * ([6] * TMath::Exp([4] * x) + (1. - [6]) * TMath::Exp([5] * x)))', 0., 400.)
for iP in range(4):
    ratio.SetParameter(iP, obsFit.GetParameter(iP))
    ratio.SetParameter(iP + 4, bkgFit.GetParameter(iP))

ratio.SetLineColor(ROOT.kBlack)
ratio.Draw()

ratio.GetHistogram().SetTitle('Data / MC ratio;P_{T} (GeV)')
ratio.GetYaxis().SetRangeUser(0., 2.)

canvas.Print('/afs/cern.ch/user/y/yiiyama/www/plots/dimuonPt/overlay.pdf')

outputFile.cd()
ratio.Write()

outputFile.Close()
