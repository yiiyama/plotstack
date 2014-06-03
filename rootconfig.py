import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = 2000
ROOT.gStyle.SetOptStat(0)

rootlogon = ROOT.gEnv.GetValue("Rint.Logon", "")
if rootlogon:
    ROOT.gROOT.Macro(rootlogon)

ROOT.gROOT.LoadMacro("ROOT/EventProcessor.cc+")
ROOT.gROOT.LoadMacro('ROOT/histogram.h+')
ROOT.gROOT.LoadMacro("ROOT/PlotMaker.cc+")
