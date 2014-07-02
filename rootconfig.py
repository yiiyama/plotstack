import os
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

rootlogon = ROOT.gEnv.GetValue("Rint.Logon", "")
if rootlogon:
    ROOT.gROOT.Macro(rootlogon)

thisdir = os.path.dirname(os.path.abspath(__file__))

ROOT.gROOT.LoadMacro(thisdir + '/ROOT/EventProcessor.cc+')
ROOT.gROOT.LoadMacro(thisdir + '/ROOT/histogram.h+')
ROOT.gROOT.LoadMacro(thisdir + '/ROOT/PlotMaker.cc+')
