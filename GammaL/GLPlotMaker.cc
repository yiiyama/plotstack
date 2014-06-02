#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "../../CommonCode/ObjectTree.h"
#include "../../CommonCode/Utilities.h"
#include "../ROOT/PlotMaker.h"
#include "plots.h"

#include <cmath>
#include <iostream>

class GLPlotMaker : public PlotMaker {
public:
  int leptonFlavor;
  int matchPhoton;
  int matchLepton;
  bool vetoZ;

  GLPlotMaker(int f) : PlotMaker(nPlots), leptonFlavor(f), matchPhoton(0), matchLepton(0), vetoZ(leptonFlavor == 0)
  {
  }

  void matchTruePhoton() { matchPhoton = 22; }
  void vetoTruePhoton() { matchPhoton = -22; }
  void matchElePhoton() { matchPhoton = 11; }
  void vetoElePhoton() { matchPhoton = -11; }
  void matchFakePhoton() { matchPhoton = -1; }
  void matchTrueLepton() { matchLepton = 1; }
  void vetoTrueLepton() { matchLepton = -1; }
  void useZ() { vetoZ = false; }

  void run(double& _sumW, double& _sumWE2)
  {
    _sumW = 0.;
    _sumWE2 = 0.;

    /////////////////////////
    //// OPEN EVENT LIST ////
    /////////////////////////

    if(!eventList){
      std::cerr << "Input not set" << std::endl;
      return;
    }

    float met(0.);
    float metPhi(0.);
    float mt(0.);
    float puWeight(0.);
    double eventSigma(0.);
    double sigmaErr(0.);
    susy::PhotonVarsArray photons;
    bool photon_matchGen[susy::NMAX];
    bool photon_matchGenE[susy::NMAX];
    susy::ElectronVarsArray electrons;
    susy::MuonVarsArray muons;
    bool lepton_matchGen[susy::NMAX];
    susy::JetVarsArray jets;
    float mass2(0.);
    float mass3(0.);
    eventList->SetBranchAddress("met", &met);
    eventList->SetBranchAddress("metPhi", &metPhi);
    eventList->SetBranchAddress("mt", &mt);
    eventList->SetBranchAddress("puWeight", &puWeight);
    eventList->SetBranchAddress("eventSigma", &eventSigma);
    eventList->SetBranchAddress("sigmaErr", &sigmaErr);
    photons.setAddress(*eventList);
    eventList->SetBranchAddress("photon.matchGen", photon_matchGen);
    eventList->SetBranchAddress("photon.matchGenE", photon_matchGenE);
    if(leptonFlavor == 0){
      electrons.setAddress(*eventList);
      eventList->SetBranchAddress("electron.matchGen", lepton_matchGen);
    }
    else{
      muons.setAddress(*eventList);
      eventList->SetBranchAddress("muon.matchGen", lepton_matchGen);
    }
    jets.setAddress(*eventList);
    eventList->SetBranchAddress("mass2", &mass2);
    eventList->SetBranchAddress("mass3", &mass3);

#ifdef WGME
    WGammaIntegral integrator("/afs/cern.ch/user/y/yiiyama/src/GammaL/wgme/MG5/Cards/param_card.dat");
#endif

    ////////////////////
    //// FILL PLOTS ////
    ////////////////////

    try{

      double sumDW(0.);

      long iEntry(0);
      while(eventList->GetEntry(iEntry++)){
        if(vetoZ && leptonFlavor == 0 && mass2 > 86. && mass2 < 96.) continue;
        
        switch(matchPhoton){
        case 22:
          if(!photon_matchGen[0]) continue;
          break;
        case -22:
          if(!photon_matchGen[0]) continue;
          break;
        case 11:
          if(!photon_matchGenE[0]) continue;
          break;
        case -11:
          if(photon_matchGenE[0]) continue;
          break;
        case -1:
          if(photon_matchGen[0] || photon_matchGenE[0]) continue;
          break;
        }

        switch(matchLepton){
        case 1:
          if(!lepton_matchGen[0]) continue;
          break;
        case -1:
          if(lepton_matchGen[0]) continue;
          break;
        }

        double wgt(eventSigma * Lnorm * puWeight);
        double wgtErr(sigmaErr * Lnorm * puWeight);

        unsigned leptonSize(0);
        double leptonPt(0.);
        double leptonEta(0.);
        double leptonPhi(0.);
        if(leptonFlavor == 0){
          leptonSize = electrons.size;
          leptonPt = electrons.pt[0];
          leptonEta = electrons.eta[0];
          leptonPhi = electrons.phi[0];
        }
        else{
          leptonSize = muons.size;
          leptonPt = muons.pt[0];
          leptonEta = muons.eta[0];
          leptonPhi = muons.phi[0];
        }

        histograms[kNPhoton]->fill(photons.size, wgt, wgtErr);
        histograms[kNLepton]->fill(leptonSize, wgt, wgtErr);
        histograms[kNPhotonNLepton]->fill(photons.size, leptonSize, wgt, wgtErr);
        histograms[kNJet]->fill(jets.size, wgt, wgtErr);
        histograms[kMass2]->fill(mass2, wgt, wgtErr);
        histograms[kMass2Wide]->fill(mass2, wgt, wgtErr);
        histograms[kMass3]->fill(mass3, wgt, wgtErr);
        histograms[kMass3Wide]->fill(mass3, wgt, wgtErr);
        histograms[kPhotonPt]->fill(photons.pt[0], wgt, wgtErr);
        histograms[kPhotonEta]->fill(photons.eta[0], wgt, wgtErr);
        histograms[kLeptonPt]->fill(leptonPt, wgt, wgtErr);
        histograms[kLeptonEta]->fill(leptonEta, wgt, wgtErr);
        histograms[kMet]->fill(met, wgt, wgtErr);
        histograms[kMt]->fill(mt, wgt, wgtErr);
        histograms[kMetMt]->fill(met, mt, wgt, wgtErr);

        double dEta(photons.eta[0] - leptonEta);
        double dPhi(TVector2::Phi_mpi_pi(photons.phi[0] - leptonPhi));
        histograms[kDEtaPhotonLepton]->fill(dEta, wgt, wgtErr);
        histograms[kDPhiPhotonLepton]->fill(dEta, wgt, wgtErr);

        double dR(std::sqrt(dEta * dEta + dPhi * dPhi));
        histograms[kDRPhotonLepton]->fill(dR, wgt, wgtErr);

        histograms[kDPhiPhotonMet]->fill(TVector2::Phi_mpi_pi(photons.phi[0] - metPhi), wgt, wgtErr);
        histograms[kDPhiLeptonMet]->fill(TVector2::Phi_mpi_pi(leptonPhi - metPhi), wgt, wgtErr);

        if(jets.size != 0){
          double minDRGJ(100.);
          double minDRLJ(100.);
          double ht(0.);
          for(unsigned iJ(0); iJ != jets.size; ++iJ){
            double dRGJ(susy::deltaR(jets.eta[iJ], jets.phi[iJ], photons.eta[0], photons.phi[0]));
            double dRLJ(susy::deltaR(jets.eta[iJ], jets.phi[iJ], leptonEta, leptonPhi));
            if(dRGJ < minDRGJ) minDRGJ = dRGJ;
            if(dRLJ < minDRLJ) minDRLJ = dRLJ;
            if(jets.pt[iJ] > 30.) ht += jets.pt[iJ];
          }
          histograms[kDRPhotonJet]->fill(minDRGJ, wgt, wgtErr);
          histograms[kDRLeptonJet]->fill(minDRLJ, wgt, wgtErr);
          histograms[kMetOverSqrtHtLowPhotonPt]->fill(met / std::sqrt(ht), wgt, wgtErr);
        }

        TLorentzVector pl;
        pl.SetPtEtaPhiM(leptonPt, leptonEta, leptonPhi, 0.);
        TVector2 metV;
        metV.SetMagPhi(met, metPhi);
        double mu2(80.385 * 80.385 + 2. * (pl.X() * metV.X() + pl.Y() * metV.Y()));
        double mm(mu2 * mu2 - 4. * leptonPt * leptonPt * met * met);
        if(mm >= 0.){
          double PS(pl.P() * std::sqrt(mm));
          double pw(std::min(TVector3(pl.X() + metV.X(), pl.Y() + metV.Y(), pl.Z() + (pl.Z() * mu2 - PS) / 2 / pl.Perp2()).Mag(),
                             TVector3(pl.X() + metV.X(), pl.Y() + metV.Y(), pl.Z() + (pl.Z() * mu2 + PS) / 2 / pl.Perp2()).Mag()));
          histograms[kPW]->fill(pw, wgt, wgtErr);

          if(pw < 200.) histograms[kMetLowPW]->fill(met, wgt, wgtErr);
        }
        else{
          TLorentzVector pnu;
          pnu.SetPtEtaPhiM(met, 0., metPhi, 0.);
          histograms[kMW]->fill((pl + pnu).M(), wgt, wgtErr);
        }

        if(photons.pt[0] > 180.)
          histograms[kMetHighPhotonPt]->fill(met, wgt, wgtErr);

        if(photons.pt[0] < 80.){
          histograms[kMetLowPhotonPt]->fill(met, wgt, wgtErr);
          histograms[kDRPhotonLeptonLowPhotonPt]->fill(dR, wgt, wgtErr);
        }

        if(met > 120.){
          histograms[kMtHighMet]->fill(mt, wgt, wgtErr);
          histograms[kPhotonPtHighMet]->fill(photons.pt[0], wgt, wgtErr);
          histograms[kLeptonPtHighMet]->fill(leptonPt, wgt, wgtErr);
        }

        if(met < 70.){
          histograms[kDRPhotonLeptonLowMet]->fill(dR, wgt, wgtErr);
          histograms[kPhotonPtLowMet]->fill(photons.pt[0], wgt, wgtErr);
          histograms[kLeptonPtLowMet]->fill(leptonPt, wgt, wgtErr);
        }

        if(photons.pt[0] < 80. && met < 70.)
          histograms[kMtLowMetLowPhotonPt]->fill(mt, wgt, wgtErr);

        if(mt > 140.)
          histograms[kMetHighMt]->fill(met, wgt, wgtErr);

#ifdef WGME
        if(jets.size == 0 && photons.size == 1 && leptonSize == 1){
          TLorentzVector pG;
          TLorentzVector pL;

          pG.SetPtEtaPhiM(photons.pt[0], photons.eta[0], photons.phi[0], 0.);
          pL.SetPtEtaPhiM(leptonPt, leptonEta, leptonPhi, 0.);

          double nlw(-std::log(integrator.integrate(pL, pG, metV, true)));
          histograms[kNLW]->fill(nlw, wgt, wgtErr);
          if(met > 120.)
            histograms[kNLWHighMet]->fill(nlw, wgt, wgtErr);
        }
#endif

        _sumW += wgt;
        _sumWE2 += wgt * wgt;
        sumDW += wgtErr;
      }

      _sumWE2 += sumDW * sumDW;
    }
    catch(std::exception &ex){
      std::cerr << ex.what() << std::endl;
      throw;
    }
  }

};
