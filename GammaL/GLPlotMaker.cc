#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "../../CommonCode/ObjectTree.h"
#include "../../CommonCode/Utilities.h"
#include "../ROOT/PlotMaker.h"

#include <cmath>
#include <iostream>

class GLPlotMaker : public PlotMaker {
public:
  int leptonFlavor;
  int matchPhoton;
  int matchLepton;
  bool vetoZ;
  bool dilepton;
  bool effCorrection;

  GLPlotMaker(int f) :
    PlotMaker(),
    leptonFlavor(f),
    matchPhoton(0),
    matchLepton(0),
    vetoZ(leptonFlavor == 0),
    dilepton(false),
    effCorrection(true)
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
  void dileptonOnly() { dilepton = true; }
  void noEffCorrection() { effCorrection = false; }

  void run()
  {
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
    double effScale(1.);
    double scaleErr(0.);
    susy::PhotonVarsArray photons;
    bool photon_matchGen[susy::NMAX];
    bool photon_matchGenE[susy::NMAX];
    unsigned lepton_size;
    float lepton_pt[susy::NMAX];
    float lepton_eta[susy::NMAX];
    float lepton_phi[susy::NMAX];
    float lepton_px[susy::NMAX];
    float lepton_py[susy::NMAX];
    float lepton_pz[susy::NMAX];
    float lepton_energy[susy::NMAX];
    bool lepton_matchGen[susy::NMAX];
    susy::JetVarsArray jets;
    float mass2(0.);
    float mass3(0.);
    unsigned char nVtx(0);
    eventList->SetBranchAddress("met", &met);
    eventList->SetBranchAddress("metPhi", &metPhi);
    eventList->SetBranchAddress("mt", &mt);
    eventList->SetBranchAddress("puWeight", &puWeight);
    eventList->SetBranchAddress("eventSigma", &eventSigma);
    eventList->SetBranchAddress("sigmaErr", &sigmaErr);
    eventList->SetBranchAddress("effScale", &effScale);
    eventList->SetBranchAddress("scaleErr", &scaleErr);
    photons.setAddress(*eventList);
    eventList->SetBranchAddress("photon.matchGen", photon_matchGen);
    eventList->SetBranchAddress("photon.matchGenE", photon_matchGenE);
    if(leptonFlavor == 0){
      eventList->SetBranchAddress("electron.size", &lepton_size);
      eventList->SetBranchAddress("electron.pt", lepton_pt);
      eventList->SetBranchAddress("electron.eta", lepton_eta);
      eventList->SetBranchAddress("electron.phi", lepton_phi);
      eventList->SetBranchAddress("electron.px", lepton_px);
      eventList->SetBranchAddress("electron.py", lepton_py);
      eventList->SetBranchAddress("electron.pz", lepton_pz);
      eventList->SetBranchAddress("electron.energy", lepton_energy);
      eventList->SetBranchAddress("electron.matchGen", lepton_matchGen);
    }
    else{
      eventList->SetBranchAddress("muon.size", &lepton_size);
      eventList->SetBranchAddress("muon.pt", lepton_pt);
      eventList->SetBranchAddress("muon.eta", lepton_eta);
      eventList->SetBranchAddress("muon.phi", lepton_phi);
      eventList->SetBranchAddress("muon.px", lepton_px);
      eventList->SetBranchAddress("muon.py", lepton_py);
      eventList->SetBranchAddress("muon.pz", lepton_pz);
      eventList->SetBranchAddress("muon.energy", lepton_energy);
      eventList->SetBranchAddress("muon.matchGen", lepton_matchGen);
      eventList->SetBranchAddress("muon.matchGen", lepton_matchGen);
    }
    jets.setAddress(*eventList);
    eventList->SetBranchAddress("mass2", &mass2);
    eventList->SetBranchAddress("mass3", &mass3);
    eventList->SetBranchAddress("nVtx", &nVtx);

#ifdef WGME
    WGammaIntegral integrator("/afs/cern.ch/user/y/yiiyama/src/GammaL/wgme/MG5/Cards/param_card.dat");
#endif

    ////////////////////
    //// FILL PLOTS ////
    ////////////////////

    try{

      long iEntry(0);
      while(eventList->GetEntry(iEntry++)){
        if(dilepton && lepton_size < 2) continue;

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

        if(dilepton){
          switch(matchLepton){
          case 1:
            if(!lepton_matchGen[1]) continue;
            break;
          case -1:
            if(lepton_matchGen[1]) continue;
            break;
          }
        }

        eventWeight = eventSigma * Lnorm * puWeight;
        double wgtRelErr2(sigmaErr * sigmaErr / eventSigma / eventSigma);
        if(effCorrection){
          eventWeight *= effScale;
          wgtRelErr2 += scaleErr * scaleErr / effScale / effScale;
        }
        eventWeightErr = eventWeight * std::sqrt(wgtRelErr2);

        fill("Mass2", mass2);
        fill("Mass2Wide", mass2);

        if(vetoZ && leptonFlavor == 0 && mass2 > 86. && mass2 < 96.) continue;

        countEvent();

        fill("NPhoton", photons.size);
        fill("NLepton", lepton_size);
        fill("NPhotonNLepton", photons.size, lepton_size);
        fill("Mass3", mass3);
        fill("Mass3Wide", mass3);
        fill("PhotonPt", photons.pt[0]);
        fill("PhotonEta", photons.eta[0]);
        fill("LeptonPt", lepton_pt[0]);
        fill("LeptonEta", lepton_eta[0]);
        fill("Met", met);
        fill("Mt", mt);
        fill("MetMt", met, mt);
        fill("NVtx", nVtx);

        double mll(0.);
        if(lepton_size == 2){
          mll = (TLorentzVector(lepton_px[0], lepton_py[0], lepton_pz[0], lepton_energy[0]) +
                 TLorentzVector(lepton_px[1], lepton_py[1], lepton_pz[1], lepton_energy[1])).M();
          fill("Mll", mll);
        }

        double dEtaGL(photons.eta[0] - lepton_eta[0]);
        fill("DEtaPhotonLepton", dEtaGL);

        double dPhiGL(TVector2::Phi_mpi_pi(photons.phi[0] - lepton_phi[0]));
        fill("DPhiPhotonLepton", dPhiGL);

        if(mll > 81. && mll < 101.)
          fill("DPhiPhotonLeptonOnZ", dPhiGL);

        double dRGL(std::sqrt(dEtaGL * dEtaGL + dPhiGL * dPhiGL));
        fill("DRPhotonLepton", dRGL);

        double dPhiGM(TVector2::Phi_mpi_pi(photons.phi[0] - metPhi));
        fill("DPhiPhotonMet", dPhiGM);

        double dPhiLM(TVector2::Phi_mpi_pi(lepton_phi[0] - metPhi));
        fill("DPhiLeptonMet", dPhiLM);

        unsigned nJet(0);
	double ht(0.);
        if(jets.size != 0){
          double minDRGJ(100.);
          int iMinGJ(-1);
          double minDRLJ(100.);
          int iMinLJ(-1);
          for(unsigned iJ(0); iJ != jets.size; ++iJ){
            double dRGJ(susy::deltaR(jets.eta[iJ], jets.phi[iJ], photons.eta[0], photons.phi[0]));
            double dRLJ(susy::deltaR(jets.eta[iJ], jets.phi[iJ], lepton_eta[0], lepton_phi[0]));
            if(dRGJ < minDRGJ){
              minDRGJ = dRGJ;
              iMinGJ = iJ;
            }
            if(dRLJ < minDRLJ){
              minDRLJ = dRLJ;
              iMinLJ = iJ;
            }

            if(jets.pt[iJ] > 30.){
              ht += jets.pt[iJ];
              ++nJet;
            }
          }

          if(iMinGJ != -1){
            fill("DPhiPhotonJet", TVector2::Phi_mpi_pi(photons.phi[0] - jets.phi[iMinGJ]));
            fill("DEtaPhotonJet", photons.eta[0] - jets.eta[iMinGJ]);
            fill("DRPhotonJet", minDRGJ);
          }
          if(iMinLJ != -1){
            fill("DPhiLeptonJet", TVector2::Phi_mpi_pi(lepton_phi[0] - jets.phi[iMinLJ]));
            fill("DEtaLeptonJet", lepton_eta[0] - jets.eta[iMinLJ]);
            fill("DRLeptonJet", minDRLJ);
          }
        }

        fill("NJet", nJet);
	fill("Ht", ht);

        TLorentzVector pl;
        pl.SetPtEtaPhiM(lepton_pt[0], lepton_eta[0], lepton_phi[0], 0.);
        TVector2 metV;
        metV.SetMagPhi(met, metPhi);
        double mu2(80.385 * 80.385 + 2. * (pl.X() * metV.X() + pl.Y() * metV.Y()));
        double mm(mu2 * mu2 - 4. * lepton_pt[0] * lepton_pt[0] * met * met);
        if(mm >= 0.){
          double PS(pl.P() * std::sqrt(mm));
          double pw(std::min(TVector3(pl.X() + metV.X(), pl.Y() + metV.Y(), pl.Z() + (pl.Z() * mu2 - PS) / 2 / pl.Perp2()).Mag(),
                             TVector3(pl.X() + metV.X(), pl.Y() + metV.Y(), pl.Z() + (pl.Z() * mu2 + PS) / 2 / pl.Perp2()).Mag()));
          fill("PW", pw);

          if(pw < 200.) fill("MetLowPW", met);
        }
        else{
          TLorentzVector pnu;
          pnu.SetPtEtaPhiM(met, 0., metPhi, 0.);
          fill("MW", (pl + pnu).M());
        }

        if(met > 120.){
          fill("PhotonPtHighMet", photons.pt[0]);
          fill("LeptonPtHighMet", lepton_pt[0]);
        }

        if(met < 70.){
          fill("PhotonPtLowMet", photons.pt[0]);
          fill("PhotonPtZoomLowMet", photons.pt[0]);
          fill("LeptonPtLowMet", lepton_pt[0]);
          fill("NJetLowMet", jets.size);

          if(lepton_size == 1){
            fill("DRPhotonLeptonLowMet1L", dRGL);
            fill("DPhiLeptonMetLowMet1L", dPhiLM);
          }
          else if(mll > 81. && mll < 101.)
            fill("PhotonPtLowMetOnZ", photons.pt[0]);

          if(dilepton){
            fill("Mass3LowMet", mass3);
            fill("MuonPtLowMet", lepton_pt[0]);
            fill("MuonPtLowMet", lepton_pt[1]);
          }
        }

        if(photons.pt[0] < 80.){
          if(mt > 100.) fill("MetHighMtLowPhotonPt", met);
          if(mt < 100.) fill("MetLowMtLowPhotonPt", met);
	  if(met < 70.){
	    fill("MtLowMetLowPhotonPt", mt);
	    fill("HtLowMetLowPhotonPt", ht);
	  }
        }
        else{
          if(mt > 100.) fill("MetHighMtHighPhotonPt", met);
	  if(met > 120.){
	    fill("MtHighMetHighPhotonPt", mt);
	    if(mt > 100.)
	      fill("HtHighMtHighMetHighPhotonPt", ht);
	  }
        }

#ifdef WGME
        if(jets.size == 0 && photons.size == 1 && leptonSize == 1){
          TLorentzVector pG;
          TLorentzVector pL;

          pG.SetPtEtaPhiM(photons.pt[0], photons.eta[0], photons.phi[0], 0.);
          pL.SetPtEtaPhiM(lepton_pt[0], lepton_eta[0], lepton_phi[0], 0.);

          double nlw(-std::log(integrator.integrate(pL, pG, metV, true)));
          fill("NLW", nlw);
          if(met > 120.)
            fill("NLWHighMet", nlw);
        }
#endif
      }
    }
    catch(std::exception &ex){
      std::cerr << ex.what() << std::endl;
      throw;
    }
  }

};
