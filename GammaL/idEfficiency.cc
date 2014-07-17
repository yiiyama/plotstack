#include "../../CommonCode/ObjectVars.h"
#include "../../CommonCode/SimpleEventProducer.h"
#include "../../CommonCode/Utilities.h"

#include "TString.h"
#include "TChain.h"
#include "TVector2.h"
#include "TFile.h"
#include "TH2D.h"

#include <cmath>
#include <iostream>
#include <algorithm>

bool BATCHMODE(true);

class EfficiencyCalculator {
public:
  EfficiencyCalculator();
  ~EfficiencyCalculator();
  bool initialize(char const*);
  void addInput(char const*);
  bool run();
  void clearInput();
  bool finalize();

private:
  TChain eventTree;
  TChain objTree;

  TTree* output;
  short pdgId;
  float pt;
  float eta;
  float phi;
  unsigned char nVtx;
  unsigned char recoId;
  float weight;
  float genPt;
  float genIso;
};

EfficiencyCalculator::EfficiencyCalculator() :
  eventTree("eventVars"),
  objTree("allObjects"),
  output(0)
{
}

EfficiencyCalculator::~EfficiencyCalculator()
{
  delete output;
}

bool
EfficiencyCalculator::initialize(char const* _outputName)
{
  TString outputName(_outputName);
  if(!outputName.Contains(".root")) outputName += "idEfficiency.root";

  TFile* outputFile(TFile::Open(outputName, "recreate"));
  if(!outputFile) return false;

  output = new TTree("idTree", "ID efficiency tree");
  output->Branch("pdgId", &pdgId, "pdgId/S");
  output->Branch("pt", &pt, "pt/F");
  output->Branch("eta", &eta, "eta/F");
  output->Branch("phi", &phi, "phi/F");
  output->Branch("nVtx", &nVtx, "nVtx/b");
  output->Branch("recoId", &recoId, "recoId/b");
  output->Branch("weight", &weight, "weight/F");
  output->Branch("genPt", &genPt, "genPt/F");
  output->Branch("genIso", &genIso, "genIso/F");

  return true;
}

void
EfficiencyCalculator::addInput(char const* _sourceName)
{
  eventTree.Add(_sourceName);
  objTree.Add(_sourceName);
}

bool
EfficiencyCalculator::run()
{
  eventTree.SetBranchStatus("*", 0);
  eventTree.SetBranchStatus("puWeight", 1);
  eventTree.SetBranchStatus("gen.*", 1);
  eventTree.SetBranchStatus("pf.size", 1);
  eventTree.SetBranchStatus("pf.pdgId", 1);
  eventTree.SetBranchStatus("pf.pt", 1);
  eventTree.SetBranchStatus("pf.eta", 1);
  eventTree.SetBranchStatus("pf.phi", 1);
  eventTree.SetBranchStatus("pf.vx", 1);
  eventTree.SetBranchStatus("pf.vy", 1);
  eventTree.SetBranchStatus("pf.vz", 1);

  objTree.SetBranchStatus("*", 0);
  objTree.SetBranchStatus("photon.size", 1);
  objTree.SetBranchStatus("photon.iSubdet", 1);
  objTree.SetBranchStatus("photon.isLoosePix", 1);
  objTree.SetBranchStatus("photon.isMediumPix", 1);
  objTree.SetBranchStatus("photon.isTightPix", 1);
  objTree.SetBranchStatus("photon.superClusterIndex", 1);
  objTree.SetBranchStatus("photon.pt", 1);
  objTree.SetBranchStatus("photon.eta", 1);
  objTree.SetBranchStatus("photon.phi", 1);
  objTree.SetBranchStatus("photon.calo*", 1);
  objTree.SetBranchStatus("electron.size", 1);
  objTree.SetBranchStatus("electron.iSubdet", 1);
  objTree.SetBranchStatus("electron.isLoose", 1);
  objTree.SetBranchStatus("electron.isMedium", 1);
  objTree.SetBranchStatus("electron.isTight", 1);
  objTree.SetBranchStatus("electron.superClusterIndex", 1);
  objTree.SetBranchStatus("electron.pt", 1);
  objTree.SetBranchStatus("electron.eta", 1);
  objTree.SetBranchStatus("electron.phi", 1);
  objTree.SetBranchStatus("electron.calo*", 1);
  objTree.SetBranchStatus("muon.size", 1);
  objTree.SetBranchStatus("muon.iSubdet", 1);
  objTree.SetBranchStatus("muon.isLoose", 1);
  objTree.SetBranchStatus("muon.isTight", 1);
  objTree.SetBranchStatus("muon.pt", 1);
  objTree.SetBranchStatus("muon.eta", 1);
  objTree.SetBranchStatus("muon.phi", 1);
  objTree.SetBranchStatus("vertex.*", 1);

  susy::SimpleEventProducer::EventVars eventVars;
  susy::PhotonVarsArray photons;
  susy::ElectronVarsArray electrons;
  susy::MuonVarsArray muons;
  susy::VertexVarsArray vertices;

  eventVars.setAddress(eventTree);
  photons.setAddress(objTree);
  electrons.setAddress(objTree);
  muons.setAddress(objTree);
  vertices.setAddress(objTree);

  long iEntry(0);
  while(eventTree.GetEntry(iEntry++)){
    if(iEntry == 500000) break;

    objTree.GetEntry(iEntry - 1);
    if(iEntry % 100000 == 1){
      if(BATCHMODE) std::cout << iEntry << std::endl;
      else (std::cout << "\r" << iEntry).flush();
    }

    unsigned ph_match[susy::NMAX];
    unsigned el_match[susy::NMAX];
    unsigned mu_match[susy::NMAX];
    double genIsoVals[susy::NMAXGEN];

    susy::genMatch(eventVars, &photons, &electrons, &muons, ph_match, el_match, mu_match, genIsoVals);

    weight = eventVars.puWeight;

    nVtx = 0;
    for(unsigned iV(0); iV != vertices.size; ++iV){
      if(vertices.isGood[iV]){
	if(nVtx == 0){
	  unsigned iG(0);
	  while(eventVars.gen_pdgId[iG] == 2212) ++iG;
	  if((TVector2(vertices.x[iV], vertices.y[iV]) - TVector2(eventVars.gen_vx[iG], eventVars.gen_vy[iG])).Mod() > 0.1) break;
	  if(std::abs(vertices.z[iV] - eventVars.gen_vz[iG]) > 0.2) break;
	}
	nVtx += 1;
      }
    }
    if(nVtx == 0) continue;

    for(unsigned iP(0); iP != photons.size; ++iP){
      unsigned match(ph_match[iP]);
      if(match > eventVars.gen_size) continue;
      if(eventVars.gen_pdgId[match] != 22) continue;
      if(photons.iSubdet[iP] != 0) continue;
      if(photons.pt[iP] < 40.) break;

      susy::PhotonVars photon(photons.at(iP));
      TVector3 caloPosition(photon.caloX, photon.caloY, photon.caloZ);

      recoId = 0;

      if(photon.isLoosePix){
        unsigned iM(0);
        for(; iM != muons.size; ++iM){
          if(muons.pt[iM] < 2.) continue;
          if(susy::deltaR(muons.eta[iM], muons.phi[iM], photon.eta, photon.phi) < 0.3) break;
        }
        if(iM == muons.size){
          unsigned iE(0);
          for(; iE != electrons.size; ++iE){
            if(electrons.pt[iE] < 2.) continue;
            if(electrons.superClusterIndex[iE] == photon.superClusterIndex) break;
            if(TVector3(electrons.caloX[iE], electrons.caloY[iE], electrons.caloZ[iE]).DeltaR(caloPosition) < 0.02) break;
          }
          if(iE == electrons.size){
            unsigned iPF(0);
            for(; iPF != eventVars.pf_size; ++iPF){
              if(eventVars.pf_pt[iPF] < 3.) continue;
              if(std::abs(eventVars.pf_pdgId[iPF]) != 211) continue;

              TVector3 dir(caloPosition);
              dir -= TVector3(eventVars.pf_vx[iPF], eventVars.pf_vy[iPF], eventVars.pf_vz[iPF]);
              if(std::abs(dir.Eta() - eventVars.pf_eta[iPF]) < 0.005 &&
                 std::abs(TVector2::Phi_mpi_pi(dir.Phi() - eventVars.pf_phi[iPF])) < 0.02) break;
            }
            if(iPF == eventVars.pf_size){
              if(photon.isTightPix) recoId = 3;
              else if(photon.isMediumPix) recoId = 2;
              else recoId = 1;
            }
          }
        }
      }

      pdgId = eventVars.gen_pdgId[match];
      pt = photon.pt;
      eta = caloPosition.Eta();
      phi = caloPosition.Phi();
      genPt = eventVars.gen_pt[match];
      genIso = genIsoVals[match];

      output->Fill();
    }

    for(unsigned iE(0); iE != electrons.size; ++iE){
      unsigned match(el_match[iE]);
      if(match > eventVars.gen_size) continue;
      if(electrons.iSubdet[iE] == -1) continue;
      if(electrons.pt[iE] < 25.) break;

      susy::ElectronVars electron(electrons.at(iE));

      TVector3 caloPosition(electron.caloX, electron.caloY, electron.caloZ);

      recoId = 0;
      if(electron.isTight) recoId = 3;
      else if(electron.isMedium) recoId = 2;
      else if(electron.isLoose) recoId = 1;

      pdgId = eventVars.gen_pdgId[match];
      pt = electron.pt;
      eta = caloPosition.Eta();
      phi = caloPosition.Phi();
      genPt = eventVars.gen_pt[match];
      genIso = genIsoVals[match];

      output->Fill();
    }

    for(unsigned iM(0); iM != muons.size; ++iM){
      unsigned match(mu_match[iM]);
      if(match > eventVars.gen_size) continue;
      if(muons.iSubdet[iM] == -1) continue;
      if(muons.pt[iM] < 25.) break;

      susy::MuonVars muon(muons.at(iM));

      recoId = 0;
      if(muon.isTight) recoId = 2;
      else if(muon.isLoose) recoId = 1;

      pdgId = eventVars.gen_pdgId[match];
      pt = muon.pt;
      eta = muon.eta;
      phi = muon.phi;
      genPt = eventVars.gen_pt[match];
      genIso = genIsoVals[match];

      output->Fill();
    }
  }

  if(!BATCHMODE) std::cout << std::endl;

  return true;
}

void
EfficiencyCalculator::clearInput()
{
  eventTree.Reset();
  objTree.Reset();
}

bool
EfficiencyCalculator::finalize()
{
  TFile* outputFile(output->GetCurrentFile());
  outputFile->cd();
  outputFile->Write();
  delete outputFile;

  output = 0;

  return true;
}

void
idEfficiency(TString const& _sourceName, TString const& _outputName)
{
  BATCHMODE = false;

  EfficiencyCalculator calc;
  if(!calc.initialize(_outputName)) return;
  calc.addInput(_sourceName);
  if(!calc.run()) return;
  calc.clearInput();
  calc.finalize();
}
