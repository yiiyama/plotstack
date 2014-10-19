#include "../../CommonCode/ObjectTree.h"
#include "../../CommonCode/Utilities.h"
#include "../../CommonCode/SimpleEventProducer.h"
#include "../../CommonCode/ObjectSelector.h"

#include "TFile.h"
#include "TTree.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TString.h"
#include "TObjArray.h"

#include <cmath>
#include <iostream>
#include <utility>
#include <algorithm>
#include <stdexcept>

#include "../ROOT/EventProcessor.h"

enum OutputTypes {
  oPhotonAndElectron,
  oPhotonAndMuon,
  oElePhotonAndElectron,
  oElePhotonAndMuon,
  oFakePhotonAndElectron,
  oFakePhotonAndMuon,
  oPhotonAndFakeElectron,
  oPhotonAndFakeMuon,
  nOutputTypes
};

TString outputNames[] = {
  "SoftPhotonAndElectron",
  "SoftPhotonAndMuon",
  "SoftElePhotonAndElectron",
  "SoftElePhotonAndMuon",
  "SoftFakePhotonAndElectron",
  "SoftFakePhotonAndMuon",
  "SoftPhotonAndFakeElectron",
  "SoftPhotonAndFakeMuon"
};

// Assumptions:
// . Treat candidates independently
// . Ignore double-fakes
// We are weighting the proxy (0 fake) sample and matching it to the 1 fake sample.
// weight = N_{1 fake} / N_{0 fake}
//        = Sum_{i}[f_i * Prod_{j!=i}[1-f_j]] / Prod_{i}[1-f_i]
//        = Sum_{i}[f_i / (1-f_i)]
// where f_i is the fake rate f = P(Signal ID|Fakeable Obj) / [P(Proxy ID|FO) + P(Signal ID|FO)].
// Therefore alternatively
// weight = Sum_{i} [P_i(SID|FO) / P_i(PID|FO)]
// Since the overall weight splits to independent terms, an event with multiple FOs will be assigned multiple entries in the combo list with the corresponding weight.

class WGEventWeight : public EventWeight {
public:
  WGEventWeight() : EventWeight("WGEventWeight", 1., 0.) {}
  WGEventWeight(TString const& _name) : EventWeight(_name, 1., 0.) {}
  WGEventWeight(TString const& _name, double _w, double _r) : EventWeight(_name, _w, _r) {}
  ~WGEventWeight() {}
  virtual void setPhoton(susy::PhotonVars const&, susy::VertexVarsArray const&, susy::SimpleEventProducer::EventVars const&) {}
  virtual void setElectron(susy::ElectronVars const&, susy::VertexVarsArray const&, susy::SimpleEventProducer::EventVars const&) {}
  virtual void setMuon(susy::MuonVars const&, susy::VertexVarsArray const&, susy::SimpleEventProducer::EventVars const&) {}
};

class WGammaProcessor : public EventProcessor {
public:
  WGammaProcessor(Dataset const*, char const*);
  ~WGammaProcessor() {} // ROOT closes files in case of crashes; deleting event list causes double free

  void process();
  void write();

private:
  bool preprocess();
  void processPhotonAndLepton(unsigned);
  void processElePhotonAndLepton(unsigned);
  void processFakePhotonAndLepton(unsigned);
  void processPhotonAndFakeLepton(unsigned);
  void processFakePhotonAndFakeLepton(unsigned);

  template<typename LeptonVarsArray> void calculateKinematics(susy::PhotonVars const& _photon, LeptonVarsArray const& _leptons)
  {
    double energy(_photon.energy + _leptons.energy[0]);
    double px(_photon.px + _leptons.px[0]);
    double py(_photon.py + _leptons.py[0]);
    double pz(_photon.pz + _leptons.pz[0]);

    mass2 = std::sqrt(energy * energy - px * px - py * py - pz * pz);

    if(_leptons.size > 1){
      energy += _leptons.energy[1];
      px += _leptons.px[1];
      py += _leptons.py[1];
      pz += _leptons.pz[1];
      mass3 = std::sqrt(energy * energy - px * px - py * py - pz * pz);
    }
    else
      mass3 = -1.;

    mt = std::sqrt(2. * eventVars.met * _leptons.pt[0] * (1. - std::cos(eventVars.metPhi - _leptons.phi[0])));

    mtJESUp = std::sqrt(2. * metJESUp * _leptons.pt[0] * (1. - std::cos(metJESUpV.Phi() - _leptons.phi[0])));
    mtJESDown = std::sqrt(2. * metJESDown * _leptons.pt[0] * (1. - std::cos(metJESDownV.Phi() - _leptons.phi[0])));
  }

  void setJets(unsigned _nObj = 0, double* _eta = 0, double* _phi = 0);

  /* CONFIGURATION */

  bool useElectronFilter;
  bool useMuonFilter;

  /* INPUT / OUTPUT */

  susy::SimpleEventProducer::EventVars eventVars;
  susy::ElectronVarsArray electrons;
  susy::MuonVarsArray muons;
  susy::PhotonVarsArray photons;
  susy::JetVarsArray jets;
  susy::VertexVarsArray vertices;
  bool electron_isCand[susy::NMAX];
  bool electron_isFake[susy::NMAX];
  bool muon_isCand[susy::NMAX];
  bool muon_isFake[susy::NMAX];
  bool photon_isCand[susy::NMAX];
  bool photon_isFake[susy::NMAX];
  bool photon_isEle[susy::NMAX];
  bool jet_isCand[susy::NMAX];

  /* OUTPUT */

  float metJESUp;
  float metJESDown;
  float metPhiJESUp;
  float metPhiJESDown;
  float mt;
  float mtJESUp;
  float mtJESDown;
  float mass2;
  float mass3;
  float ht;
  float htJESUp;
  float htJESDown;
  float genBoost;
  unsigned char nVtx;
  unsigned char nJet;
  unsigned char nJetJESUp;
  unsigned char nJetJESDown;
  susy::PhotonVarsArray photonsOut;
  bool photonsOut_matchGen[susy::NMAX];
  bool photonsOut_matchGenE[susy::NMAX];
  float photonsOut_genIso[susy::NMAX];
  susy::ElectronVarsArray electronsOut;
  bool electronsOut_matchGen[susy::NMAX];
  susy::MuonVarsArray muonsOut;
  bool muonsOut_matchGen[susy::NMAX];
  susy::JetVarsArray jetsOut;

  TH1D* cutflowE_;
  TH1D* cutflowM_;

  /* UTIL */

  bool photon_matchGen[susy::NMAX];
  bool photon_matchGenE[susy::NMAX];
  float photon_genIso[susy::NMAX];
  bool electron_matchGen[susy::NMAX];
  bool electron_matchHLT[susy::NMAX];
  bool muon_matchGen[susy::NMAX];
  bool muon_match21HLT[susy::NMAX];
  bool muon_matchHLT[susy::NMAX];

  TVector2 metJESUpV;
  TVector2 metJESDownV;

  TH1* idsfTable[3]; // photon, electron, muon
  TH1* ideffTable[3];
  TH1* hltsfTable[4]; // photon_e, photon_mu, electron, muon
  TH1* hlteffTable[4];
};

WGammaProcessor::WGammaProcessor(Dataset const* _dataset, char const* _outputDir) :
  EventProcessor(nOutputTypes, _dataset, _outputDir),
  useElectronFilter(false),
  useMuonFilter(false),
  eventVars(),
  electrons(),
  muons(),
  photons(),
  jets(),
  vertices(),
  metJESUp(0.),
  metJESDown(0.),
  metPhiJESUp(0.),
  metPhiJESDown(0.),
  mt(0.),
  mtJESUp(0.),
  mtJESDown(0.),
  mass2(0.),
  mass3(0.),
  ht(0.),
  htJESUp(0.),
  htJESDown(0.),
  genBoost(0.),
  nVtx(0),
  nJet(0),
  nJetJESUp(0.),
  nJetJESDown(0.),
  photonsOut(),
  electronsOut(),
  muonsOut(),
  jetsOut(),
  cutflowE_(0),
  cutflowM_(0)
{
  std::fill_n(electron_isCand, susy::NMAX, false);
  std::fill_n(electron_isFake, susy::NMAX, false);
  std::fill_n(muon_isCand, susy::NMAX, false);
  std::fill_n(muon_isFake, susy::NMAX, false);
  std::fill_n(photon_isCand, susy::NMAX, false);
  std::fill_n(photon_isFake, susy::NMAX, false);
  std::fill_n(photon_isEle, susy::NMAX, false);
  std::fill_n(jet_isCand, susy::NMAX, false);

  std::fill_n(photonsOut_matchGen, susy::NMAX, false);
  std::fill_n(photonsOut_matchGenE, susy::NMAX, false);
  std::fill_n(photonsOut_genIso, susy::NMAX, 0.);
  std::fill_n(electronsOut_matchGen, susy::NMAX, false);
  std::fill_n(muonsOut_matchGen, susy::NMAX, false);

  std::fill_n(photon_matchGen, susy::NMAX, false);
  std::fill_n(photon_matchGenE, susy::NMAX, false);
  std::fill_n(photon_genIso, susy::NMAX, 0.);
  std::fill_n(electron_matchGen, susy::NMAX, false);
  std::fill_n(electron_matchHLT, susy::NMAX, false);
  std::fill_n(muon_matchGen, susy::NMAX, false);
  std::fill_n(muon_match21HLT, susy::NMAX, false);
  std::fill_n(muon_matchHLT, susy::NMAX, false);

  for(unsigned iH(0); iH != 3; ++iH)
    idsfTable[iH] = ideffTable[iH] = 0;

  for(unsigned iH(0); iH != 4; ++iH)
    hltsfTable[iH] = hlteffTable[iH] = 0;

  for(unsigned iF(0); iF != nOutputTypes; ++iF)
    outputIndices[outputNames[iF]] = iF;
}

void
WGammaProcessor::process()
{
  useElectronFilter =
    produceOutput[oPhotonAndElectron] ||
    produceOutput[oElePhotonAndElectron] ||
    produceOutput[oFakePhotonAndElectron] ||
    produceOutput[oPhotonAndFakeElectron];
  useMuonFilter =
    produceOutput[oPhotonAndMuon] ||
    produceOutput[oElePhotonAndMuon] ||
    produceOutput[oFakePhotonAndMuon] ||
    produceOutput[oPhotonAndFakeMuon];

  ///////////////
  //// INPUT ////
  ///////////////

  TChain filterTree("eventVars");
  TChain eventTree("eventVars");
  TChain objectTree("allObjects");
  TChain cutTree("cutTree");

  for(unsigned iP(0); iP != inputPaths.size(); ++iP){
    filterTree.Add(inputPaths[iP]);
    eventTree.Add(inputPaths[iP]);
    objectTree.Add(inputPaths[iP]);
    cutTree.Add(inputPaths[iP]);
  }

  if(filterTree.GetEntries() == 0){
    std::cerr << "Empty input" << std::endl;
    return;
  }

  filterTree.SetBranchStatus("*", 0);
  if(dataset.prescale != 1) filterTree.SetBranchStatus("eventNumber", 1);
  if(useElectronFilter){
    filterTree.SetBranchStatus("SoftPhotonAndElectron", 1);
    filterTree.SetBranchStatus("SoftElePhotonAndElectron", 1);
    filterTree.SetBranchStatus("SoftFakePhotonAndElectron", 1);
    filterTree.SetBranchStatus("SoftPhotonAndFakeElectron", 1);

    filterTree.SetBranchStatus("HLT_Ele27_WP80", 1);
  }
  filterTree.SetBranchStatus("SoftPhotonAndMuon", 1);
  filterTree.SetBranchStatus("SoftElePhotonAndMuon", 1);
  filterTree.SetBranchStatus("SoftFakePhotonAndMuon", 1);
  filterTree.SetBranchStatus("SoftPhotonAndFakeMuon", 1);
  filterTree.SetBranchStatus("HLT_IsoMu24", 1);
  filterTree.SetBranchStatus("HLT_IsoMu24_eta2p1", 1);

  unsigned eventNumber(0);
  
  bool PhotonAndElectron(false);
  bool ElePhotonAndElectron(false);
  bool FakePhotonAndElectron(false);
  bool PhotonAndFakeElectron(false);
  bool PhotonAndMuon(false);
  bool ElePhotonAndMuon(false);
  bool FakePhotonAndMuon(false);
  bool PhotonAndFakeMuon(false);

  bool Ele27WP80(false);
  bool IsoMu24(false);
  bool IsoMu24eta2p1(false);

  if(dataset.prescale != 1) filterTree.SetBranchAddress("eventNumber", &eventNumber);
  if(useElectronFilter){
    filterTree.SetBranchAddress("SoftPhotonAndElectron", &PhotonAndElectron);
    filterTree.SetBranchAddress("SoftElePhotonAndElectron", &ElePhotonAndElectron);
    filterTree.SetBranchAddress("SoftFakePhotonAndElectron", &FakePhotonAndElectron);
    filterTree.SetBranchAddress("SoftPhotonAndFakeElectron", &PhotonAndFakeElectron);

    filterTree.SetBranchAddress("HLT_Ele27_WP80", &Ele27WP80);
  }
  filterTree.SetBranchAddress("SoftPhotonAndMuon", &PhotonAndMuon);
  filterTree.SetBranchAddress("SoftElePhotonAndMuon", &ElePhotonAndMuon);
  filterTree.SetBranchAddress("SoftFakePhotonAndMuon", &FakePhotonAndMuon);
  filterTree.SetBranchAddress("SoftPhotonAndFakeMuon", &PhotonAndFakeMuon);
  filterTree.SetBranchAddress("HLT_IsoMu24", &IsoMu24);
  filterTree.SetBranchAddress("HLT_IsoMu24_eta2p1", &IsoMu24eta2p1);

  eventTree.SetBranchStatus("*", 0);
  eventTree.SetBranchStatus("runNumber", 1);
  eventTree.SetBranchStatus("lumiNumber", 1);
  eventTree.SetBranchStatus("eventNumber", 1);
  if(eventTree.GetBranch("puWeight")) eventTree.SetBranchStatus("puWeight", 1);
  else eventVars.puWeight = 1.;
  eventTree.SetBranchStatus("met", 1);
  eventTree.SetBranchStatus("metPhi", 1);
  eventTree.SetBranchStatus("rho", 1);
  if(eventTree.GetBranch("gen.size")) eventTree.SetBranchStatus("gen*", 1);

  eventVars.setAddress(eventTree);

  objectTree.SetBranchStatus("*", 0);
  objectTree.SetBranchStatus("electron*", 1);
  objectTree.SetBranchStatus("muon*", 1);
  objectTree.SetBranchStatus("photon*", 1);
  objectTree.SetBranchStatus("vertex*", 1);
  objectTree.SetBranchStatus("jet*", 1);

  electrons.setAddress(objectTree);
  objectTree.SetBranchAddress("electron.isCand", electron_isCand);
  objectTree.SetBranchAddress("electron.isFake", electron_isFake);
  objectTree.SetBranchAddress("electron.hltEle27WP80TrackIsoFilter", electron_matchHLT);
  muons.setAddress(objectTree);
  objectTree.SetBranchAddress("muon.isCand", muon_isCand);
  objectTree.SetBranchAddress("muon.isFake", muon_isFake);
  objectTree.SetBranchAddress("muon.hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15", muon_match21HLT);
  objectTree.SetBranchAddress("muon.hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10", muon_match21HLT);
  objectTree.SetBranchAddress("muon.hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15", muon_matchHLT);
  photons.setAddress(objectTree);
  objectTree.SetBranchAddress("photon.isSoftCand", photon_isCand);
  objectTree.SetBranchAddress("photon.isSoftFake", photon_isFake);
  objectTree.SetBranchAddress("photon.isSoftEle", photon_isEle);
  jets.setAddress(objectTree);
  objectTree.SetBranchAddress("jet.isCand", &jet_isCand);
  vertices.setAddress(objectTree);

  ////////////////
  //// OUTPUT ////
  ////////////////

  for(unsigned iF(0); iF != nOutputTypes; ++iF){
    if(!eventList[iF]) continue;

    eventList[iF]->Branch("run", &eventVars.runNumber, "run/i");
    eventList[iF]->Branch("lumi", &eventVars.lumiNumber, "lumi/i");
    eventList[iF]->Branch("event", &eventVars.eventNumber, "event/i");
    eventList[iF]->Branch("puWeight", &eventVars.puWeight, "puWeight/F");
    eventList[iF]->Branch("rho", &eventVars.rho, "rho/F");
    eventList[iF]->Branch("met", &eventVars.met, "met/F");
    eventList[iF]->Branch("metPhi", &eventVars.metPhi, "metPhi/F");
    eventList[iF]->Branch("metJESUp", &metJESUp, "metJESUp/F");
    eventList[iF]->Branch("metJESDown", &metJESDown, "metJESDown/F");
    eventList[iF]->Branch("metPhiJESUp", &metPhiJESUp, "metPhiJESUp/F");
    eventList[iF]->Branch("metPhiJESDown", &metPhiJESDown, "metPhiJESDown/F");
    eventList[iF]->Branch("mt", &mt, "mt/F");
    eventList[iF]->Branch("mtJESUp", &mtJESUp, "mtJESUp/F");
    eventList[iF]->Branch("mtJESDown", &mtJESDown, "mtJESDown/F");
    eventList[iF]->Branch("genBoost", &genBoost, "genBoost/F");
    eventList[iF]->Branch("mass2", &mass2, "mass2/F");
    eventList[iF]->Branch("mass3", &mass3, "mass3/F");
    eventList[iF]->Branch("ht", &ht, "ht/F");
    eventList[iF]->Branch("htJESUp", &htJESUp, "htJESUp/F");
    eventList[iF]->Branch("htJESDown", &htJESDown, "htJESDown/F");
    eventList[iF]->Branch("nVtx", &nVtx, "nVtx/b");
    eventList[iF]->Branch("nJet", &nJet, "nJet/b");
    eventList[iF]->Branch("nJetJESUp", &nJetJESUp, "nJetJESUp/b");
    eventList[iF]->Branch("nJetJESDown", &nJetJESDown, "nJetJESDown/b");
    photonsOut.setBranches(*eventList[iF]);
    eventList[iF]->Branch("photon.matchGen", photonsOut_matchGen, "matchGen[photon.size]/O");
    eventList[iF]->Branch("photon.matchGenE", photonsOut_matchGenE, "matchGenE[photon.size]/O");
    eventList[iF]->Branch("photon.genIso", photonsOut_genIso, "genIso[photon.size]/F");
    if(iF == oPhotonAndElectron || iF == oElePhotonAndElectron || iF == oFakePhotonAndElectron || iF == oPhotonAndFakeElectron){
      electronsOut.setBranches(*eventList[iF]);
      eventList[iF]->Branch("electron.matchGen", electronsOut_matchGen, "matchGen[electron.size]/O");
    }
    if(iF == oPhotonAndMuon || iF == oElePhotonAndMuon || iF == oFakePhotonAndMuon || iF == oPhotonAndFakeMuon){
      muonsOut.setBranches(*eventList[iF]);
      eventList[iF]->Branch("muon.matchGen", muonsOut_matchGen, "matchGen[muon.size]/O");
    }
    jetsOut.setBranches(*eventList[iF]);
  }

  /////////////////
  //// CUTFLOW ////
  /////////////////

  if(produceOutput[oPhotonAndElectron] || produceOutput[oPhotonAndMuon]){
    TString cuts[] = {
      "AllEvents",
      "HLT",
      "GoodLumi",
      "MetFilter",
      "GoodVertex",
      "RadiationVeto",
      "GoodPhoton",
      "GoodLepton",
      "FSRVeto"
    };
    unsigned nCuts(sizeof(cuts) / sizeof(TString));

    if(produceOutput[oPhotonAndElectron]){
      eventList[oPhotonAndElectron]->GetCurrentFile()->cd();
      cutflowE_ = new TH1D("cutflowE", "Cutflow (Electron Channel)", nCuts, 0., nCuts);
      for(unsigned iC(0); iC != nCuts; ++iC){
        cutflowE_->GetXaxis()->SetBinLabel(iC + 1, cuts[iC]);
        if(iC != nCuts - 1) cutflowE_->SetBinContent(iC + 1, cutTree.GetEntries(TString::Format("cutflowE >= %d", iC)));
      }
    }
    if(produceOutput[oPhotonAndMuon]){
      eventList[oPhotonAndMuon]->GetCurrentFile()->cd();
      cutflowM_ = new TH1D("cutflowM", "Cutflow (Muon Channel)", nCuts, 0., nCuts);
      for(unsigned iC(0); iC != nCuts; ++iC){
        cutflowM_->GetXaxis()->SetBinLabel(iC + 1, cuts[iC]);
        if(iC != nCuts - 1) cutflowM_->SetBinContent(iC + 1, cutTree.GetEntries(TString::Format("cutflowM >= %d", iC)));
      }
    }
  }

  ////////////////////
  //// START LOOP ////
  ////////////////////

  std::cout << "Processing " << dataset.name << std::endl;

  long iEntry(0);
  while(filterTree.GetEntry(iEntry++) > 0){
    try{
      if(iEntry % 1000 == 0) (std::cout << "\r" << iEntry).flush();
      if(eventNumber % dataset.prescale != 0) continue;

      if(!(Ele27WP80 &&
           ((PhotonAndElectron && produceOutput[oPhotonAndElectron]) ||
            (ElePhotonAndElectron && produceOutput[oElePhotonAndElectron]) ||
            (FakePhotonAndElectron && produceOutput[oFakePhotonAndElectron]) ||
            (PhotonAndFakeElectron && produceOutput[oPhotonAndFakeElectron])
            )) &&
         !((IsoMu24 || IsoMu24eta2p1) &&
           ((PhotonAndMuon && produceOutput[oPhotonAndMuon]) ||
            (ElePhotonAndMuon && produceOutput[oElePhotonAndMuon]) ||
            (FakePhotonAndMuon && produceOutput[oFakePhotonAndMuon]) ||
            (PhotonAndFakeMuon && produceOutput[oPhotonAndFakeMuon])
            ))) continue;

      eventTree.GetEntry(iEntry - 1);
      objectTree.GetEntry(iEntry - 1);

      if(!preprocess()) continue;

      if(PhotonAndMuon){
        if(produceOutput[oPhotonAndMuon]) processPhotonAndLepton(oPhotonAndMuon);
        continue;
      }
      else if(PhotonAndElectron){
        if(produceOutput[oPhotonAndElectron]) processPhotonAndLepton(oPhotonAndElectron);
        continue;
      }

      ///////////////////
      /// FAKE EVENTS ///
      ///////////////////
      /* Since we only consider single fakes at the moment, we never have the "second candidate". */
      /* If e.g. we had a fake photon and a candidate photon, that is a candidate event. */
      /* Consequently, the array of the faking object will always have only one element. */

      if(ElePhotonAndElectron && !ElePhotonAndMuon && produceOutput[oElePhotonAndElectron])
        processElePhotonAndLepton(oElePhotonAndElectron);

      if(ElePhotonAndMuon && produceOutput[oElePhotonAndMuon])
        processElePhotonAndLepton(oElePhotonAndMuon);
      
      if(FakePhotonAndElectron && produceOutput[oFakePhotonAndElectron])
        processFakePhotonAndLepton(oFakePhotonAndElectron);

      if(FakePhotonAndMuon && produceOutput[oFakePhotonAndMuon])
        processFakePhotonAndLepton(oFakePhotonAndMuon);

      if(PhotonAndFakeElectron && produceOutput[oPhotonAndFakeElectron])
        processPhotonAndFakeLepton(oPhotonAndFakeElectron);

      if(PhotonAndFakeMuon && produceOutput[oPhotonAndFakeMuon])
        processPhotonAndFakeLepton(oPhotonAndFakeMuon);
    }
    catch(std::exception& _ex){
      std::cerr << std::endl << _ex.what() << std::endl;
    }
  }
  std::cout << std::endl;
}

void
WGammaProcessor::write()
{
  if(cutflowE_){
    eventList[oPhotonAndElectron]->GetCurrentFile()->cd();
    cutflowE_->Write();
    delete cutflowE_;
    cutflowE_ = 0;
  }
  if(cutflowM_){
    eventList[oPhotonAndMuon]->GetCurrentFile()->cd();
    cutflowM_->Write();
    delete cutflowM_;
    cutflowM_ = 0;
  }

  EventProcessor::write();
}

bool
WGammaProcessor::preprocess()
{
  if(dataset.dataType != Dataset::kRealData){
    unsigned photonIndex[susy::NMAX];
    unsigned electronIndex[susy::NMAX];
    unsigned muonIndex[susy::NMAX];
    double genIsoVals[susy::NMAXGEN];
    susy::genMatch(eventVars, &photons, useElectronFilter ? &electrons : 0, useMuonFilter ? &muons : 0, photonIndex, useElectronFilter ? electronIndex : 0, useMuonFilter ? muonIndex : 0, genIsoVals);

    for(unsigned iP(0); iP != photons.size; ++iP){
      photon_matchGen[iP] = photon_matchGenE[iP] = false;
      photon_genIso[iP] = -1.;
      if(photonIndex[iP] < eventVars.gen_size){
        photon_genIso[iP] = genIsoVals[photonIndex[iP]];
        if(eventVars.gen_pdgId[photonIndex[iP]] == 22) photon_matchGen[iP] = true;
        else photon_matchGenE[iP] = true;
      }
    }

    if(useElectronFilter){
      for(unsigned iL(0); iL != electrons.size; ++iL)
        electron_matchGen[iL] = (electronIndex[iL] < eventVars.gen_size);
    }

    if(useMuonFilter){
      for(unsigned iL(0); iL != muons.size; ++iL)
        muon_matchGen[iL] = (muonIndex[iL] < eventVars.gen_size);
    }

    std::set<unsigned> status3; // list of status 3 particles that decayed to non-status 3
    for(unsigned iG(0); iG != eventVars.gen_size; ++iG){
      if(eventVars.gen_status[iG] != 3) continue;
      std::set<unsigned>::iterator sItr(status3.find(eventVars.gen_motherIndex[iG]));
      if(sItr != status3.end()) status3.erase(sItr);
      status3.insert(iG);
    }

    std::set<short> poi;
    for(std::set<unsigned>::iterator gItr(status3.begin()); gItr != status3.end(); ++gItr){
      short idx(*gItr);
      while(idx != -1 && (eventVars.gen_pdgId[idx] == 21 || eventVars.gen_pdgId[idx] == 2212 || std::abs(eventVars.gen_pdgId[idx]) < 6)) idx = eventVars.gen_motherIndex[idx];
      if(idx == -1) continue; // nothing of interest was in the decay chain

      // go up to initial quark / gluon splitting
      while(eventVars.gen_motherIndex[idx] != -1 && !(eventVars.gen_pdgId[eventVars.gen_motherIndex[idx]] == 21 || eventVars.gen_pdgId[eventVars.gen_motherIndex[idx]] == 2212 || std::abs(eventVars.gen_pdgId[eventVars.gen_motherIndex[idx]]) < 6)) idx = eventVars.gen_motherIndex[idx];
      if(eventVars.gen_motherIndex[idx] == -1)
        throw std::logic_error("Gen decay list hit index -1 before id 2212");

      poi.insert(idx);
    }

    TVector2 genPt;
    for(std::set<short>::iterator pItr(poi.begin()); pItr != poi.end(); ++pItr)
      genPt += TVector2(eventVars.gen_px[*pItr], eventVars.gen_py[*pItr]);

    genBoost = genPt.Mod();
  }

  nVtx = 0;
  for(unsigned iV(0); iV != vertices.size; ++iV)
    if(vertices.isGood[iV]) ++nVtx;

  // corrMet' = rawMet - sum[corrJet' - rawJet]t
  // corrMet  = rawMet - sum[corrJet  - rawJet]t
  //    diff  = sum[corrJet - corrJet']
  TVector2 jetDiff(0., 0.);
  for(unsigned iJ(0); iJ != jets.size; ++iJ){
    if(jets.pt[iJ] < 10.) continue;
    // (1 - (jec+delta)/jec)
    double diffFactor(-jets.jecUncert[iJ] / jets.jecScale[iJ]);
    jetDiff += TVector2(jets.px[iJ] * diffFactor, jets.py[iJ] * diffFactor);
  }
  TVector2 metV;
  metV.SetMagPhi(eventVars.met, eventVars.metPhi);
  metJESUpV = metV + jetDiff;
  metJESDownV = metV - jetDiff;

  metJESUp = metJESUpV.Mod();
  metJESDown = metJESDownV.Mod();

  metPhiJESUp = TVector2::Phi_mpi_pi(metJESUpV.Phi());
  metPhiJESDown = TVector2::Phi_mpi_pi(metJESDownV.Phi());

  return true;
}

void
WGammaProcessor::processPhotonAndLepton(unsigned _outputType)
{
  if(_outputType != oPhotonAndElectron && _outputType != oPhotonAndMuon)
    throw std::runtime_error("Incorrect filter type");

  effScale = 1.;
  double photonRelErr(0.);
  double electronRelErr(0.);
  double muonRelErr(0.);
  double scaleRelErr2(0.);

  photonsOut.clear();

  for(unsigned iP(0); iP != photons.size; ++iP){
    if(photons.iSubdet[iP] == -1) continue;
    if(photons.pt[iP] < 15.) break;
    if(!photon_isCand[iP]) continue;

    photonsOut_matchGen[photonsOut.size] = photon_matchGen[iP];
    photonsOut_matchGenE[photonsOut.size] = photon_matchGenE[iP];
    photonsOut_genIso[photonsOut.size] = photon_genIso[iP];

    photonsOut.push_back(photons.at(iP));
  }

  // uncomment when switching to liso trees
  if(photonsOut.size == 0) return;

  if(_outputType == oPhotonAndElectron){
    electronsOut.clear();

    for(unsigned iL(0); iL != electrons.size; ++iL){
      if(electrons.iSubdet[iL] == -1) continue;
      if(electrons.pt[iL] < 30.) break;
      if(!electron_matchHLT[iL]) continue;
      if(!electron_isCand[iL]) continue;

      electronsOut_matchGen[electronsOut.size] = electron_matchGen[iL];
      electronsOut.push_back(electrons.at(iL));
    }

    if(electronsOut.size == 0) return;
  }
  else if(_outputType == oPhotonAndMuon){
    muonsOut.clear();

    unsigned nMuons(0);

    for(unsigned iL(0); iL != muons.size; ++iL){
      if(muons.iSubdet[iL] == -1) continue;
      if(muons.pt[iL] > 10.) ++nMuons;
      if(muons.pt[iL] < 26.) continue;
      if(!muon_matchHLT[iL] && !muon_match21HLT[iL]) continue;
      if(!muon_isCand[iL]) continue;

      muonsOut_matchGen[muonsOut.size] = muon_matchGen[iL];
      muonsOut.push_back(muons.at(iL));
    }

    if(nMuons != 1 || muonsOut.size == 0) return;
  }

  scaleRelErr2 += photonRelErr * photonRelErr + electronRelErr * electronRelErr + muonRelErr * muonRelErr;
  scaleErr = effScale * std::sqrt(scaleRelErr2);

  susy::PhotonVars photon(photonsOut.at(0));

  if(_outputType == oPhotonAndElectron){
    if(susy::deltaR(photonsOut.eta[0], photonsOut.phi[0], electronsOut.eta[0], electronsOut.phi[0]) < 0.7) return;
    calculateKinematics(photon, electronsOut);

    if(cutflowE_) cutflowE_->Fill(cutflowE_->GetNbinsX() - 0.5);
  }
  else if(_outputType == oPhotonAndMuon){
    if(susy::deltaR(photonsOut.eta[0], photonsOut.phi[0], muonsOut.eta[0], muonsOut.phi[0]) < 0.7) return;
    calculateKinematics(photon, muonsOut);

    if(cutflowM_) cutflowM_->Fill(cutflowM_->GetNbinsX() - 0.5);
  }

  setJets();
  
  if(weightCalc[_outputType]){
    WGEventWeight* calc(static_cast<WGEventWeight*>(weightCalc[_outputType]));
    calc->setPhoton(photon, vertices, eventVars);
    if(_outputType == oPhotonAndElectron)
      calc->setElectron(electronsOut.at(0), vertices, eventVars);
    if(_outputType == oPhotonAndMuon)
      calc->setMuon(muonsOut.at(0), vertices, eventVars);
  }

  fill(_outputType);
}

void
WGammaProcessor::processElePhotonAndLepton(unsigned _outputType)
{
  if(_outputType == oElePhotonAndMuon){
    muonsOut.clear();

    unsigned nMuons(0);

    for(unsigned iL(0); iL != muons.size; ++iL){
      if(muons.iSubdet[iL] == -1) continue;
      if(muons.pt[iL] > 10.) ++nMuons;
      if(muons.pt[iL] < 26.) continue;
      if(!muon_matchHLT[iL] && !muon_match21HLT[iL]) continue;
      if(!muon_isCand[iL]) continue;
      muonsOut_matchGen[muonsOut.size] = muon_matchGen[iL];
      muonsOut.push_back(muons.at(iL));
    }

    if(nMuons != 1 || muonsOut.size == 0) return;
  }
  else if(_outputType != oElePhotonAndElectron)
    throw std::runtime_error("Incorrect filter type");

  for(unsigned iEP(0); iEP != photons.size; ++iEP){
    if(photons.iSubdet[iEP] == -1) continue;
    if(photons.pt[iEP] < 15.) break;
    if(!photon_isEle[iEP]) continue;

    susy::PhotonVars elePhoton(photons.at(iEP));

    TVector3 caloPosition(elePhoton.caloX, elePhoton.caloY, elePhoton.caloZ);

    if(_outputType == oElePhotonAndElectron){
      electronsOut.clear();

      for(unsigned iL(0); iL != electrons.size; ++iL){
        if(electrons.iSubdet[iL] == -1) continue;
        if(electrons.pt[iL] < 30.) break;
        if(!electron_matchHLT[iL]) continue;
        if(!electron_isCand[iL]) continue;
        if(electrons.superClusterIndex[iL] == elePhoton.superClusterIndex ||
           TVector3(electrons.caloX[iL], electrons.caloY[iL], electrons.caloZ[iL]).DeltaR(caloPosition) < 0.02)
          continue;
        electronsOut_matchGen[electronsOut.size] = electron_matchGen[iL];
        electronsOut.push_back(electrons.at(iL));
      }

      if(electronsOut.size == 0) continue;
    }

    if(_outputType == oElePhotonAndElectron){
      if(susy::deltaR(elePhoton.eta, elePhoton.phi, electronsOut.eta[0], electronsOut.phi[0]) < 0.7) continue;
      calculateKinematics(elePhoton, electronsOut);
    }
    else if(_outputType == oElePhotonAndMuon){
      if(susy::deltaR(elePhoton.eta, elePhoton.phi, muonsOut.eta[0], muonsOut.phi[0]) < 0.7) continue;
      calculateKinematics(elePhoton, muonsOut);
    }

    photonsOut.clear();

    photonsOut_matchGen[photonsOut.size] = photon_matchGen[iEP];
    photonsOut_matchGenE[photonsOut.size] = photon_matchGenE[iEP];
    photonsOut_genIso[photonsOut.size] = photon_genIso[iEP];
    photonsOut.push_back(elePhoton);

    setJets();

    if(weightCalc[_outputType])
      static_cast<WGEventWeight*>(weightCalc[_outputType])->setPhoton(elePhoton, vertices, eventVars);

    fill(_outputType);
  }
}

void
WGammaProcessor::processFakePhotonAndLepton(unsigned _outputType)
{
  if(_outputType == oFakePhotonAndElectron){
    electronsOut.clear();

    for(unsigned iL(0); iL != electrons.size; ++iL){
      if(electrons.iSubdet[iL] == -1) continue;
      if(electrons.pt[iL] < 30.) break;
      if(!electron_matchHLT[iL]) continue;
      if(!electron_isCand[iL]) continue;
      electronsOut_matchGen[electronsOut.size] = electron_matchGen[iL];
      electronsOut.push_back(electrons.at(iL));
    }

    if(electronsOut.size == 0) return;
  }
  else if(_outputType == oFakePhotonAndMuon){
    muonsOut.clear();

    unsigned nMuons(0);

    for(unsigned iL(0); iL != muons.size; ++iL){
      if(muons.iSubdet[iL] == -1) continue;
      if(muons.pt[iL] > 10.) ++nMuons;
      if(muons.pt[iL] < 26.) continue;
      if(!muon_matchHLT[iL] && !muon_match21HLT[iL]) continue;
      if(!muon_isCand[iL]) continue;
      muonsOut_matchGen[muonsOut.size] = muon_matchGen[iL];
      muonsOut.push_back(muons.at(iL));
    }

    if(nMuons != 1 || muonsOut.size == 0) return;
  }
  else
    throw std::runtime_error("Incorrect filter type");

  for(unsigned iFP(0); iFP != photons.size; ++iFP){
    if(photons.iSubdet[iFP] == -1) continue;
    if(photons.pt[iFP] < 15.) break;
    if(!photon_isFake[iFP]) continue;

    if(photons.hOverE[iFP] > 0.05 || photons.sigmaIetaIeta[iFP] > 0.014) continue;
    if(_outputType == oFakePhotonAndMuon && (photons.chargedHadronIso[iFP] > 15. || photons.neutralHadronIso[iFP] > 3.5 || photons.photonIso[iFP] > 1.3)) continue;

    susy::PhotonVars fakePhoton(photons.at(iFP));

    if(_outputType == oFakePhotonAndElectron){
      if(susy::deltaR(fakePhoton.eta, fakePhoton.phi, electronsOut.eta[0], electronsOut.phi[0]) < 0.7) continue;
      calculateKinematics(fakePhoton, electronsOut);
    }
    else if(_outputType == oFakePhotonAndMuon){
      if(susy::deltaR(fakePhoton.eta, fakePhoton.phi, muonsOut.eta[0], muonsOut.phi[0]) < 0.7) continue;
      calculateKinematics(fakePhoton, muonsOut);
    }

    photonsOut.clear();

    photonsOut_matchGen[photonsOut.size] = photon_matchGen[iFP];
    photonsOut_matchGenE[photonsOut.size] = photon_matchGenE[iFP];
    photonsOut_genIso[photonsOut.size] = photon_genIso[iFP];
    photonsOut.push_back(fakePhoton);

    double etaVeto[1] = {fakePhoton.eta};
    double phiVeto[1] = {fakePhoton.phi};
    setJets(1, etaVeto, phiVeto);

    if(weightCalc[_outputType])
      static_cast<WGEventWeight*>(weightCalc[_outputType])->setPhoton(fakePhoton, vertices, eventVars);

    fill(_outputType);
  }
}

void
WGammaProcessor::processPhotonAndFakeLepton(unsigned _outputType)
{
  if(_outputType != oPhotonAndFakeElectron && _outputType != oPhotonAndFakeMuon)
    throw std::runtime_error("Incorrect filter type");

  photonsOut.clear();

  for(unsigned iP(0); iP != photons.size; ++iP){
    if(photons.iSubdet[iP] == -1) continue;
    if(photons.pt[iP] < 15.) break;
    if(!photon_isCand[iP]) continue;

    photonsOut_matchGen[photonsOut.size] = photon_matchGen[iP];
    photonsOut_matchGenE[photonsOut.size] = photon_matchGenE[iP];
    photonsOut_genIso[photonsOut.size] = photon_genIso[iP];
    photonsOut.push_back(photons.at(iP));
  }

  if(photonsOut.size == 0) return;

  susy::PhotonVars photon(photonsOut.at(0));

  unsigned size(0);
  if(_outputType == oPhotonAndFakeElectron) size = electrons.size;
  else if(_outputType == oPhotonAndFakeMuon) size = muons.size;

  WGEventWeight* calc(static_cast<WGEventWeight*>(weightCalc[_outputType]));

  std::bitset<susy::nElectronCriteria> elIdResults;
  std::bitset<susy::nElectronCriteria> elBaseline(susy::ObjectSelector::elReferences[susy::ElMedium12]);
  elBaseline.reset(susy::ElCombIso);
  elBaseline.reset(susy::ElDeltaEta);
  elBaseline.reset(susy::ElDeltaPhi);
  std::bitset<susy::nMuonCriteria> muIdResults;
  std::bitset<susy::nMuonCriteria> muBaseline(susy::ObjectSelector::muReferences[susy::MuTight12]);
  muBaseline.reset(susy::MuCombIso);

  if(_outputType == oPhotonAndFakeMuon){
    unsigned nMuons(0);
    for(unsigned iFL(0); iFL != size; ++iFL)
      if(muons.iSubdet[iFL] != -1 && muons.pt[iFL] > 10.) ++nMuons;
    if(nMuons != 1) return;
  }

  for(unsigned iFL(0); iFL != size; ++iFL){
    double lEta(0.);
    double lPhi(0.);

    if(_outputType == oPhotonAndFakeElectron){
      if(electrons.iSubdet[iFL] == -1) continue;
      if(electrons.pt[iFL] < 30.) break;
      if(!electron_matchHLT[iFL]) continue;
      if(!electron_isFake[iFL]) continue;

      if(electrons.combRelIso[iFL] * electrons.pt[iFL] < 10.) continue;

      susy::ElectronVars fakeElectron(electrons.at(iFL));

      susy::ObjectSelector::isGoodElectron(fakeElectron, susy::ElMedium12, &elIdResults);
      if((elIdResults & elBaseline) != elBaseline) continue;
      if(elIdResults[susy::ElDeltaEta] && elIdResults[susy::ElDeltaPhi]) continue;

      if(susy::deltaR(photon, fakeElectron) < 0.7) continue;

      electronsOut.clear();

      electronsOut_matchGen[electronsOut.size] = electron_matchGen[iFL];
      electronsOut.push_back(fakeElectron);

      calculateKinematics(photon, electronsOut);

      lEta = fakeElectron.eta;
      lPhi = fakeElectron.phi;

      if(calc) calc->setElectron(fakeElectron, vertices, eventVars);
    }
    else if(_outputType == oPhotonAndFakeMuon){
      if(muons.iSubdet[iFL] == -1) continue;
      if(muons.pt[iFL] < 26.) break;
      if(!muon_matchHLT[iFL] && !muon_match21HLT[iFL]) continue;
      if(!muon_isFake[iFL]) continue;

      if(muons.combRelIso[iFL] < 0.15 || muons.combRelIso[iFL] > 0.6) continue;

      susy::MuonVars fakeMuon(muons.at(iFL));

      susy::ObjectSelector::isGoodMuon(fakeMuon, susy::MuTight12, &muIdResults);
      if((muIdResults & muBaseline) != muBaseline) continue;

      if(susy::deltaR(photon, fakeMuon) < 0.7) continue;

      muonsOut.clear();

      muonsOut_matchGen[muonsOut.size] = muon_matchGen[iFL];
      muonsOut.push_back(fakeMuon);

      calculateKinematics(photon, muonsOut);

      lEta = fakeMuon.eta;
      lPhi = fakeMuon.phi;

      if(calc) calc->setMuon(fakeMuon, vertices, eventVars);
    }

    double etaVeto[1] = {lEta};
    double phiVeto[1] = {lPhi};
    setJets(1, etaVeto, phiVeto);

    fill(_outputType);
  }
}

void
WGammaProcessor::setJets(unsigned _nObj/* = 0*/, double* _eta/* = 0*/, double* _phi/* = 0*/)
{
  jetsOut.clear();
  ht = 0.;
  htJESUp = 0.;
  htJESDown = 0.;
  nJet = 0;
  nJetJESUp = 0;
  nJetJESDown = 0;
  for(unsigned iJ(0); iJ != jets.size; ++iJ){
    if(!jet_isCand[iJ]) continue;

    unsigned iObj(0);
    for(; iObj != _nObj; ++iObj)
      if(susy::deltaR(_eta[iObj], _phi[iObj], jets.eta[iJ], jets.phi[iJ]) < 0.5) break;
    if(iObj != _nObj) continue;

    jetsOut.push_back(jets.at(iJ));
    if(std::abs(jets.eta[iJ]) > 2.5) continue;
    if(jets.pt[iJ] > 30.){
      ht += jets.pt[iJ];
      ++nJet;
    }
    double jetPtUp(jets.pt[iJ] * (jets.jecScale[iJ] + jets.jecUncert[iJ]) / jets.jecScale[iJ]);
    if(jetPtUp > 30.){
      htJESUp += jetPtUp;
      ++nJetJESUp;
    }
    double jetPtDown(jets.pt[iJ] * (jets.jecScale[iJ] - jets.jecUncert[iJ]) / jets.jecScale[iJ]);
    if(jetPtDown > 30.){
      htJESDown += jetPtDown;
      ++nJetJESDown;
    }
  }
}


class MCJetPhotonWeight : public WGEventWeight {
public:
  MCJetPhotonWeight(int _lepton) :
    WGEventWeight("MCJetPhotonHLTIsoWeight"),
    tfact_(0)
  {
    TFile* source(0);
    if(_lepton == 0)
      source = TFile::Open("/afs/cern.ch/user/y/yiiyama/output/GammaL/jetGammaFake/jetGammaMC_el.root");
    else
      source = TFile::Open("/afs/cern.ch/user/y/yiiyama/output/GammaL/jetGammaFake/jetGammaMC_mu.root");
    if(!source)
      throw std::runtime_error("MCJetPhotonWeight source not found");
    tfact_ = static_cast<TH1D*>(source->Get("tfactTrue"));
    tfact_->SetDirectory(0);
    delete source;
  }
  ~MCJetPhotonWeight()
  {
    delete tfact_;
  }

  void setPhoton(susy::PhotonVars const& _photon, susy::VertexVarsArray const&, susy::SimpleEventProducer::EventVars const&)
  {
    int iBin(tfact_->FindFixBin(_photon.pt));
    if(iBin == 0) iBin = 1;
    else if(iBin > tfact_->GetNbinsX()) iBin = tfact_->GetNbinsX();
    
    weight = tfact_->GetBinContent(iBin);
    relErr = tfact_->GetBinError(iBin) / weight;
  }

private:
  TH1D* tfact_;
};

class JetPhotonHLTIsoWeight : public WGEventWeight {
public:
  JetPhotonHLTIsoWeight() : WGEventWeight("JetPhotonHLTIsoWeight") {}
  ~JetPhotonHLTIsoWeight() {}
  void setPhoton(susy::PhotonVars const& _photon, susy::VertexVarsArray const&, susy::SimpleEventProducer::EventVars const&)
  {
    weight = (0.0415 + 5.22 * std::exp(-_photon.pt / 23.7)) / (1. + 30.5 * std::exp(-_photon.pt / 16.4));
    relErr = (0.010275140445 * _photon.pt + -0.212842251317) * 1.0666;
  }
};

class JetPhotonWeight : public WGEventWeight {
public:
  JetPhotonWeight() : WGEventWeight("JetPhotonWeight") {}
  ~JetPhotonWeight() {}
  void setPhoton(susy::PhotonVars const& _photon, susy::VertexVarsArray const&, susy::SimpleEventProducer::EventVars const&)
  {
    weight = (0.461 + 28.0 * std::exp(-_photon.pt / 14.0)) / (1. + 42.1 * std::exp(-_photon.pt / 13.7));
    relErr = (0.00968391095319 * _photon.pt + -0.326302767227) * 1.2035;
  }
};

class MCElePhotonFunctionalWeight : public WGEventWeight {
public:
  MCElePhotonFunctionalWeight(bool _muCorr = false) : WGEventWeight("MCElePhotonFunctionalWeight"), muCorr_(_muCorr) {}
  ~MCElePhotonFunctionalWeight() {}
  void setPhoton(susy::PhotonVars const& _photon, susy::VertexVarsArray const& _vertices, susy::SimpleEventProducer::EventVars const&)
  {
    unsigned nVtx(0);
    int iPV(-1);
    for(unsigned iV(0); iV != _vertices.size; ++iV){
      if(!_vertices.isGood[iV]) continue;
      if(iPV < 0) iPV = iV;
      ++nVtx;
    }
    double ineff(0.00143056);
    if(muCorr_) ineff = 0.0017645;

    double fakerate(1. - (1. - ineff) * (1. - std::pow(1.051e+00 * _photon.pt + 1., -1.613e+00)) * (1. - 1.903e-01 * std::exp(-3.940e-01 * _vertices.nTracks[iPV])) * (1. - 1.312e-04 * nVtx));
    weight = fakerate / (1. - fakerate);
    relErr = 0.050;
  }

private:
  bool muCorr_;
};

class ElePhotonFunctionalWeight : public WGEventWeight {
public:
  ElePhotonFunctionalWeight(bool _muCorr = false) : WGEventWeight("ElePhotonFunctionalWeight"), muCorr_(_muCorr) {}
  ~ElePhotonFunctionalWeight() {}
  void setPhoton(susy::PhotonVars const& _photon, susy::VertexVarsArray const& _vertices, susy::SimpleEventProducer::EventVars const&)
  {
    unsigned nVtx(0);
    int iPV(-1);
    for(unsigned iV(0); iV != _vertices.size; ++iV){
      if(!_vertices.isGood[iV]) continue;
      if(iPV < 0) iPV = iV;
      ++nVtx;
    }
    double ineff(0.0019387);
    if(muCorr_) ineff = 0.002344;

    double fakerate(1. - (1. - ineff) * (1. - std::pow(7.094e-02 * _photon.pt + 1., -4.945e+00)) * (1. - 1.427e-01 * std::exp(-2.961e-01 * _vertices.nTracks[iPV])) * (1. - 3.149e-04 * nVtx));
    weight = fakerate / (1. - fakerate);
    relErr = 0.14;
  }

private:
  bool muCorr_;
};

class PhotonGenIsoWeight : public WGEventWeight {
public:
  PhotonGenIsoWeight() :
    WGEventWeight("PhotonGenIsoWeight"),
    table_(0)
  {
    TFile* ttgSource(TFile::Open("/afs/cern.ch/user/y/yiiyama/work/TTGJets.root"));
    TTree* ttgTree(static_cast<TTree*>(ttgSource->Get("photonTree")));
    TFile* whizSource(TFile::Open("/afs/cern.ch/user/y/yiiyama/work/WHIZARD.root"));
    TTree* whizTree(static_cast<TTree*>(whizSource->Get("photonTree")));

    double binning[100];
    for(unsigned i(0); i != 50; ++i) binning[i] = i;
    for(unsigned i(50); i != 75; ++i) binning[i] = 50. + 2. * (i - 50);
    for(unsigned i(75); i != 80; ++i) binning[i] = 100. + 10. * (i - 75);
    binning[80] = 150.;
    binning[81] = 200.;
    binning[82] = 300.;
    TH1D* ttg(new TH1D("ttg", "ttg", 82, binning));
    TH1D* whiz(new TH1D("whiz", "whiz", 82, binning));
    ttg->Sumw2();
    whiz->Sumw2();
    ttgTree->Draw("iso>>+ttg", "", "goff e");
    whizTree->Draw("iso>>+whiz", "", "goff e");

    table_ = static_cast<TH1*>(whiz->Clone("PhotonGenIsoWeight"));
    table_->SetDirectory(0);
    table_->Divide(ttg);

    delete ttg;
    delete whiz;
    delete ttgSource;
    delete whizSource;
  }
  ~PhotonGenIsoWeight()
  {
    delete table_;
  }
  void setPhoton(susy::PhotonVars const& _photon, susy::VertexVarsArray const&, susy::SimpleEventProducer::EventVars const& _eventVars)
  {
    TVector3 caloPosition(_photon.caloX, _photon.caloY, _photon.caloZ);

    unsigned iPh(0);
    for(; iPh != _eventVars.gen_size; ++iPh){
      if(_eventVars.gen_status[iPh] != 1 || _eventVars.gen_pdgId[iPh] != 22) continue;
      TVector3 dir(-_eventVars.gen_vx[iPh], -_eventVars.gen_vy[iPh], -_eventVars.gen_vz[iPh]);
      dir += caloPosition;
      double dEta(dir.Eta() - _eventVars.gen_eta[iPh]);
      double dPhi(TVector2::Phi_mpi_pi(dir.Phi() - _eventVars.gen_phi[iPh]));
      if(dEta * dEta + dPhi * dPhi < 0.01) break;
    }
    if(iPh == _eventVars.gen_size){
      weight = 0.;
      relErr = 0.;
      return;
    }

    double iso(0.);
    for(unsigned iH(0); iH != _eventVars.gen_size; ++iH){
      if(iH == iPh) continue;
      if(_eventVars.gen_status[iH] != 1) continue;
      double dEta(_eventVars.gen_eta[iH] - _eventVars.gen_eta[iPh]);
      double dPhi(TVector2::Phi_mpi_pi(_eventVars.gen_phi[iH] - _eventVars.gen_phi[iPh]));
      if(dEta * dEta + dPhi * dPhi < 0.01) iso += _eventVars.gen_pt[iH];
    }

    int iBin(table_->FindFixBin(iso));
    if(iBin > table_->GetNbinsX()){
      weight = 0.;
      relErr = 0.;
    }
    else{
      weight = table_->GetBinContent(iBin);
      relErr = table_->GetBinError(iBin) / weight;
    }
  }

private:
  TH1* table_;
};
