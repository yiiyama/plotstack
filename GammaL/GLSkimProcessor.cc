#include "../../CommonCode/ObjectTree.h"
#include "../../CommonCode/Utilities.h"
#include "../../CommonCode/SimpleEventProducer.h"
#include "../../CommonCode/ObjectSelector.h"

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TChain.h"
#include "TString.h"

#include <cmath>
#include <iostream>
#include <utility>
#include <algorithm>
#include <stdexcept>

#include "../ROOT/EventProcessor.h"

enum FilterTypes {
  kPhotonAndElectron,
  kPhotonAndMuon,
  kElePhotonAndElectron,
  kElePhotonAndMuon,
  kFakePhotonAndElectron,
  kFakePhotonAndMuon,
  kPhotonAndFakeElectron,
  kPhotonAndFakeMuon,
  nFilterTypes
};

TString filterNames[] = {
  "PhotonAndElectron",
  "PhotonAndMuon",
  "ElePhotonAndElectron",
  "ElePhotonAndMuon",
  "FakePhotonAndElectron",
  "FakePhotonAndMuon",
  "PhotonAndFakeElectron",
  "PhotonAndFakeMuon"
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

class GLEventWeight : public EventWeight {
public:
  GLEventWeight() : EventWeight() {}
  GLEventWeight(double _w, double _r) : EventWeight(_w, _r) {}
  ~GLEventWeight() {}
  virtual void setPhoton(susy::PhotonVars const&, susy::VertexVarsArray const&, susy::SimpleEventProducer::EventVars const&) {}
  virtual void setElectron(susy::ElectronVars const&, susy::VertexVarsArray const&, susy::SimpleEventProducer::EventVars const&) {}
  virtual void setMuon(susy::MuonVars const&, susy::VertexVarsArray const&, susy::SimpleEventProducer::EventVars const&) {}
};

class GLSkimProcessor : public EventProcessor {
public:
  GLSkimProcessor(char const*, double, double, char const*);
  ~GLSkimProcessor() {} // ROOT closes files in case of crashes; deleting event list causes double free

  void process();

private:
  bool preprocess();
  void processPhotonAndLepton(unsigned);
  void processElePhotonAndLepton(unsigned);
  void processFakePhotonAndLepton(unsigned);
  void processPhotonAndFakeLepton(unsigned);

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
  }

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

  bool hltBit;
  float mt;
  float mass2;
  float mass3;
  susy::PhotonVarsArray photonsOut;
  bool photonsOut_matchGen[susy::NMAX];
  bool photonsOut_matchGenE[susy::NMAX];
  bool photonsOut_matchHLT[susy::NMAX];
  susy::ElectronVarsArray electronsOut;
  bool electronsOut_matchGen[susy::NMAX];
  bool electronsOut_matchHLT[susy::NMAX];
  susy::MuonVarsArray muonsOut;
  bool muonsOut_matchGen[susy::NMAX];
  bool muonsOut_matchHLT[susy::NMAX];
  susy::JetVarsArray jetsOut;

  /* UTIL */

  bool photon_muonIso[susy::NMAX];
  bool photon_matchGen[susy::NMAX];
  bool photon_matchGenE[susy::NMAX];
  bool photon_matchEHLT[2][susy::NMAX];
  bool photon_matchMHLT[susy::NMAX];
  bool electron_matchGen[susy::NMAX];
  bool electron_matchHLT[susy::NMAX];
  bool muon_matchGen[susy::NMAX];
  bool muon_matchHLT[susy::NMAX];
  bool electronHLT;
  bool muonHLT;
};

GLSkimProcessor::GLSkimProcessor(char const* _datasetName, double _Leff, double _sigmaRelErr2, char const* _outputDir) :
  EventProcessor(nFilterTypes, _datasetName, _Leff, _sigmaRelErr2, _outputDir),
  useElectronFilter(false),
  useMuonFilter(false),
  eventVars(),
  electrons(),
  muons(),
  photons(),
  jets(),
  vertices(),
  hltBit(false),
  mt(0.),
  mass2(0.),
  mass3(0.),
  photonsOut(),
  electronsOut(),
  muonsOut(),
  jetsOut(),
  electronHLT(false),
  muonHLT(false)
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
  std::fill_n(photonsOut_matchHLT, susy::NMAX, false);
  std::fill_n(electronsOut_matchGen, susy::NMAX, false);
  std::fill_n(electronsOut_matchHLT, susy::NMAX, false);
  std::fill_n(muonsOut_matchGen, susy::NMAX, false);
  std::fill_n(muonsOut_matchHLT, susy::NMAX, false);

  std::fill_n(photon_muonIso, susy::NMAX, true);
  std::fill_n(photon_matchGen, susy::NMAX, false);
  std::fill_n(photon_matchGenE, susy::NMAX, false);
  std::fill_n(photon_matchEHLT[0], susy::NMAX, false);
  std::fill_n(photon_matchEHLT[1], susy::NMAX, false);
  std::fill_n(photon_matchMHLT, susy::NMAX, false);
  std::fill_n(electron_matchGen, susy::NMAX, false);
  std::fill_n(electron_matchHLT, susy::NMAX, false);
  std::fill_n(muon_matchGen, susy::NMAX, false);
  std::fill_n(muon_matchHLT, susy::NMAX, false);

  for(unsigned iF(0); iF != nFilterTypes; ++iF)
    filterIndices[filterNames[iF]] = iF;
}

void
GLSkimProcessor::process()
{
  for(unsigned iF(0); iF != nFilterTypes; ++iF)
    if(useFilter[iF] && !weightCalc[iF])
      throw runtime_error("Event weight calculator not set");

  useElectronFilter = useFilter[kPhotonAndElectron] || useFilter[kElePhotonAndElectron] || useFilter[kFakePhotonAndElectron] || useFilter[kPhotonAndFakeElectron];
  useMuonFilter = useFilter[kPhotonAndMuon] || useFilter[kElePhotonAndMuon] || useFilter[kFakePhotonAndMuon] || useFilter[kPhotonAndFakeMuon];

  TChain filterTree("eventVars");
  TChain eventTree("eventVars");
  TChain objectTree("allObjects");

  for(unsigned iP(0); iP != inputPaths.size(); ++iP){
    filterTree.Add(inputPaths[iP]);
    eventTree.Add(inputPaths[iP]);
    objectTree.Add(inputPaths[iP]);
  }

  filterTree.SetBranchStatus("*", 0);
  if(useElectronFilter){
    filterTree.SetBranchStatus("PhotonAndElectron", 1);
    filterTree.SetBranchStatus("ElePhotonAndElectron", 1);
    filterTree.SetBranchStatus("FakePhotonAndElectron", 1);
    filterTree.SetBranchStatus("PhotonAndFakeElectron", 1);
  }
  filterTree.SetBranchStatus("PhotonAndMuon", 1);
  filterTree.SetBranchStatus("ElePhotonAndMuon", 1);
  filterTree.SetBranchStatus("FakePhotonAndMuon", 1);
  filterTree.SetBranchStatus("PhotonAndFakeMuon", 1);

  bool filterBits[nFilterTypes];
  std::fill_n(filterBits, unsigned(nFilterTypes), false);

  if(useElectronFilter){
    filterTree.SetBranchAddress("PhotonAndElectron", filterBits + kPhotonAndElectron);
    filterTree.SetBranchAddress("ElePhotonAndElectron", filterBits + kElePhotonAndElectron);
    filterTree.SetBranchAddress("FakePhotonAndElectron", filterBits + kFakePhotonAndElectron);
    filterTree.SetBranchAddress("PhotonAndFakeElectron", filterBits + kPhotonAndFakeElectron);
  }
  filterTree.SetBranchAddress("PhotonAndMuon", filterBits + kPhotonAndMuon);
  filterTree.SetBranchAddress("ElePhotonAndMuon", filterBits + kElePhotonAndMuon);
  filterTree.SetBranchAddress("FakePhotonAndMuon", filterBits + kFakePhotonAndMuon);
  filterTree.SetBranchAddress("PhotonAndFakeMuon", filterBits + kPhotonAndFakeMuon);


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
  if(useElectronFilter) eventTree.SetBranchStatus("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50", 1);
  if(useMuonFilter) eventTree.SetBranchStatus("HLT_Mu22_Photon22_CaloIdL", 1);

  eventVars.setAddress(eventTree);
  if(useElectronFilter) eventTree.SetBranchAddress("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50", &electronHLT);
  if(useMuonFilter) eventTree.SetBranchAddress("HLT_Mu22_Photon22_CaloIdL", &muonHLT);


  objectTree.SetBranchStatus("*", 0);
  if(useElectronFilter) objectTree.SetBranchStatus("electron*", 1);
  if(useMuonFilter) objectTree.SetBranchStatus("muon*", 1);
  objectTree.SetBranchStatus("photon*", 1);
  objectTree.SetBranchStatus("vertex*", 1);
  objectTree.SetBranchStatus("jet*", 1);

  if(useElectronFilter){
    electrons.setAddress(objectTree);
    objectTree.SetBranchAddress("electron.isCand", electron_isCand);
    objectTree.SetBranchAddress("electron.isFake", electron_isFake);
    objectTree.SetBranchAddress("electron.hltEG22CaloId10Iso50TrackIsoDoubleLastFilterUnseeded", electron_matchHLT);
  }
  if(useMuonFilter){
    muons.setAddress(objectTree);
    objectTree.SetBranchAddress("muon.isCand", muon_isCand);
    objectTree.SetBranchAddress("muon.isFake", muon_isFake);
    objectTree.SetBranchAddress("muon.hltL1Mu3p5EG12L3Filtered22", muon_matchHLT);
  }
  photons.setAddress(objectTree);
  if(useElectronFilter){
    objectTree.SetBranchAddress("photon.hltEG36CaloId10Iso50HcalIsoLastFilter", photon_matchEHLT[0]);
    objectTree.SetBranchAddress("photon.hltEG22CaloId10Iso50TrackIsoDoubleLastFilterUnseeded", photon_matchEHLT[1]);
  }
  if(useMuonFilter)
    objectTree.SetBranchAddress("photon.hltMu22Photon22CaloIdLHEFilter", photon_matchMHLT);
  objectTree.SetBranchAddress("photon.isCand", photon_isCand);
  objectTree.SetBranchAddress("photon.isFake", photon_isFake);
  objectTree.SetBranchAddress("photon.isEle", photon_isEle);
  jets.setAddress(objectTree);
  objectTree.SetBranchAddress("jet.isCand", &jet_isCand);
  vertices.setAddress(objectTree);

  for(unsigned iF(0); iF != nFilterTypes; ++iF){
    if(!eventList[iF]) continue;

    eventList[iF]->Branch("run", &eventVars.runNumber, "run/i");
    eventList[iF]->Branch("lumi", &eventVars.lumiNumber, "lumi/i");
    eventList[iF]->Branch("event", &eventVars.eventNumber, "event/i");
    eventList[iF]->Branch("hltBit", &hltBit, "hlt/O");
    eventList[iF]->Branch("puWeight", &eventVars.puWeight, "puWeight/F");
    eventList[iF]->Branch("met", &eventVars.met, "met/F");
    eventList[iF]->Branch("metPhi", &eventVars.metPhi, "metPhi/F");
    eventList[iF]->Branch("rho", &eventVars.rho, "rho/F");
    eventList[iF]->Branch("mt", &mt, "mt/F");
    eventList[iF]->Branch("mass2", &mass2, "mass2/F");
    eventList[iF]->Branch("mass3", &mass3, "mass3/F");
    photonsOut.setBranches(*eventList[iF]);
    eventList[iF]->Branch("photon.matchGen", photonsOut_matchGen, "matchGen[photon.size]/O");
    eventList[iF]->Branch("photon.matchGenE", photonsOut_matchGenE, "matchGenE[photon.size]/O");
    eventList[iF]->Branch("photon.matchHLT", photonsOut_matchHLT, "matchHLT[photon.size]/O");
    if(iF == kPhotonAndElectron || iF == kElePhotonAndElectron || iF == kFakePhotonAndElectron || iF == kPhotonAndFakeElectron){
      electronsOut.setBranches(*eventList[iF]);
      eventList[iF]->Branch("electron.matchGen", electronsOut_matchGen, "matchGen[electron.size]/O");
      eventList[iF]->Branch("electron.matchHLT", electronsOut_matchHLT, "matchHLT[electron.size]/O");
    }
    if(iF == kPhotonAndMuon || iF == kElePhotonAndMuon || iF == kFakePhotonAndMuon || iF == kPhotonAndFakeMuon){
      muonsOut.setBranches(*eventList[iF]);
      eventList[iF]->Branch("muon.matchGen", muonsOut_matchGen, "matchGen[muon.size]/O");
      eventList[iF]->Branch("muon.matchHLT", muonsOut_matchHLT, "matchHLT[muon.size]/O");
    }
    jetsOut.setBranches(*eventList[iF]);
  }

  ////////////////////
  //// START LOOP ////
  ////////////////////

  std::cout << "Processing " << datasetName << std::endl;

  try{

    long iEntry(0);
    while(filterTree.GetEntry(iEntry++) > 0){
      if(iEntry % 1000 == 1) (std::cout << "\r" << iEntry).flush();

      unsigned iF(0);
      for(; iF != nFilterTypes; ++iF)
        if(useFilter[iF] && filterBits[iF]) break;
      if(iF == nFilterTypes) continue;

      eventTree.GetEntry(iEntry - 1);
      objectTree.GetEntry(iEntry - 1);

      if(!preprocess()) continue;

      if(filterBits[kPhotonAndElectron] && !filterBits[kPhotonAndMuon] && useFilter[kPhotonAndElectron]){
        processPhotonAndLepton(kPhotonAndElectron);
        continue;
      }

      if(filterBits[kPhotonAndMuon] && useFilter[kPhotonAndMuon]){
        processPhotonAndLepton(kPhotonAndMuon);
        continue;
      }

      ///////////////////
      /// FAKE EVENTS ///
      ///////////////////
      /* Since we only consider single fakes at the moment, we never have the "second candidate". */
      /* If e.g. we had a fake photon and a candidate photon, that is a candidate event. */
      /* Consequently, the array of the faking object will always have only one element. */

      if(filterBits[kElePhotonAndElectron] && !filterBits[kElePhotonAndMuon] && useFilter[kElePhotonAndElectron])
        processElePhotonAndLepton(kElePhotonAndElectron);

      if(filterBits[kElePhotonAndMuon] && useFilter[kElePhotonAndMuon])
        processElePhotonAndLepton(kElePhotonAndMuon);

      if(filterBits[kFakePhotonAndElectron] && !filterBits[kFakePhotonAndMuon] && useFilter[kFakePhotonAndElectron])
        processFakePhotonAndLepton(kFakePhotonAndElectron);

      if(filterBits[kFakePhotonAndMuon] && useFilter[kFakePhotonAndMuon])
        processFakePhotonAndLepton(kFakePhotonAndMuon);

      if(filterBits[kPhotonAndFakeElectron] && !filterBits[kPhotonAndMuon] && useFilter[kPhotonAndFakeElectron])
        processPhotonAndFakeLepton(kPhotonAndFakeElectron);

      if(filterBits[kPhotonAndFakeMuon] && useFilter[kPhotonAndFakeMuon])
        processPhotonAndFakeLepton(kPhotonAndFakeMuon);
    }

    std::cout << std::endl;

  }
  catch(std::exception& _ex){
    std::cerr << std::endl << _ex.what() << std::endl;
    throw;
  }

}

bool
GLSkimProcessor::preprocess()
{
  for(unsigned iP(0); iP != photons.size; ++iP){
    unsigned iL(0);
    for(; iL != muons.size; ++iL)
      if(susy::deltaR(muons.eta[iL], muons.phi[iL], photons.eta[iP], photons.phi[iP]) < 0.3) break;
    photon_muonIso[iP] = iL == muons.size;

    photon_matchGen[iP] = false;
    photon_matchGenE[iP] = false;
    for(unsigned iG(0); iG != eventVars.gen_size; ++iG){
      if(eventVars.gen_status[iG] != 1) continue;
      unsigned absId(std::abs(eventVars.gen_pdgId[iG]));
      if(absId != 22 && absId != 11) continue;
      short mIdx(eventVars.gen_motherIndex[iG]);
      if(mIdx < 0) continue;
      if(std::abs(eventVars.gen_pdgId[mIdx]) > 99) continue;

      TVector3 dir(photons.caloX[iP] - eventVars.gen_vx[iG], photons.caloY[iP] - eventVars.gen_vy[iG], photons.caloZ[iP] - eventVars.gen_vz[iG]);
      if(susy::deltaR(eventVars.gen_eta[iG], eventVars.gen_phi[iG], dir.Eta(), dir.Phi()) < 0.1){
        if(absId == 11){
          photon_matchGenE[iP] = true;
          break;
        }
        else if(absId == 22)
          photon_matchGen[iP] = true;
      }
    }
  }

  unsigned iP(0);
  for(; iP != photons.size; ++iP)
    if(photon_muonIso[iP]) break;

  if(iP == photons.size) return false;

  if(useElectronFilter){
    for(unsigned iE(0); iE != electrons.size; ++iE){
      electron_matchGen[iE] = false;
      for(unsigned iG(0); iG != eventVars.gen_size; ++iG){
        if(eventVars.gen_status[iG] != 1 || std::abs(eventVars.gen_pdgId[iG]) != 11) continue;

        TVector3 dir(electrons.caloX[iE] - eventVars.gen_vx[iG], electrons.caloY[iE] - eventVars.gen_vy[iG], electrons.caloZ[iE] - eventVars.gen_vz[iG]);
        if(susy::deltaR(eventVars.gen_eta[iG], eventVars.gen_phi[iG], dir.Eta(), dir.Phi()) < 0.1){
          electron_matchGen[iE] = true;
          break;
        }
      }
    }
  }

  if(useMuonFilter){
    for(unsigned iM(0); iM != muons.size; ++iM){
      muon_matchGen[iM] = false;
      for(unsigned iG(0); iG != eventVars.gen_size; ++iG){
        if(eventVars.gen_status[iG] != 1 || std::abs(eventVars.gen_pdgId[iG]) != 13) continue;

        if(susy::deltaR(eventVars.gen_eta[iG], eventVars.gen_phi[iG], muons.eta[iM], muons.phi[iM]) < 0.1){
          muon_matchGen[iM] = true;
          break;
        }
      }
    }
  }

  return true;
}

void
GLSkimProcessor::processPhotonAndLepton(unsigned _filterType)
{
  if(_filterType != kPhotonAndElectron && _filterType != kPhotonAndMuon)
    throw runtime_error("Incorrect filter type");

  if(_filterType == kPhotonAndElectron)
    hltBit = electronHLT;
  else
    hltBit = muonHLT;

  photonsOut.clear();

  for(unsigned iP(0); iP != photons.size; ++iP){
    if(!photon_isCand[iP]) continue;
    if(!photon_muonIso[iP]) continue;
    photonsOut_matchGen[photonsOut.size] = photon_matchGen[iP];
    photonsOut_matchGenE[photonsOut.size] = photon_matchGenE[iP];
    if(_filterType == kPhotonAndElectron)
      photonsOut_matchHLT[photonsOut.size] = photon_matchEHLT[0][iP] && photon_matchEHLT[1][iP];
    else
      photonsOut_matchHLT[photonsOut.size] = photon_matchMHLT[iP];
    photonsOut.push_back(photons.at(iP));
  }

  if(photonsOut.size == 0) return;

  if(_filterType == kPhotonAndElectron){
    electronsOut.clear();

    for(unsigned iL(0); iL != electrons.size; ++iL){
      if(!electron_isCand[iL]) continue;
      electronsOut_matchGen[electronsOut.size] = electron_matchGen[iL];
      electronsOut_matchHLT[electronsOut.size] = electron_matchHLT[iL];
      electronsOut.push_back(electrons.at(iL));
    }

    if(electronsOut.size == 0) return;
  }
  else if(_filterType == kPhotonAndMuon){
    muonsOut.clear();

    for(unsigned iL(0); iL != muons.size; ++iL){
      if(!muon_isCand[iL]) continue;
      muonsOut_matchGen[muonsOut.size] = muon_matchGen[iL];
      muonsOut_matchHLT[muonsOut.size] = muon_matchHLT[iL];
      muonsOut.push_back(muons.at(iL));
    }

    if(muonsOut.size == 0) return;
  }

  if(_filterType == kPhotonAndElectron){
    if(susy::deltaR(photonsOut.eta[0], photonsOut.phi[0], electronsOut.eta[0], electronsOut.phi[0]) < 0.8) return;
    calculateKinematics(photonsOut.at(0), electronsOut);
  }
  else if(_filterType == kPhotonAndMuon){
    if(susy::deltaR(photonsOut.eta[0], photonsOut.phi[0], muonsOut.eta[0], muonsOut.phi[0]) < 0.8) return;
    calculateKinematics(photonsOut.at(0), muonsOut);
  }

  jetsOut.clear();
  for(unsigned iJ(0); iJ != jets.size; ++iJ)
    if(jet_isCand[iJ]) jetsOut.push_back(jets.at(iJ));

  fill(_filterType);
}

void
GLSkimProcessor::processElePhotonAndLepton(unsigned _filterType)
{
  if(_filterType == kElePhotonAndMuon){
    muonsOut.clear();

    for(unsigned iL(0); iL != muons.size; ++iL){
      if(!muon_isCand[iL]) continue;
      muonsOut_matchGen[muonsOut.size] = muon_matchGen[iL];
      muonsOut_matchHLT[muonsOut.size] = muon_matchHLT[iL];
      muonsOut.push_back(muons.at(iL));
    }

    if(muonsOut.size == 0) return;
  }
  else if(_filterType != kElePhotonAndElectron)
    throw runtime_error("Incorrect filter type");

  if(_filterType == kElePhotonAndElectron)
    hltBit = electronHLT;
  else
    hltBit = muonHLT;

  for(unsigned iEP(0); iEP != photons.size; ++iEP){
    if(!photon_isEle[iEP]) continue;
    if(!photon_muonIso[iEP]) continue;

    susy::PhotonVars elePhoton(photons.at(iEP));

    TVector3 caloPosition(elePhoton.caloX, elePhoton.caloY, elePhoton.caloZ);

    if(_filterType == kElePhotonAndElectron){
      electronsOut.clear();

      for(unsigned iL(0); iL != electrons.size; ++iL){
        if(!electron_isCand[iL]) continue;
        if(electrons.superClusterIndex[iL] == elePhoton.superClusterIndex ||
           TVector3(electrons.caloX[iL], electrons.caloY[iL], electrons.caloZ[iL]).DeltaR(caloPosition) < 0.04)
          continue;
        electronsOut_matchGen[electronsOut.size] = electron_matchGen[iL];
        electronsOut_matchHLT[electronsOut.size] = electron_matchHLT[iL];
        electronsOut.push_back(electrons.at(iL));
      }

      if(electronsOut.size == 0) continue;
    }

    if(_filterType == kElePhotonAndElectron){
      if(susy::deltaR(elePhoton.eta, elePhoton.phi, electronsOut.eta[0], electronsOut.phi[0]) < 0.8) continue;
      calculateKinematics(elePhoton, electronsOut);
    }
    else if(_filterType == kElePhotonAndMuon){
      if(susy::deltaR(elePhoton.eta, elePhoton.phi, muonsOut.eta[0], muonsOut.phi[0]) < 0.8) continue;
      calculateKinematics(elePhoton, muonsOut);
    }

    photonsOut.clear();

    photonsOut_matchGen[photonsOut.size] = photon_matchGen[iEP];
    photonsOut_matchGenE[photonsOut.size] = photon_matchGenE[iEP];
    if(_filterType == kElePhotonAndElectron)
      photonsOut_matchHLT[photonsOut.size] = photon_matchEHLT[0][iEP] && photon_matchEHLT[1][iEP];
    else
      photonsOut_matchHLT[photonsOut.size] = photon_matchMHLT[iEP];
    photonsOut.push_back(elePhoton);

    jetsOut.clear();
    for(unsigned iJ(0); iJ != jets.size; ++iJ)
      if(jet_isCand[iJ]) jetsOut.push_back(jets.at(iJ));

    static_cast<GLEventWeight*>(weightCalc[_filterType])->setPhoton(elePhoton, vertices, eventVars);

    fill(_filterType);
  }
}

void
GLSkimProcessor::processFakePhotonAndLepton(unsigned _filterType)
{
  if(_filterType == kFakePhotonAndElectron){
    electronsOut.clear();

    for(unsigned iL(0); iL != electrons.size; ++iL){
      if(!electron_isCand[iL]) continue;
      electronsOut_matchGen[electronsOut.size] = electron_matchGen[iL];
      electronsOut_matchHLT[electronsOut.size] = electron_matchHLT[iL];
      electronsOut.push_back(electrons.at(iL));
    }

    if(electronsOut.size == 0) return;

    hltBit = electronHLT;
  }
  else if(_filterType == kFakePhotonAndMuon){
    muonsOut.clear();

    for(unsigned iL(0); iL != muons.size; ++iL){
      if(!muon_isCand[iL]) continue;
      muonsOut_matchGen[muonsOut.size] = muon_matchGen[iL];
      muonsOut_matchHLT[muonsOut.size] = muon_matchHLT[iL];
      muonsOut.push_back(muons.at(iL));
    }

    if(muonsOut.size == 0) return;

    hltBit = muonHLT;
  }
  else
    throw runtime_error("Incorrect filter type");

  for(unsigned iFP(0); iFP != photons.size; ++iFP){
    if(!photon_isFake[iFP]) continue;
    if(!photon_muonIso[iFP]) continue;
    if(photons.hOverE[iFP] > 0.05 || photons.sigmaIetaIeta[iFP] > 0.014) continue;
    if(_filterType == kFakePhotonAndMuon && (photons.chargedHadronIso[iFP] > 15. || photons.neutralHadronIso[iFP] > 3.5 || photons.photonIso[iFP] > 1.3)) continue;

    susy::PhotonVars fakePhoton(photons.at(iFP));

    if(_filterType == kFakePhotonAndElectron){
      if(susy::deltaR(fakePhoton.eta, fakePhoton.phi, electronsOut.eta[0], electronsOut.phi[0]) < 0.8) continue;
      calculateKinematics(fakePhoton, electronsOut);
    }
    else if(_filterType == kFakePhotonAndMuon){
      if(susy::deltaR(fakePhoton.eta, fakePhoton.phi, muonsOut.eta[0], muonsOut.phi[0]) < 0.8) continue;
      calculateKinematics(fakePhoton, muonsOut);
    }

    photonsOut.clear();

    photonsOut_matchGen[photonsOut.size] = photon_matchGen[iFP];
    photonsOut_matchGenE[photonsOut.size] = photon_matchGenE[iFP];
    if(_filterType == kFakePhotonAndElectron)
      photonsOut_matchHLT[photonsOut.size] = photon_matchEHLT[0][iFP] && photon_matchEHLT[1][iFP];
    else
      photonsOut_matchHLT[photonsOut.size] = photon_matchMHLT[iFP];
    photonsOut.push_back(fakePhoton);

    jetsOut.clear();
    for(unsigned iJ(0); iJ != jets.size; ++iJ)
      if(jet_isCand[iJ] && susy::deltaR(fakePhoton.eta, fakePhoton.phi, jets.eta[iJ], jets.phi[iJ]) > 0.5)
        jetsOut.push_back(jets.at(iJ));

    static_cast<GLEventWeight*>(weightCalc[_filterType])->setPhoton(fakePhoton, vertices, eventVars);

    fill(_filterType);
  }
}

void
GLSkimProcessor::processPhotonAndFakeLepton(unsigned _filterType)
{
  if(_filterType != kPhotonAndFakeElectron && _filterType != kPhotonAndFakeMuon)
    throw runtime_error("Incorrect filter type");

  if(_filterType == kPhotonAndFakeElectron)
    hltBit = electronHLT;
  else
    hltBit = muonHLT;

  photonsOut.clear();

  for(unsigned iP(0); iP != photons.size; ++iP){
    if(!photon_isCand[iP]) continue;
    if(!photon_muonIso[iP]) continue;
    photonsOut_matchGen[photonsOut.size] = photon_matchGen[iP];
    photonsOut_matchGenE[photonsOut.size] = photon_matchGenE[iP];
    if(_filterType == kPhotonAndFakeElectron)
      photonsOut_matchHLT[photonsOut.size] = photon_matchEHLT[0][iP] && photon_matchEHLT[1][iP];
    else
      photonsOut_matchHLT[photonsOut.size] = photon_matchMHLT[iP];
    photonsOut.push_back(photons.at(iP));
  }

  if(photonsOut.size == 0) return;

  susy::PhotonVars photon(photonsOut.at(0));

  unsigned size(0);
  if(_filterType == kPhotonAndFakeElectron) size = electrons.size;
  else if(_filterType == kPhotonAndFakeMuon) size = muons.size;

  GLEventWeight* calc(static_cast<GLEventWeight*>(weightCalc[_filterType]));

  std::bitset<susy::nElectronCriteria> elIdResults;
  std::bitset<susy::nElectronCriteria> elBaseline(susy::ObjectSelector::elReferences[susy::ElMedium12]);
  elBaseline.reset(susy::ElCombIso);
  elBaseline.reset(susy::ElDeltaEta);
  elBaseline.reset(susy::ElDeltaPhi);
  std::bitset<susy::nMuonCriteria> muIdResults;
  std::bitset<susy::nMuonCriteria> muBaseline(susy::ObjectSelector::muReferences[susy::MuTight12]);
  muBaseline.reset(susy::MuCombIso);

  for(unsigned iFL(0); iFL != size; ++iFL){
    double lEta(0.);
    double lPhi(0.);

    if(_filterType == kPhotonAndFakeElectron){
      if(!electron_isFake[iFL]) continue;

      if(electrons.combRelIso[iFL] * electrons.pt[iFL] < 10.) continue;

      susy::ElectronVars fakeElectron(electrons.at(iFL));

      susy::ObjectSelector::isGoodElectron(fakeElectron, susy::ElMedium12, &elIdResults);
      if((elIdResults & elBaseline) != elBaseline) continue;
      if(elIdResults[susy::ElDeltaEta] && elIdResults[susy::ElDeltaPhi]) continue;

      if(susy::deltaR(photon, fakeElectron) < 0.8) continue;

      electronsOut.clear();

      electronsOut_matchGen[electronsOut.size] = electron_matchGen[iFL];
      electronsOut_matchHLT[electronsOut.size] = electron_matchHLT[iFL];
      electronsOut.push_back(fakeElectron);

      calculateKinematics(photon, electronsOut);

      lEta = fakeElectron.eta;
      lPhi = fakeElectron.phi;

      calc->setElectron(fakeElectron, vertices, eventVars);
    }
    else if(_filterType == kPhotonAndFakeMuon){
      if(!muon_isFake[iFL]) continue;

      if(muons.combRelIso[iFL] < 0.15 || muons.combRelIso[iFL] > 0.6) continue;

      susy::MuonVars fakeMuon(muons.at(iFL));

      susy::ObjectSelector::isGoodMuon(fakeMuon, susy::MuTight12, &muIdResults);
      if((muIdResults & muBaseline) != muBaseline) continue;

      if(susy::deltaR(photon, fakeMuon) < 0.8) continue;

      muonsOut.clear();

      muonsOut_matchGen[muonsOut.size] = muon_matchGen[iFL];
      muonsOut_matchHLT[muonsOut.size] = muon_matchHLT[iFL];
      muonsOut.push_back(fakeMuon);

      calculateKinematics(photon, muonsOut);

      lEta = fakeMuon.eta;
      lPhi = fakeMuon.phi;

      calc->setMuon(fakeMuon, vertices, eventVars);
    }

    jetsOut.clear();
    for(unsigned iJ(0); iJ != jets.size; ++iJ)
      if(jet_isCand[iJ] && susy::deltaR(lEta, lPhi, jets.eta[iJ], jets.phi[iJ]) > 0.5)
        jetsOut.push_back(jets.at(iJ));

    fill(_filterType);
  }
}


class MCJetPhotonHLTIsoWeight : public GLEventWeight {
public:
  MCJetPhotonHLTIsoWeight() : GLEventWeight(0.443,
                                          std::sqrt(std::pow(0.030 / 0.443, 2.) +
                                                    std::pow(0.0016 / 0.0529, 2.) +
                                                    std::pow((0.443 - 0.381) / 0.443, 2.))
                                          ) {} // stat + fake fraction nonclosure + HLT emulation - match discrepancy
  ~MCJetPhotonHLTIsoWeight() {}
};

class JetPhotonHLTIsoWeight : public GLEventWeight {
public:
  JetPhotonHLTIsoWeight() : GLEventWeight() {}
  ~JetPhotonHLTIsoWeight() {}
  void setPhoton(susy::PhotonVars const& _photon, susy::VertexVarsArray const&, susy::SimpleEventProducer::EventVars const&)
  {
    weight = 0.72 + 482. * std::pow(_photon.pt, -1.73) - 13.5 * std::pow(_photon.pt, -0.666);
    relErr = 0.17; // nonclosure + emul-match + 10% on data measurement
  }
};

class MCJetPhotonWeight : public GLEventWeight {
public:
  MCJetPhotonWeight() : GLEventWeight(0.425,
                                    std::sqrt(std::pow(0.026 / 0.425, 2.) +
                                              std::pow(0.0025 / 0.0499, 2.))
                                    ) {} // stat + fake fraction nonclosure
  ~MCJetPhotonWeight() {}
};

class JetPhotonWeight : public GLEventWeight {
public:
  JetPhotonWeight() : GLEventWeight() {}
  ~JetPhotonWeight() {}
  void setPhoton(susy::PhotonVars const& _photon, susy::VertexVarsArray const&, susy::SimpleEventProducer::EventVars const&)
  {
//     weight = 0.450 + 3.82e9 * std::pow(_photon.pt, -6.37);
//     relErr = 0.11; // nonclosure + 10% on data measurement
    weight = 1.734 * std::pow(_photon.pt, -0.363);
    relErr = 0.25;
  }
};

class MCElePhotonFunctionalWeight : public GLEventWeight {
public:
  MCElePhotonFunctionalWeight() : GLEventWeight() {}
  ~MCElePhotonFunctionalWeight() {}
  void setPhoton(susy::PhotonVars const& _photon, susy::VertexVarsArray const& _vertices, susy::SimpleEventProducer::EventVars const&)
  {
    unsigned nV(0);
    int iPV(-1);
    for(unsigned iV(0); iV != _vertices.size; ++iV){
      if(!_vertices.isGood[iV]) continue;
      if(iPV < 0) iPV = iV;
      ++nV;
    }
    double fakerate(1. - 0.997839 * (1. - std::pow(1.153e-01 * _photon.pt + 1., -3.930e+00)) * (1. - 1.981e-01 * std::exp(-3.995e-01 * _vertices.nTracks[iPV])) * (1. - 1.268e-04 * nV));
    weight = fakerate / (1. - fakerate);
    relErr = 0.092;
  }
};

class ElePhotonFunctionalWeight : public GLEventWeight {
public:
  ElePhotonFunctionalWeight() : GLEventWeight() {}
  ~ElePhotonFunctionalWeight() {}
  void setPhoton(susy::PhotonVars const& _photon, susy::VertexVarsArray const& _vertices, susy::SimpleEventProducer::EventVars const&)
  {
    unsigned nV(0);
    int iPV(-1);
    for(unsigned iV(0); iV != _vertices.size; ++iV){
      if(!_vertices.isGood[iV]) continue;
      if(iPV < 0) iPV = iV;
      ++nV;
    }

    double fakerate(1. - 0.998456 * (1. - std::pow(1.110e-01 * _photon.pt + 1., -3.769e+00)) * (1. - 1.424e-01 * std::exp(-2.967e-01 * _vertices.nTracks[iPV])) * (1. - 3.154e-04 * nV));
    weight = fakerate / (1. - fakerate);
    relErr = 0.106;
  }
};

class ElePhotonConstantWeight : public GLEventWeight {
public:
  ElePhotonConstantWeight() : GLEventWeight()
  {
    double fakerate(8.49759907050168760e-03);
    weight = fakerate / (1. - fakerate);
    relErr = 0.106;
  }
  ~ElePhotonConstantWeight() {}
};

class JetElectronBinnedWeight : public GLEventWeight {
public:
  JetElectronBinnedWeight() :
    GLEventWeight()
  {
    source_ = TFile::Open("/afs/cern.ch/user/y/yiiyama/output/GammaL/jetEFake/scTrackMatchingFit.root");
    if(!source_ || source_->IsZombie()){
      delete source_;
      throw std::runtime_error("PhotonAndFakeElectron sources not found");
    }

    source_->GetObject("h_scales_pt_eta0", h_[0]);
    source_->GetObject("h_scales_pt_eta1", h_[1]);
  }
  ~JetElectronBinnedWeight()
  {
    delete source_;
  }
  void setElectron(susy::ElectronVars const& _electron, susy::VertexVarsArray const&, susy::SimpleEventProducer::EventVars const&)
  {
    // P(SID) / P(PID) binned in Pt

    int bin(h_[_electron.iSubdet]->FindFixBin(_electron.pt));
    if(bin > h_[_electron.iSubdet]->GetNbinsX()) bin = h_[_electron.iSubdet]->GetNbinsX();

    weight = h_[_electron.iSubdet]->GetBinContent(bin);
    relErr = h_[_electron.iSubdet]->GetBinError(bin) / weight;
  }

private:
  TFile* source_;
  TH1* h_[2];
};
