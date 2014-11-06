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
#include "TPRegexp.h"
#include "TEntryListFromFile.h"

#include <cmath>
#include <iostream>
#include <utility>
#include <algorithm>
#include <stdexcept>

#include "../ROOT/EventProcessor.h"

bool const ALTFAKE(false);

enum OutputTypes {
  oPhotonAndElectron,
  oPhotonAndMuon,
  oElePhotonAndElectron,
  oElePhotonAndMuon,
  oFakePhotonAndElectron,
  oFakePhotonAndMuon,
  oPhotonAndFakeElectron,
  oPhotonAndFakeMuon,
  oPhotonAndDimuon,
  oElePhotonAndDimuon,
  oFakePhotonAndDimuon,
  oElePhotonAndFakeElectron,
  oElePhotonAndFakeMuon,
  oFakePhotonAndFakeElectron,
  oFakePhotonAndFakeMuon,
  nOutputTypes
};

TString outputNames[] = {
  "PhotonAndElectron",
  "PhotonAndMuon",
  "ElePhotonAndElectron",
  "ElePhotonAndMuon",
  "FakePhotonAndElectron",
  "FakePhotonAndMuon",
  "PhotonAndFakeElectron",
  "PhotonAndFakeMuon",
  "PhotonAndDimuon",
  "ElePhotonAndDimuon",
  "FakePhotonAndDimuon",
  "ElePhotonAndFakeElectron",
  "ElePhotonAndFakeMuon",
  "FakePhotonAndFakeElectron",
  "FakePhotonAndFakeMuon"
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
  GLEventWeight() : EventWeight("GLEventWeight", 1., 0.) {}
  GLEventWeight(TString const& _name) : EventWeight(_name, 1., 0.) {}
  GLEventWeight(TString const& _name, double _w, double _r) : EventWeight(_name, _w, _r) {}
  ~GLEventWeight() {}
  virtual void setPhoton(susy::PhotonVars const&, susy::VertexVarsArray const&, susy::SimpleEventProducer::EventVars const&) {}
  virtual void setElectron(susy::ElectronVars const&, susy::VertexVarsArray const&, susy::SimpleEventProducer::EventVars const&) {}
  virtual void setMuon(susy::MuonVars const&, susy::VertexVarsArray const&, susy::SimpleEventProducer::EventVars const&) {}
};

class GLSkimProcessor : public EventProcessor {
public:
  GLSkimProcessor(Dataset const*, char const*);
  ~GLSkimProcessor() {} // ROOT closes files in case of crashes; deleting event list causes double free

  void process();
  void write();

private:
  bool preprocess();
  void processPhotonAndLepton(unsigned);
  void processElePhotonAndLepton(unsigned);
  void processFakePhotonAndLepton(unsigned);
  void processPhotonAndFakeLepton(unsigned);
  void processFakePhotonAndFakeLepton(unsigned);
  void processPhotonAndDimuon();
  void processElePhotonAndDimuon();
  void processFakePhotonAndDimuon();

  void photonEfficiencyScaleFactor(unsigned, unsigned, double&, double&);
  void electronEfficiencyScaleFactor(unsigned, double&, double&);
  void muonEfficiencyScaleFactor(unsigned, double&, double&);
  void photonInefficiencyScaleFactor(unsigned, unsigned, double&, double&);
  void electronInefficiencyScaleFactor(unsigned, double&, double&);
  void muonInefficiencyScaleFactor(unsigned, double&, double&);
  void getValueFromTable(double&, double&, TH1 const&, double, double = 0.);

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

  double maxDR;

  bool photon_muonIso[susy::NMAX];
  bool photon_electronIso[susy::NMAX];
  bool photon_matchGen[susy::NMAX];
  bool photon_matchGenE[susy::NMAX];
  float photon_genIso[susy::NMAX];
  bool photon_matchEHLT[2][susy::NMAX];
  bool photon_matchMHLT[susy::NMAX];
  bool electron_matchGen[susy::NMAX];
  bool electron_matchHLT[susy::NMAX];
  bool muon_matchGen[susy::NMAX];
  bool muon_matchHLT[susy::NMAX];

  TVector2 metJESUpV;
  TVector2 metJESDownV;

  TH1* idsfTable[3]; // photon, electron, muon
  TH1* ideffTable[3];
  TH1* hltsfTable[4]; // photon_e, photon_mu, electron, muon
  TH1* hlteffTable[4];
};

GLSkimProcessor::GLSkimProcessor(Dataset const* _dataset, char const* _outputDir) :
  EventProcessor(nOutputTypes, "eventVars", _dataset, _outputDir),
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
  cutflowM_(0),
  maxDR(0.8)
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

  std::fill_n(photon_muonIso, susy::NMAX, true);
  std::fill_n(photon_electronIso, susy::NMAX, true);
  std::fill_n(photon_matchGen, susy::NMAX, false);
  std::fill_n(photon_matchGenE, susy::NMAX, false);
  std::fill_n(photon_genIso, susy::NMAX, 0.);
  std::fill_n(photon_matchEHLT[0], susy::NMAX, false);
  std::fill_n(photon_matchEHLT[1], susy::NMAX, false);
  std::fill_n(photon_matchMHLT, susy::NMAX, false);
  std::fill_n(electron_matchGen, susy::NMAX, false);
  std::fill_n(electron_matchHLT, susy::NMAX, false);
  std::fill_n(muon_matchGen, susy::NMAX, false);
  std::fill_n(muon_matchHLT, susy::NMAX, false);

  for(unsigned iH(0); iH != 3; ++iH)
    idsfTable[iH] = ideffTable[iH] = 0;

  for(unsigned iH(0); iH != 4; ++iH)
    hltsfTable[iH] = hlteffTable[iH] = 0;

  for(unsigned iF(0); iF != nOutputTypes; ++iF)
    outputIndices[outputNames[iF]] = iF;
}

void
GLSkimProcessor::process()
{
  std::cout << "Processing " << dataset.name << std::endl;

  long nRows(inputChain.GetEntryList() ? static_cast<TEntryListFromFile*>(inputChain.GetEntryList())->GetEntries() : inputChain.GetEntries());
  if(nRows == 0){
    std::cerr << "Empty input" << std::endl;
    return;
  }

  std::cout << nRows << " entries" << std::endl;

  useElectronFilter =
    produceOutput[oPhotonAndElectron] ||
    produceOutput[oElePhotonAndElectron] ||
    produceOutput[oFakePhotonAndElectron] ||
    produceOutput[oPhotonAndFakeElectron] ||
    produceOutput[oElePhotonAndFakeElectron] ||
    produceOutput[oFakePhotonAndFakeElectron];
  useMuonFilter =
    produceOutput[oPhotonAndMuon] ||
    produceOutput[oElePhotonAndMuon] ||
    produceOutput[oFakePhotonAndMuon] ||
    produceOutput[oPhotonAndFakeMuon] ||
    produceOutput[oPhotonAndDimuon] ||
    produceOutput[oElePhotonAndDimuon] ||
    produceOutput[oFakePhotonAndDimuon] ||
    produceOutput[oElePhotonAndFakeMuon] ||
    produceOutput[oFakePhotonAndFakeMuon];

  ///////////////
  //// INPUT ////
  ///////////////

  TChain eventTree("eventVars");
  TChain objectTree("allObjects");

  TObjArray* fileNames(inputChain.GetListOfFiles());
  for(int iF(0); iF != fileNames->GetEntries(); ++iF){
    eventTree.Add(fileNames->At(iF)->GetTitle());
    objectTree.Add(fileNames->At(iF)->GetTitle());
  }

  inputChain.SetBranchStatus("*", 0);
  if(useElectronFilter){
    inputChain.SetBranchStatus("PhotonAndElectron", 1);
    inputChain.SetBranchStatus("ElePhotonAndElectron", 1);
    inputChain.SetBranchStatus("FakePhotonAndElectron", 1);
    inputChain.SetBranchStatus("PhotonAndFakeElectron", 1);
    if(inputChain.GetBranch("ElePhotonAndFakeElectron")) inputChain.SetBranchStatus("ElePhotonAndFakeElectron", 1);
    if(inputChain.GetBranch("FakePhotonAndFakeElectron")) inputChain.SetBranchStatus("FakePhotonAndFakeElectron", 1);
    inputChain.SetBranchStatus("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50", 1);
  }
  inputChain.SetBranchStatus("PhotonAndMuon", 1);
  inputChain.SetBranchStatus("ElePhotonAndMuon", 1);
  inputChain.SetBranchStatus("FakePhotonAndMuon", 1);
  inputChain.SetBranchStatus("PhotonAndFakeMuon", 1);
  if(inputChain.GetBranch("ElePhotonAndFAkeMuon")) inputChain.SetBranchStatus("ElePhotonAndFakeMuon", 1);
  if(inputChain.GetBranch("FakePhotonAndFakeMuon")) inputChain.SetBranchStatus("FakePhotonAndFakeMuon", 1);
  inputChain.SetBranchStatus("HLT_Mu22_Photon22_CaloIdL", 1);

  bool PhotonAndElectron(false);
  bool ElePhotonAndElectron(false);
  bool FakePhotonAndElectron(false);
  bool PhotonAndFakeElectron(false);
  bool ElePhotonAndFakeElectron(false);
  bool FakePhotonAndFakeElectron(false);
  bool PhotonAndMuon(false);
  bool ElePhotonAndMuon(false);
  bool FakePhotonAndMuon(false);
  bool PhotonAndFakeMuon(false);
  bool ElePhotonAndFakeMuon(false);
  bool FakePhotonAndFakeMuon(false);
  bool elHLT(false);
  bool muHLT(false);

  if(useElectronFilter){
    inputChain.SetBranchAddress("PhotonAndElectron", &PhotonAndElectron);
    inputChain.SetBranchAddress("ElePhotonAndElectron", &ElePhotonAndElectron);
    inputChain.SetBranchAddress("FakePhotonAndElectron", &FakePhotonAndElectron);
    inputChain.SetBranchAddress("PhotonAndFakeElectron", &PhotonAndFakeElectron);
    if(inputChain.GetBranch("ElePhotonAndFakeElectron")) inputChain.SetBranchAddress("ElePhotonAndFakeElectron", &ElePhotonAndFakeElectron);
    if(inputChain.GetBranch("FakePhotonAndFakeElectron")) inputChain.SetBranchAddress("FakePhotonAndFakeElectron", &FakePhotonAndFakeElectron);
    inputChain.SetBranchAddress("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50", &elHLT);
  }
  inputChain.SetBranchAddress("PhotonAndMuon", &PhotonAndMuon);
  inputChain.SetBranchAddress("ElePhotonAndMuon", &ElePhotonAndMuon);
  inputChain.SetBranchAddress("FakePhotonAndMuon", &FakePhotonAndMuon);
  inputChain.SetBranchAddress("PhotonAndFakeMuon", &PhotonAndFakeMuon);
  if(inputChain.GetBranch("ElePhotonAndFakeMuon")) inputChain.SetBranchAddress("ElePhotonAndFakeMuon", &ElePhotonAndFakeMuon);
  if(inputChain.GetBranch("FakePhotonAndFakeMuon")) inputChain.SetBranchAddress("FakePhotonAndFakeMuon", &FakePhotonAndFakeMuon);
  inputChain.SetBranchAddress("HLT_Mu22_Photon22_CaloIdL", &muHLT);

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
  objectTree.SetBranchAddress("electron.hltEG22CaloId10Iso50TrackIsoDoubleLastFilterUnseeded", electron_matchHLT);
  muons.setAddress(objectTree);
  objectTree.SetBranchAddress("muon.isCand", muon_isCand);
  objectTree.SetBranchAddress("muon.isFake", muon_isFake);
  objectTree.SetBranchAddress("muon.hltL1Mu3p5EG12L3Filtered22", muon_matchHLT);
  photons.setAddress(objectTree);
  objectTree.SetBranchAddress("photon.hltEG36CaloId10Iso50HcalIsoLastFilter", photon_matchEHLT[0]);
  objectTree.SetBranchAddress("photon.hltEG22CaloId10Iso50TrackIsoDoubleLastFilterUnseeded", photon_matchEHLT[1]);
  objectTree.SetBranchAddress("photon.hltMu22Photon22CaloIdLHEFilter", photon_matchMHLT);
  objectTree.SetBranchAddress("photon.isCand", photon_isCand);
  objectTree.SetBranchAddress("photon.isFake", photon_isFake);
  objectTree.SetBranchAddress("photon.isEle", photon_isEle);
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
    if(iF == oPhotonAndMuon || iF == oElePhotonAndMuon || iF == oFakePhotonAndMuon || iF == oPhotonAndFakeMuon ||
       iF == oPhotonAndDimuon || iF == oElePhotonAndDimuon || iF == oFakePhotonAndDimuon){
      muonsOut.setBranches(*eventList[iF]);
      eventList[iF]->Branch("muon.matchGen", muonsOut_matchGen, "matchGen[muon.size]/O");
    }
    jetsOut.setBranches(*eventList[iF]);
  }

  /////////////////
  //// CUTFLOW ////
  /////////////////

  if(produceOutput[oPhotonAndElectron] || produceOutput[oPhotonAndMuon]){
    TChain cutTree("cutTree");
    for(int iF(0); iF != fileNames->GetEntries(); ++iF)
      cutTree.Add(fileNames->At(iF)->GetTitle());

    enum Flag {
      kGoodLumi,
      kHLTE,
      kHLTM,
      kMetFilter,
      kGoodVertex,
      kGoodPhotonE,
      kGoodPhotonM,
      kGoodElectron,
      kGoodMuon,
      nFlags
    };
    TString flagNames[] = {
      "GoodLumi",
      "HLTE",
      "HLTM",
      "MetFilter",
      "GoodVertex",
      "GoodPhotonE",
      "GoodPhotonM",
      "GoodElectron",
      "GoodMuon"
    };
    unsigned const nCounterBins(1 << nFlags);

    TString countExpr;
    for(unsigned iC(0); iC != nFlags; ++iC){
      countExpr += TString::Format("%d * ", iC + 1);
      countExpr += flagNames[iC];
      if(iC != nFlags - 1) countExpr += " + ";
    }
    
    cutTree.Draw(countExpr + TString::Format(">>counter(%d,0,%d)", nCounterBins, nCounterBins), "", "goff");
    TH1* counter(static_cast<TH1*>(gDirectory->Get("counter")));

    TString signalCuts[] = {
      "AllEvents",
      "GoodLumi",
      "HLT",
      "MetFilter",
      "GoodVertex",
      "GoodPhoton",
      "GoodLepton",
      "FSRVeto"
    };
    unsigned const nSignalCuts(sizeof(signalCuts) / sizeof(TString));

    if(produceOutput[oPhotonAndElectron]){
      eventList[oPhotonAndElectron]->GetCurrentFile()->cd();
      cutflowE_ = new TH1D("cutflowE", "Cutflow (Electron Channel)", nSignalCuts, 0., nSignalCuts);
      for(unsigned iC(0); iC != nSignalCuts; ++iC)
        cutflowE_->GetXaxis()->SetBinLabel(iC + 1, signalCuts[iC]);

      unsigned iFlags[nSignalCuts] = {nFlags, kGoodLumi, kHLTE, kMetFilter, kGoodVertex, kGoodPhotonE, kGoodElectron, nFlags};
      double entries[nSignalCuts];
      std::fill_n(entries, nSignalCuts, 0);
      for(unsigned iX(1); iX <= nCounterBins; ++iX){
        double cont(counter->GetBinContent(iX));
        for(unsigned iS(0); iS != nSignalCuts; ++iS)
          if(((iX >> iFlags[iS]) & 1) == 1) entries[iS] += cont;
      }

      cutflowE_->SetBinContent(1, counter->GetEntries());
      for(unsigned iS(1); iS != nSignalCuts - 1; ++iS)
        cutflowE_->SetBinContent(iS + 1, entries[iS]);
    }
    if(produceOutput[oPhotonAndMuon]){
      eventList[oPhotonAndMuon]->GetCurrentFile()->cd();
      cutflowM_ = new TH1D("cutflowM", "Cutflow (Muon Channel)", nSignalCuts, 0., nSignalCuts);
      for(unsigned iC(0); iC != nSignalCuts; ++iC)
        cutflowM_->GetXaxis()->SetBinLabel(iC + 1, signalCuts[iC]);

      unsigned iFlags[nSignalCuts] = {nFlags, kGoodLumi, kHLTM, kMetFilter, kGoodVertex, kGoodPhotonM, kGoodMuon};
      double entries[nSignalCuts];
      std::fill_n(entries, nSignalCuts, 0);
      for(unsigned iX(1); iX <= nCounterBins; ++iX){
        double cont(counter->GetBinContent(iX));
        for(unsigned iS(0); iS != nSignalCuts; ++iS)
          if(((iX >> iFlags[iS]) & 1) == 1) entries[iS] += cont;
      }

      cutflowM_->SetBinContent(1, counter->GetEntries());
      for(unsigned iS(1); iS != nSignalCuts - 1; ++iS)
        cutflowM_->SetBinContent(iS + 1, entries[iS]);
    }

    delete counter;
  }

  ////////////////////////////////
  //// EFFICIENCY CORRECTIONS ////
  ////////////////////////////////

  if(produceOutput[oPhotonAndElectron] || produceOutput[oPhotonAndMuon] || produceOutput[oPhotonAndDimuon]){

    TFile* simIdSFSource(0);
    if(dataset.dataType == Dataset::kFastSim){
      simIdSFSource = TFile::Open("/afs/cern.ch/user/y/yiiyama/output/GammaL/main/idEfficiencies/full_to_fast.root");
      if(!simIdSFSource)
        throw std::runtime_error("ID Full/FastSim SF source not available");
    }

    if(dataset.dataType != Dataset::kRealData){
      TFile* idSFSource(TFile::Open("/afs/cern.ch/user/y/yiiyama/output/GammaL/main/idEfficiencies/scalefactors.root"));
      if(!idSFSource)
        throw std::runtime_error("ID SF source not available");

      idsfTable[0] = static_cast<TH1*>(idSFSource->Get("photon"));
      idsfTable[0]->SetDirectory(0);
      if(simIdSFSource) idsfTable[0]->Multiply(static_cast<TH1*>(simIdSFSource->Get("photon_sf")));
      if(useElectronFilter){
        idsfTable[1] = static_cast<TH1*>(idSFSource->Get("electron"));
        idsfTable[1]->SetDirectory(0);
        if(simIdSFSource) idsfTable[1]->Multiply(static_cast<TH1*>(simIdSFSource->Get("electron_sf")));
      }
      if(useMuonFilter){
        idsfTable[2] = static_cast<TH1*>(idSFSource->Get("muon"));
        idsfTable[2]->SetDirectory(0);
        if(simIdSFSource) idsfTable[2]->Multiply(static_cast<TH1*>(simIdSFSource->Get("muon_sf")));
      }

      delete idSFSource;
      delete simIdSFSource;

      TFile* hltSFSource(TFile::Open("/afs/cern.ch/user/y/yiiyama/output/GammaL/main/hltEfficiencies/scalefactors.root"));
      if(!hltSFSource)
        throw std::runtime_error("HLT SF source not available");

      if(useElectronFilter){
        hltsfTable[0] = static_cast<TH1*>(hltSFSource->Get("photon_e"));
        hltsfTable[0]->SetDirectory(0);
        hltsfTable[2] = static_cast<TH1*>(hltSFSource->Get("electron"));
        hltsfTable[2]->SetDirectory(0);
      }
      if(useMuonFilter){
        TH1* mueg(static_cast<TH1*>(hltSFSource->Get("mueg")));
        hltsfTable[1] = static_cast<TH1*>(hltSFSource->Get("photon_mu"));
        hltsfTable[1]->Multiply(mueg);
        hltsfTable[1]->SetDirectory(0);
        hltsfTable[3] = static_cast<TH1*>(hltSFSource->Get("muon"));
        hltsfTable[3]->SetDirectory(0);
        delete mueg;
      }

      delete hltSFSource;

      TString inputName;
      if(dataset.name.Contains("TChiwg") || dataset.name.Contains("T5wg") || dataset.name.Contains("Spectra_gW")){
        TPRegexp pat("^([^_]+(?:|_[a-zA-Z]+))_([0-9]+(?:|_[0-9]+))$");
        TObjArray* matches(pat.MatchS(dataset.name));
        if(matches->GetEntries() == 0)
          throw std::runtime_error(("Invalid signal dataset name " + dataset.name).Data());

        inputName = matches->At(1)->GetName();
        inputName += "/";
        inputName += matches->At(2)->GetName();
      }
      else
        inputName = dataset.name;

      TFile* idEffSource(TFile::Open("/afs/cern.ch/user/y/yiiyama/output/GammaL/main/idEfficiencies/" + inputName + ".root"));
      if(!idEffSource)
        throw std::runtime_error("ID Eff source not available");

      ideffTable[0] = static_cast<TH1*>(idEffSource->Get("photon_eff"));
      ideffTable[0]->SetDirectory(0);
      if(useElectronFilter){
        ideffTable[1] = static_cast<TH1*>(idEffSource->Get("electron_eff"));
        ideffTable[1]->SetDirectory(0);
      }
      if(useMuonFilter){
        ideffTable[2] = static_cast<TH1*>(idEffSource->Get("muon_eff"));
        ideffTable[2]->SetDirectory(0);
      }

      delete idEffSource;

      TFile* hltEffSource(TFile::Open("/afs/cern.ch/user/y/yiiyama/output/GammaL/main/hltEfficiencies/" + inputName + ".root"));
      if(!hltEffSource)
        throw std::runtime_error("HLT Eff source not available");

      if(useElectronFilter){
        hlteffTable[0] = static_cast<TH1*>(hltEffSource->Get("photon_e_eff"));
        hlteffTable[0]->SetDirectory(0);
        hlteffTable[2] = static_cast<TH1*>(hltEffSource->Get("electron_eff"));
        hlteffTable[2]->SetDirectory(0);
      }
      if(useMuonFilter){
        hlteffTable[1] = static_cast<TH1*>(hltEffSource->Get("photon_mu_eff"));
        hlteffTable[1]->SetDirectory(0);
        hlteffTable[3] = static_cast<TH1*>(hltEffSource->Get("muon_eff"));
        hlteffTable[3]->SetDirectory(0);
      }

      delete hltEffSource;
    }
  }

  ////////////////////
  //// START LOOP ////
  ////////////////////
  
  long iRow(0);
  long nIncrements(0);
  long iEntry;
  while(iRow != nRows){
    try{
      iEntry = inputChain.GetEntryNumber(iRow);
      if(iEntry < 0) break;

      inputChain.GetEntry(iEntry);
      if(nIncrements++ % 1000 == 0) (std::cout << "\r" << iRow).flush();

      iRow += dataset.prescale;

      if(!(elHLT && ((PhotonAndElectron && produceOutput[oPhotonAndElectron]) ||
                     (ElePhotonAndElectron && produceOutput[oElePhotonAndElectron]) ||
                     (FakePhotonAndElectron && produceOutput[oFakePhotonAndElectron]) ||
                     (PhotonAndFakeElectron && produceOutput[oPhotonAndFakeElectron]) ||
                     (ElePhotonAndFakeElectron && produceOutput[oElePhotonAndFakeElectron]) ||
                     (FakePhotonAndFakeElectron && produceOutput[oFakePhotonAndFakeElectron])
                     )) &&
         !(muHLT && ((PhotonAndMuon && (produceOutput[oPhotonAndMuon] || produceOutput[oPhotonAndDimuon])) ||
                     (ElePhotonAndMuon && produceOutput[oElePhotonAndMuon]) ||
                     (FakePhotonAndMuon && produceOutput[oFakePhotonAndMuon]) ||
                     (PhotonAndFakeMuon && produceOutput[oPhotonAndFakeMuon]) ||
                     (ElePhotonAndFakeMuon && produceOutput[oElePhotonAndFakeMuon]) ||
                     (FakePhotonAndFakeMuon && produceOutput[oFakePhotonAndFakeMuon])))) continue;


      eventTree.GetEntry(iEntry);
      objectTree.GetEntry(iEntry);

      if(!preprocess()) continue;

      if(PhotonAndMuon){
        if(produceOutput[oPhotonAndMuon]) processPhotonAndLepton(oPhotonAndMuon);
        if(produceOutput[oPhotonAndDimuon]) processPhotonAndDimuon();
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
      
      if(ElePhotonAndMuon && produceOutput[oElePhotonAndDimuon])
        processElePhotonAndDimuon();

      if(FakePhotonAndElectron && produceOutput[oFakePhotonAndElectron])
        processFakePhotonAndLepton(oFakePhotonAndElectron);

      if(FakePhotonAndMuon && produceOutput[oFakePhotonAndMuon])
        processFakePhotonAndLepton(oFakePhotonAndMuon);

      if(FakePhotonAndMuon && produceOutput[oFakePhotonAndDimuon])
        processFakePhotonAndDimuon();

      if(PhotonAndFakeElectron && produceOutput[oPhotonAndFakeElectron])
        processPhotonAndFakeLepton(oPhotonAndFakeElectron);

      if(PhotonAndFakeMuon && produceOutput[oPhotonAndFakeMuon])
        processPhotonAndFakeLepton(oPhotonAndFakeMuon);

      if(ElePhotonAndFakeElectron && produceOutput[oElePhotonAndFakeElectron])
        processFakePhotonAndFakeLepton(oElePhotonAndFakeElectron);

      if(ElePhotonAndFakeMuon && produceOutput[oElePhotonAndFakeMuon])
        processFakePhotonAndFakeLepton(oElePhotonAndFakeMuon);

      if(FakePhotonAndFakeElectron && produceOutput[oFakePhotonAndFakeElectron])
        processFakePhotonAndFakeLepton(oFakePhotonAndFakeElectron);

      if(FakePhotonAndFakeMuon && produceOutput[oFakePhotonAndFakeMuon])
        processFakePhotonAndFakeLepton(oFakePhotonAndFakeMuon);
    }
    catch(std::exception& _ex){
      std::cerr << std::endl << _ex.what() << std::endl;
    }
  }
  std::cout << std::endl;
}

void
GLSkimProcessor::write()
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
GLSkimProcessor::preprocess()
{
  for(unsigned iP(0); iP != photons.size; ++iP){
    unsigned iL(0);
    for(; iL != muons.size; ++iL)
      if(susy::deltaR(muons.eta[iL], muons.phi[iL], photons.eta[iP], photons.phi[iP]) < 0.3) break;
    photon_muonIso[iP] = iL == muons.size;
    TVector3 caloPosition(photons.caloX[iP], photons.caloY[iP], photons.caloZ[iP]);
    iL = 0;
    for(; iL != electrons.size; ++iL){
      double dR(TVector3(electrons.caloX[iL], electrons.caloY[iL], electrons.caloZ[iL]).DeltaR(caloPosition));
      if(dR > 0.02 && dR < 0.3) break;
    }
    photon_electronIso[iP] = iL == electrons.size;
  }

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
GLSkimProcessor::processPhotonAndLepton(unsigned _outputType)
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
    if(photons.iSubdet[iP] != 0) continue;
    if(photons.pt[iP] < 40.) break;
    if(!photon_muonIso[iP]) continue;
    if(!photon_electronIso[iP]) continue;
    if(_outputType == oPhotonAndElectron){
      if(!photon_matchEHLT[0][iP] || !photon_matchEHLT[1][iP]) continue;
    }
    else{
      if(!photon_matchMHLT[iP]) continue;
    }
    if(!photon_isCand[iP]){
      if(photon_matchGen[iP])
        photonInefficiencyScaleFactor(iP, _outputType == oPhotonAndElectron ? 0 : 1, effScale, photonRelErr);
      continue;
    }
    if(photon_matchGen[iP])
      photonEfficiencyScaleFactor(iP, _outputType == oPhotonAndElectron ? 0 : 1, effScale, photonRelErr);

    photonsOut_matchGen[photonsOut.size] = photon_matchGen[iP];
    photonsOut_matchGenE[photonsOut.size] = photon_matchGenE[iP];
    photonsOut_genIso[photonsOut.size] = photon_genIso[iP];

    photonsOut.push_back(photons.at(iP));
  }

  if(photonsOut.size == 0) return;

  if(_outputType == oPhotonAndElectron){
    electronsOut.clear();

    for(unsigned iL(0); iL != electrons.size; ++iL){
      if(electrons.iSubdet[iL] == -1) continue;
      if(electrons.pt[iL] < 25.) break;
      if(!electron_matchHLT[iL]) continue;
      if(!electron_isCand[iL]){
        if(electron_matchGen[iL])
          electronInefficiencyScaleFactor(iL, effScale, electronRelErr);
        continue;
      }
      if(electron_matchGen[iL])
        electronEfficiencyScaleFactor(iL, effScale, electronRelErr);

      electronsOut_matchGen[electronsOut.size] = electron_matchGen[iL];
      electronsOut.push_back(electrons.at(iL));
    }

    if(electronsOut.size == 0) return;
  }
  else if(_outputType == oPhotonAndMuon){
    muonsOut.clear();

    for(unsigned iL(0); iL != muons.size; ++iL){
      if(muons.iSubdet[iL] == -1) continue;
      if(muons.pt[iL] < 25.) break;
      if(!muon_matchHLT[iL]) continue;
      if(!muon_isCand[iL]){
        if(muon_matchGen[iL])
          muonInefficiencyScaleFactor(iL, effScale, muonRelErr);
        continue;
      }
      if(muon_matchGen[iL])
        muonEfficiencyScaleFactor(iL, effScale, muonRelErr);

      muonsOut_matchGen[muonsOut.size] = muon_matchGen[iL];
      muonsOut.push_back(muons.at(iL));
    }

    if(muonsOut.size == 0) return;
  }

  scaleRelErr2 += photonRelErr * photonRelErr + electronRelErr * electronRelErr + muonRelErr * muonRelErr;
  scaleErr = effScale * std::sqrt(scaleRelErr2);

  susy::PhotonVars photon(photonsOut.at(0));

  if(_outputType == oPhotonAndElectron){
    if(susy::deltaR(photonsOut.eta[0], photonsOut.phi[0], electronsOut.eta[0], electronsOut.phi[0]) < maxDR) return;
    calculateKinematics(photon, electronsOut);

    if(cutflowE_) cutflowE_->Fill(cutflowE_->GetNbinsX() - 0.5);
  }
  else if(_outputType == oPhotonAndMuon){
    if(susy::deltaR(photonsOut.eta[0], photonsOut.phi[0], muonsOut.eta[0], muonsOut.phi[0]) < maxDR) return;
    calculateKinematics(photon, muonsOut);

    if(cutflowM_) cutflowM_->Fill(cutflowM_->GetNbinsX() - 0.5);
  }

  setJets();
  
  if(weightCalc[_outputType]){
    GLEventWeight* calc(static_cast<GLEventWeight*>(weightCalc[_outputType]));
    calc->setPhoton(photon, vertices, eventVars);
    if(_outputType == oPhotonAndElectron)
      calc->setElectron(electronsOut.at(0), vertices, eventVars);
    if(_outputType == oPhotonAndMuon)
      calc->setMuon(muonsOut.at(0), vertices, eventVars);
  }

  fill(_outputType);
}

void
GLSkimProcessor::processElePhotonAndLepton(unsigned _outputType)
{
  if(_outputType == oElePhotonAndMuon){
    muonsOut.clear();

    for(unsigned iL(0); iL != muons.size; ++iL){
      if(muons.iSubdet[iL] == -1) continue;
      if(muons.pt[iL] < 25.) break;
      if(!muon_matchHLT[iL]) continue;
      if(!muon_isCand[iL]) continue;
      muonsOut_matchGen[muonsOut.size] = muon_matchGen[iL];
      muonsOut.push_back(muons.at(iL));
    }

    if(muonsOut.size == 0) return;
  }
  else if(_outputType != oElePhotonAndElectron)
    throw std::runtime_error("Incorrect filter type");

  for(unsigned iEP(0); iEP != photons.size; ++iEP){
    if(photons.iSubdet[iEP] != 0) continue;
    if(photons.pt[iEP] < 40.) break;
    if(!photon_isEle[iEP]) continue;
    if(!photon_muonIso[iEP]) continue;
    if(!photon_electronIso[iEP]) continue;
    if(_outputType == oElePhotonAndElectron){
      if(!photon_matchEHLT[0][iEP] || !photon_matchEHLT[1][iEP]) continue;
    }
    else{
      if(!photon_matchMHLT[iEP]) continue;
    }

    susy::PhotonVars elePhoton(photons.at(iEP));

    TVector3 caloPosition(elePhoton.caloX, elePhoton.caloY, elePhoton.caloZ);

    if(_outputType == oElePhotonAndElectron){
      electronsOut.clear();

      for(unsigned iL(0); iL != electrons.size; ++iL){
        if(electrons.iSubdet[iL] == -1) continue;
        if(electrons.pt[iL] < 25.) break;
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
      if(susy::deltaR(elePhoton.eta, elePhoton.phi, electronsOut.eta[0], electronsOut.phi[0]) < maxDR) continue;
      calculateKinematics(elePhoton, electronsOut);
    }
    else if(_outputType == oElePhotonAndMuon){
      if(susy::deltaR(elePhoton.eta, elePhoton.phi, muonsOut.eta[0], muonsOut.phi[0]) < maxDR) continue;
      calculateKinematics(elePhoton, muonsOut);
    }

    photonsOut.clear();

    photonsOut_matchGen[photonsOut.size] = photon_matchGen[iEP];
    photonsOut_matchGenE[photonsOut.size] = photon_matchGenE[iEP];
    photonsOut_genIso[photonsOut.size] = photon_genIso[iEP];
    photonsOut.push_back(elePhoton);

    setJets();

    if(weightCalc[_outputType])
      static_cast<GLEventWeight*>(weightCalc[_outputType])->setPhoton(elePhoton, vertices, eventVars);

    fill(_outputType);
  }
}

void
GLSkimProcessor::processFakePhotonAndLepton(unsigned _outputType)
{
  if(_outputType == oFakePhotonAndElectron){
    electronsOut.clear();

    for(unsigned iL(0); iL != electrons.size; ++iL){
      if(electrons.iSubdet[iL] == -1) continue;
      if(electrons.pt[iL] < 25.) break;
      if(!electron_matchHLT[iL]) continue;
      if(!electron_isCand[iL]) continue;
      electronsOut_matchGen[electronsOut.size] = electron_matchGen[iL];
      electronsOut.push_back(electrons.at(iL));
    }

    if(electronsOut.size == 0) return;
  }
  else if(_outputType == oFakePhotonAndMuon){
    muonsOut.clear();

    for(unsigned iL(0); iL != muons.size; ++iL){
      if(muons.iSubdet[iL] == -1) continue;
      if(muons.pt[iL] < 25.) break;
      if(!muon_matchHLT[iL]) continue;
      if(!muon_isCand[iL]) continue;
      muonsOut_matchGen[muonsOut.size] = muon_matchGen[iL];
      muonsOut.push_back(muons.at(iL));
    }

    if(muonsOut.size == 0) return;
  }
  else
    throw std::runtime_error("Incorrect filter type");

  for(unsigned iFP(0); iFP != photons.size; ++iFP){
    if(photons.iSubdet[iFP] != 0) continue;
    if(photons.pt[iFP] < 40.) break;
    if(!photon_isFake[iFP]) continue;
    if(!photon_muonIso[iFP]) continue;
    if(!photon_electronIso[iFP]) continue;
    if(_outputType == oFakePhotonAndElectron){
      if(!photon_matchEHLT[0][iFP] || !photon_matchEHLT[1][iFP]) continue;
    }
    else{
      if(!photon_matchMHLT[iFP]) continue;
    }
    if(!ALTFAKE){
      if(photons.hOverE[iFP] > 0.05 || photons.sigmaIetaIeta[iFP] > 0.014) continue;
      if(_outputType == oFakePhotonAndMuon && (photons.chargedHadronIso[iFP] > 15. || photons.neutralHadronIso[iFP] > 3.5 || photons.photonIso[iFP] > 1.3)) continue;
    }

    susy::PhotonVars fakePhoton(photons.at(iFP));

    if(_outputType == oFakePhotonAndElectron){
      if(susy::deltaR(fakePhoton.eta, fakePhoton.phi, electronsOut.eta[0], electronsOut.phi[0]) < maxDR) continue;
      calculateKinematics(fakePhoton, electronsOut);
    }
    else if(_outputType == oFakePhotonAndMuon){
      if(susy::deltaR(fakePhoton.eta, fakePhoton.phi, muonsOut.eta[0], muonsOut.phi[0]) < maxDR) continue;
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
      static_cast<GLEventWeight*>(weightCalc[_outputType])->setPhoton(fakePhoton, vertices, eventVars);

    fill(_outputType);
  }
}

void
GLSkimProcessor::processPhotonAndFakeLepton(unsigned _outputType)
{
  if(_outputType != oPhotonAndFakeElectron && _outputType != oPhotonAndFakeMuon)
    throw std::runtime_error("Incorrect filter type");

  photonsOut.clear();

  for(unsigned iP(0); iP != photons.size; ++iP){
    if(photons.iSubdet[iP] != 0) continue;
    if(photons.pt[iP] < 40.) break;
    if(!photon_isCand[iP]) continue;
    if(!photon_muonIso[iP]) continue;
    if(!photon_electronIso[iP]) continue;
    if(_outputType == oPhotonAndFakeElectron){
      if(!photon_matchEHLT[0][iP] || !photon_matchEHLT[1][iP]) continue;
    }
    else{
      if(!photon_matchMHLT[iP]) continue;
    }
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

  GLEventWeight* calc(static_cast<GLEventWeight*>(weightCalc[_outputType]));

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

    if(_outputType == oPhotonAndFakeElectron){
      if(electrons.iSubdet[iFL] == -1) continue;
      if(electrons.pt[iFL] < 25.) break;
      if(!electron_matchHLT[iFL]) continue;
      if(!electron_isFake[iFL]) continue;

      if(electrons.combRelIso[iFL] * electrons.pt[iFL] < 10.) continue;

      susy::ElectronVars fakeElectron(electrons.at(iFL));

      if(!ALTFAKE){
        susy::ObjectSelector::isGoodElectron(fakeElectron, susy::ElMedium12, &elIdResults);
        if((elIdResults & elBaseline) != elBaseline) continue;
        if(elIdResults[susy::ElDeltaEta] && elIdResults[susy::ElDeltaPhi]) continue;
      }

      if(susy::deltaR(photon, fakeElectron) < maxDR) continue;

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
      if(muons.pt[iFL] < 25.) break;
      if(!muon_matchHLT[iFL]) continue;
      if(!muon_isFake[iFL]) continue;

      if(muons.combRelIso[iFL] < 0.15) continue;
      if(!ALTFAKE && muons.combRelIso[iFL] > 0.6) continue;

      susy::MuonVars fakeMuon(muons.at(iFL));

      if(!ALTFAKE){
        susy::ObjectSelector::isGoodMuon(fakeMuon, susy::MuTight12, &muIdResults);
        if((muIdResults & muBaseline) != muBaseline) continue;
      }

      if(susy::deltaR(photon, fakeMuon) < maxDR) continue;

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
GLSkimProcessor::processFakePhotonAndFakeLepton(unsigned _outputType)
{
  bool elePhotonEvent(false);
  bool electronEvent(false);
  switch(_outputType){
  case oElePhotonAndFakeElectron:
    elePhotonEvent = true;
    //fallthrough
  case oFakePhotonAndFakeElectron:
    electronEvent = true;
    break;
  case oElePhotonAndFakeMuon:
    elePhotonEvent = true;
    //fallthrough
  case oFakePhotonAndFakeMuon:
    break;
  default:
    throw std::runtime_error("Incorrect filter type");
  }

  std::bitset<susy::nElectronCriteria> elIdResults;
  std::bitset<susy::nElectronCriteria> elBaseline(susy::ObjectSelector::elReferences[susy::ElMedium12]);
  elBaseline.reset(susy::ElCombIso);
  elBaseline.reset(susy::ElDeltaEta);
  elBaseline.reset(susy::ElDeltaPhi);
  std::bitset<susy::nMuonCriteria> muIdResults;
  std::bitset<susy::nMuonCriteria> muBaseline(susy::ObjectSelector::muReferences[susy::MuTight12]);
  muBaseline.reset(susy::MuCombIso);

  GLEventWeight* calc(static_cast<GLEventWeight*>(weightCalc[_outputType]));

  for(unsigned iFP(0); iFP != photons.size; ++iFP){
    if(photons.iSubdet[iFP] != 0) continue;
    if(photons.pt[iFP] < 40.) break;
    if(elePhotonEvent){
      if(!photon_isEle[iFP]) continue;
    }
    else{
      if(!photon_isFake[iFP]) continue;
      if(photons.hOverE[iFP] > 0.05 || photons.sigmaIetaIeta[iFP] > 0.014) continue;
      if(!electronEvent && (photons.chargedHadronIso[iFP] > 15. || photons.neutralHadronIso[iFP] > 3.5 || photons.photonIso[iFP] > 1.3)) continue;
    }
    if(!photon_muonIso[iFP]) continue;
    if(!photon_electronIso[iFP]) continue;
    if(electronEvent){
      if(!photon_matchEHLT[0][iFP] || !photon_matchEHLT[1][iFP]) continue;
    }
    else{
      if(!photon_matchMHLT[iFP]) continue;
    }

    photonsOut_matchGen[0] = photon_matchGen[iFP];
    photonsOut_matchGenE[0] = photon_matchGenE[iFP];
    photonsOut_genIso[0] = photon_genIso[iFP];

    susy::PhotonVars fakePhoton(photons.at(iFP));

    TVector3 caloPosition(fakePhoton.caloX, fakePhoton.caloY, fakePhoton.caloZ);

    if(electronEvent){
      for(unsigned iFL(0); iFL != electrons.size; ++iFL){
        if(electrons.iSubdet[iFL] == -1) continue;
        if(electrons.pt[iFL] < 25.) break;
        if(!electron_matchHLT[iFL]) continue;
        if(!electron_isFake[iFL]) continue;

        if(electrons.combRelIso[iFL] * electrons.pt[iFL] < 10.) continue;

        if(elePhotonEvent){
          if(electrons.superClusterIndex[iFL] == fakePhoton.superClusterIndex ||
             TVector3(electrons.caloX[iFL], electrons.caloY[iFL], electrons.caloZ[iFL]).DeltaR(caloPosition) < 0.02)
            continue;
        }

        susy::ElectronVars fakeElectron(electrons.at(iFL));

        susy::ObjectSelector::isGoodElectron(fakeElectron, susy::ElMedium12, &elIdResults);
        if((elIdResults & elBaseline) != elBaseline) continue;
        if(elIdResults[susy::ElDeltaEta] && elIdResults[susy::ElDeltaPhi]) continue;

        if(susy::deltaR(fakePhoton, fakeElectron) < maxDR) continue;

        photonsOut.clear();
        electronsOut.clear();

        electronsOut_matchGen[0] = electron_matchGen[iFL];

        photonsOut.push_back(fakePhoton);
        electronsOut.push_back(fakeElectron);

        calculateKinematics(fakePhoton, electronsOut);

        if(calc) calc->setElectron(fakeElectron, vertices, eventVars);

        // elePhoton is already vetoed in filter.cd
        unsigned nVeto(1);
        double etaVeto[2] = {fakeElectron.eta};
        double phiVeto[2] = {fakeElectron.phi};
        if(!elePhotonEvent){
          nVeto = 2;
          etaVeto[1] = fakePhoton.eta;
          phiVeto[1] = fakePhoton.phi;
        }
        setJets(nVeto, etaVeto, phiVeto);

        fill(_outputType);
      }
    }
    else{
      for(unsigned iFL(0); iFL != muons.size; ++iFL){
        if(muons.iSubdet[iFL] == -1) continue;
        if(muons.pt[iFL] < 25.) break;
        if(!muon_matchHLT[iFL]) continue;
        if(!muon_isFake[iFL]) continue;

        if(muons.combRelIso[iFL] < 0.15 || muons.combRelIso[iFL] > 0.6) continue;

        susy::MuonVars fakeMuon(muons.at(iFL));

        susy::ObjectSelector::isGoodMuon(fakeMuon, susy::MuTight12, &muIdResults);
        if((muIdResults & muBaseline) != muBaseline) continue;

        if(susy::deltaR(fakePhoton, fakeMuon) < maxDR) continue;

        photonsOut.clear();
        muonsOut.clear();

        muonsOut_matchGen[0] = muon_matchGen[iFL];

        photonsOut.push_back(fakePhoton);
        muonsOut.push_back(fakeMuon);

        calculateKinematics(fakePhoton, muonsOut);

        if(calc) calc->setMuon(fakeMuon, vertices, eventVars);

        unsigned nVeto(1);
        double etaVeto[2] = {fakeMuon.eta};
        double phiVeto[2] = {fakeMuon.phi};
        if(!elePhotonEvent){
          nVeto = 2;
          etaVeto[1] = fakePhoton.eta;
          phiVeto[1] = fakePhoton.phi;
        }
        setJets(nVeto, etaVeto, phiVeto);

        fill(_outputType);
      }
    }
  }
}

void
GLSkimProcessor::processPhotonAndDimuon()
{
  effScale = 1.;
  double photonRelErr(0.);
  double muonRelErr(0.);
  double scaleRelErr2(0.);

  photonsOut.clear();

  for(unsigned iP(0); iP != photons.size; ++iP){
    if(photons.iSubdet[iP] != 0) continue;
    if(photons.pt[iP] < 40.) break;
    if(!photon_muonIso[iP]) continue;
    if(!photon_electronIso[iP]) continue;
    if(!photon_matchMHLT[iP]) continue;
    if(!photon_isCand[iP]){
      if(photon_matchGen[iP])
        photonInefficiencyScaleFactor(iP, 1, effScale, photonRelErr);
      continue;
    }
    if(photon_matchGen[iP])
      photonEfficiencyScaleFactor(iP, 1, effScale, photonRelErr);

    photonsOut_matchGen[photonsOut.size] = photon_matchGen[iP];
    photonsOut_matchGenE[photonsOut.size] = photon_matchGenE[iP];
    photonsOut_genIso[photonsOut.size] = photon_genIso[iP];

    photonsOut.push_back(photons.at(iP));
  }

  if(photonsOut.size == 0) return;

  muonsOut.clear();

  // efficiency correction is not done right here - need SFs for loose muons to be precise

  int iCand(-1);
  for(unsigned iL(0); iL != muons.size; ++iL){
    if(muons.iSubdet[iL] == -1) continue;
    if(muons.pt[iL] < 25.) break;
    if(!muon_matchHLT[iL]) continue;
    if(!muon_isCand[iL]){
      if(muon_matchGen[iL])
        muonInefficiencyScaleFactor(iL, effScale, muonRelErr);
      continue;
    }
    if(muon_matchGen[iL])
      muonEfficiencyScaleFactor(iL, effScale, muonRelErr);

    muonsOut_matchGen[0] = muon_matchGen[iL];
    muonsOut.push_back(muons.at(iL));

    iCand = iL;
    break;
  }

  if(iCand == -1) return;

  for(unsigned iL(0); iL != muons.size; ++iL){
    if(muons.iSubdet[iL] == -1) continue;
    if(iL == unsigned(iCand)) continue;
    if(!muons.isLoose[iL]) continue;

    muonsOut_matchGen[muonsOut.size] = muon_matchGen[iL];
    muonsOut.push_back(muons.at(iL));
  }  

  if(muonsOut.size != 2) return;

  scaleRelErr2 += photonRelErr * photonRelErr + muonRelErr * muonRelErr;
  scaleErr = effScale * std::sqrt(scaleRelErr2);

  calculateKinematics(photonsOut.at(0), muonsOut);

  setJets();

  fill(oPhotonAndDimuon);
}

void
GLSkimProcessor::processElePhotonAndDimuon()
{
  muonsOut.clear();

  int iCand(-1);
  for(unsigned iL(0); iL != muons.size; ++iL){
    if(muons.iSubdet[iL] == -1) continue;
    if(muons.pt[iL] < 25.) break;
    if(!muon_matchHLT[iL]) continue;
    if(!muon_isCand[iL]) continue;
    muonsOut_matchGen[0] = muon_matchGen[iL];
    muonsOut.push_back(muons.at(iL));

    iCand = iL;
    break;
  }

  if(iCand == -1) return;

  for(unsigned iL(0); iL != muons.size; ++iL){
    if(muons.iSubdet[iL] == -1) continue;
    if(iL == unsigned(iCand)) continue;
    if(!muons.isLoose[iL]) continue;
    
    muonsOut_matchGen[muonsOut.size] = muon_matchGen[iL];
    muonsOut.push_back(muons.at(iL));
  }
  
  if(muonsOut.size != 2) return;

  for(unsigned iEP(0); iEP != photons.size; ++iEP){
    if(photons.iSubdet[iEP] != 0) continue;
    if(photons.pt[iEP] < 40.) break;
    if(!photon_isEle[iEP]) continue;
    if(!photon_muonIso[iEP]) continue;
    if(!photon_matchMHLT[iEP]) continue;

    susy::PhotonVars elePhoton(photons.at(iEP));

    calculateKinematics(elePhoton, muonsOut);

    photonsOut.clear();

    photonsOut_matchGen[photonsOut.size] = photon_matchGen[iEP];
    photonsOut_matchGenE[photonsOut.size] = photon_matchGenE[iEP];
    photonsOut_genIso[photonsOut.size] = photon_genIso[iEP];
    photonsOut.push_back(elePhoton);

    setJets();

    if(weightCalc[oElePhotonAndDimuon])
      static_cast<GLEventWeight*>(weightCalc[oElePhotonAndDimuon])->setPhoton(elePhoton, vertices, eventVars);

    fill(oElePhotonAndDimuon);
  }
}

void
GLSkimProcessor::processFakePhotonAndDimuon()
{
  muonsOut.clear();

  int iCand(-1);
  for(unsigned iL(0); iL != muons.size; ++iL){
    if(muons.iSubdet[iL] == -1) continue;
    if(muons.pt[iL] < 25.) break;
    if(!muon_matchHLT[iL]) continue;
    if(!muon_isCand[iL]) continue;
    muonsOut_matchGen[0] = muon_matchGen[iL];
    muonsOut.push_back(muons.at(iL));

    iCand = iL;
    break;
  }

  if(iCand == -1) return;

  for(unsigned iL(0); iL != muons.size; ++iL){
    if(muons.iSubdet[iL] == -1) continue;
    if(iL == unsigned(iCand)) continue;
    if(!muons.isLoose[iL]) continue;
    
    muonsOut_matchGen[muonsOut.size] = muon_matchGen[iL];
    muonsOut.push_back(muons.at(iL));
  }

  if(muonsOut.size != 2) return;

  for(unsigned iFP(0); iFP != photons.size; ++iFP){
    if(photons.iSubdet[iFP] != 0) continue;
    if(photons.pt[iFP] < 40.) break;
    if(!photon_isFake[iFP]) continue;
    if(!photon_muonIso[iFP]) continue;
    if(!photon_electronIso[iFP]) continue;
    if(!photon_matchMHLT[iFP]) continue;
    if(photons.hOverE[iFP] > 0.05 || photons.sigmaIetaIeta[iFP] > 0.014) continue;
    if(photons.chargedHadronIso[iFP] > 15. || photons.neutralHadronIso[iFP] > 3.5 || photons.photonIso[iFP] > 1.3) continue;

    susy::PhotonVars fakePhoton(photons.at(iFP));

    calculateKinematics(fakePhoton, muonsOut);

    photonsOut.clear();

    photonsOut_matchGen[photonsOut.size] = photon_matchGen[iFP];
    photonsOut_matchGenE[photonsOut.size] = photon_matchGenE[iFP];
    photonsOut_genIso[photonsOut.size] = photon_genIso[iFP];
    photonsOut.push_back(fakePhoton);

    double etaVeto[1] = {fakePhoton.eta};
    double phiVeto[1] = {fakePhoton.phi};
    setJets(1, etaVeto, phiVeto);
    
    if(weightCalc[oFakePhotonAndDimuon])
      static_cast<GLEventWeight*>(weightCalc[oFakePhotonAndDimuon])->setPhoton(fakePhoton, vertices, eventVars);

    fill(oFakePhotonAndDimuon);
  }
}

void
GLSkimProcessor::photonEfficiencyScaleFactor(unsigned _iP, unsigned _channel, double& _effScale, double& _relErr)
{
  double eta(TVector3(photons.caloX[_iP], photons.caloY[_iP], photons.caloZ[_iP]).Eta());

  double idsf(1.);
  double idsfErr(0.);
  getValueFromTable(idsf, idsfErr, *idsfTable[0], photons.pt[_iP], eta);

  double hltsf(1.);
  double hltsfErr(0.);
  if(_channel == 0){
    getValueFromTable(hltsf, hltsfErr, *hltsfTable[0], eta, nVtx);
    //   getValueFromTable(hltsf, hltsfErr, *hltsfTable[0], photons.pt[_iP], eta);
    if(dataset.dataType == Dataset::kFastSim) hltsf *= 0.96405; // Photon36IdIso
  }
  else{
    getValueFromTable(hltsf, hltsfErr, *hltsfTable[1], 0.);
    if(dataset.dataType == Dataset::kFastSim) hltsf *= 1.01549; // Mu3p5EG12
  }

  double sf(idsf * hltsf);

  _effScale *= sf;
  _relErr += std::sqrt(std::pow(idsfErr / idsf, 2.) + std::pow(hltsfErr / hltsf, 2.));
}

void
GLSkimProcessor::electronEfficiencyScaleFactor(unsigned _iE, double& _effScale, double& _relErr)
{
  double eta(TVector3(electrons.caloX[_iE], electrons.caloY[_iE], electrons.caloZ[_iE]).Eta());

  double idsf(1.);
  double idsfErr(0.);
  getValueFromTable(idsf, idsfErr, *idsfTable[1], electrons.pt[_iE], eta);

  double hltsf(1.);
  double hltsfErr(0.);
  getValueFromTable(hltsf, hltsfErr, *hltsfTable[2], eta, nVtx);
  //  getValueFromTable(hltsf, hltsfErr, *hltsfTable[2], electrons.pt[_iE], eta);
  if(dataset.dataType == Dataset::kFastSim) hltsf *= 0.96147; // Photon22IdIso

  double sf(idsf * hltsf);

  _effScale *= sf;
  _relErr += std::sqrt(std::pow(idsfErr / idsf, 2.) + std::pow(hltsfErr / hltsf, 2.));
}

void
GLSkimProcessor::muonEfficiencyScaleFactor(unsigned _iM, double& _effScale, double& _relErr)
{
  double idsf(1.);
  double idsfErr(0.);
  getValueFromTable(idsf, idsfErr, *idsfTable[2], muons.pt[_iM], muons.eta[_iM]);

  double hltsf(1.);
  double hltsfErr(0.);
  getValueFromTable(hltsf, hltsfErr, *hltsfTable[3], muons.eta[_iM]);
  if(dataset.dataType == Dataset::kFastSim) hltsf *= 0.96957; // Mu22

  double sf(idsf * hltsf);

  _effScale *= sf;
  _relErr += std::sqrt(std::pow(idsfErr / idsf, 2.) + std::pow(hltsfErr / hltsf, 2.));
}

void
GLSkimProcessor::photonInefficiencyScaleFactor(unsigned _iP, unsigned _channel, double& _effScale, double& _relErr)
{
  double eta(TVector3(photons.caloX[_iP], photons.caloY[_iP], photons.caloZ[_iP]).Eta());

  double idsf(1.);
  double idsfErr(0.);
  double ideff(1.);
  double ideffErr(0.);
  getValueFromTable(idsf, idsfErr, *idsfTable[0], photons.pt[_iP], eta);
  getValueFromTable(ideff, ideffErr, *ideffTable[0], photons.pt[_iP], eta);

  double hltsf(1.);
  double hltsfErr(0.);
  double hlteff(1.);
  double hlteffErr(0.);
  if(_channel == 0){
    getValueFromTable(hltsf, hltsfErr, *hltsfTable[0], eta, nVtx);
    getValueFromTable(hlteff, hlteffErr, *hlteffTable[0], eta, nVtx);
  //   getValueFromTable(hltsf, hltsfErr, *hltsfTable[0], photons.pt[_iP], eta);
  //   getValueFromTable(hlteff, hlteffErr, *hlteffTable[0], photons.pt[_iP], eta);

    if(dataset.dataType == Dataset::kFastSim){
      hltsf *= 0.96405;
      hlteff *= 0.96405;
    }
  }
  else{
    getValueFromTable(hltsf, hltsfErr, *hltsfTable[1], 0.);
    getValueFromTable(hlteff, hlteffErr, *hlteffTable[1], 0.);

    if(dataset.dataType == Dataset::kFastSim){
      hltsf *= 1.01549; // Mu3p5EG12
      hlteff *= 1.01549; // Mu3p5EG12
    }
  }

  double sf(idsf * hltsf);
  double eff(ideff * hlteff);
  double sfErr(sf * std::sqrt(std::pow(idsfErr / idsf, 2.) + std::pow(hltsfErr / hltsf, 2.)));
  double effErr(eff * std::sqrt(std::pow(ideffErr / ideff, 2.) + std::pow(hlteffErr / hlteff, 2.)));

  if(eff >= 1.)
    throw std::runtime_error("MC efficiency >= 1");

  _effScale *= (1. - sf * eff) / (1. - eff);
  _relErr -= 1. / (1. - sf * eff) * std::sqrt(std::pow(eff * sfErr, 2.) + std::pow((1 - sf) / (1. - eff) * effErr, 2.));
}

void
GLSkimProcessor::electronInefficiencyScaleFactor(unsigned _iE, double& _effScale, double& _relErr)
{
  double eta(TVector3(electrons.caloX[_iE], electrons.caloY[_iE], electrons.caloZ[_iE]).Eta());

  double idsf(1.);
  double idsfErr(0.);
  double ideff(1.);
  double ideffErr(0.);
  getValueFromTable(idsf, idsfErr, *idsfTable[1], electrons.pt[_iE], eta);
  getValueFromTable(ideff, ideffErr, *ideffTable[1], electrons.pt[_iE], eta);

  double hltsf(1.);
  double hltsfErr(0.);
  double hlteff(1.);
  double hlteffErr(0.);
  getValueFromTable(hltsf, hltsfErr, *hltsfTable[2], eta, nVtx);
  getValueFromTable(hlteff, hlteffErr, *hlteffTable[2], eta, nVtx);
  // getValueFromTable(hltsf, hltsfErr, *hltsfTable[2], electrons.pt[_iE], eta);
  // getValueFromTable(hlteff, hlteffErr, *hlteffTable[2], electrons.pt[_iE], eta);

  if(dataset.dataType == Dataset::kFastSim){
    hltsf *= 0.96147; // Photon22IdIso
    hlteff *= 0.96147; // Photon22IdIso
  }

  double sf(idsf * hltsf);
  double eff(ideff * hlteff);
  double sfErr(sf * std::sqrt(std::pow(idsfErr / idsf, 2.) + std::pow(hltsfErr / hltsf, 2.)));
  double effErr(eff * std::sqrt(std::pow(ideffErr / ideff, 2.) + std::pow(hlteffErr / hlteff, 2.)));

  if(eff >= 1.)
    throw std::runtime_error("MC efficiency >= 1");

  _effScale *= (1. - sf * eff) / (1. - eff);
  _relErr -= 1. / (1. - sf * eff) * std::sqrt(std::pow(eff * sfErr, 2.) + std::pow((1 - sf) / (1. - eff) * effErr, 2.));
}

void
GLSkimProcessor::muonInefficiencyScaleFactor(unsigned _iM, double& _effScale, double& _relErr)
{
  double idsf(1.);
  double idsfErr(0.);
  double ideff(1.);
  double ideffErr(0.);
  getValueFromTable(idsf, idsfErr, *idsfTable[2], muons.pt[_iM], muons.eta[_iM]);
  getValueFromTable(ideff, ideffErr, *ideffTable[2], muons.pt[_iM], muons.eta[_iM]);

  double hltsf(1.);
  double hltsfErr(0.);
  double hlteff(1.);
  double hlteffErr(0.);
  getValueFromTable(hltsf, hltsfErr, *hltsfTable[3], muons.eta[_iM]);
  getValueFromTable(hlteff, hlteffErr, *hlteffTable[3], muons.eta[_iM]);

  if(dataset.dataType == Dataset::kFastSim){
    hltsf *= 0.96957; // Mu22
    hlteff *= 0.96957; // Mu22
  }

  double sf(idsf * hltsf);
  double eff(ideff * hlteff);
  double sfErr(sf * std::sqrt(std::pow(idsfErr / idsf, 2.) + std::pow(hltsfErr / hltsf, 2.)));
  double effErr(eff * std::sqrt(std::pow(ideffErr / ideff, 2.) + std::pow(hlteffErr / hlteff, 2.)));

  if(eff >= 1.)
    throw std::runtime_error("MC efficiency >= 1");

  _effScale *= (1. - sf * eff) / (1. - eff);
  _relErr -= 1. / (1. - sf * eff) * std::sqrt(std::pow(eff * sfErr, 2.) + std::pow((1 - sf) / (1. - eff) * effErr, 2.));
}

void
GLSkimProcessor::getValueFromTable(double& _val, double& _err, TH1 const& _table, double _x, double _y/* = 0.*/)
{
  int xbin(0);
  if(_table.GetXaxis()->GetXmin() == 0.)
    xbin = _table.GetXaxis()->FindFixBin(std::abs(_x));
  else
    xbin = _table.GetXaxis()->FindFixBin(_x);

  if(xbin == 0 || xbin > _table.GetNbinsX()) throw std::runtime_error(TString::Format("Invalid x value %.2f while looking up %s", _x, _table.GetName()).Data());

  int bin(xbin);

  if(_table.GetDimension() == 2){
    int ybin(0);
    if(_table.GetYaxis()->GetXmin() == 0.)
      ybin = _table.GetYaxis()->FindFixBin(std::abs(_y));
    else
      ybin = _table.GetYaxis()->FindFixBin(_y);

    if(ybin == 0 || ybin > _table.GetNbinsY()) throw std::runtime_error(TString::Format("Invalid y value %.3f while looking up %s", _y, _table.GetName()).Data());

    bin = _table.GetBin(xbin, ybin);
  }

  _val = _table.GetBinContent(bin);
  _err = _table.GetBinError(bin);
}

void
GLSkimProcessor::setJets(unsigned _nObj/* = 0*/, double* _eta/* = 0*/, double* _phi/* = 0*/)
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

class MCJetPhotonWeight : public GLEventWeight {
public:
  MCJetPhotonWeight(int _lepton) :
    GLEventWeight("MCJetPhotonHLTIsoWeight"),
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

class JetPhotonHLTIsoWeight : public GLEventWeight {
public:
  JetPhotonHLTIsoWeight() : GLEventWeight("JetPhotonHLTIsoWeight") {}
  ~JetPhotonHLTIsoWeight() {}
  void setPhoton(susy::PhotonVars const& _photon, susy::VertexVarsArray const&, susy::SimpleEventProducer::EventVars const&)
  {
    weight = (0.19 + 6.8 * TMath::Exp(-_photon.pt * 0.0526)) / (1. + 30. * TMath::Exp(-_photon.pt * 0.0625));
    relErr = (0.010275140445 * _photon.pt + -0.212842251317) * 1.25;
  }
};

class JetPhotonWeight : public GLEventWeight {
public:
  JetPhotonWeight() : GLEventWeight("JetPhotonWeight") {}
  ~JetPhotonWeight() {}
  void setPhoton(susy::PhotonVars const& _photon, susy::VertexVarsArray const&, susy::SimpleEventProducer::EventVars const&)
  {
    double exp(TMath::Exp(-_photon.pt * 0.07143));
    weight = (0.39 + 27 * exp) / (1. + 42. * exp);
    relErr = (0.00968391095319 * _photon.pt + -0.326302767227) * 1.33;
  }
};

class MCElePhotonFunctionalWeight : public GLEventWeight {
public:
  MCElePhotonFunctionalWeight(bool _muCorr = false) : GLEventWeight("MCElePhotonFunctionalWeight"), muCorr_(_muCorr) {}
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

class ElePhotonFunctionalWeight : public GLEventWeight {
public:
  ElePhotonFunctionalWeight(bool _muCorr = false) : GLEventWeight("ElePhotonFunctionalWeight"), muCorr_(_muCorr) {}
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

class PhotonGenIsoWeight : public GLEventWeight {
public:
  PhotonGenIsoWeight() :
    GLEventWeight("PhotonGenIsoWeight"),
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
