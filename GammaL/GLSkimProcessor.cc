#include "../../CommonCode/ObjectTree.h"
#include "../../CommonCode/Utilities.h"
#include "../../CommonCode/SimpleEventProducer.h"
#include "../../CommonCode/ObjectSelector.h"

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH1D.h"
#include "TGraph.h"
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
  oPhotonAndDimuon,
  oElePhotonAndDimuon,
  oFakePhotonAndDimuon,
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
  "FakePhotonAndDimuon"
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
  GLSkimProcessor(char const*, unsigned, double, double, char const*);
  ~GLSkimProcessor() {} // ROOT closes files in case of crashes; deleting event list causes double free

  void process();
  void write();

private:
  bool preprocess();
  void processPhotonAndLepton(unsigned);
  void processElePhotonAndLepton(unsigned);
  void processFakePhotonAndLepton(unsigned);
  void processPhotonAndFakeLepton(unsigned);
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

  float mt;
  float metJESUp;
  float metJESDown;
  float mtJESUp;
  float mtJESDown;
  float mass2;
  float mass3;
  unsigned char nVtx;
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

  bool photon_muonIso[susy::NMAX];
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

GLSkimProcessor::GLSkimProcessor(char const* _datasetName, unsigned _dataType, double _Leff, double _sigmaRelErr2, char const* _outputDir) :
  EventProcessor(nOutputTypes, _datasetName, _dataType, _Leff, _sigmaRelErr2, _outputDir),
  useElectronFilter(false),
  useMuonFilter(false),
  eventVars(),
  electrons(),
  muons(),
  photons(),
  jets(),
  vertices(),
  mt(0.),
  mass2(0.),
  mass3(0.),
  nVtx(0),
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

  std::fill_n(photon_muonIso, susy::NMAX, true);
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
  for(unsigned iF(0); iF != nOutputTypes; ++iF)
    if(produceOutput[iF] && !weightCalc[iF])
      throw std::runtime_error("Event weight calculator not set");

  useElectronFilter = produceOutput[oPhotonAndElectron] || produceOutput[oElePhotonAndElectron] || produceOutput[oFakePhotonAndElectron] || produceOutput[oPhotonAndFakeElectron];
  useMuonFilter = produceOutput[oPhotonAndMuon] || produceOutput[oElePhotonAndMuon] || produceOutput[oFakePhotonAndMuon] || produceOutput[oPhotonAndFakeMuon] ||
    produceOutput[oPhotonAndDimuon] || produceOutput[oElePhotonAndDimuon] || produceOutput[oFakePhotonAndDimuon];

  ///////////////
  //// INPUT ////
  ///////////////

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
    filterTree.SetBranchStatus("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50", 1);
  }
  filterTree.SetBranchStatus("PhotonAndMuon", 1);
  filterTree.SetBranchStatus("ElePhotonAndMuon", 1);
  filterTree.SetBranchStatus("FakePhotonAndMuon", 1);
  filterTree.SetBranchStatus("PhotonAndFakeMuon", 1);
  filterTree.SetBranchStatus("HLT_Mu22_Photon22_CaloIdL", 1);

  bool PhotonAndElectron(false);
  bool ElePhotonAndElectron(false);
  bool FakePhotonAndElectron(false);
  bool PhotonAndFakeElectron(false);
  bool PhotonAndMuon(false);
  bool ElePhotonAndMuon(false);
  bool FakePhotonAndMuon(false);
  bool PhotonAndFakeMuon(false);
  bool elHLT(false);
  bool muHLT(false);

  if(useElectronFilter){
    filterTree.SetBranchAddress("PhotonAndElectron", &PhotonAndElectron);
    filterTree.SetBranchAddress("ElePhotonAndElectron", &ElePhotonAndElectron);
    filterTree.SetBranchAddress("FakePhotonAndElectron", &FakePhotonAndElectron);
    filterTree.SetBranchAddress("PhotonAndFakeElectron", &PhotonAndFakeElectron);
    filterTree.SetBranchAddress("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50", &elHLT);
  }
  filterTree.SetBranchAddress("PhotonAndMuon", &PhotonAndMuon);
  filterTree.SetBranchAddress("ElePhotonAndMuon", &ElePhotonAndMuon);
  filterTree.SetBranchAddress("FakePhotonAndMuon", &FakePhotonAndMuon);
  filterTree.SetBranchAddress("PhotonAndFakeMuon", &PhotonAndFakeMuon);
  filterTree.SetBranchAddress("HLT_Mu22_Photon22_CaloIdL", &muHLT);

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

  ////////////////
  //// OUTPUT ////
  ////////////////

  for(unsigned iF(0); iF != nOutputTypes; ++iF){
    if(!eventList[iF]) continue;

    eventList[iF]->Branch("run", &eventVars.runNumber, "run/i");
    eventList[iF]->Branch("lumi", &eventVars.lumiNumber, "lumi/i");
    eventList[iF]->Branch("event", &eventVars.eventNumber, "event/i");
    eventList[iF]->Branch("puWeight", &eventVars.puWeight, "puWeight/F");
    eventList[iF]->Branch("met", &eventVars.met, "met/F");
    eventList[iF]->Branch("metPhi", &eventVars.metPhi, "metPhi/F");
    eventList[iF]->Branch("rho", &eventVars.rho, "rho/F");
    eventList[iF]->Branch("mt", &mt, "mt/F");
    eventList[iF]->Branch("metJESUp", &metJESUp, "metJESUp/F");
    eventList[iF]->Branch("metJESDown", &metJESDown, "metJESDown/F");
    eventList[iF]->Branch("mtJESUp", &mtJESUp, "mtJESUp/F");
    eventList[iF]->Branch("mtJESDown", &mtJESDown, "mtJESDown/F");
    eventList[iF]->Branch("mass2", &mass2, "mass2/F");
    eventList[iF]->Branch("mass3", &mass3, "mass3/F");
    eventList[iF]->Branch("nVtx", &nVtx, "nVtx/b");
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

    TObjArray* fileNames(filterTree.GetListOfFiles());
    for(int iF(0); iF != fileNames->GetEntries(); ++iF){
      TString fileName(static_cast<TChainElement*>(fileNames->At(iF))->GetTitle());
      TFile* file(TFile::Open(fileName));
      if(produceOutput[oPhotonAndElectron]){
        if(!cutflowE_) cutflowE_ = static_cast<TH1D*>(file->Get("cutflowE"));
        else cutflowE_->Add(static_cast<TH1D*>(file->Get("cutflowE")));
      }
      if(produceOutput[oPhotonAndMuon]){
        if(!cutflowM_) cutflowM_ = static_cast<TH1D*>(file->Get("cutflowM"));
        else cutflowM_->Add(static_cast<TH1D*>(file->Get("cutflowM")));
      }
    }

    if(cutflowE_){
      cutflowE_->SetDirectory(eventList[oPhotonAndElectron]->GetCurrentFile());
      cutflowE_->SetBins(cutflowE_->GetNbinsX() + 1, cutflowE_->GetXaxis()->GetXmin(), cutflowE_->GetXaxis()->GetXmax() + 1.);
      cutflowE_->GetXaxis()->SetBinLabel(cutflowE_->GetNbinsX(), "FSRVeto");
    }
    if(cutflowM_){
      cutflowM_->SetDirectory(eventList[oPhotonAndMuon]->GetCurrentFile());
      cutflowM_->SetBins(cutflowM_->GetNbinsX() + 1, cutflowM_->GetXaxis()->GetXmin(), cutflowM_->GetXaxis()->GetXmax() + 1.);
      cutflowM_->GetXaxis()->SetBinLabel(cutflowM_->GetNbinsX(), "FSRVeto");
    }

  }

  ////////////////////////////////
  //// EFFICIENCY CORRECTIONS ////
  ////////////////////////////////

  if(produceOutput[oPhotonAndElectron] || produceOutput[oPhotonAndMuon] || produceOutput[oPhotonAndDimuon]){

    if(dataType != kRealData){
      TFile* idSFSource(TFile::Open("/afs/cern.ch/user/y/yiiyama/output/GammaL/main/idEfficiency/scalefactors.root"));
      if(!idSFSource)
        throw std::runtime_error("ID SF source not available");

      idsfTable[0] = static_cast<TH1*>(idSFSource->Get("photon"));
      idsfTable[0]->SetDirectory(0);
      if(useElectronFilter){
        idsfTable[1] = static_cast<TH1*>(idSFSource->Get("electron"));
        idsfTable[1]->SetDirectory(0);
      }
      if(useMuonFilter){
        idsfTable[2] = static_cast<TH1*>(idSFSource->Get("muon"));
        idsfTable[2]->SetDirectory(0);
      }

      delete idSFSource;

      TFile* hltSFSource(TFile::Open("/afs/cern.ch/user/y/yiiyama/output/GammaL/main/hltEfficiency/scalefactors.root"));
      if(!hltSFSource)
        throw std::runtime_error("HLT SF source not available");

      if(useElectronFilter){
        hltsfTable[0] = static_cast<TH1*>(hltSFSource->Get("photon_e"));
        hltsfTable[0]->SetDirectory(0);
        hltsfTable[2] = static_cast<TH1*>(hltSFSource->Get("electron"));
        hltsfTable[2]->SetDirectory(0);
      }
      if(useMuonFilter){
        TGraph* mueg(static_cast<TGraph*>(hltSFSource->Get("mueg")));
        TGraph* photon_mu(static_cast<TGraph*>(hltSFSource->Get("photon_mu")));
	hltsfTable[1] = new TH1D("photon_mu", "HLT SF", 1, 0., 1.);
	hltsfTable[1]->SetDirectory(0);
	hltsfTable[1]->SetBinContent(1, mueg->GetY()[0] * photon_mu->GetY()[0]);
	hltsfTable[1]->SetBinError(1, std::sqrt(std::pow(mueg->GetErrorY(0) * photon_mu->GetY()[0], 2.) + std::pow(mueg->GetY()[0] * photon_mu->GetErrorY(0), 2.)));
        delete mueg;
        delete photon_mu;

        hltsfTable[3] = static_cast<TH1*>(hltSFSource->Get("muon"));
	hltsfTable[3]->SetDirectory(0);
      }

      delete hltSFSource;

      TString effSourceName;
      if(datasetName.Contains("WG") || datasetName.Contains("WW")) effSourceName = "WGToLNuG";
      else if(datasetName.Contains("ZG") || datasetName.Contains("DY")) effSourceName = "ZGToLLG";
      else if(datasetName.Contains("TT") || datasetName.Contains("T5wg")) effSourceName = "TTGJets";
      else throw std::logic_error("Efficiency table missing");

      std::cout << "Using efficiency source " << effSourceName << std::endl;

      TFile* idEffSource(TFile::Open("/afs/cern.ch/user/y/yiiyama/output/GammaL/main/idEfficiency/" + effSourceName + ".root"));
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

      TFile* hltEffSource(TFile::Open("/afs/cern.ch/user/y/yiiyama/output/GammaL/main/hltEfficiency/" + effSourceName + ".root"));
      if(!hltEffSource)
        throw std::runtime_error("HLT Eff source not available");

      if(useElectronFilter){
        hlteffTable[0] = static_cast<TH1*>(hltEffSource->Get("photon_e_eff"));
        hlteffTable[0]->SetDirectory(0);
        hlteffTable[2] = static_cast<TH1*>(hltEffSource->Get("electron_eff"));
        hlteffTable[2]->SetDirectory(0);
      }
      if(useMuonFilter){
        TGraph* photon_mu(static_cast<TGraph*>(hltEffSource->Get("photon_mu_eff")));
	hlteffTable[1] = new TH1D("photon_mu_eff", "HLT EFficiency", 1, 0., 1.);
	hlteffTable[1]->SetDirectory(0);
	hlteffTable[1]->SetBinContent(1, photon_mu->GetY()[0]);
	hlteffTable[1]->SetBinError(1, photon_mu->GetErrorY(0));
	delete photon_mu;

	hlteffTable[3] = static_cast<TH1*>(hltEffSource->Get("muon_eff"));
	hlteffTable[3]->SetDirectory(0);
      }

      delete hltEffSource;
    }
  }

  ////////////////////
  //// START LOOP ////
  ////////////////////

  std::cout << "Processing " << datasetName << std::endl;

  long iEntry(0);
  while(filterTree.GetEntry(iEntry++) > 0){
    try{
      if(iEntry % 1000 == 1) (std::cout << "\r" << iEntry).flush();

      if(!(elHLT && ((PhotonAndElectron && produceOutput[oPhotonAndElectron]) ||
                     (ElePhotonAndElectron && produceOutput[oElePhotonAndElectron]) ||
                     (FakePhotonAndElectron && produceOutput[oFakePhotonAndElectron]) ||
                     (PhotonAndFakeElectron && produceOutput[oPhotonAndFakeElectron]))) &&
         !(muHLT && ((PhotonAndMuon && (produceOutput[oPhotonAndMuon] || produceOutput[oPhotonAndDimuon])) ||
                     (ElePhotonAndMuon && produceOutput[oElePhotonAndMuon]) ||
                     (FakePhotonAndMuon && produceOutput[oFakePhotonAndMuon]) ||
                     (PhotonAndFakeMuon && produceOutput[oPhotonAndFakeMuon])))) continue;


      eventTree.GetEntry(iEntry - 1);
      objectTree.GetEntry(iEntry - 1);

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

      if(FakePhotonAndElectron && !FakePhotonAndMuon && produceOutput[oFakePhotonAndElectron])
        processFakePhotonAndLepton(oFakePhotonAndElectron);

      if(FakePhotonAndMuon && produceOutput[oFakePhotonAndMuon])
        processFakePhotonAndLepton(oFakePhotonAndMuon);

      if(FakePhotonAndMuon && produceOutput[oFakePhotonAndDimuon])
        processFakePhotonAndDimuon();

      if(PhotonAndFakeElectron && !PhotonAndMuon && produceOutput[oPhotonAndFakeElectron])
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
  }

  if(dataType != kRealData){
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

  if(_outputType == oPhotonAndElectron){
    if(susy::deltaR(photonsOut.eta[0], photonsOut.phi[0], electronsOut.eta[0], electronsOut.phi[0]) < 0.8) return;
    calculateKinematics(photonsOut.at(0), electronsOut);

    if(cutflowE_) cutflowE_->Fill(cutflowE_->GetNbinsX() - 1.);
  }
  else if(_outputType == oPhotonAndMuon){
    if(susy::deltaR(photonsOut.eta[0], photonsOut.phi[0], muonsOut.eta[0], muonsOut.phi[0]) < 0.8) return;
    calculateKinematics(photonsOut.at(0), muonsOut);

    if(cutflowM_) cutflowM_->Fill(cutflowM_->GetNbinsX() - 1.);
  }

  jetsOut.clear();
  for(unsigned iJ(0); iJ != jets.size; ++iJ)
    if(jet_isCand[iJ]) jetsOut.push_back(jets.at(iJ));

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
      if(susy::deltaR(elePhoton.eta, elePhoton.phi, electronsOut.eta[0], electronsOut.phi[0]) < 0.8) continue;
      calculateKinematics(elePhoton, electronsOut);
    }
    else if(_outputType == oElePhotonAndMuon){
      if(susy::deltaR(elePhoton.eta, elePhoton.phi, muonsOut.eta[0], muonsOut.phi[0]) < 0.8) continue;
      calculateKinematics(elePhoton, muonsOut);
    }

    photonsOut.clear();

    photonsOut_matchGen[photonsOut.size] = photon_matchGen[iEP];
    photonsOut_matchGenE[photonsOut.size] = photon_matchGenE[iEP];
    photonsOut_genIso[photonsOut.size] = photon_genIso[iEP];
    photonsOut.push_back(elePhoton);

    jetsOut.clear();
    for(unsigned iJ(0); iJ != jets.size; ++iJ)
      if(jet_isCand[iJ]) jetsOut.push_back(jets.at(iJ));

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
    if(_outputType == oFakePhotonAndElectron){
      if(!photon_matchEHLT[0][iFP] || !photon_matchEHLT[1][iFP]) continue;
    }
    else{
      if(!photon_matchMHLT[iFP]) continue;
    }
    if(photons.hOverE[iFP] > 0.05 || photons.sigmaIetaIeta[iFP] > 0.014) continue;
    if(_outputType == oFakePhotonAndMuon && (photons.chargedHadronIso[iFP] > 15. || photons.neutralHadronIso[iFP] > 3.5 || photons.photonIso[iFP] > 1.3)) continue;

    susy::PhotonVars fakePhoton(photons.at(iFP));

    if(_outputType == oFakePhotonAndElectron){
      if(susy::deltaR(fakePhoton.eta, fakePhoton.phi, electronsOut.eta[0], electronsOut.phi[0]) < 0.8) continue;
      calculateKinematics(fakePhoton, electronsOut);
    }
    else if(_outputType == oFakePhotonAndMuon){
      if(susy::deltaR(fakePhoton.eta, fakePhoton.phi, muonsOut.eta[0], muonsOut.phi[0]) < 0.8) continue;
      calculateKinematics(fakePhoton, muonsOut);
    }

    photonsOut.clear();

    photonsOut_matchGen[photonsOut.size] = photon_matchGen[iFP];
    photonsOut_matchGenE[photonsOut.size] = photon_matchGenE[iFP];
    photonsOut_genIso[photonsOut.size] = photon_genIso[iFP];
    photonsOut.push_back(fakePhoton);

    jetsOut.clear();
    for(unsigned iJ(0); iJ != jets.size; ++iJ)
      if(jet_isCand[iJ] && susy::deltaR(fakePhoton.eta, fakePhoton.phi, jets.eta[iJ], jets.phi[iJ]) > 0.5)
        jetsOut.push_back(jets.at(iJ));

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

      susy::ObjectSelector::isGoodElectron(fakeElectron, susy::ElMedium12, &elIdResults);
      if((elIdResults & elBaseline) != elBaseline) continue;
      if(elIdResults[susy::ElDeltaEta] && elIdResults[susy::ElDeltaPhi]) continue;

      if(susy::deltaR(photon, fakeElectron) < 0.8) continue;

      electronsOut.clear();

      electronsOut_matchGen[electronsOut.size] = electron_matchGen[iFL];
      electronsOut.push_back(fakeElectron);

      calculateKinematics(photon, electronsOut);

      lEta = fakeElectron.eta;
      lPhi = fakeElectron.phi;

      calc->setElectron(fakeElectron, vertices, eventVars);
    }
    else if(_outputType == oPhotonAndFakeMuon){
      if(muons.iSubdet[iFL] == -1) continue;
      if(muons.pt[iFL] < 25.) break;
      if(!muon_matchHLT[iFL]) continue;
      if(!muon_isFake[iFL]) continue;

      if(muons.combRelIso[iFL] < 0.15 || muons.combRelIso[iFL] > 0.6) continue;

      susy::MuonVars fakeMuon(muons.at(iFL));

      susy::ObjectSelector::isGoodMuon(fakeMuon, susy::MuTight12, &muIdResults);
      if((muIdResults & muBaseline) != muBaseline) continue;

      if(susy::deltaR(photon, fakeMuon) < 0.8) continue;

      muonsOut.clear();

      muonsOut_matchGen[muonsOut.size] = muon_matchGen[iFL];
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

    fill(_outputType);
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

  jetsOut.clear();
  for(unsigned iJ(0); iJ != jets.size; ++iJ)
    if(jet_isCand[iJ]) jetsOut.push_back(jets.at(iJ));

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

    jetsOut.clear();
    for(unsigned iJ(0); iJ != jets.size; ++iJ)
      if(jet_isCand[iJ]) jetsOut.push_back(jets.at(iJ));

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

    jetsOut.clear();
    for(unsigned iJ(0); iJ != jets.size; ++iJ)
      if(jet_isCand[iJ] && susy::deltaR(fakePhoton.eta, fakePhoton.phi, jets.eta[iJ], jets.phi[iJ]) > 0.5)
        jetsOut.push_back(jets.at(iJ));

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
  // if(_channel == 0)
  //   getValueFromTable(hltsf, hltsfErr, *hltsfTable[0], eta, nVtx);
  // else
  //   getValueFromTable(hltsf, hltsfErr, *hltsfTable[1], eta, nVtx);
  if(_channel == 0)
    getValueFromTable(hltsf, hltsfErr, *hltsfTable[0], photons.pt[_iP], eta);
  else
    getValueFromTable(hltsf, hltsfErr, *hltsfTable[1], 0.);

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
  //  getValueFromTable(hltsf, hltsfErr, *hltsfTable[2], eta, nVtx);
  getValueFromTable(hltsf, hltsfErr, *hltsfTable[2], electrons.pt[_iE], eta);

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
  getValueFromTable(hltsf, hltsfErr, *hltsfTable[3], muons.pt[_iM], muons.eta[_iM]);

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
  // if(_channel == 0){
  //   getValueFromTable(hltsf, hltsfErr, *hltsfTable[0], eta, nVtx);
  //   getValueFromTable(hlteff, hlteffErr, *hlteffTable[0], eta, nVtx);
  // }
  // else{
  //   getValueFromTable(hltsf, hltsfErr, *hltsfTable[1], eta, nVtx);
  //   getValueFromTable(hlteff, hlteffErr, *hlteffTable[1], eta, nVtx);
  // }
  if(_channel == 0){
    getValueFromTable(hltsf, hltsfErr, *hltsfTable[0], photons.pt[_iP], eta);
    getValueFromTable(hlteff, hlteffErr, *hlteffTable[0], photons.pt[_iP], eta);
  }
  else{
    getValueFromTable(hltsf, hltsfErr, *hltsfTable[1], 0.);
    getValueFromTable(hlteff, hlteffErr, *hlteffTable[1], 0.);
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
  // getValueFromTable(hltsf, hltsfErr, *hltsfTable[2], eta, nVtx);
  // getValueFromTable(hlteff, hlteffErr, *hlteffTable[2], eta, nVtx);
  getValueFromTable(hltsf, hltsfErr, *hltsfTable[2], electrons.pt[_iE], eta);
  getValueFromTable(hlteff, hlteffErr, *hlteffTable[2], electrons.pt[_iE], eta);

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
  getValueFromTable(hltsf, hltsfErr, *hltsfTable[3], muons.pt[_iM], muons.eta[_iM]);
  getValueFromTable(hlteff, hlteffErr, *hlteffTable[3], muons.pt[_iM], muons.eta[_iM]);

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

class MCJetPhotonHLTIsoWeight : public GLEventWeight {
public:
  MCJetPhotonHLTIsoWeight() : GLEventWeight("MCJetPhotonHLTIsoWeight")
  {
    TFile* source(TFile::Open("/afs/cern.ch/user/y/yiiyama/output/GammaL/jetGammaFake/jetGammaMC_el.root"));
    if(!source)
      throw std::runtime_error("MCJetPhotonHLTIsoWeight source not found");
    TH1D* hFake(static_cast<TH1D*>(source->Get("hFakePt")));
    TH1D* hProxy(static_cast<TH1D*>(source->Get("hProxyPt")));
    int low(hFake->FindFixBin(40.));
    double fakeErr;
    double fake(hFake->IntegralAndError(low, hFake->GetNbinsX(), fakeErr));
    double proxyErr;
    double proxy(hProxy->IntegralAndError(low, hProxy->GetNbinsX(), proxyErr));

    delete source;

    weight = fake / proxy;
    relErr = std::sqrt(std::pow(fakeErr / fake, 2.) + std::pow(proxyErr / proxy, 2.));
  }
  ~MCJetPhotonHLTIsoWeight() {}
};

class JetPhotonHLTIsoWeight : public GLEventWeight {
public:
  JetPhotonHLTIsoWeight() : GLEventWeight("JetPhotonHLTIsoWeight") {}
  ~JetPhotonHLTIsoWeight() {}
  void setPhoton(susy::PhotonVars const& _photon, susy::VertexVarsArray const&, susy::SimpleEventProducer::EventVars const&)
  {
    weight = (8.39 * std::exp(-_photon.pt / 21.1) + 0.0498 * std::exp(_photon.pt / 87.8)) / (1. + 30.1 * std::exp(-_photon.pt / 16.6));
    relErr = 0.25; // nonclosure + emul-match + 10% on data measurement
  }
};

class MCJetPhotonWeight : public GLEventWeight {
public:
  MCJetPhotonWeight() :
    GLEventWeight("MCJetPhotonWeight")
  {
    TFile* source(TFile::Open("/afs/cern.ch/user/y/yiiyama/output/GammaL/jetGammaFake/jetGammaMC_mu.root"));
    if(!source)
      throw std::runtime_error("MCJetPhotonWeight source not found");
    TH1D* hFake(static_cast<TH1D*>(source->Get("hFakePt")));
    TH1D* hProxy(static_cast<TH1D*>(source->Get("hProxyPt")));
    int low(hFake->FindFixBin(40.));
    double fakeErr;
    double fake(hFake->IntegralAndError(low, hFake->GetNbinsX(), fakeErr));
    double proxyErr;
    double proxy(hProxy->IntegralAndError(low, hProxy->GetNbinsX(), proxyErr));

    delete source;

    weight = fake / proxy;
    relErr = std::sqrt(std::pow(fakeErr / fake, 2.) + std::pow(proxyErr / proxy, 2.));
  }
  ~MCJetPhotonWeight() {}
};

class JetPhotonWeight : public GLEventWeight {
public:
  JetPhotonWeight() : GLEventWeight("JetPhotonWeight") {}
  ~JetPhotonWeight() {}
  void setPhoton(susy::PhotonVars const& _photon, susy::VertexVarsArray const&, susy::SimpleEventProducer::EventVars const&)
  {
    weight = (28.0 * std::exp(-_photon.pt / 14.0) + 0.461 * std::exp(-_photon.pt / 681.)) / (1. + 42.1 * std::exp(-_photon.pt / 13.7));
    relErr = 0.11; // nonclosure + 10% on data measurement
  }
};

class MCElePhotonFunctionalWeight : public GLEventWeight {
public:
  MCElePhotonFunctionalWeight(bool _muCorr = false) : GLEventWeight("MCElePhotonFunctionalWeight"), muCorr_(_muCorr) {}
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
    double ineff(0.00133962);
    if(muCorr_) ineff *= 1.;

    double fakerate(1. - (1. - ineff) * (1. - std::pow(1.053e+00 * _photon.pt + 1., -1.612e+00)) * (1. - 1.903e-01 * std::exp(-3.939e-01 * _vertices.nTracks[iPV])) * (1. - 1.311e-04 * nV));
    weight = fakerate / (1. - fakerate);
    if(muCorr_)
      relErr = 0.;
    else
      relErr = 0.047;
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
    unsigned nV(0);
    int iPV(-1);
    for(unsigned iV(0); iV != _vertices.size; ++iV){
      if(!_vertices.isGood[iV]) continue;
      if(iPV < 0) iPV = iV;
      ++nV;
    }
    double ineff(0.00167252);
    if(muCorr_) ineff *= 1.;

    double fakerate(1. - (1. - ineff) * (1. - std::pow(9.635e-02 * _photon.pt + 1., -4.150e+00)) * (1. - 1.552e-01 * std::exp(-3.073e-01 * _vertices.nTracks[iPV])) * (1. - 3.140e-04 * nV));
    weight = fakerate / (1. - fakerate);
    if(muCorr_)
      relErr = 0.;
    else
      relErr = 0.091;
  }

private:
  bool muCorr_;
};

class JetElectronBinnedWeight : public GLEventWeight {
public:
  JetElectronBinnedWeight() :
    GLEventWeight("JetElectronBinnedWeight")
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
