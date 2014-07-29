#include "EventProcessor.h"

#include "TSystem.h"
#include "TPRegexp.h"

#include <cmath>

EventProcessor::EventProcessor(unsigned _nOutputs, Dataset const* _dataset, char const* _outputDir) :
  dataset(*_dataset),
  outputDir(_outputDir),
  produceOutput(_nOutputs, false),
  outputIndices(),
  eventList(_nOutputs, 0),
  eventSigma(0.),
  sigmaErr(0.),
  effScale(1.),
  scaleErr(0.),
  weightCalc(_nOutputs, 0)
{
}

void
EventProcessor::addInput(char const* _inputPath)
{
  // inputPath can contain regular expression in the base name

  TString dirName(gSystem->DirName(_inputPath));
  void* dir(gSystem->OpenDirectory(dirName));
  if(!dir) return;
  TString entry;
  TPRegexp pat(TString("^") + gSystem->BaseName(_inputPath) + "$");
  while(true){
    entry = gSystem->GetDirEntry(dir);
    if(entry.Length() == 0) break;
    if(entry == "." || entry == "..") continue;
    if(!pat.MatchB(entry)) continue;

    inputPaths.push_back(dirName + "/" + entry);
  }
  gSystem->FreeDirectory(dir);
}

void
EventProcessor::setOutput(char const* _outputName, EventWeight* _calc/* = 0*/)
{
  std::map<TString, unsigned>::iterator idxItr(outputIndices.find(_outputName));
  if(idxItr == outputIndices.end()) return;
    
  produceOutput[idxItr->second] = true;
  weightCalc[idxItr->second] = _calc;
}

void
EventProcessor::book()
{
  for(std::map<TString, unsigned>::iterator fItr(outputIndices.begin()); fItr != outputIndices.end(); ++fItr){
    if(!produceOutput[fItr->second]) continue;

    TString outputName(outputDir + "/" + dataset.name + "_" + fItr->first);
    if(dataset.prescale != 1) outputName += TString::Format("_ps%d", dataset.prescale);
    outputName += ".root";
    TFile::Open(outputName, "recreate");
    if(weightCalc[fItr->second])
      TObjString(weightCalc[fItr->second]->name).Write();
    
    eventList[fItr->second] = new TTree("eventList", "Event List");
    eventList[fItr->second]->SetAutoSave(10000000);

    eventList[fItr->second]->Branch("eventSigma", &eventSigma, "eventSigma/D");
    eventList[fItr->second]->Branch("sigmaErr", &sigmaErr, "sigmaErr/D");
    eventList[fItr->second]->Branch("effScale", &effScale, "effScale/D");
    eventList[fItr->second]->Branch("scaleErr", &scaleErr, "scaleErr/D");
  }
}

void
EventProcessor::write()
{
  for(unsigned iF(0); iF != eventList.size(); ++iF){
    if(eventList[iF]){
      TFile* file(eventList[iF]->GetCurrentFile());
      file->cd();
      file->Write();
      delete file;
    }
    eventList[iF] = 0;
  }
}

void
EventProcessor::fill(unsigned _iF)
{
  if(weightCalc[_iF]){
    eventSigma = weightCalc[_iF]->weight / dataset.Leff;
    sigmaErr = eventSigma * std::sqrt(weightCalc[_iF]->relErr * weightCalc[_iF]->relErr + dataset.sigmaRelErr2);
  }
  else{
    eventSigma = 1. / dataset.Leff;
    sigmaErr = eventSigma * std::sqrt(dataset.sigmaRelErr2);
  }

  eventList[_iF]->Fill();
}
