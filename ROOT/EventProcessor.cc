#include "EventProcessor.h"

#include "TSystem.h"
#include "TPRegexp.h"

#include <cmath>

EventProcessor::EventProcessor(unsigned _nOutputs, char const* _datasetName, unsigned _dataType, double _Leff, double _sigmaRelErr2, char const* _outputDir) :
  datasetName(_datasetName),
  dataType(_dataType),
  Leff(_Leff),
  sigmaRelErr2(_sigmaRelErr2),
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
EventProcessor::setOutput(char const* _outputName, EventWeight* _calc)
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

    TFile::Open(outputDir + "/" + datasetName + "_" + fItr->first + ".root", "recreate");
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
  eventSigma = weightCalc[_iF]->weight / Leff;
  sigmaErr = eventSigma * std::sqrt(weightCalc[_iF]->relErr * weightCalc[_iF]->relErr + sigmaRelErr2);

  eventList[_iF]->Fill();
}
