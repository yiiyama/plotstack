#include "EventProcessor.h"

#include "TSystem.h"
#include "TPRegexp.h"

#include <cmath>

EventProcessor::EventProcessor(int _nFilters, char const* _datasetName, double _Leff, double _sigmaRelErr2, char const* _outputDir) :
  nFilters(_nFilters),
  datasetName(_datasetName),
  Leff(_Leff),
  sigmaRelErr2(_sigmaRelErr2),
  outputDir(_outputDir),
  useFilter(nFilters, false),
  filterIndices(),
  eventList(nFilters, 0),
  weightCalc(nFilters, 0)
{}

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
EventProcessor::setFilter(char const* _filterName, EventWeight* _calc)
{
  std::map<TString, unsigned>::iterator idxItr(filterIndices.find(_filterName));
  if(idxItr == filterIndices.end()) return;
    
  useFilter[idxItr->second] = true;
  weightCalc[idxItr->second] = _calc;
}

void
EventProcessor::book()
{
  for(std::map<TString, unsigned>::iterator fItr(filterIndices.begin()); fItr != filterIndices.end(); ++fItr){
    if(!useFilter[fItr->second]) continue;

    TFile::Open(outputDir + "/" + datasetName + "_" + fItr->first + ".root", "recreate");
    
    eventList[fItr->second] = new TTree("eventList", "Event List");
    eventList[fItr->second]->SetAutoSave(10000000);

    eventList[fItr->second]->Branch("eventSigma", &eventSigma, "eventSigma/D");
    eventList[fItr->second]->Branch("sigmaErr", &sigmaErr, "sigmaErr/D");
  }
}

void
EventProcessor::write()
{
  for(unsigned iF(0); iF != nFilters; ++iF){
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
