#ifndef EventProcessor_h
#define EventProcessor_h

#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TEntryList.h"

#include "Dataset.h"

#include <vector>
#include <map>

class EventWeight {
public:
  TString name;
  double weight;
  double relErr;

  EventWeight() : name(""), weight(1.), relErr(0.) {}
  EventWeight(TString const& _name, double _w, double _r) : name(_name), weight(_w), relErr(_r) {}
  virtual ~EventWeight() {}
};

class EventProcessor {
public:
  EventProcessor(unsigned, Dataset const*, char const*);
  virtual ~EventProcessor();

  virtual void addInput(char const*, char const* = 0);
  virtual void setOutput(char const*, EventWeight* = 0);
  virtual void book();
  virtual void write();
  virtual void fill(unsigned);

  virtual void process() {}

  /* CONFIGURATION */

  Dataset dataset;
  TString outputDir;
  std::vector<bool> produceOutput;
  std::map<TString, unsigned> outputIndices;

  /* INPUT / OUTPUT */

  TString inputDir;
  std::vector<TString> inputPaths;
  TEntryList* entrylist;

  /* OUTPUT */

  std::vector<TTree*> eventList;
  double eventSigma; // effective cross section per event
  double sigmaErr;
  double effScale;
  double scaleErr;

  /* UTIL */

  std::vector<EventWeight*> weightCalc;
};

#endif
