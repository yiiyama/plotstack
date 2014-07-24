#ifndef EventProcessor_h
#define EventProcessor_h

#include "TString.h"
#include "TTree.h"
#include "TFile.h"

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
  virtual ~EventProcessor() {} // ROOT closes files in case of crashes; deleting event list causes double free

  virtual void addInput(char const*);
  virtual void setOutput(char const*, EventWeight*);
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
