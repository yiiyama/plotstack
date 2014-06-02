#ifndef EventProcessor_h
#define EventProcessor_h

#include "TString.h"
#include "TTree.h"
#include "TFile.h"

#include <vector>
#include <map>

class EventWeight {
public:
  double weight;
  double relErr;

  EventWeight() : weight(1.), relErr(0.) {}
  EventWeight(double _w, double _r) : weight(_w), relErr(_r) {}
  virtual ~EventWeight() {}
};

class EventProcessor {
public:
  EventProcessor(int, char const*, double, double, char const*);
  virtual ~EventProcessor() {} // ROOT closes files in case of crashes; deleting event list causes double free

  void addInput(char const*);
  void setFilter(char const*, EventWeight*);
  void book();
  void write();
  void fill(unsigned);

  virtual void process() {}

  /* DATASET */

  unsigned nFilters;
  TString datasetName;
  double Leff;
  double sigmaRelErr2;

  /* CONFIGURATION */

  TString outputDir;
  std::vector<bool> useFilter;
  std::map<TString, unsigned> filterIndices;

  /* INPUT / OUTPUT */

  TString inputDir;
  std::vector<TString> inputPaths;

  /* OUTPUT */

  std::vector<TTree*> eventList;
  double eventSigma; // effective cross section per event
  double sigmaErr;

  /* UTIL */

  std::vector<EventWeight*> weightCalc;
};

#endif
