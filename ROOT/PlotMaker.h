#ifndef PlotMaker_h
#define PlotMaker_h

#include "TTree.h"
#include "TH1.h"

#include "histogram.h"

#include <map>

class PlotMaker {
 public:
  std::map<TString, Histogram*> histograms;
  TTree* eventList;
  TH1* counter;
  double Lnorm;
  double eventWeight;
  double eventWeightErr;

  PlotMaker();
  virtual ~PlotMaker() {}
  void setCounter(TH1* _counter) { counter = _counter; }
  void setHistogram(char const* _name, Histogram* _h) { histograms[_name] = _h; }
  void countEvent();
  void fill(TString const&, double);
  void fill(TString const&, double, double);
  virtual void run() = 0;
};

#endif
