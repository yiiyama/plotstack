#ifndef PlotMaker_h
#define PlotMaker_h

#include "TTree.h"

#include "histogram.h"

#include <vector>

class PlotMaker {
 public:
  std::vector<Histogram*> histograms;
  TTree* eventList;
  double Lnorm;

  PlotMaker(unsigned);
  virtual ~PlotMaker() {}
  void setHistogram(unsigned, Histogram*);
  virtual void run(double&, double&) {}
};

#endif
