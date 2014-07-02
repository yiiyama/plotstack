#include "PlotMaker.h"

PlotMaker::PlotMaker() :
  histograms(),
  eventList(0),
  counter(0),
  Lnorm(0.),
  eventWeight(1.),
  eventWeightErr(0.)
{
}

void
PlotMaker::countEvent()
{
  if(!counter) return;
  counter->Fill(0.5, eventWeight);
  counter->Fill(1.5, eventWeightErr);
  counter->Fill(2.5);
}

void
PlotMaker::fill(TString const& _name, double _x)
{
  std::map<TString, Histogram*>::iterator hItr(histograms.find(_name));
  if(hItr == histograms.end()) return;
  
  hItr->second->fill(_x, eventWeight, eventWeightErr);
}

void
PlotMaker::fill(TString const& _name, double _x, double _y)
{
  std::map<TString, Histogram*>::iterator hItr(histograms.find(_name));
  if(hItr == histograms.end()) return;
  
  hItr->second->fill(_x, _y, eventWeight, eventWeightErr);
}
