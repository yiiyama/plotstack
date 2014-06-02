#include "PlotMaker.h"

PlotMaker::PlotMaker(unsigned _nP) :
  histograms(_nP, 0),
  eventList(),
  Lnorm(0.)
{
}

void
PlotMaker::setHistogram(unsigned _iH, Histogram* _h)
{
  if(_iH >= histograms.size()) return;
  histograms[_iH] = _h;
}
