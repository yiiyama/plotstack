#include "TH1.h"
#include "TH2.h"

#include <stdexcept>

class Histogram {
public:
  TH1* hWeighted;
  TH1* hScaleUp;
  TH1* hScaleDown;
  TH1* hRaw;
  bool overflowable;

  Histogram() :
    hWeighted(0), hScaleUp(0), hScaleDown(0), hRaw(0)
  {
  }

  Histogram(TH1* _hWeighted, TH1* _hScaleUp, TH1* _hScaleDown, TH1* _hRaw) :
    hWeighted(_hWeighted), hScaleUp(_hScaleUp), hScaleDown(_hScaleDown), hRaw(_hRaw)
  {
    if(hWeighted->GetSumw2N() == 0) hWeighted->Sumw2();
    if(hScaleUp->GetSumw2N() == 0) hScaleUp->Sumw2();
    if(hScaleDown->GetSumw2N() == 0) hScaleDown->Sumw2();
  }

  ~Histogram()
  {
  }

  void
  fill(double _x, double _weight, double _weightErr)
  {
    if(hWeighted->GetDimension() != 1)
      throw std::runtime_error("1D fill called for 2D histogram");

    double width(1.);

    int maxBin(overflowable ? hWeighted->GetNbinsX() - 1 : hWeighted->GetNbinsX());

    if(overflowable && _x >= hWeighted->GetXaxis()->GetBinUpEdge(maxBin))
      _x = hWeighted->GetXaxis()->GetBinCenter(maxBin + 1);
    else{
      int bin(hWeighted->FindFixBin(_x));
      if(bin > 0 && bin <= maxBin)
        width = hWeighted->GetXaxis()->GetBinWidth(bin);
    }

    hWeighted->Fill(_x, _weight / width);
    hScaleUp->Fill(_x, (_weight + _weightErr) / width);
    hScaleDown->Fill(_x, (_weight - _weightErr) / width);
    hRaw->Fill(_x);
  }

  void
  fill(double _x, double _y, double _weight, double _weightErr)
  {
    if(hWeighted->GetDimension() != 2)
      throw std::runtime_error("2D fill called for 1D histogram");

    double area(1.);

    int maxBinX(overflowable ? hWeighted->GetNbinsX() - 1 : hWeighted->GetNbinsX());
    int maxBinY(overflowable ? hWeighted->GetNbinsY() - 1 : hWeighted->GetNbinsY());

    if(overflowable && _x >= hWeighted->GetXaxis()->GetBinUpEdge(maxBinX))
      _x = hWeighted->GetXaxis()->GetBinCenter(maxBinX + 1);
    if(overflowable && _y >= hWeighted->GetYaxis()->GetBinUpEdge(maxBinY))
      _y = hWeighted->GetYaxis()->GetBinCenter(maxBinY + 1);

    int xbin(hWeighted->GetXaxis()->FindFixBin(_x));
    int ybin(hWeighted->GetYaxis()->FindFixBin(_y));
    if(xbin > 0 && xbin <= maxBinX && ybin > 0 && ybin <= maxBinY)
      area = hWeighted->GetXaxis()->GetBinWidth(xbin) * hWeighted->GetYaxis()->GetBinWidth(ybin);

    static_cast<TH2*>(hWeighted)->Fill(_x, _y, _weight / area);
    static_cast<TH2*>(hScaleUp)->Fill(_x, _y, (_weight + _weightErr) / area);
    static_cast<TH2*>(hScaleDown)->Fill(_x, _y, (_weight - _weightErr) / area);
    static_cast<TH2*>(hRaw)->Fill(_x, _y);
  }
};
