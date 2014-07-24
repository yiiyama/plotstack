#ifndef Dataset_h
#define Dataset_h

#include "TString.h"

struct Dataset {
  enum DataType {
    kRealData,
    kFullSim,
    kFastSim,
    nDataTypes
  };

  Dataset() {}
Dataset(char const* _name, unsigned _dataType, double _Leff, double _sigmaRelErr, unsigned _prescale) :
    name(_name),
    dataType(DataType(_dataType)),
    Leff(_Leff),
    sigmaRelErr2(_sigmaRelErr * _sigmaRelErr),
    prescale(_prescale)
  {}
  Dataset(Dataset const& _orig):
    name(_orig.name),
    dataType(_orig.dataType),
    Leff(_orig.Leff),
    sigmaRelErr2(_orig.sigmaRelErr2),
    prescale(_orig.prescale)
  {}  

  TString name;
  DataType dataType;
  double Leff;
  double sigmaRelErr2;
  unsigned prescale;
};

#endif
