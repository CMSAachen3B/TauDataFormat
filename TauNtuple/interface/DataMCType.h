#ifndef DataMCType_h
#define DataMCType_h

#include "TString.h"

class DataMCType{
 public:
  enum Type {Data=1, H_tautau=10, Hpm_taunu=15, ttbar=20, W_lnu=20, DY_ll=30, ZZ=50, WW=51, WZ=52, QCD=60, unknown=999};

  DataMCType();
  ~DataMCType();

  unsigned int GetType(TString name);

};
#endif
