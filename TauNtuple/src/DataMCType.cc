#include "TauDataFormat/TauNtuple/interface/DataMCType.h"


DataMCType::DataMCType(){
}

DataMCType::~DataMCType(){
}

unsigned int DataMCType::GetType(TString name){
  name.ToLower();
  if(name=="data")      return Data;
  if(name=="h_tautau")  return H_tautau;
  if(name=="hpm_taunu") return Hpm_taunu;
  if(name=="ttbar")     return ttbar;
  if(name=="w_taunu")   return W_lnu;
  if(name=="dy_ee")     return DY_ll;
  if(name=="ZZ")        return ZZ;
  if(name=="WW")        return WW;
  if(name=="WZ")        return WZ;
  if(name=="qcd")       return QCD;
  return unknown;
}
