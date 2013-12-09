#include "TauDataFormat/TauNtuple/interface/DataMCType.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "TauDataFormat/TauNtuple/interface/TauDecay.h"
#include <iostream>
#include <cstdlib>

TString DataMCType::type="unknown";


DataMCType::DataMCType(){
}

DataMCType::~DataMCType(){
}

unsigned int DataMCType::GetType(TString name){
  name.ToLower();
  StoreType(name);
  if(name=="data")      return Data;
  if(name=="h_tautau")  return H_tautau;
  if(name=="hpm_taunu") return Hpm_taunu;
  if(name=="ttbar")     return ttbar;
  if(name=="w_lnu")     return W_lnu;
  if(name=="w_enu")     return W_enu;
  if(name=="w_munu")    return W_munu;
  if(name=="w_taunu")   return W_taunu;
  if(name=="dy_ll")     return DY_ll;
  if(name=="dy_ee")     return DY_ee;
  if(name=="dy_mumu")   return DY_mumu;
  if(name=="dy_tautau") return DY_tautau;
  if(name=="dy_emu_embedded") return DY_emu_embedded;
  if(name=="zz")        return ZZ;
  if(name=="ww")        return WW;
  if(name=="wz")        return WZ;
  if(name=="qcd")       return QCD;
  if(name=="dy_emu")    return DY_emu;
  if(name=="dy_mue")    return DY_emu;
  if(name=="dy_mutau")  return DY_mutau;
  if(name=="dy_etau")   return DY_etau;
  std::cout << "ERROR: Data/MC Type UNKNOWN!!!! " << std::endl;
  return unknown;
}

unsigned int DataMCType::SignalCode(unsigned int type,unsigned int JAK_ID1, unsigned int nprong1,unsigned int JAK_ID2, unsigned int nprong2){
  if(type==Data)return type;
//   if(JAK_ID1==TauDecay::JAK_A1_3PI && nprong1==3 && (nprong2==1 || (JAK_ID2==TauDecay::JAK_A1_3PI && nprong2==3))) return type+(JAK_ID1+100*nprong1)*100+(JAK_ID2+nprong2*100)*1000*100;
//   if(JAK_ID2==TauDecay::JAK_A1_3PI && nprong2==3 && (nprong1==1 || (JAK_ID2==TauDecay::JAK_A1_3PI && nprong1==3))) return type+(JAK_ID2+100*nprong2)*100+(JAK_ID1+nprong1*100)*1000*100;
  if(JAK_ID1==TauDecay::JAK_PION && nprong1==1 && (nprong2==1 || (JAK_ID2==TauDecay::JAK_PION && nprong2==1))) return type+(JAK_ID1+100*nprong1)*100+(JAK_ID2+nprong2*100)*1000*100;
  if(JAK_ID2==TauDecay::JAK_PION && nprong2==1 && (nprong1==1 || (JAK_ID1==TauDecay::JAK_PION && nprong1==1))) return type+(JAK_ID2+100*nprong2)*100+(JAK_ID1+nprong1*100)*1000*100;

  if(JAK_ID1==TauDecay::JAK_RHO_PIPI0 && nprong1==1 && (nprong2==1 || (JAK_ID2==TauDecay::JAK_RHO_PIPI0 && nprong2==1))) return type+(JAK_ID1+100*nprong1)*100+(JAK_ID2+nprong2*100)*1000*100;
  if(JAK_ID2==TauDecay::JAK_RHO_PIPI0 && nprong2==1 && (nprong1==1 || (JAK_ID1==TauDecay::JAK_RHO_PIPI0 && nprong1==1))) return type+(JAK_ID2+100*nprong2)*100+(JAK_ID1+nprong1*100)*1000*100;

  if(JAK_ID1==TauDecay::JAK_A1_3PI && nprong1==3 && (nprong2==1 || (JAK_ID2==TauDecay::JAK_A1_3PI && nprong2==3))) return type+(JAK_ID1+100*nprong1)*100+(JAK_ID2+nprong2*100)*1000*100;
  if(JAK_ID2==TauDecay::JAK_A1_3PI && nprong2==3 && (nprong1==1 || (JAK_ID1==TauDecay::JAK_A1_3PI && nprong1==3))) return type+(JAK_ID2+100*nprong2)*100+(JAK_ID1+nprong1*100)*1000*100;
                                                                                                                  
  if(JAK_ID1==TauDecay::JAK_3PIPI0 && nprong1==3 && (JAK_ID2==TauDecay::JAK_MUON && nprong2==1)) return type+(JAK_ID1+100*nprong1)*100+(JAK_ID2+nprong2*100)*1000*100;
  if(JAK_ID2==TauDecay::JAK_3PIPI0 && nprong2==3 && (JAK_ID1==TauDecay::JAK_MUON && nprong1==1)) return type+(JAK_ID2+100*nprong2)*100+(JAK_ID1+nprong1*100)*1000*100;

  if(JAK_ID1==TauDecay::JAK_KPIK && nprong1==3 && (JAK_ID2==TauDecay::JAK_MUON && nprong2==1))   return type+(JAK_ID1+100*nprong1)*100+(JAK_ID2+nprong2*100)*1000*100;
  if(JAK_ID2==TauDecay::JAK_KPIK && nprong2==3 && (JAK_ID1==TauDecay::JAK_MUON && nprong1==1))   return type+(JAK_ID2+100*nprong2)*100+(JAK_ID1+nprong1*100)*1000*100;

  if(JAK_ID1==TauDecay::JAK_KPIPI && nprong1==3 && (JAK_ID2==TauDecay::JAK_MUON && nprong2==1)) return type+(JAK_ID1+100*nprong1)*100+(JAK_ID2+nprong2*100)*1000*100;
  if(JAK_ID2==TauDecay::JAK_KPIPI && nprong2==3 && (JAK_ID1==TauDecay::JAK_MUON && nprong1==1)) return type+(JAK_ID2+100*nprong2)*100+(JAK_ID1+nprong1*100)*1000*100;

  return type;
}

void DataMCType::DecodeSignal(unsigned int code,unsigned int &type,unsigned int &JAK_ID1, unsigned int &nprong1,unsigned int &JAK_ID2, unsigned int &nprong2){
  type=code%100;
  JAK_ID1=(code/100)%100;
  nprong1=(code/100*100)%10;
  JAK_ID2=(code/100*100*10)%100;
  nprong2=(code/100*100*10*100)%10;
}

bool DataMCType::isSignalParticle(int pdg_id){
  unsigned int pdgid=abs(pdg_id);
  if(pdgid==PDGInfo::Z0 || pdgid==PDGInfo::W_plus || pdgid==PDGInfo::Higgs0 || pdgid==PDGInfo::Higgs_plus){
    return true;
  }
  return false;
}

