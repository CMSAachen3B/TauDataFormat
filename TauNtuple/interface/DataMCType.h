#ifndef DataMCType_h
#define DataMCType_h

#include "TString.h"

class DataMCType{
 public:
  enum Type {Data=1,
	  	 H_tautau=10,
	     H_tautau_ggF=11,
	     H_tautau_VBF=12,
	     H_tautau_WHZHTTH=13,
	     Hpm_taunu=15,
	     HplusBWB=16,
	     W_lnu=20,
	     W_enu=21,
	     W_munu=22,
	     W_taunu=23, 
	     DY_ll=30, 
	     DY_ee=31,
	     DY_mumu=32,
	     DY_tautau=33,
	     DY_emu_embedded=34,
	     ZZ=50, 
	     WW=51, 
	     WZ=52,
	     ZZ_4l=53,
	     ZZ_2l2nu=54,
	     ZZ_2l2q=55,
	     WW_2l2nu=56,
	     WZ_2l2q=57,
	     WZ_3l1nu=58,
	     QCD=60, 
	     ttbar=70,
	     tw=71,
	     tbarw=72,
	     unknown=999,
	     DY_ll_Signal=10230530, //Extra Flags for Ntuple Analysis  
	     DY_tautau_Signal=10230533, 
	     Signal=998,
	     DY_ll_tautau=97,
	     DY_ll_mmee=96,
	     DY_emu=40,
	     DY_mutau=41,
	     DY_etau=42
  };

  DataMCType();
  ~DataMCType();

  unsigned int GetType(TString name);
  unsigned int SignalCode(unsigned int type,unsigned int JAK_ID1, unsigned int nprong1,unsigned int JAK_ID2, unsigned int nprong2);
  void DecodeSignal(unsigned int code,unsigned int &type,unsigned int &JAK_ID1, unsigned int &nprong1,unsigned int &JAK_ID2, unsigned int &nprong2);
  bool isSignalParticle(int pdg_id);
  void StoreType(TString t){type=t;}
  unsigned int GetType(){return GetType(type);}

 private:
  static TString type;

};
#endif
