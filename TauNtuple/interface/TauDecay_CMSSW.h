// -*- C++ -*-
//
// Package:    TauNtuple
// Class:      TauDecay_CMSSW
// 
/**\class TauDecay TauDecay_CMSSW.cc TauDataFormat/TauNtuple/src/TauDecay_CMSSW.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ian Nugent  
//         Created:  Fri Nov 18 13:49:02 CET 2011
// $Id: TauDecay_CMSSW.h,v 1.2 2011/12/15 15:59:20 inugent Exp $
//
//
#ifndef TauDecay_CMSSW_h
#define TauDecay_CMSSW_h

#include "Validation/EventGenerator/interface/TauDecay.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include <SimDataFormats/GeneratorProducts/interface/HepMCProduct.h>

//
// class declaration
//
class TauDecay_CMSSW : public TauDecay {
public:
  TauDecay_CMSSW();
  ~TauDecay_CMSSW();

  //Function to analyze the tau
  bool AnalyzeTau(const reco::GenParticle *Tau,unsigned int &JAK_ID,unsigned int &TauBitMask);
  // Functions to get results
  std::vector<const reco::GenParticle* > Get_TauDecayProducts(){return TauDecayProducts;}
  std::vector<unsigned int> Get_MotherIdx(){return MotherIdx;}
  void CheckForSignal(unsigned int &type,edm::Handle<reco::GenParticleCollection> &genParticles);

private:
  // recursive function to loop through tau decay products
  void Analyze(const reco::GenParticle *Particle,unsigned int midx);
  void AddPi0Info(const reco::GenParticle *Particle,unsigned int midx);
  //varibles
  std::vector<const reco::GenParticle*> TauDecayProducts;
  std::vector<unsigned int> MotherIdx;
  unsigned int JAK_ID, TauBitMask;

};
#endif
