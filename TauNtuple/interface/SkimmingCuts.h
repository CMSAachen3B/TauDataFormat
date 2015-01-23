// -*- C++ -*-
//
// Package:    SkimmingCuts
// Class:      SkimmingCuts
// 
/**\class SkimmingCuts SkimmingCuts.cc SkimmingTools/SkimmingCuts/src/SkimmingCuts.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Vladimir Cherepanov
//         Created:  Thu Feb 23 18:53:30 CET 2012
// $Id: SkimmingCuts.h,v 1.1 2013/07/02 09:35:44 inugent Exp $
//
//
#ifndef SkimmingCuts_h
#define SkimmingCuts_h


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/MuonReco/interface/MuonIsolation.h>
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include <DataFormats/Candidate/interface/Candidate.h>
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "TauDataFormat/TauNtuple/interface/TauNtuple.h"
#include <TMath.h>
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TVector3.h"

#include "TVectorT.h"
#include "TH1.h"
#include "TLorentzVector.h"


//
// class declaration
//

class SkimmingCuts : public edm::EDFilter {
public:
  explicit SkimmingCuts(const edm::ParameterSet&);
  ~SkimmingCuts();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob();
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  
  virtual bool beginRun(edm::Run&, edm::EventSetup const&){return true;}
  virtual bool endRun(edm::Run&, edm::EventSetup const&){return true;}
  virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&){return true;}
  virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&){return true;}
  
  // ----------member data ---------------------------
  bool MuonCuts(edm::Event& iEvent, const edm::EventSetup& iSetup);
  bool ElectronCuts(edm::Event& iEvent, const edm::EventSetup& iSetup);
  bool PFTauCuts(edm::Event& iEvent, const edm::EventSetup& iSetup);
  bool DoubleMu(edm::Event& iEvent, const edm::EventSetup& iSetup);
  bool DoubleEle(edm::Event& iEvent, const edm::EventSetup& iSetup);
  
  edm::Event * iEvent_;
  std::string preselection_;

  int cnt_;
  int cntFound_;

};

#endif
