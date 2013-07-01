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
// $Id: SkimmingCuts.cc,v 1.8 2013/06/15 10:32:26 inugent Exp $
//
//


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
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
  bool MuonsCuts(edm::Event& iEvent, const edm::EventSetup& iSetup);
  bool KFitTausCuts(edm::Event& iEvent, const edm::EventSetup& iSetup);
  bool ElectronCuts(edm::Event& iEvent, const edm::EventSetup& iSetup);
  bool PFTausCuts(edm::Event& iEvent, const edm::EventSetup& iSetup);

  edm::Event * iEvent_;
  edm::InputTag hpsTauProducer_;

  double MuonPtCut_;
  bool MuonIsGlo_;
  double NMuons_;
  double MuonEtaCut_;
  double PFTauPtCut_;
  double PFTauEtaCut_;

  double ElectronPtCut_;
  double ElectronEtaCut_;
  edm::InputTag KinFitAdvanced_;
  bool doMuonOnly_;

  int nMuon_;
  int nMuonPass_;
  int nPFTaus_;
  int nPFTausPass_;


  int cnt_;
  int cntFound_;

};

SkimmingCuts::SkimmingCuts(const edm::ParameterSet& iConfig):
  hpsTauProducer_( iConfig.getParameter<edm::InputTag>( "hpsTauProducer" ) ),
  MuonPtCut_( iConfig.getParameter<double>("MuonPtCut") ),
  MuonIsGlo_( iConfig.getParameter<bool>("MuonIsGlobal") ),
  NMuons_( iConfig.getParameter<double>("NMuons") ),
  MuonEtaCut_( iConfig.getParameter<double>("MuonEtaCut") ),
  PFTauPtCut_( iConfig.getParameter<double>("PFTauPtCut") ),
  PFTauEtaCut_( iConfig.getParameter<double>("PFTauEtaCut") ),
  ElectronPtCut_( iConfig.getParameter<double>("ElectronPtCut") ),
  ElectronEtaCut_( iConfig.getParameter<double>("ElectronEtaCut") ),
  doMuonOnly_( iConfig.getParameter<bool>("doMuonOnly") )
{
}


SkimmingCuts::~SkimmingCuts(){}

bool SkimmingCuts::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){
  cnt_++;
  bool pass = false;
  iEvent_=&iEvent;
  bool AcceptMuon = MuonsCuts(iEvent, iSetup);
  if(doMuonOnly_){
    if(AcceptMuon){
      pass = true;
      cntFound_++;
    }
    return pass;
  }
  //bool AcceptKFTau = KFitTausCuts(iEvent, iSetup);
  bool AcceptPFTau = PFTausCuts(iEvent, iSetup);
  //bool AcceptElectron = ElectronCuts(iEvent, iSetup);

    //----------- This Blos is for private production to analyse muon tau decay with KF
  if(AcceptMuon && (/*AcceptKFTau || AcceptElectron ||*/ AcceptPFTau)){
    pass = true;
    cntFound_++;
  }
  return pass;
}

bool SkimmingCuts::MuonsCuts(edm::Event& iEvent, const edm::EventSetup& iSetup){
  edm::Handle< reco::MuonCollection > muonCollection;
  iEvent_->getByLabel("muons",  muonCollection);
  
  for(unsigned int iMuon = 0; iMuon< muonCollection->size(); iMuon++){
    reco::MuonRef RefMuon(muonCollection, iMuon);
    if(RefMuon.isNonnull()){
      if(RefMuon->p4().Pt() > MuonPtCut_ && fabs(RefMuon->p4().Eta()) < MuonEtaCut_){
	if(RefMuon->isGlobalMuon() == MuonIsGlo_ && RefMuon->isPFMuon()){
	  return true;
	}
      }
    }
  }
  return false;
}

bool SkimmingCuts::ElectronCuts(edm::Event& iEvent, const edm::EventSetup& iSetup){
  edm::Handle< reco::GsfElectronCollection > electronCollection;
  iEvent_->getByLabel("gsfElectrons", electronCollection);
  
  for(unsigned int iElectron=0; iElectron<electronCollection->size(); iElectron++){
    reco::GsfElectronRef RefElectron(electronCollection, iElectron);
    if(RefElectron.isNonnull()){
      if(RefElectron->p4().Pt()>ElectronPtCut_ && fabs(RefElectron->p4().Eta())<ElectronEtaCut_){
	return true;
      }
    }
  }
  return false;
}



bool SkimmingCuts::PFTausCuts(edm::Event& iEvent, const edm::EventSetup& iSetup){
  edm::Handle<std::vector<reco::PFTau> > PFTaus;
  iEvent.getByLabel("hpsPFTauProducer", PFTaus);

  edm::Handle<reco::PFTauDiscriminator> HPSAgainstElectronsTight;
  iEvent.getByLabel("hpsPFTauDiscriminationByLooseElectronRejection", HPSAgainstElectronsTight);
  
  edm::Handle<reco::PFTauDiscriminator> HPSAgainstMuonTight;
  iEvent.getByLabel("hpsPFTauDiscriminationByLooseMuonRejection", HPSAgainstMuonTight);

  edm::Handle<reco::PFTauDiscriminator> HPSByDecayModeFinding;
  iEvent.getByLabel("hpsPFTauDiscriminationByDecayModeFinding", HPSByDecayModeFinding);

  //edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByLooseIsolationMVA2;
  //iEvent.getByLabel("hpsPFTauDiscriminationByLooseIsolationMVA2", HPSPFTauDiscriminationByLooseIsolationMVA2);

  //edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByMediumIsolationMVA2;
  //iEvent.getByLabel("hpsPFTauDiscriminationByMediumIsolationMVA2", HPSPFTauDiscriminationByMediumIsolationMVA2);
  
  edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByTightIsolationMVA2;
  iEvent.getByLabel("hpsPFTauDiscriminationByTightIsolationMVA2", HPSPFTauDiscriminationByTightIsolationMVA2);

  //================== KinematicFit Info ===================
  for ( unsigned int iPFTau = 0; iPFTau < PFTaus->size(); iPFTau++ ) {
    reco::PFTauRef PFTauCand(PFTaus, iPFTau);
    if(PFTauCand.isNonnull()){
      if(PFTauCand->p4().Pt() >PFTauPtCut_ && fabs(PFTauCand->p4().eta()) <PFTauEtaCut_){
	if((*HPSByDecayModeFinding)[PFTauCand] && PFTauCand->decayMode()==10 &&  (*HPSAgainstElectronsTight)[PFTauCand] && (*HPSAgainstMuonTight)[PFTauCand] && (*HPSPFTauDiscriminationByTightIsolationMVA2)[PFTauCand]){
	  return true;
	}
      }
    }
  }
  return false;
}


// ------------ method called once each job just before starting event loop  ------------
void 
SkimmingCuts::beginJob()
{
  nMuon_=0;
  nMuonPass_=0;
  nPFTaus_=0;
  nPFTausPass_=0;
  cnt_ =0;
  cntFound_ =0;


  //   std::cout<<"Starting preselection ...  " <<std::endl;
 }

// ------------ method called once each job just after ending the event loop  ------------
void 
SkimmingCuts::endJob() {

//   float ratioMuon = 0.0;
//   if(nMuon_!=0) ratioMuon=(float)nMuonPass_/nMuon_;

//   float ratioPFTau = 0.0;
//   if(nPFTaus_!=0) ratioPFTau=(float)nPFTausPass_/nPFTaus_;


  float ratio = 0.0;
  if(cnt_!=0) ratio=(float)cntFound_/cnt_;

  //  std::cout<<"NMuons: " << nMuon_ <<"   NMuonsPass: "<< nMuonPass_ <<"   Efficiency: "<< ratioMuon*100.0 <<"%"<<std::endl;
  // std::cout<<"NPFtaus: " << nPFTaus_ <<"   NPFtausPass: "<< nPFTausPass_ <<"   Efficiency: "<< ratioPFTau*100.0 <<"%"<<std::endl;


    std::cout<<"[SkimmingCuts]-->  "<<"NEvents: " << cnt_ <<"   NEventsPass: "<< cntFound_ <<"   Efficiency: "<< ratio*100.0 <<"%"<<std::endl;

}

// ------------ method called when starting to processes a run  ------------
bool 
SkimmingCuts::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
SkimmingCuts::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
SkimmingCuts::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
SkimmingCuts::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SkimmingCuts::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(SkimmingCuts);
