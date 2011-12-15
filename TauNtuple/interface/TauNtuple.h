// -*- C++ -*-
//
// Package:    TauNtuple
// Class:      TauNtuple
// 
/**\class TauNtuple TauNtuple.cc TauDataFormat/TauNtuple/src/TauNtuple.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ian Nugent  and  Vladimir Cherepanov
//         Created:  Mon Nov 14 13:49:02 CET 2011
// $Id: TauNtuple.h,v 1.6 2011/12/12 13:32:41 cherepan Exp $
//
//
#ifndef TauNtuple_h
#define TauNtuple_h


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


//  ROOT 
#include "TTree.h"
#include "TFile.h"


#include <utility>

// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"

#include <TMath.h>
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TVector3.h"

#include "TVectorT.h"
#include "TH1.h"
#include "TLorentzVector.h"


#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include <DataFormats/Candidate/interface/Candidate.h>
#include <SimDataFormats/GeneratorProducts/interface/HepMCProduct.h>
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"


//  Jet  stuff
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include <DataFormats/Common/interface/Ptr.h>

// MET stuff
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"


#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/MuonReco/interface/MuonIsolation.h>


#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"//VertexCollection
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include <DataFormats/TrackReco/interface/TrackBase.h>
#include "TVector3.h"

#include "Math/SMatrix.h"
#include "Math/Vector4D.h"
//#include "LorentzVector.h"
#include "TMatrixT.h"
#include "DataFormats/Math/interface/Error.h"
#include <DataFormats/Candidate/interface/Candidate.h>
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
//  to be checked. Not all of the includes are required
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include <DataFormats/Candidate/interface/Candidate.h>
#include <SimDataFormats/GeneratorProducts/interface/HepMCProduct.h>
#include "DataFormats/KinematicFit/interface/SelectedKinematicDecay.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"


// PU
#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/Lumi3DReWeighting.h"
//




//  Trigger 
#include "DataFormats/Common/interface/TriggerResults.h"
#include <FWCore/Common/interface/TriggerNames.h>
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

//
//
// class declaration
//

class TauNtuple : public edm::EDProducer {
   public:
      explicit TauNtuple(const edm::ParameterSet&);
      ~TauNtuple();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
  void fillMCTruth(edm::Event& iEvent, const edm::EventSetup& iSetup);
  void fillPrimeVertex(edm::Event& iEvent, const edm::EventSetup& iSetup,edm::Handle< std::vector<reco::Track>  > &trackCollection);
  void fillMuons(edm::Event& iEvent, const edm::EventSetup& iSetup,edm::Handle< std::vector<reco::Track>  > &trackCollection);
  void fillTracks(edm::Handle< std::vector<reco::Track>  > &trackCollection);
  void fillPFTaus(edm::Event& iEvent, const edm::EventSetup& iSetup,edm::Handle< std::vector<reco::Track>  > &trackCollection);
  void fillKinFitTaus(edm::Event& iEvent, const edm::EventSetup& iSetup,edm::Handle< std::vector<reco::Track>  > &trackCollection);
  void fillPFJets(edm::Event& iEvent, const edm::EventSetup& iSetup,edm::Handle< std::vector<reco::Track>  > &trackCollection);
  void fillElectrons(edm::Event& iEvent, const edm::EventSetup& iSetup,edm::Handle< std::vector<reco::Track>  > &trackCollection);
  void fillMET(edm::Event& iEvent, const edm::EventSetup& iSetup);

  void fillEventInfo(edm::Event& iEvent, const edm::EventSetup& iSetup);
  std::vector<bool> CheckTauDiscriminators(std::vector<edm::Handle<reco::PFTauDiscriminator> > tauDiscriminators, reco::PFTauRef tauRef);
  reco::PFTauRef getMatchedHPSTau(edm::Handle<std::vector<reco::PFTau> > & HPStaus,   std::vector<float>  &UnmodifiedTau, unsigned int &match);
  reco::PFTauRef getHPSTauMatchedToJet(edm::Handle<std::vector<reco::PFTau> > & HPStaus,   std::vector<float>  &Jet, unsigned int &match);

  bool getTrackMatch(edm::Handle< std::vector<reco::Track>  > &trackCollection, reco::TrackRef &refTrack, int &match);
  bool getTrackMatch(edm::Handle< std::vector<reco::Track>  > &trackCollection, reco::TrackRefVector &refTracks, std::vector<int> &matches,std::vector<bool> &found); 	
  double DeltaPhi(double phi1, double phi2);
  void ClearEvent();


  edm::Event * iEvent_;
  edm::InputTag primVtxTag_;
  edm::InputTag muonsTag_;
  edm::InputTag hpsTauProducer_;
  edm::InputTag hpsPFTauDiscriminationByTightIsolation_;
  edm::InputTag hpsPFTauDiscriminationByMediumIsolation_;
  edm::InputTag hpsPFTauDiscriminationByLooseIsolation_;
  edm::InputTag pfMETTag_;
  edm::InputTag kinTausTag_;
  edm::InputTag KinFitAdvanced_;
  edm::InputTag pfjetsTag_;
  edm::InputTag generalTracks_;
  edm::InputTag gensrc_;
  edm::InputTag GenEventInfo_;
  std::vector<std::string> discriminators_;
  std::string DataMC_Type_;
  unsigned int DataMC_Type_idx;

  // PU
  std::string ScaleFactor_;
  std::string PUInputFile_;
  std::string PUInputHistoMC_;
  std::string PUInputHistoData_;


  edm::Lumi3DReWeighting LumiWeights_;
  //


  TFile *output;
  TTree *output_tree;

  int cnt_;

  //=======  Vertex ===
  std::vector<float> Vtx_chi2;
  std::vector<float> Vtx_nTrk;
  std::vector<float> Vtx_ndof;
  std::vector<float> Vtx_y;
  std::vector<float> Vtx_x;
  std::vector<float> Vtx_z;
  std::vector<std::vector<std::vector<float> > >  Vtx_Cov;
  std::vector<std::vector<int> >    Vtx_Track_idx;
  std::vector<float> Vtx_isFake;

  //=======  Muons ===
  std::vector<std::vector<float> > Muon_p4;
  std::vector<std::vector<double > > Muon_Poca;
  std::vector<bool> Muon_isGlobalMuon;
  std::vector<bool> Muon_isStandAloneMuon;
  std::vector<bool> Muon_isTrackerMuon;
  std::vector<bool> Muon_isCaloMuon;
  std::vector<bool> Muon_isIsolationValid;
  std::vector<bool> Muon_isQualityValid;
  std::vector<bool> Muon_isTimeValid;

  std::vector<float>	Muon_emEt03          ;
  std::vector<float>	Muon_emVetoEt03      ;
  std::vector<float>	Muon_hadEt03         ;
  std::vector<float>	Muon_hadVetoEt03     ;
  std::vector<float>	Muon_nJets03         ;
  std::vector<float>	Muon_nTracks03       ;
  std::vector<float>	Muon_sumPt03         ;
  std::vector<float>	Muon_trackerVetoPt03 ;
  
  std::vector<float>	Muon_emEt05          ;
  std::vector<float>	Muon_emVetoEt05      ;
  std::vector<float>	Muon_hadEt05         ;
  std::vector<float>	Muon_hadVetoEt05     ;
  std::vector<float>	Muon_nJets05         ;
  std::vector<float>	Muon_nTracks05       ;
  std::vector<float>	Muon_sumPt05         ;
  std::vector<float>	Muon_trackerVetoPt05 ;

  std::vector<int>      Muon_numberOfChambers;
  std::vector<int>      Muon_Charge;
  std::vector<unsigned int>  Muon_Track_idx;

  std::vector<float> Muon_hitPattern_pixelLayerwithMeas;
  std::vector<float> Muon_numberOfMatchedStations;
  std::vector<float> Muon_normChi2;
  std::vector<float> Muon_hitPattern_numberOfValidMuonHits;
  std::vector<float> Muon_innerTrack_numberofValidHits;
  std::vector<float> Muon_numberOfMatches;


  //======= PFTaus ===

  std::vector<std::vector<float> > PFTau_p4;
  //std::vector<std::vector<double > > PFTau_vertex;  // preserve for any vtx info
  std::vector<bool> PFTau_isTightIsolation;
  std::vector<bool> PFTau_isMediumIsolation;
  std::vector<bool> PFTau_isLooseIsolation;
  std::vector<int>  PFTau_hpsDecayMode;   
  std::vector<int>  PFTau_Charge;
  std::vector<std::vector<int> > PFTau_Track_idx;

  // to include: HPS discriminators against muon and electron
  //======= KinFitTaus ===

  std::vector<bool> KFTau_discriminatorByKFit;
  std::vector<bool> KFTau_discriminatorByQC;
  int  KFTau_nKinTaus;
  std::vector<std::vector<float> > KFTau_TauVis_p4;
  std::vector<std::vector<float> > KFTau_TauFit_p4;
  std::vector<std::vector<float> > KFTau_Neutrino_p4;
  std::vector<unsigned int> KFTau_MatchedHPS_idx;
  std::vector<std::vector<int> > KFTau_Track_idx;
  std::vector<int> KFTau_indexOfFitInfo;

  std::vector<float>	KFTau_Fit_chi2;
  std::vector<float>	KFTau_Fit_ndf;
  std::vector<int>	KFTau_Fit_ambiguity;
  std::vector<int>	KFTau_Fit_charge;
  std::vector<int>	KFTau_Fit_csum;     
  std::vector<int>      KFTau_Fit_iterations;
  std::vector<int> KFTau_Fit_IndexToPrimVertexVector;
  std::vector<std::vector<float> > KFTau_Fit_TauPrimVtx;
  //=======  Electrons ===




  //=======  PFJets ===
  std::vector<std::vector<float> > PFJet_p4;
  std::vector<float> PFJet_chargedEmEnergy;
  std::vector<float> PFJet_chargedHadronEnergy;
  std::vector<float> PFJet_chargedHadronMultiplicity;
  std::vector<float> PFJet_chargedMuEnergy;
  std::vector<float> PFJet_chargedMultiplicity;
  std::vector<float> PFJet_electronEnergy;
  std::vector<float> PFJet_electronMultiplicity;
  std::vector<float> PFJet_HFEMEnergy;
  std::vector<float> PFJet_HFEMMultiplicity;
  std::vector<float> PFJet_HFHadronEnergy;
  std::vector<float> PFJet_HFHadronMultiplicity;
  std::vector<float> PFJet_muonEnergy;
  std::vector<float> PFJet_muonMultiplicity;
  std::vector<float> PFJet_neutralEmEnergy;
  std::vector<float> PFJet_neutralHadronEnergy;
  std::vector<float> PFJet_neutralHadronMultiplicity;
  std::vector<float> PFJet_photonEnergy;
  std::vector<float> PFJet_photonMultiplicity;
  std::vector<float> PFJet_jetArea; 
  std::vector<float> PFJet_maxDistance;
  std::vector<int>   PFJet_nConstituents;
  std::vector<float> PFJet_pileup;  
  std::vector<float> PFJet_etaetaMoment;
  std::vector<float> PFJet_etaphiMoment;
  std::vector<std::vector<int> > PFJet_Track_idx;
  std::vector<unsigned int> PFJet_MatchedHPS_idx;

  std::vector<int>   PFJet_numberOfDaughters;
  std::vector<float> PFJet_chargedEmEnergyFraction;
  std::vector<float> PFJet_HFHadronEnergyFraction;
  std::vector<float> PFJet_neutralHadronEnergyFraction;
  std::vector<float> PFJet_PFJet_neutralEmEnergyFraction;

  //=======  MET ===
  // now only PFMET
  double MET_et;
  double MET_phi;
  double MET_sumET;

  //=======  Event ===

  int EventNumber;

  unsigned int  Event_EventNumber;	 
  unsigned int  Event_RunNumber;	 
  int  Event_bunchCrossing;	 
  int  Event_orbitNumber;	 
  unsigned int  Event_luminosityBlock;	 
  bool  Event_isRealData;        
  int PileupInfo_NumInteractions_nm1;
  int PileupInfo_NumInteractions_n0;
  int PileupInfo_NumInteractions_np1;
  double EvtWeight3D;


  //====== Tracks ======= 
  std::vector<std::vector<float> > Track_p4;
  std::vector<std::vector<double > > Track_Poca;
  std::vector<int> Track_charge;
  std::vector<float> Track_chi2;
  std::vector<float> Track_ndof;
  std::vector<unsigned short> Track_numberOfLostHits;
  std::vector<unsigned short> Track_numberOfValidHits;
  std::vector<unsigned int> Track_qualityMask;
  std::vector<std::vector<float> > Track_par;
  std::vector<std::vector<std::vector<float> > > Track_parCov;


  //====== Trigger ======= 

  //====== MCTruth =======
  // verbose variables
  std::vector<double> GenEventInfoProduct_weights;
  float GenEventInfoProduct_signalProcessID;
  float GenEventInfoProduct_weight;
  float GenEventInfoProduct_qScale;
  float GenEventInfoProduct_alphaQED;
  float GenEventInfoProduct_alphaQCD;

  //do All
  bool do_MCComplete_;
  std::vector<std::vector<float> > MC_p4;
  std::vector<int> MC_pdgid;
  std::vector<int> MC_charge;
  std::vector<unsigned int> MC_midx;

  // Do MC Signal Summary
  bool do_MCSummary_;

  // Signal particles Z, W, H0, Hpm
  std::vector<std::vector<float> > MCSignalParticle_p4;
  std::vector<int> MCSignalParticle_pdgid;
  std::vector<int> MCSignalParticle_charge;
  std::vector<std::vector<double > > MCSignalParticle_Poca;
  std::vector<std::vector<unsigned int> > MCSignalParticle_Tauidx;

  // MC Tau Info
  std::vector<std::vector<std::vector<float> > > MCTauandProd_p4;
  std::vector<std::vector<int> > MCTauandProd_pdgid;
  std::vector<std::vector<unsigned int> > MCTauandProd_midx;
  std::vector<std::vector<int> > MCTauandProd_charge;
  std::vector<unsigned int>      MCTau_JAK;
  std::vector<unsigned int>      MCTau_DecayBitMask;


  // Object Truth
  /*std::vector<std::vector<float> > MCJet_p4;
  double MET_et;
  double MET_phi;
  double MET_sumET;
  double MET_et;
  double MET_phi;
  double MET_sumET;
  */



  std::string processName_; // process name of (HLT) process for which to get HLT configuration
  HLTConfigProvider hltConfig_;/// The instance of the HLTConfigProvider as a data member


};
//define this as a plug-in
DEFINE_FWK_MODULE(TauNtuple);
#endif
