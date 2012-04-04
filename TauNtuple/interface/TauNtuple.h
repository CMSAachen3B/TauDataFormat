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
// $Id: TauNtuple.h,v 1.18 2012/04/03 00:05:39 inugent Exp $
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
#include <stdint.h>

//  ROOT 
#include "TTree.h"
#include "TFile.h"


#include <utility>

// system include files
#include <memory>
#include <iostream>
#include <stdint.h>
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
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include <FWCore/Common/interface/TriggerNames.h>
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTanalyzers/interface/JetUtil.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

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
  void fillTriggerInfo(edm::Event& iEvent, const edm::EventSetup& iSetup);
  template <class T>
  void TriggerMatch(edm::Handle<trigger::TriggerEvent> &triggerEvent,unsigned int triggerIndex,T obj,double drmax,std::vector<float> &match);
  void fillEventInfo(edm::Event& iEvent, const edm::EventSetup& iSetup);
  std::vector<bool> CheckTauDiscriminators(std::vector<edm::Handle<reco::PFTauDiscriminator> > tauDiscriminators, const reco::PFTauRef tauRef);
  reco::PFTauRef getMatchedHPSTau(edm::Handle<std::vector<reco::PFTau> > & HPStaus,   std::vector<float>  &UnmodifiedTau, unsigned int &match);
  reco::PFTauRef getHPSTauMatchedToJet(edm::Handle<std::vector<reco::PFTau> > & HPStaus,   std::vector<float>  &Jet, unsigned int &match);

  bool getTrackMatch(edm::Handle< std::vector<reco::Track>  > &trackCollection, reco::TrackRef &refTrack, int &match);
  double DeltaPhi(double phi1, double phi2);
  void ClearEvent();



  edm::InputTag primVtxTag_;
  edm::InputTag muonsTag_;
  edm::InputTag hpsTauProducer_;
  edm::InputTag hpsPFTauDiscriminationByTightIsolation_;
  edm::InputTag hpsPFTauDiscriminationByMediumIsolation_;
  edm::InputTag hpsPFTauDiscriminationByLooseIsolation_;


  edm::InputTag hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr_;
  edm::InputTag hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr_; 
  edm::InputTag hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr_;
  edm::InputTag hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr_; 


  edm::InputTag hpsPFTauDiscriminationAgainstElectronLoose_; 
  edm::InputTag hpsPFTauDiscriminationAgainstElectronMedium_;
  edm::InputTag hpsPFTauDiscriminationAgainstElectronTight_; 
  edm::InputTag hpsPFTauDiscriminationAgainstMuonLoose_;     
  edm::InputTag hpsPFTauDiscriminationAgainstMuonTight_;     
  edm::InputTag hpsPFTauDiscriminationByDecayModeFinding_;  




  edm::InputTag pfMETTag_;
  edm::InputTag kinTausTag_;
  edm::InputTag KinFitAdvanced_;
  edm::InputTag tauPrimaryVtx_;
  edm::InputTag pfjetsTag_;
  edm::InputTag PFElectronTag_;
  edm::InputTag generalTracks_;
  edm::InputTag gensrc_;
  edm::InputTag GenEventInfo_;
  std::vector<std::string> discriminators_;

  // PU
  std::string ScaleFactor_;
  std::string PUInputFile_;
  std::string PUInputHistoMC_;
  std::string PUInputHistoData_;
  std::string PUOutputFile_;

  edm::Lumi3DReWeighting LumiWeights_;

  // MC Signal
  bool do_MCSummary_;
  bool do_MCComplete_;

  // Trigger
  std::string processName_; 
  HLTConfigProvider hltConfig_; 
  L1GtUtils m_l1GtUtils_;
  edm::InputTag TriggerInfoName_;
  edm::InputTag TriggerEvent_;
  edm::InputTag TriggerResults_;
  edm::InputTag l1GtTriggerMenuLite_;
  float TriggerJetMatchingdr_;
  float TriggerMuonMatchingdr_;
  float TriggerElectronMatchingdr_;
  float TriggerTauMatchingdr_;
  float TriggerMETMatchingdr_;
  edm::Handle< L1GtTriggerMenuLite > triggerMenuLite_;


  // Control flags
  bool doBJets_;
  bool doPFJets_;
  bool doMuons_;
  bool doElectrons_;
  bool doPFTaus_;
  bool doTracks_;
  bool doKinFitTaus_;
  bool doTrigger_;
  bool doPrimeVertex_;
  bool doMET_;
  bool doMC_;
  bool doPatJets_;
  bool doPatElectrons_;
  bool doPatMuons_;
  bool doPatMET_;

  // Pat objects
  std::string srcPatJets_;
  std::string PatJetScale_;
  std::string BTagAlgorithim_;
  std::string srcPatMET_;

  // outputfile
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
  std::vector<std::vector<float > > Muon_Poca;
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
  std::vector<std::vector<float > > PFTau_Poca;
  std::vector<bool> PFTau_isTightIsolation;
  std::vector<bool> PFTau_isMediumIsolation;
  std::vector<bool> PFTau_isLooseIsolation;

  std::vector<bool> PFTau_isTightIsolationDBSumPtCorr;
  std::vector<bool> PFTau_isMediumIsolationDBSumPtCorr;
  std::vector<bool> PFTau_isLooseIsolationDBSumPtCorr;
  std::vector<bool> PFTau_isVLooseIsolationDBSumPtCorr;

  std::vector<bool> PFTau_isHPSAgainstElectronsLoose;
  std::vector<bool> PFTau_isHPSAgainstElectronsMedium;
  std::vector<bool> PFTau_isHPSAgainstElectronsTight;
  std::vector<bool> PFTau_isHPSAgainstMuonLoose;
  std::vector<bool> PFTau_isHPSAgainstMuonTight;
  std::vector<bool> PFTau_isHPSByDecayModeFinding;     



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
  std::vector<std::vector<float> > KFTau_Fit_TauPrimVtx;

  std::vector<float> KFTau_Fit_TauEnergyFraction;
  std::vector<float> KFTau_Fit_RefitVisibleMass;
  std::vector<float> KFTau_Fit_Chi2;
  std::vector<float> KFTau_Fit_PV_PV_significance;
  std::vector<float> KFTau_Fit_SV_PV_significance;

  std::vector<std::vector<int> > KFTau_Daughter_pdgid;
  std::vector<std::vector<int> > KFTau_Daughter_charge;
  std::vector<std::vector<float> > KFTau_Daughter_ambiguity;

  std::vector<std::vector<std::vector<float> > > KFTau_Daughter_par;
  std::vector<std::vector<std::vector<float> > > KFTau_Daughter_parCov;
  std::vector<std::vector<std::vector<float> > > KFTau_Daughter_inputpar;
  std::vector<std::vector<std::vector<float> > > KFTau_Daughter_inputparCov;

  std::vector<float> ReducedVtx_chi2;
  std::vector<float> ReducedVtx_nTrk;
  std::vector<float> ReducedVtx_ndof;
  std::vector<float> ReducedVtx_y;
  std::vector<float> ReducedVtx_x;
  std::vector<float> ReducedVtx_z;
  std::vector<std::vector<std::vector<float> > >  ReducedVtx_Cov;
  std::vector<std::vector<int> >    ReducedVtx_Track_idx;
  std::vector<float> ReducedVtx_isFake;

  //=======  Electrons ===
  std::vector<std::vector<float> > Electron_p4;
  std::vector<std::vector<float > > Electron_Poca;
  std::vector<float> Electron_Charge;
  std::vector<float> Electron_Gsf_deltaEtaEleClusterTrackAtCalo;
  std::vector<float> Electron_Gsf_deltaEtaSeedClusterTrackAtCalo;
  std::vector<float> Electron_Gsf_deltaEtaSuperClusterTrackAtVtx;
  std::vector<float> Electron_Gsf_deltaPhiEleClusterTrackAtCalo;
  std::vector<float> Electron_Gsf_deltaPhiSeedClusterTrackAtCalo;
  std::vector<float> Electron_Gsf_deltaPhiSuperClusterTrackAtVtx;
  std::vector<float> Electron_Gsf_dr03EcalRecHitSumE;
  std::vector<float> Electron_Gsf_dr03HcalDepth1TowerSumEt;
  std::vector<float> Electron_Gsf_dr03HcalDepth1TowerSumEtBc;
  std::vector<float> Electron_Gsf_dr03HcalDepth2TowerSumEt;
  std::vector<float> Electron_Gsf_dr03HcalDepth2TowerSumEtBc;
  std::vector<float> Electron_Gsf_dr03HcalTowerSumEt;
  std::vector<float> Electron_Gsf_dr03HcalTowerSumEtBc;
  std::vector<float> Electron_Gsf_dr03TkSumPt;
  std::vector<bool>  Electron_Gsf_passingCutBasedPreselection;
  std::vector<bool>  Electron_Gsf_passingMvaPreselection;
  std::vector<int>   Electron_gsftrack_trackerExpectedHitsInner_numberOfLostHits;
  std::vector<float> Electron_supercluster_e;
  std::vector<float> Electron_supercluster_phi;
  std::vector<float> Electron_supercluster_eta;
  std::vector<float> Electron_supercluster_centroid_x;
  std::vector<float> Electron_supercluster_centroid_y;
  std::vector<float> Electron_supercluster_centroid_z;
  std::vector<unsigned int> Electron_Track_idx;


  //=======  PFJets ===
  std::vector<std::vector<float> > PFJet_p4;
  std::vector<std::vector<float > > PFJet_Poca;
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
  std::vector<float> PFJet_chargedHadronEnergyFraction;
  std::vector<float> PFJet_neutralHadronEnergyFraction;
  std::vector<float> PFJet_neutralEmEnergyFraction;

  std::vector<float>   PFJet_partonFlavour;
  std::vector<float>   PFJet_bDiscriminator;
  std::vector<std::vector<float> >   PFJet_BTagWeight;



  //=======  MET ===
  // now only PFMET
  float MET_et;
  float MET_phi;
  float MET_sumET;
  float MET_metSignificance;
  float MET_MuonEtFraction;
  float MET_NeutralEMFraction;
  float MET_NeutralHadEtFraction;
  float MET_Type6EtFraction;
  float MET_Type7EtFraction;


  //=======  Event ===
  int EventNumber;
  unsigned int DataMC_Type_idx;
  unsigned int  Event_EventNumber;	 
  unsigned int  Event_RunNumber;	 
  int  Event_bunchCrossing;	 
  int  Event_orbitNumber;	 
  unsigned int  Event_luminosityBlock;	 
  bool  Event_isRealData;        
  int PileupInfo_NumInteractions_nm1;
  int PileupInfo_NumInteractions_n0;
  int PileupInfo_NumInteractions_np1;
  float EvtWeight3D;


  //====== Tracks ======= 
  std::vector<std::vector<float> > Track_p4;
  std::vector<std::vector<float > > Track_Poca;
  std::vector<int> Track_charge;
  std::vector<float> Track_chi2;
  std::vector<float> Track_ndof;
  std::vector<unsigned short> Track_numberOfLostHits;
  std::vector<unsigned short> Track_numberOfValidHits;
  std::vector<unsigned int> Track_qualityMask;
  std::vector<std::vector<float> > Track_par;
  std::vector<std::vector<std::vector<float> > > Track_parCov;

  //====== Trigger ======= 
  std::vector<std::string>  HTLTriggerName;
  std::vector<unsigned int> HLTPrescale;
  std::vector<unsigned int> NHLTL1GTSeeds;
  std::vector<unsigned int> L1SEEDPrescale;
  std::vector<bool>         L1SEEDInvalidPrescale;
  std::vector<bool>         TriggerAccept;
  std::vector<bool>         TriggerWasRun;
  std::vector<bool>         TriggerError;
  std::vector<std::vector<float> >        MuonTriggerMatch;
  std::vector<std::vector<float> >        ElectronTriggerMatch;
  std::vector<std::vector<float> >        JetTriggerMatch;
  std::vector<std::vector<float> >        TauTriggerMatch;
  bool TriggerOK;

  //====== MCTruth =======
  // verbose variables
  std::vector<double> GenEventInfoProduct_weights;
  float GenEventInfoProduct_signalProcessID;
  float GenEventInfoProduct_weight;
  float GenEventInfoProduct_qScale;
  float GenEventInfoProduct_alphaQED;
  float GenEventInfoProduct_alphaQCD;

  //do All
  std::vector<std::vector<float> > MC_p4;
  std::vector<int> MC_pdgid;
  std::vector<int> MC_charge;
  std::vector<unsigned int> MC_midx;


  // Signal particles Z, W, H0, Hpm
  std::vector<std::vector<float> > MCSignalParticle_p4;
  std::vector<int> MCSignalParticle_pdgid;
  std::vector<int> MCSignalParticle_charge;
  std::vector<std::vector<float > > MCSignalParticle_Poca;
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
  float MET_et;
  float MET_phi;
  float MET_sumET;
  float MET_et;
  float MET_phi;
  float MET_sumET;
  */




};
//define this as a plug-in
DEFINE_FWK_MODULE(TauNtuple);
#endif
