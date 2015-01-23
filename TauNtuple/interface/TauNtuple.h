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
// $Id: TauNtuple.h,v 1.39 2013/07/04 17:15:33 bkargoll Exp $
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
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"

// MET stuff
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

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
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "RecoTauTag/RecoTau/interface/RecoTauPiZeroPlugins.h"

// PU
#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/Lumi3DReWeighting.h"
//
#include <DataFormats/EgammaCandidates/interface/GsfElectron.h>
#include <RecoEgamma/EgammaTools/interface/ConversionTools.h>

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
#include "TauDataFormat/TauNtuple/interface/SkimmingCuts.h"
#include "TauDataFormat/TauNtuple/interface/MultiTriggerFilter.h"
#include <DataFormats/PatCandidates/interface/Jet.h>
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

// Electron MVA ID
#include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimator.h"

// embedded samples
#include "SimDataFormats/GeneratorProducts/interface/GenFilterInfo.h"

// sorting struct
struct sortIdxByValue {
    bool operator()(const std::pair<int,double> &left, const std::pair<int,double> &right) {
        return left.second > right.second;
    }
};

//
//
// class declaration
//

class TauNtuple: public edm::EDProducer {
	friend class SkimmingCuts;
	friend class MultiTriggerFilter;
public:
	explicit TauNtuple(const edm::ParameterSet&);
	~TauNtuple();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
	virtual void beginJob();
	virtual void produce(edm::Event&, const edm::EventSetup&);
	virtual void endJob();

	virtual void beginRun(edm::Run&, edm::EventSetup const&);
	virtual void endRun(edm::Run&, edm::EventSetup const&);
	virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
	virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

	// ----------member data ---------------------------
	void fillMCTruth(edm::Event& iEvent, const edm::EventSetup& iSetup);
	void fillPrimeVertex(edm::Event& iEvent, const edm::EventSetup& iSetup, edm::Handle<std::vector<reco::Track> > &trackCollection);
	void fillMuons(edm::Event& iEvent, const edm::EventSetup& iSetup, edm::Handle<std::vector<reco::Track> > &trackCollection);
	void fillTracks(edm::Handle<std::vector<reco::Track> > &trackCollection, const edm::EventSetup& iSetup);
	void fillPFTaus(edm::Event& iEvent, const edm::EventSetup& iSetup, edm::Handle<std::vector<reco::Track> > &trackCollection);
	void fillKinFitTaus(edm::Event& iEvent, const edm::EventSetup& iSetup, edm::Handle<std::vector<reco::Track> > &trackCollection);
	void fillPFJets(edm::Event& iEvent, const edm::EventSetup& iSetup, edm::Handle<std::vector<reco::Track> > &trackCollection);
	void fillElectrons(edm::Event& iEvent, const edm::EventSetup& iSetup, edm::Handle<std::vector<reco::Track> > &trackCollection);
	void fillMET(edm::Event& iEvent, const edm::EventSetup& iSetup);
	void fillTriggerInfo(edm::Event& iEvent, const edm::EventSetup& iSetup);
	template<class T>
	void TriggerMatch(edm::Handle<trigger::TriggerEvent> &triggerEvent, unsigned int triggerIndex, T obj, double drmax, std::vector<float> &match);
	void fillEventInfo(edm::Event& iEvent, const edm::EventSetup& iSetup);
	std::vector<bool> CheckTauDiscriminators(std::vector<edm::Handle<reco::PFTauDiscriminator> > tauDiscriminators, const reco::PFTauRef tauRef);
	reco::PFTauRef getMatchedHPSTau(edm::Handle<std::vector<reco::PFTau> > & HPStaus, std::vector<double> &UnmodifiedTau, int &match);
	reco::PFTauRef getHPSTauMatchedToJet(edm::Handle<std::vector<reco::PFTau> > & HPStaus, std::vector<double> &Jet, int &match);
	reco::PFJetRef getJetIndexMatchedToGivenHPSTauCandidate(edm::Handle<std::vector<reco::PFJet> > & PFJets, std::vector<double> &Tau, unsigned int &match);

	std::vector<reco::PFCandidatePtr> pfCandidates(const reco::PFJet& jet, int particleId, bool sort = true);
	std::vector<reco::PFCandidatePtr> pfCandidates(const reco::PFJet& jet, const std::vector<int>& particleIds, bool sort);
	std::vector<reco::PFCandidatePtr> pfphotons(const reco::PFJet& jet, bool sort = true);

	bool getTrackMatch(edm::Handle<std::vector<reco::Track> > &trackCollection, reco::TrackRef &refTrack, int &match);
	bool getTrackMatch(edm::Handle<std::vector<reco::Track> > &trackCollection, reco::GsfElectronRef &refElectron, int &match);
	double DeltaPhi(double phi1, double phi2);
	void ClearEvent();

	static bool isGoodMuon(reco::MuonRef &RefMuon);
	static bool isGoodTau(reco::PFTauRef &RefTau, edm::Handle<reco::PFTauDiscriminator> &Dis);
	static bool isGoodElectron(reco::GsfElectronRef &RefElectron);
	static bool isGoodVertex(const reco::Vertex &pv);
	static bool isGoodJet(reco::PFJetRef &RefJet);
	static bool isGoodJet(pat::JetRef &RefJet);
	static bool isGoodGenJet(reco::GenJetRef &RefGenJet);

	EGammaMvaEleEstimator* myMVATrigNoIP2012;
	EGammaMvaEleEstimator* myMVATrig2012;
	EGammaMvaEleEstimator* myMVANonTrig2012;
	std::vector<std::string> myManualCatWeightsTrigNoIP2012;
	std::vector<std::string> myManualCatWeightsTrig2012;
	std::vector<std::string> myManualCatWeightsNonTrig2012;

	static double MuonPtCut_;
	static double MuonEtaCut_;
	static double TauPtCut_;
	static double TauEtaCut_;
	static double ElectronPtCut_;
	static double ElectronEtaCut_;
	static double JetPtCut_;
	static double JetEtaCut_;

	bool RemoveMuonTracks_;
	bool RemoveElectronTracks_;
	edm::InputTag beamSpotTag_;
	bool useBeamSpot_;

	static edm::InputTag primVtxTag_;
	static edm::InputTag muonsTag_;
	static edm::InputTag hpsTauProducer_;
	static edm::InputTag PFElectronTag_;
	static std::vector<std::string> MyTriggerInfoNames;

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
	edm::InputTag hpsPFTauDiscriminationAgainstMuonMedium_;
	edm::InputTag hpsPFTauDiscriminationAgainstMuonTight_;
	edm::InputTag hpsPFTauDiscriminationByDecayModeFinding_;

	edm::InputTag pfMETCorrT0rt_;
	edm::InputTag pfMETCorrT0rtT1_;
	edm::InputTag pfMETCorrT0pc_;
	edm::InputTag pfMETCorrT0pcT1_;
	edm::InputTag pfMETCorrT0rtTxy_;
	edm::InputTag pfMETCorrT0rtT1Txy_;
	edm::InputTag pfMETCorrT0pcTxy_;
	edm::InputTag pfMETCorrT0pcT1Txy_;
	edm::InputTag pfMETCorrT1_;
	edm::InputTag pfMETCorrT1Txy_;
	edm::InputTag caloMETCorrT1_;
	edm::InputTag caloMETCorrT1T2_;
	edm::InputTag pfMETCorrMVA_;
	edm::InputTag muonsForPfMetCorrMVA_;
	edm::InputTag tausForPfMetCorrMVA_;
	edm::InputTag elecsForPfMetCorrMVA_;
	edm::InputTag pfMETCorrMVAMuTau_;
	edm::InputTag muonsForPfMetCorrMVAMuTau_;
	edm::InputTag tausForPfMetCorrMVAMuTau_;
	edm::InputTag pfMETUncorr_;
	edm::InputTag patMETCorrT0rt_;
	edm::InputTag patMETCorrT0rtT1_;
	edm::InputTag patMETCorrT0pc_;
	edm::InputTag patMETCorrT0pcT1_;
	edm::InputTag patMETCorrT0rtTxy_;
	edm::InputTag patMETCorrT0rtT1Txy_;
	edm::InputTag patMETCorrT0pcTxy_;
	edm::InputTag patMETCorrT0pcT1Txy_;
	edm::InputTag patMETCorrT1_;
	edm::InputTag patMETCorrT1Txy_;
	edm::InputTag patCaloMETCorrT1_;
	edm::InputTag patCaloMETCorrT1T2_;
	edm::InputTag patMETCorrMVA_;
	edm::InputTag patMETCorrMVAMuTau_;
	edm::InputTag patMETUncorr_;
	//input tags for MET systematics
	edm::InputTag patMETType1Corr_;
	edm::InputTag patMETType1CorrEleEnUp_;
	edm::InputTag patMETType1CorrEleEnDown_;
	edm::InputTag patMETType1CorrMuEnUp_;
	edm::InputTag patMETType1CorrMuEnDown_;
	edm::InputTag patMETType1CorrTauEnUp_;
	edm::InputTag patMETType1CorrTauEnDown_;
	edm::InputTag patMETType1CorrJetResUp_;
	edm::InputTag patMETType1CorrJetResDown_;
	edm::InputTag patMETType1CorrJetEnUp_;
	edm::InputTag patMETType1CorrJetEnDown_;
	edm::InputTag patMETType1CorrUnclusteredUp_;
	edm::InputTag patMETType1CorrUnclusteredDown_;
	edm::InputTag patMETType1p2Corr_;
	edm::InputTag patMETType1p2CorrEleEnUp_;
	edm::InputTag patMETType1p2CorrEleEnDown_;
	edm::InputTag patMETType1p2CorrMuEnUp_;
	edm::InputTag patMETType1p2CorrMuEnDown_;
	edm::InputTag patMETType1p2CorrTauEnUp_;
	edm::InputTag patMETType1p2CorrTauEnDown_;
	edm::InputTag patMETType1p2CorrJetResUp_;
	edm::InputTag patMETType1p2CorrJetResDown_;
	edm::InputTag patMETType1p2CorrJetEnUp_;
	edm::InputTag patMETType1p2CorrJetEnDown_;
	edm::InputTag patMETType1p2CorrUnclusteredUp_;
	edm::InputTag patMETType1p2CorrUnclusteredDown_;

	edm::InputTag pfjetsTag_;
	edm::InputTag genjetsTag_;
	edm::InputTag genjetsNoNuTag_;
	edm::InputTag rhoIsolAllInputTag_;
	edm::InputTag generalTracks_;
	edm::InputTag gensrc_;
	edm::InputTag GenEventInfo_;
	bool Embedded_; //embedding
	//edm::InputTag reducedEBRecHitCollection_;
	//edm::InputTag reducedEERecHitCollection_;
	// Electron MVA ID
	std::string ElectronMVATrigWeights1_;
	std::string ElectronMVATrigWeights2_;
	std::string ElectronMVATrigWeights3_;
	std::string ElectronMVATrigWeights4_;
	std::string ElectronMVATrigWeights5_;
	std::string ElectronMVATrigWeights6_;
	std::string ElectronMVATrigNoIPWeights1_;
	std::string ElectronMVATrigNoIPWeights2_;
	std::string ElectronMVATrigNoIPWeights3_;
	std::string ElectronMVATrigNoIPWeights4_;
	std::string ElectronMVATrigNoIPWeights5_;
	std::string ElectronMVATrigNoIPWeights6_;
	std::string ElectronMVANonTrigWeights1_;
	std::string ElectronMVANonTrigWeights2_;
	std::string ElectronMVANonTrigWeights3_;
	std::string ElectronMVANonTrigWeights4_;
	std::string ElectronMVANonTrigWeights5_;
	std::string ElectronMVANonTrigWeights6_;

	std::string JECuncData_;
	std::string JECuncMC_;

	// PU
	std::string ScaleFactor_;
	std::string PUInputFile_;
	std::string PUInputHistoMC_;
	std::string PUInputHistoData_;
	std::string PUInputHistoData_p5_;
	std::string PUInputHistoData_m5_;
	std::string PUInputHistoMCFineBins_;
	std::string PUInputHistoDataFineBins_;
	std::string PUOutputFile_;

	edm::LumiReWeighting LumiWeights_;
	edm::LumiReWeighting LumiWeights_p5_;
	edm::LumiReWeighting LumiWeights_m5_;
	edm::LumiReWeighting LumiWeightsFineBinning_;
	edm::Lumi3DReWeighting LumiWeights3D_;
	edm::Lumi3DReWeighting LumiWeights3D_p5_;
	edm::Lumi3DReWeighting LumiWeights3D_m5_;

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
	edm::Handle<L1GtTriggerMenuLite> triggerMenuLite_;
	bool doL1Triggers_;
	std::vector<std::string> l1TriggerNames_;
	std::vector<std::string> useFilterModules_;

	float TriggerJetMatchingdr_;
	float TriggerMuonMatchingdr_;
	float TriggerElectronMatchingdr_;
	float TriggerTauMatchingdr_;
	float TriggerMETMatchingdr_;

	// Control flags
	bool doBJets_;
	bool doPFJets_;
	bool doMuons_;
	bool doElectrons_;
	bool doPFTaus_;
	bool doTracks_;
	bool doTrigger_;
	bool doPrimeVertex_;
	bool doMET_;
	bool doMC_;
	bool doPatJets_;
	bool doPatElectrons_;
	bool doPatMuons_;
	bool doPatMET_;
	bool doMVAMET_;

	// Pat objects
	std::string srcPatJets_;
	std::string PatJetScale_;

	// b-tagging
	std::string BTagAlgorithm_;
	edm::InputTag jetFlavourTag_;

	// PUJetID
	edm::InputTag PUJetIdDisc_;
	edm::InputTag PUJetIdFlag_;

	// outputfile
	TFile *output;
	TTree *output_tree;

	int cnt_;

	//=======  Vertex ===
	std::vector<double> Vtx_chi2;
	std::vector<unsigned> Vtx_nTrk;
	std::vector<float> Vtx_ndof;
	std::vector<double> Vtx_y;
	std::vector<double> Vtx_x;
	std::vector<double> Vtx_z;
	std::vector<std::vector<std::vector<double> > > Vtx_Cov;
	std::vector<std::vector<int> > Vtx_Track_idx;
	std::vector<std::vector<float> > Vtx_Track_Weights;
	std::vector<bool> Vtx_isFake;

	std::vector<std::vector<std::vector<double> > > Vtx_TracksP4;

	//=======  Muons ===
	std::vector<std::vector<double> > Muon_p4;
	std::vector<std::vector<double> > Muon_Poca;
	std::vector<bool> Muon_isGlobalMuon;
	std::vector<bool> Muon_isStandAloneMuon;
	std::vector<bool> Muon_isTrackerMuon;
	std::vector<bool> Muon_isCaloMuon;
	std::vector<bool> Muon_isIsolationValid;
	std::vector<bool> Muon_isQualityValid;
	std::vector<bool> Muon_isTimeValid;
	std::vector<bool> Muon_isPFMuon;

	std::vector<float> Muon_emEt03;
	std::vector<float> Muon_emVetoEt03;
	std::vector<float> Muon_hadEt03;
	std::vector<float> Muon_hadVetoEt03;
	std::vector<int> Muon_nJets03;
	std::vector<int> Muon_nTracks03;
	std::vector<float> Muon_sumPt03;
	std::vector<float> Muon_trackerVetoPt03;

	std::vector<float> Muon_emEt05;
	std::vector<float> Muon_emVetoEt05;
	std::vector<float> Muon_hadEt05;
	std::vector<float> Muon_hadVetoEt05;
	std::vector<int> Muon_nJets05;
	std::vector<int> Muon_nTracks05;
	std::vector<float> Muon_sumPt05;
	std::vector<float> Muon_trackerVetoPt05;

	std::vector<float> Muon_sumChargedHadronPt03;              // sum-pt of charged Hadron
	std::vector<float> Muon_sumChargedParticlePt03;            // sum-pt of charged Particles(inludes e/mu)
	std::vector<float> Muon_sumNeutralHadronEt03;              // sum pt of neutral hadrons
	std::vector<float> Muon_sumNeutralHadronEtHighThreshold03; // sum pt of neutral hadrons with a higher threshold
	std::vector<float> Muon_sumPhotonEt03;                     // sum pt of PF photons
	std::vector<float> Muon_sumPhotonEtHighThreshold03;        // sum pt of PF photons with a higher threshold
	std::vector<float> Muon_sumPUPt03;                         // sum pt of charged Particles not from PV (for Pu corrections)

	std::vector<float> Muon_sumChargedHadronPt04;              // sum-pt of charged Hadron
	std::vector<float> Muon_sumChargedParticlePt04;            // sum-pt of charged Particles(inludes e/mu)
	std::vector<float> Muon_sumNeutralHadronEt04;              // sum pt of neutral hadrons
	std::vector<float> Muon_sumNeutralHadronEtHighThreshold04; // sum pt of neutral hadrons with a higher threshold
	std::vector<float> Muon_sumPhotonEt04;                     // sum pt of PF photons
	std::vector<float> Muon_sumPhotonEtHighThreshold04;        // sum pt of PF photons with a higher threshold
	std::vector<float> Muon_sumPUPt04;                         // sum pt of charged Particles not from PV (for Pu corrections)

	std::vector<int> Muon_numberOfChambers;
	std::vector<unsigned int> Muon_Track_idx;

	std::vector<int> Muon_charge;
	std::vector<int> Muon_trackCharge;
	std::vector<int> Muon_pdgid;
	std::vector<double> Muon_B;
	std::vector<double> Muon_M;
	std::vector<std::vector<double> > Muon_par;
	std::vector<std::vector<double> > Muon_cov;

	std::vector<int> Muon_hitPattern_pixelLayerwithMeas;
	std::vector<int> Muon_numberOfMatchedStations;
	std::vector<float> Muon_normChi2;
	std::vector<int> Muon_hitPattern_numberOfValidMuonHits;
	std::vector<int> Muon_innerTrack_numberofValidHits;
	std::vector<int> Muon_numberofValidPixelHits;
	std::vector<int> Muon_numberOfMatches;
	std::vector<int> Muon_trackerLayersWithMeasurement;

	//======= PFTaus ===
	std::vector<std::vector<double> > PFTau_p4;
	std::vector<std::vector<double> > PFTau_Poca;
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
	std::vector<bool> PFTau_isHPSAgainstMuonMedium;
	std::vector<bool> PFTau_isHPSAgainstMuonTight;
	std::vector<bool> PFTau_isHPSAgainstMuonLoose2;
	std::vector<bool> PFTau_isHPSAgainstMuonMedium2;
	std::vector<bool> PFTau_isHPSAgainstMuonTight2;
	std::vector<bool> PFTau_isHPSByDecayModeFinding;

	//  std::vector<bool> PFTau_HPSPFTauDiscriminationByMVA3rawElectronRejection;
	std::vector<bool> PFTau_HPSPFTauDiscriminationByMVA3LooseElectronRejection;
	std::vector<bool> PFTau_HPSPFTauDiscriminationByMVA3MediumElectronRejection;
	std::vector<bool> PFTau_HPSPFTauDiscriminationByMVA3TightElectronRejection;
	std::vector<bool> PFTau_HPSPFTauDiscriminationByMVA3VTightElectronRejection;
	//  std::vector<bool> PFTau_HPSPFTauDiscriminationByDeadECALElectronRejection;
	std::vector<bool> PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits;
	std::vector<bool> PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits;
	std::vector<bool> PFTau_HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits;
	std::vector<float> PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits;
	std::vector<bool> PFTau_HPSPFTauDiscriminationByLooseIsolationMVA;
	std::vector<bool> PFTau_HPSPFTauDiscriminationByMediumIsolationMVA;
	std::vector<bool> PFTau_HPSPFTauDiscriminationByTightIsolationMVA;
	std::vector<bool> PFTau_HPSPFTauDiscriminationByLooseIsolationMVA2;
	std::vector<bool> PFTau_HPSPFTauDiscriminationByMediumIsolationMVA2;
	std::vector<bool> PFTau_HPSPFTauDiscriminationByTightIsolationMVA2;
	std::vector<int> PFTau_hpsDecayMode;
	std::vector<int> PFTau_Charge;
	std::vector<std::vector<int> > PFTau_Track_idx;

	std::vector<std::vector<double> > PFTau_TIP_primaryVertex_pos;
	std::vector<std::vector<double> > PFTau_TIP_primaryVertex_cov;
	std::vector<std::vector<double> > PFTau_TIP_secondaryVertex_pos;
	std::vector<std::vector<double> > PFTau_TIP_secondaryVertex_cov;
	std::vector<std::vector<double> > PFTau_TIP_secondaryVertex_vtxchi2;
	std::vector<std::vector<double> > PFTau_TIP_secondaryVertex_vtxndof;
	std::vector<std::vector<double> > PFTau_TIP_primaryVertex_vtxchi2;
	std::vector<std::vector<double> > PFTau_TIP_primaryVertex_vtxndof;

	std::vector<std::vector<double> > PFTau_a1_lvp;
	std::vector<std::vector<double> > PFTau_a1_cov;
	std::vector<std::vector<int> > PFTau_a1_charge;
	std::vector<std::vector<int> > PFTau_a1_pdgid;
	std::vector<std::vector<double> > PFTau_a1_B;
	std::vector<std::vector<double> > PFTau_a1_M;

	std::vector<std::vector<std::vector<double> > > PFTau_daughterTracks;
	std::vector<std::vector<std::vector<double> > > PFTau_daughterTracks_cov;
	std::vector<std::vector<int> > PFTau_daughterTracks_charge;
	std::vector<std::vector<int> > PFTau_daughterTracks_pdgid;
	std::vector<std::vector<double> > PFTau_daughterTracks_B;
	std::vector<std::vector<double> > PFTau_daughterTracks_M;
	std::vector<std::vector<std::vector<double> > > PFTau_daughterTracks_poca;

	std::vector<std::vector<std::vector<double> > > PFTau_PionsP4;
	std::vector<std::vector<int> >  PFTau_PionsCharge;


	std::vector<std::vector<double> > PFTau_3PS_A1_LV;
	std::vector<std::vector<double> > PFTau_3PS_M_A1;
	std::vector<std::vector<double> > PFTau_3PS_M_12;
	std::vector<std::vector<double> > PFTau_3PS_M_13;
	std::vector<std::vector<double> > PFTau_3PS_M_23;
	std::vector<std::vector<int> > PFTau_3PS_Tau_Charge;
	std::vector<std::vector<float> > PFTau_3PS_LCchi2;
	std::vector<std::vector<int> > PFTau_3PS_has3ProngSolution;
	std::vector<std::vector<std::vector<double> > > PFTau_3PS_Tau_LV;

	std::vector<std::vector<double> > PFTau_TIP_flightLength;
	std::vector<std::vector<double> > PFTau_TIP_flightLengthSig;

	std::vector<std::vector<std::vector<double> > > PFTau_PiZeroP4;
	std::vector<std::vector<int> > PFTau_PiZeroNumOfPhotons;
	std::vector<std::vector<int> > PFTau_PiZeroNumOfElectrons;
	std::vector<std::vector<std::vector<double> > > PFTau_ChargedHadronsP4;
	std::vector<std::vector<std::vector<int> > > PFTau_ChargedHadronsCharge;
	std::vector<std::vector<std::vector<double> > > PFTau_GammaP4;
	std::vector<std::vector<std::vector<double> > > PFTau_Photons_p4_inDR05;
	//-------- Gamma information ---------

	std::vector<std::vector<double> > PFTau_MatchedPFJetP4;
	std::vector<std::vector<std::vector<double> > > PFTau_MatchedPFJetGammasP4;
	std::vector<std::vector<std::vector<double> > > PFTau_MatchedPFJetSCVariables;
	std::vector<std::vector<std::vector<double> > > PFTau_MatchedPFJetPhotonVariables;

	std::vector<double> PFTau_PhotonEnergyFraction;

        std::vector<std::vector<int> > PFTau_photon_hasPixelSeed;
        std::vector<std::vector<float> > PFTau_photon_hadronicOverEm;
        std::vector<std::vector<float> > PFTau_photon_sigmaIetaIeta;
        std::vector<std::vector<float> > PFTau_photon_trkSumPtHollowConeDR04;
        std::vector<std::vector<float> > PFTau_photon_ecalRecHitSumEtConeDR04;
        std::vector<std::vector<float> > PFTau_photon_hcalTowerSumEtConeDR04;
        std::vector<float > PFTau_photon_rho;


//   std::vector<std::vector<int> > PFTau_hasSC;  
//   std::vector<std::vector<int> > PFTau_hasPhoton;

	//=======  Electrons ===
	double RhoIsolationAllInputTags;
	std::vector<std::vector<double> > Electron_p4;
	std::vector<std::vector<double> > Electron_Poca;
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
	std::vector<bool> Electron_Gsf_passingCutBasedPreselection;
	std::vector<bool> Electron_Gsf_passingMvaPreselection;
	std::vector<int> Electron_gsftrack_trackerExpectedHitsInner_numberOfLostHits;
	std::vector<double> Electron_supercluster_e;
	std::vector<double> Electron_supercluster_phi;
	std::vector<double> Electron_supercluster_eta;
	std::vector<float> Electron_supercluster_centroid_x;
	std::vector<float> Electron_supercluster_centroid_y;
	std::vector<float> Electron_supercluster_centroid_z;
	std::vector<unsigned int> Electron_Track_idx;

	std::vector<float> Electron_ecalRecHitSumEt03;
	std::vector<float> Electron_hcalDepth1TowerSumEt03;
	std::vector<float> Electron_hcalDepth1TowerSumEtBc03;
	std::vector<float> Electron_hcalDepth2TowerSumEt03;
	std::vector<float> Electron_hcalDepth2TowerSumEtBc03;
	std::vector<float> Electron_tkSumPt03;
	std::vector<float> Electron_ecalRecHitSumEt04;
	std::vector<float> Electron_hcalDepth1TowerSumEt04;
	std::vector<float> Electron_hcalDepth1TowerSumEtBc04;
	std::vector<float> Electron_hcalDepth2TowerSumEt04;
	std::vector<float> Electron_hcalDepth2TowerSumEtBc04;
	std::vector<float> Electron_tkSumPt04;

	std::vector<float> Electron_chargedHadronIso;
	std::vector<float> Electron_neutralHadronIso;
	std::vector<float> Electron_photonIso;
	std::vector<double> Electron_isoDeposits_chargedHadronIso04;
	std::vector<double> Electron_isoDeposits_neutralHadronIso04;
	std::vector<double> Electron_isoDeposits_photonIso04;
	std::vector<double> Electron_isoDeposits_chargedHadronIso03;
	std::vector<double> Electron_isoDeposits_neutralHadronIso03;
	std::vector<double> Electron_isoDeposits_photonIso03;

	std::vector<float> Electron_sigmaIetaIeta;
	std::vector<float> Electron_hadronicOverEm;
	std::vector<float> Electron_fbrem;
	std::vector<float> Electron_eSuperClusterOverP;
	std::vector<float> Electron_ecalEnergy;
	std::vector<float> Electron_trackMomentumAtVtx;
	std::vector<int> Electron_numberOfMissedHits;  //number of missing hits conversion rejection
	std::vector<bool> Electron_HasMatchedConversions;

	// Electron energy calibration
	std::vector<double> Electron_RegEnergy;
	std::vector<double> Electron_RegEnergyError;

	// Electron MVA ID
	std::vector<double> Electron_Rho_kt6PFJets;
	std::vector<double> Electron_MVA_Trig_discriminator;
	std::vector<double> Electron_MVA_TrigNoIP_discriminator;
	std::vector<double> Electron_MVA_NonTrig_discriminator;

	std::vector<int> Electron_charge;
	std::vector<int> Electron_trackCharge;
	std::vector<int> Electron_pdgid;
	std::vector<double> Electron_B;
	std::vector<double> Electron_M;
	std::vector<std::vector<double> > Electron_par;
	std::vector<std::vector<double> > Electron_cov;

	//=======  PFJets ===
	std::vector<std::vector<double> > PFJet_p4;
	std::vector<float> PFJet_chargedEmEnergy;
	std::vector<float> PFJet_chargedHadronEnergy;
	std::vector<int> PFJet_chargedHadronMultiplicity;
	std::vector<float> PFJet_chargedMuEnergy;
	std::vector<int> PFJet_chargedMultiplicity;
	std::vector<float> PFJet_electronEnergy;
	std::vector<int> PFJet_electronMultiplicity;
	std::vector<float> PFJet_HFEMEnergy;
	std::vector<int> PFJet_HFEMMultiplicity;
	std::vector<float> PFJet_HFHadronEnergy;
	std::vector<int> PFJet_HFHadronMultiplicity;
	std::vector<float> PFJet_muonEnergy;
	std::vector<int> PFJet_muonMultiplicity;
	std::vector<float> PFJet_neutralEmEnergy;
	std::vector<float> PFJet_neutralHadronEnergy;
	std::vector<int> PFJet_neutralHadronMultiplicity;
	std::vector<float> PFJet_photonEnergy;
	std::vector<int> PFJet_photonMultiplicity;
	std::vector<float> PFJet_jetArea;
	std::vector<float> PFJet_maxDistance;
	std::vector<int> PFJet_nConstituents;
	std::vector<float> PFJet_pileup;
	std::vector<float> PFJet_etaetaMoment;
	std::vector<float> PFJet_etaphiMoment;
	std::vector<std::vector<int> > PFJet_Track_idx;
	std::vector<int> PFJet_MatchedHPS_idx;

	std::vector<float> PFJet_PUJetID_discr;
	std::vector<bool> PFJet_PUJetID_looseWP;
	std::vector<bool> PFJet_PUJetID_mediumWP;
	std::vector<bool> PFJet_PUJetID_tightWP;

	std::vector<std::vector<std::vector<double> > > PFJet_TracksP4;
	std::vector<int> PFJet_nTrk;

	std::vector<int> PFJet_numberOfDaughters;
	std::vector<float> PFJet_chargedEmEnergyFraction;
	std::vector<float> PFJet_chargedHadronEnergyFraction;
	std::vector<float> PFJet_neutralHadronEnergyFraction;
	std::vector<float> PFJet_neutralEmEnergyFraction;

	std::vector<int> PFJet_partonFlavour;
	std::vector<float> PFJet_bDiscriminator;
	std::vector<std::vector<float> > PFJet_BTagWeight;
	//std::vector<std::string> PFJet_bTagAlgorithmName;
	//std::vector<float>   PFJet_bTagAlgorithmValue;

	std::vector<float> PFJet_JECuncertainty;

	std::vector<std::vector<double> > PFJet_GenJet_p4;
	std::vector<std::vector<std::vector<double> > > PFJet_GenJet_Constituents_p4;
	std::vector<std::vector<double> > PFJet_GenJetNoNu_p4;
	std::vector<std::vector<std::vector<double> > > PFJet_GenJetNoNu_Constituents_p4;

	//=======  MET ===
	// now only PFMET
	double MET_Uncorr_et;
	float MET_Uncorr_pt;
	double MET_Uncorr_phi;
	float MET_Uncorr_sumET;
	double MET_Uncorr_significance;
	double MET_Uncorr_significance_xx;
	double MET_Uncorr_significance_xy;
	double MET_Uncorr_significance_yy;
	float MET_Uncorr_MuonEtFraction;
	float MET_Uncorr_NeutralEMFraction;
	float MET_Uncorr_NeutralHadEtFraction;
	float MET_Uncorr_Type6EtFraction;
	float MET_Uncorr_Type7EtFraction;

	double MET_CorrT0rt_et;
	float MET_CorrT0rt_pt;
	double MET_CorrT0rt_phi;
	float MET_CorrT0rt_sumET;
	float MET_CorrT0rt_MuonEtFraction;
	float MET_CorrT0rt_NeutralEMFraction;
	float MET_CorrT0rt_NeutralHadEtFraction;
	float MET_CorrT0rt_Type6EtFraction;
	float MET_CorrT0rt_Type7EtFraction;

	double MET_CorrT0rtT1_et;
	float MET_CorrT0rtT1_pt;
	double MET_CorrT0rtT1_phi;
	float MET_CorrT0rtT1_sumET;
	float MET_CorrT0rtT1_MuonEtFraction;
	float MET_CorrT0rtT1_NeutralEMFraction;
	float MET_CorrT0rtT1_NeutralHadEtFraction;
	float MET_CorrT0rtT1_Type6EtFraction;
	float MET_CorrT0rtT1_Type7EtFraction;

	double MET_CorrT0pc_et;
	float MET_CorrT0pc_pt;
	double MET_CorrT0pc_phi;
	float MET_CorrT0pc_sumET;
	float MET_CorrT0pc_MuonEtFraction;
	float MET_CorrT0pc_NeutralEMFraction;
	float MET_CorrT0pc_NeutralHadEtFraction;
	float MET_CorrT0pc_Type6EtFraction;
	float MET_CorrT0pc_Type7EtFraction;

	double MET_CorrT0pcT1_et;
	float MET_CorrT0pcT1_pt;
	double MET_CorrT0pcT1_phi;
	float MET_CorrT0pcT1_sumET;
	float MET_CorrT0pcT1_MuonEtFraction;
	float MET_CorrT0pcT1_NeutralEMFraction;
	float MET_CorrT0pcT1_NeutralHadEtFraction;
	float MET_CorrT0pcT1_Type6EtFraction;
	float MET_CorrT0pcT1_Type7EtFraction;

	double MET_CorrT0rtTxy_et;
	float MET_CorrT0rtTxy_pt;
	double MET_CorrT0rtTxy_phi;
	float MET_CorrT0rtTxy_sumET;
	float MET_CorrT0rtTxy_MuonEtFraction;
	float MET_CorrT0rtTxy_NeutralEMFraction;
	float MET_CorrT0rtTxy_NeutralHadEtFraction;
	float MET_CorrT0rtTxy_Type6EtFraction;
	float MET_CorrT0rtTxy_Type7EtFraction;

	double MET_CorrT0rtT1Txy_et;
	float MET_CorrT0rtT1Txy_pt;
	double MET_CorrT0rtT1Txy_phi;
	float MET_CorrT0rtT1Txy_sumET;
	float MET_CorrT0rtT1Txy_MuonEtFraction;
	float MET_CorrT0rtT1Txy_NeutralEMFraction;
	float MET_CorrT0rtT1Txy_NeutralHadEtFraction;
	float MET_CorrT0rtT1Txy_Type6EtFraction;
	float MET_CorrT0rtT1Txy_Type7EtFraction;

	double MET_CorrT0pcTxy_et;
	float MET_CorrT0pcTxy_pt;
	double MET_CorrT0pcTxy_phi;
	float MET_CorrT0pcTxy_sumET;
	float MET_CorrT0pcTxy_MuonEtFraction;
	float MET_CorrT0pcTxy_NeutralEMFraction;
	float MET_CorrT0pcTxy_NeutralHadEtFraction;
	float MET_CorrT0pcTxy_Type6EtFraction;
	float MET_CorrT0pcTxy_Type7EtFraction;

	double MET_CorrT0pcT1Txy_et;
	float MET_CorrT0pcT1Txy_pt;
	double MET_CorrT0pcT1Txy_phi;
	float MET_CorrT0pcT1Txy_sumET;
	float MET_CorrT0pcT1Txy_MuonEtFraction;
	float MET_CorrT0pcT1Txy_NeutralEMFraction;
	float MET_CorrT0pcT1Txy_NeutralHadEtFraction;
	float MET_CorrT0pcT1Txy_Type6EtFraction;
	float MET_CorrT0pcT1Txy_Type7EtFraction;

	double MET_CorrT1_et;
	float MET_CorrT1_pt;
	double MET_CorrT1_phi;
	float MET_CorrT1_sumET;
	float MET_CorrT1_MuonEtFraction;
	float MET_CorrT1_NeutralEMFraction;
	float MET_CorrT1_NeutralHadEtFraction;
	float MET_CorrT1_Type6EtFraction;
	float MET_CorrT1_Type7EtFraction;

	double MET_CorrT1Txy_et;
	float MET_CorrT1Txy_pt;
	double MET_CorrT1Txy_phi;
	float MET_CorrT1Txy_sumET;
	float MET_CorrT1Txy_MuonEtFraction;
	float MET_CorrT1Txy_NeutralEMFraction;
	float MET_CorrT1Txy_NeutralHadEtFraction;
	float MET_CorrT1Txy_Type6EtFraction;
	float MET_CorrT1Txy_Type7EtFraction;

	double MET_CorrCaloT1_et;
	float MET_CorrCaloT1_pt;
	double MET_CorrCaloT1_phi;
	float MET_CorrCaloT1_sumET;

	double MET_CorrCaloT1T2_et;
	float MET_CorrCaloT1T2_pt;
	double MET_CorrCaloT1T2_phi;
	float MET_CorrCaloT1T2_sumET;

	double MET_CorrMVA_et;
	float MET_CorrMVA_pt;
	double MET_CorrMVA_phi;
	float MET_CorrMVA_sumET;
	double MET_CorrMVA_significance;
	double MET_CorrMVA_significance_xx;
	double MET_CorrMVA_significance_xy;
	double MET_CorrMVA_significance_yy;
	float MET_CorrMVA_MuonEtFraction;
	float MET_CorrMVA_NeutralEMFraction;
	float MET_CorrMVA_NeutralHadEtFraction;
	float MET_CorrMVA_Type6EtFraction;
	float MET_CorrMVA_Type7EtFraction;
	std::vector<std::vector<double> > MET_CorrMVA_srcMuon_p4;
	std::vector<std::vector<double> > MET_CorrMVA_srcElectron_p4;
	std::vector<std::vector<double> > MET_CorrMVA_srcTau_p4;

	double MET_CorrMVAMuTau_et;
	float MET_CorrMVAMuTau_pt;
	double MET_CorrMVAMuTau_phi;
	float MET_CorrMVAMuTau_sumET;
	double MET_CorrMVAMuTau_significance;
	double MET_CorrMVAMuTau_significance_xx;
	double MET_CorrMVAMuTau_significance_xy;
	double MET_CorrMVAMuTau_significance_yy;
	float MET_CorrMVAMuTau_MuonEtFraction;
	float MET_CorrMVAMuTau_NeutralEMFraction;
	float MET_CorrMVAMuTau_NeutralHadEtFraction;
	float MET_CorrMVAMuTau_Type6EtFraction;
	float MET_CorrMVAMuTau_Type7EtFraction;
	std::vector<std::vector<double> > MET_CorrMVAMuTau_srcMuon_p4;
	std::vector<std::vector<double> > MET_CorrMVAMuTau_srcTau_p4;

	double MET_Type1Corr_et;
	float MET_Type1Corr_pt;
	double MET_Type1Corr_phi;
	float MET_Type1Corr_sumET;
	float MET_Type1Corr_MuonEtFraction;
	float MET_Type1Corr_NeutralEMFraction;
	float MET_Type1Corr_NeutralHadEtFraction;
	float MET_Type1Corr_Type6EtFraction;
	float MET_Type1Corr_Type7EtFraction;

	double MET_Type1p2Corr_et;
	float MET_Type1p2Corr_pt;
	double MET_Type1p2Corr_phi;
	float MET_Type1p2Corr_sumET;
	float MET_Type1p2Corr_MuonEtFraction;
	float MET_Type1p2Corr_NeutralEMFraction;
	float MET_Type1p2Corr_NeutralHadEtFraction;
	float MET_Type1p2Corr_Type6EtFraction;
	float MET_Type1p2Corr_Type7EtFraction;

	double MET_Type1CorrElectronUp_et;
	double MET_Type1CorrElectronDown_et;
	double MET_Type1CorrMuonUp_et;
	double MET_Type1CorrMuonDown_et;
	double MET_Type1CorrTauUp_et;
	double MET_Type1CorrTauDown_et;
	double MET_Type1CorrJetResUp_et;
	double MET_Type1CorrJetResDown_et;
	double MET_Type1CorrJetEnUp_et;
	double MET_Type1CorrJetEnDown_et;
	double MET_Type1CorrUnclusteredUp_et;
	double MET_Type1CorrUnclusteredDown_et;

	double MET_Type1p2CorrElectronUp_et;
	double MET_Type1p2CorrElectronDown_et;
	double MET_Type1p2CorrMuonUp_et;
	double MET_Type1p2CorrMuonDown_et;
	double MET_Type1p2CorrTauUp_et;
	double MET_Type1p2CorrTauDown_et;
	double MET_Type1p2CorrJetResUp_et;
	double MET_Type1p2CorrJetResDown_et;
	double MET_Type1p2CorrJetEnUp_et;
	double MET_Type1p2CorrJetEnDown_et;
	double MET_Type1p2CorrUnclusteredUp_et;
	double MET_Type1p2CorrUnclusteredDown_et;

	//=======  Event ===
	int EventNumber;
	unsigned int DataMC_Type_idx;
	unsigned int Event_EventNumber;
	unsigned int Event_RunNumber;
	int Event_bunchCrossing;
	int Event_orbitNumber;
	unsigned int Event_luminosityBlock;
	bool Event_isRealData;
	float PileupInfo_TrueNumInteractions_nm1;
	float PileupInfo_TrueNumInteractions_n0;
	float PileupInfo_TrueNumInteractions_np1;
	float PUWeight;
	float PUWeight_p5;
	float PUWeight_m5;
	float PUWeight3D;
	float PUWeight3D_p5;
	float PUWeight3D_m5;
	float PUWeightFineBins;

	std::vector<double> beamspot_par;
	std::vector<double> beamspot_cov;
	double beamspot_emittanceX;
	double beamspot_emittanceY;
	double beamspot_betaStar;

	// for embedded samples
	float TauSpinnerWeight;
	float SelEffWeight;
	float MinVisPtFilter;
	float KinWeightPt;
	float KinWeightEta;
	float KinWeightMassPt;
	float EmbeddedWeight;

	//====== Tracks =======
	std::vector<std::vector<double> > Track_p4;
	std::vector<std::vector<double> > Track_Poca;
	std::vector<double> Track_chi2;
	std::vector<double> Track_ndof;
	std::vector<unsigned short> Track_numberOfLostHits;
	std::vector<unsigned short> Track_numberOfValidHits;
	std::vector<unsigned int> Track_qualityMask;

	std::vector<int> Track_charge;
	std::vector<int> Track_pdgid;
	std::vector<double> Track_B;
	std::vector<double> Track_M;
	std::vector<std::vector<double> > Track_par;
	std::vector<std::vector<double> > Track_cov;

	//====== Trigger =======
	std::vector<std::string> HTLTriggerName;
	std::vector<unsigned int> HLTPrescale;
	std::vector<unsigned int> NHLTL1GTSeeds;
	std::vector<unsigned int> L1SEEDPrescale;
	std::vector<bool> L1SEEDisTechBit;
	std::vector<bool> L1SEEDInvalidPrescale;
	std::vector<bool> TriggerAccept;
	std::vector<bool> TriggerWasRun;
	std::vector<bool> TriggerError;
	std::vector<std::vector<float> > TriggerMatchMuon;
	std::vector<std::vector<float> > TriggerMatchElectron;
	std::vector<std::vector<float> > TriggerMatchJet;
	std::vector<std::vector<float> > TriggerMatchTau;
	std::vector<std::vector<float> > HLTTrigger_objs_Pt;
	std::vector<std::vector<float> > HLTTrigger_objs_Eta;
	std::vector<std::vector<float> > HLTTrigger_objs_Phi;
	std::vector<std::vector<float> > HLTTrigger_objs_E;
	std::vector<std::vector<int> > HLTTrigger_objs_Id;
	std::vector<std::string> HLTTrigger_objs_trigger;

	std::vector<std::string> L1TriggerName;
	std::vector<bool> L1TriggerDecision;
	std::vector<int> L1ErrorCode;
	std::vector<unsigned int> L1Prescale;

	bool TriggerOK;

	//====== MCTruth =======
	// verbose variables
	std::vector<double> GenEventInfoProduct_weights;
	unsigned GenEventInfoProduct_signalProcessID;
	float GenEventInfoProduct_weight;
	float GenEventInfoProduct_qScale;
	float GenEventInfoProduct_alphaQED;
	float GenEventInfoProduct_alphaQCD;
	int GenEventInfoProduct_id1;
	int GenEventInfoProduct_id2;
	double GenEventInfoProduct_x1;
	double GenEventInfoProduct_x2;
	double GenEventInfoProduct_scalePDF;

	//do All
	std::vector<std::vector<double> > MC_p4;
	std::vector<int> MC_pdgid;
	std::vector<std::vector<int> > MC_childpdgid;
	std::vector<int> MC_charge;
	std::vector<unsigned int> MC_midx;
	std::vector<int> MC_status;

	// Signal particles Z, W, H0, Hpm
	std::vector<std::vector<double> > MCSignalParticle_p4;
	std::vector<int> MCSignalParticle_pdgid;
	std::vector<std::vector<int> > MCSignalParticle_childpdgid;
	std::vector<int> MCSignalParticle_charge;
	std::vector<std::vector<double> > MCSignalParticle_Poca;
	std::vector<std::vector<unsigned int> > MCSignalParticle_Tauidx;

	// MC Tau Info
	std::vector<std::vector<std::vector<double> > > MCTauandProd_p4;
	std::vector<std::vector<std::vector<double> > > MCTauandProd_Vertex;
	std::vector<std::vector<int> > MCTauandProd_pdgid;
	std::vector<std::vector<unsigned int> > MCTauandProd_midx;
	std::vector<std::vector<int> > MCTauandProd_charge;
	std::vector<unsigned int> MCTau_JAK;
	std::vector<unsigned int> MCTau_DecayBitMask;

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

#endif
