#include "TauDataFormat/TauNtuple/interface/TauNtuple.h"
#include "TauDataFormat/TauNtuple/interface/TauDecay_CMSSW.h"
#include "Validation/EventGenerator/interface/PdtPdgMini.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"
#include <vector>
#include <map>
#include "TMatrixT.h"

#include "DataFormats/JetReco/interface/PileupJetIdentifier.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include <sys/types.h>
#include <dirent.h>
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/ParticleFlowReco/interface/GsfPFRecTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

#include <DataFormats/EgammaReco/interface/SuperCluster.h>

#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/PatCandidates/interface/MET.h>

// MVA electron ID
#include <cmath>
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

//Simple Fits
#include "SimpleFits/FitSoftware/interface/Particle.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "TauDataFormat/TauNtuple/interface/ParticleBuilder.h"
#include "SimpleFits/FitSoftware/interface/TauA1NuConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/Chi2VertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TLorentzVector.h"
#include "SimpleFits/FitSoftware/interface/Chi2VertexFitter.h"                                                                                                                                                                         
#include "SimpleFits/FitSoftware/interface/ChiSquareFunctionUpdator.h"                                                                                                                                                                 
#include "SimpleFits/FitSoftware/interface/TauA1NuConstrainedFitter.h"                                                                                                                                                                 
#include "SimpleFits/FitSoftware/interface/LagrangeMultipliersFitter.h"                                                                                                                                                                
#include "SimpleFits/FitSoftware/interface/TrackHelixVertexFitter.h"                                                                                                                                                                   
#include "SimpleFits/FitSoftware/interface/TrackTools.h"                                                                                                                                                                               
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"          

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

// JEC uncertainties
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include "DataFormats/Common/interface/ValueMap.h"

double TauNtuple::MuonPtCut_(-1.);
double TauNtuple::MuonEtaCut_(999);
double TauNtuple::TauPtCut_(-1.);
double TauNtuple::TauEtaCut_(999);
double TauNtuple::ElectronPtCut_(-1.);
double TauNtuple::ElectronEtaCut_(999);
double TauNtuple::JetPtCut_(-1.);
double TauNtuple::JetEtaCut_(999);

edm::InputTag TauNtuple::primVtxTag_;
edm::InputTag TauNtuple::muonsTag_;
edm::InputTag TauNtuple::hpsTauProducer_;
edm::InputTag TauNtuple::PFElectronTag_;

std::vector<std::string> TauNtuple::MyTriggerInfoNames;

TauNtuple::TauNtuple(const edm::ParameterSet& iConfig) :
		RemoveMuonTracks_(iConfig.getUntrackedParameter("RemoveMuonTracks", (bool) true)),
		RemoveElectronTracks_(iConfig.getUntrackedParameter("RemoveuonTracks", (bool) true)),
		beamSpotTag_(iConfig.getParameter<edm::InputTag>("beamSpot")),
		useBeamSpot_(iConfig.getUntrackedParameter("useBeamSpot", (bool) true)),
		hpsPFTauDiscriminationByTightIsolation_(iConfig.getParameter<edm::InputTag>("hpsPFTauDiscriminationByTightIsolation")),
		hpsPFTauDiscriminationByMediumIsolation_(iConfig.getParameter<edm::InputTag>("hpsPFTauDiscriminationByMediumIsolation")),
		hpsPFTauDiscriminationByLooseIsolation_(iConfig.getParameter<edm::InputTag>("hpsPFTauDiscriminationByLooseIsolation")),
		hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr_(iConfig.getParameter<edm::InputTag>("hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr")),
		hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr_(iConfig.getParameter<edm::InputTag>("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr")),
		hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr_(iConfig.getParameter<edm::InputTag>("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr")),
		hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr_(iConfig.getParameter<edm::InputTag>("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr")),
		hpsPFTauDiscriminationAgainstElectronLoose_(iConfig.getParameter<edm::InputTag>("hpsPFTauDiscriminationAgainstElectronLoose")),
		hpsPFTauDiscriminationAgainstElectronMedium_(iConfig.getParameter<edm::InputTag>("hpsPFTauDiscriminationAgainstElectronMedium")),
		hpsPFTauDiscriminationAgainstElectronTight_(iConfig.getParameter<edm::InputTag>("hpsPFTauDiscriminationAgainstElectronTight")),
		hpsPFTauDiscriminationAgainstMuonLoose_(iConfig.getParameter<edm::InputTag>("hpsPFTauDiscriminationAgainstMuonLoose")),
		hpsPFTauDiscriminationAgainstMuonMedium_(iConfig.getParameter<edm::InputTag>("hpsPFTauDiscriminationAgainstMuonMedium")),
		hpsPFTauDiscriminationAgainstMuonTight_(iConfig.getParameter<edm::InputTag>("hpsPFTauDiscriminationAgainstMuonTight")),
		hpsPFTauDiscriminationByDecayModeFinding_(iConfig.getParameter<edm::InputTag>("hpsPFTauDiscriminationByDecayModeFinding")),
		pfMETCorrT0rt_(iConfig.getParameter<edm::InputTag>("pfMetCorrT0rt")),
		pfMETCorrT0rtT1_(iConfig.getParameter<edm::InputTag>("pfMetCorrT0rtT1")),
		pfMETCorrT0pc_(iConfig.getParameter<edm::InputTag>("pfMetCorrT0pc")),
		pfMETCorrT0pcT1_(iConfig.getParameter<edm::InputTag>("pfMetCorrT0pcT1")),
		pfMETCorrT0rtTxy_(iConfig.getParameter<edm::InputTag>("pfMetCorrT0rtTxy")),
		pfMETCorrT0rtT1Txy_(iConfig.getParameter<edm::InputTag>("pfMetCorrT0rtT1Txy")),
		pfMETCorrT0pcTxy_(iConfig.getParameter<edm::InputTag>("pfMetCorrT0pcTxy")),
		pfMETCorrT0pcT1Txy_(iConfig.getParameter<edm::InputTag>("pfMetCorrT0pcT1Txy")),
		pfMETCorrT1_(iConfig.getParameter<edm::InputTag>("pfMetCorrT1")),
		pfMETCorrT1Txy_(iConfig.getParameter<edm::InputTag>("pfMetCorrT1Txy")),
		caloMETCorrT1_(iConfig.getParameter<edm::InputTag>("caloMetCorrT1")),
		caloMETCorrT1T2_(iConfig.getParameter<edm::InputTag>("caloMetCorrT1T2")),
		pfMETCorrMVA_(iConfig.getParameter<edm::InputTag>("pfMetCorrMVA")),
		muonsForPfMetCorrMVA_(iConfig.getParameter<edm::InputTag>("muonsForPfMetCorrMVA")),
		tausForPfMetCorrMVA_(iConfig.getParameter<edm::InputTag>("tausForPfMetCorrMVA")),
		elecsForPfMetCorrMVA_(iConfig.getParameter<edm::InputTag>("elecsForPfMetCorrMVA")),
		pfMETCorrMVAMuTau_(iConfig.getParameter<edm::InputTag>("pfMetCorrMVAMuTau")),
		muonsForPfMetCorrMVAMuTau_(iConfig.getParameter<edm::InputTag>("muonsForPfMetCorrMVAMuTau")),
		tausForPfMetCorrMVAMuTau_(iConfig.getParameter<edm::InputTag>("tausForPfMetCorrMVAMuTau")),
		pfMETUncorr_(iConfig.getParameter<edm::InputTag>("pfMetUncorr")),
		patMETCorrT0rt_(iConfig.getParameter<edm::InputTag>("patMetCorrT0rt")),
		patMETCorrT0rtT1_(iConfig.getParameter<edm::InputTag>("patMetCorrT0rtT1")),
		patMETCorrT0pc_(iConfig.getParameter<edm::InputTag>("patMetCorrT0pc")),
		patMETCorrT0pcT1_(iConfig.getParameter<edm::InputTag>("patMetCorrT0pcT1")),
		patMETCorrT0rtTxy_(iConfig.getParameter<edm::InputTag>("patMetCorrT0rtTxy")),
		patMETCorrT0rtT1Txy_(iConfig.getParameter<edm::InputTag>("patMetCorrT0rtT1Txy")),
		patMETCorrT0pcTxy_(iConfig.getParameter<edm::InputTag>("patMetCorrT0pcTxy")),
		patMETCorrT0pcT1Txy_(iConfig.getParameter<edm::InputTag>("patMetCorrT0pcT1Txy")),
		patMETCorrT1_(iConfig.getParameter<edm::InputTag>("patMetCorrT1")),
		patMETCorrT1Txy_(iConfig.getParameter<edm::InputTag>("patMetCorrT1Txy")),
		patCaloMETCorrT1_(iConfig.getParameter<edm::InputTag>("patCaloMetCorrT1")),
		patCaloMETCorrT1T2_(iConfig.getParameter<edm::InputTag>("patCaloMetCorrT1T2")),
		patMETCorrMVA_(iConfig.getParameter<edm::InputTag>("patMetCorrMVA")),
		patMETCorrMVAMuTau_(iConfig.getParameter<edm::InputTag>("patMetCorrMVAMuTau")),
		patMETUncorr_(iConfig.getParameter<edm::InputTag>("patMetUncorr")),
		patMETType1Corr_(iConfig.getParameter<edm::InputTag>("patMETType1Corr")),
		patMETType1CorrEleEnUp_(iConfig.getParameter<edm::InputTag>("patMETType1CorrEleEnUp")),
		patMETType1CorrEleEnDown_(iConfig.getParameter<edm::InputTag>("patMETType1CorrEleEndown")),
		patMETType1CorrMuEnUp_(iConfig.getParameter<edm::InputTag>("patMETType1CorrMuEnUp")),
		patMETType1CorrMuEnDown_(iConfig.getParameter<edm::InputTag>("patMETType1CorrMuEnDown")),
		patMETType1CorrTauEnUp_(iConfig.getParameter<edm::InputTag>("patMETType1CorrTauEnUp")),
		patMETType1CorrTauEnDown_(iConfig.getParameter<edm::InputTag>("patMETType1CorrTauEnDown")),
		patMETType1CorrJetResUp_(iConfig.getParameter<edm::InputTag>("patMETType1CorrJetResUp")),
		patMETType1CorrJetResDown_(iConfig.getParameter<edm::InputTag>("patMETType1CorrJetResDown")),
		patMETType1CorrJetEnUp_(iConfig.getParameter<edm::InputTag>("patMETType1CorrJetEnUp")),
		patMETType1CorrJetEnDown_(iConfig.getParameter<edm::InputTag>("patMETType1CorrJetEnDown")),
		patMETType1CorrUnclusteredUp_(iConfig.getParameter<edm::InputTag>("patMETType1CorrUnclusteredUp")),
		patMETType1CorrUnclusteredDown_(iConfig.getParameter<edm::InputTag>("patMETType1CorrUnclusteredDown")),
		patMETType1p2Corr_(iConfig.getParameter<edm::InputTag>("patMETType1p2Corr")),
		patMETType1p2CorrEleEnUp_(iConfig.getParameter<edm::InputTag>("patMETType1p2CorrEleEnUp")),
		patMETType1p2CorrEleEnDown_(iConfig.getParameter<edm::InputTag>("patMETType1p2CorrEleEndown")),
		patMETType1p2CorrMuEnUp_(iConfig.getParameter<edm::InputTag>("patMETType1p2CorrMuEnUp")),
		patMETType1p2CorrMuEnDown_(iConfig.getParameter<edm::InputTag>("patMETType1p2CorrMuEnDown")),
		patMETType1p2CorrTauEnUp_(iConfig.getParameter<edm::InputTag>("patMETType1p2CorrTauEnUp")),
		patMETType1p2CorrTauEnDown_(iConfig.getParameter<edm::InputTag>("patMETType1p2CorrTauEnDown")),
		patMETType1p2CorrJetResUp_(iConfig.getParameter<edm::InputTag>("patMETType1p2CorrJetResUp")),
		patMETType1p2CorrJetResDown_(iConfig.getParameter<edm::InputTag>("patMETType1p2CorrJetResDown")),
		patMETType1p2CorrJetEnUp_(iConfig.getParameter<edm::InputTag>("patMETType1p2CorrJetEnUp")),
		patMETType1p2CorrJetEnDown_(iConfig.getParameter<edm::InputTag>("patMETType1p2CorrJetEnDown")),
		patMETType1p2CorrUnclusteredUp_(iConfig.getParameter<edm::InputTag>("patMETType1p2CorrUnclusteredUp")),
		patMETType1p2CorrUnclusteredDown_(iConfig.getParameter<edm::InputTag>("patMETType1p2CorrUnclusteredDown")),
		pfjetsTag_(iConfig.getParameter<edm::InputTag>("pfjets")),
		genjetsTag_(iConfig.getParameter<edm::InputTag>("genjets")),
		genjetsNoNuTag_(iConfig.getParameter<edm::InputTag>("genjetsNoNu")),
		rhoIsolAllInputTag_(iConfig.getParameter<edm::InputTag>("RhoIsolAllInputTag")),
		generalTracks_(iConfig.getParameter<edm::InputTag>("generalTracks")),
		gensrc_(iConfig.getParameter<edm::InputTag>("gensrc")),
		GenEventInfo_(iConfig.getParameter<edm::InputTag>("GenEventInfo")),
		Embedded_(iConfig.getUntrackedParameter("Embedded", (bool) false)),
		//embedding
		ElectronMVATrigWeights1_(iConfig.getUntrackedParameter<std::string>("EleMVATrigWeights1")), // Electron MVA ID
		ElectronMVATrigWeights2_(iConfig.getUntrackedParameter<std::string>("EleMVATrigWeights2")), //  |
		ElectronMVATrigWeights3_(iConfig.getUntrackedParameter<std::string>("EleMVATrigWeights3")), //  |
		ElectronMVATrigWeights4_(iConfig.getUntrackedParameter<std::string>("EleMVATrigWeights4")), // \ /
		ElectronMVATrigWeights5_(iConfig.getUntrackedParameter<std::string>("EleMVATrigWeights5")), //  v
		ElectronMVATrigWeights6_(iConfig.getUntrackedParameter<std::string>("EleMVATrigWeights6")),
		ElectronMVATrigNoIPWeights1_(iConfig.getUntrackedParameter<std::string>("EleMVATrigNoIPWeights1")),
		ElectronMVATrigNoIPWeights2_(iConfig.getUntrackedParameter<std::string>("EleMVATrigNoIPWeights2")),
		ElectronMVATrigNoIPWeights3_(iConfig.getUntrackedParameter<std::string>("EleMVATrigNoIPWeights3")),
		ElectronMVATrigNoIPWeights4_(iConfig.getUntrackedParameter<std::string>("EleMVATrigNoIPWeights4")),
		ElectronMVATrigNoIPWeights5_(iConfig.getUntrackedParameter<std::string>("EleMVATrigNoIPWeights5")),
		ElectronMVATrigNoIPWeights6_(iConfig.getUntrackedParameter<std::string>("EleMVATrigNoIPWeights6")),
		ElectronMVANonTrigWeights1_(iConfig.getUntrackedParameter<std::string>("EleMVANonTrigWeights1")),
		ElectronMVANonTrigWeights2_(iConfig.getUntrackedParameter<std::string>("EleMVANonTrigWeights2")),
		ElectronMVANonTrigWeights3_(iConfig.getUntrackedParameter<std::string>("EleMVANonTrigWeights3")),
		ElectronMVANonTrigWeights4_(iConfig.getUntrackedParameter<std::string>("EleMVANonTrigWeights4")),
		ElectronMVANonTrigWeights5_(iConfig.getUntrackedParameter<std::string>("EleMVANonTrigWeights5")),
		ElectronMVANonTrigWeights6_(iConfig.getUntrackedParameter<std::string>("EleMVANonTrigWeights6")),
		JECuncData_(iConfig.getUntrackedParameter<std::string>("JECuncData")),
		JECuncMC_(iConfig.getUntrackedParameter<std::string>("JECuncMC")),
		ScaleFactor_(iConfig.getUntrackedParameter<std::string>("ScaleFactor")),
		PUInputFile_(iConfig.getUntrackedParameter<std::string>("PUInputFile")),
		PUInputHistoMC_(iConfig.getUntrackedParameter<std::string>("PUInputHistoMC")),
		PUInputHistoData_(iConfig.getUntrackedParameter<std::string>("PUInputHistoData")),
		PUInputHistoData_p5_(iConfig.getUntrackedParameter<std::string>("PUInputHistoData_p5")),
		PUInputHistoData_m5_(iConfig.getUntrackedParameter<std::string>("PUInputHistoData_m5")),
		PUInputHistoMCFineBins_(iConfig.getUntrackedParameter<std::string>("PUInputHistoMCFineBins")),
		PUInputHistoDataFineBins_(iConfig.getUntrackedParameter<std::string>("PUInputHistoDataFineBins")),
		PUOutputFile_(iConfig.getUntrackedParameter("PUOutputFile", (std::string) ("Weight3D.root"))),
		do_MCSummary_(iConfig.getUntrackedParameter("do_MCSummary", (bool) (true))),
		do_MCComplete_(iConfig.getUntrackedParameter("do_MCComplete", (bool) (true))),
		processName_(iConfig.getUntrackedParameter("TriggerProcessName", (std::string) "HLT")),
		TriggerInfoName_(iConfig.getParameter<edm::InputTag>("TriggerInfoName")),
		TriggerEvent_(iConfig.getParameter<edm::InputTag>("TriggerEvent")),
		TriggerResults_(iConfig.getParameter<edm::InputTag>("TriggerResults")),
		l1GtTriggerMenuLite_(iConfig.getParameter<edm::InputTag>("L1GtTriggerMenuLite")),
		doL1Triggers_(iConfig.getUntrackedParameter("doL1Triggers_", (bool) (false))),
		l1TriggerNames_(iConfig.getParameter<std::vector<std::string> >("l1TriggerNames")),
		useFilterModules_(iConfig.getParameter<std::vector<std::string> >("useFilterModules")),
		TriggerJetMatchingdr_(iConfig.getUntrackedParameter("TriggerJetMatchingdr", (double) 0.3)),
		TriggerMuonMatchingdr_(iConfig.getUntrackedParameter("TriggerMuonMatchingdr", (double) 0.3)),
		TriggerElectronMatchingdr_(iConfig.getUntrackedParameter("TriggerElectronMatchingdr", (double) 0.3)),
		TriggerTauMatchingdr_(iConfig.getUntrackedParameter("TriggerTauMatchingdr", (double) 0.3)),
		doBJets_(iConfig.getUntrackedParameter("doBJets", (bool) (true))),
		doPFJets_(iConfig.getUntrackedParameter("doPFJets", (bool) (true))),
		doMuons_(iConfig.getUntrackedParameter("doMuons", (bool) (true))),
		doElectrons_(iConfig.getUntrackedParameter("doElectrons", (bool) (true))),
		doPFTaus_(iConfig.getUntrackedParameter("doPFTaus", (bool) (true))),
		doTracks_(iConfig.getUntrackedParameter("doTrack", (bool) (true))),
		doTrigger_(iConfig.getUntrackedParameter("doTrigger", (bool) (true))),
		doPrimeVertex_(iConfig.getUntrackedParameter("doPrimeVertex", (bool) (true))),
		doMET_(iConfig.getUntrackedParameter("doMET", (bool) (true))),
		doMC_(iConfig.getUntrackedParameter("doMC", (bool) (true))),
		doPatJets_(iConfig.getUntrackedParameter("doPatJets", (bool) (false))),
		doPatElectrons_(iConfig.getUntrackedParameter("doPatElectrons", (bool) (false))),
		doPatMuons_(iConfig.getUntrackedParameter("doPatMuons", (bool) (false))),
		doPatMET_(iConfig.getUntrackedParameter("doPatMET", (bool) (false))),
		doMVAMET_(iConfig.getUntrackedParameter("doMVAMET", (bool) (false))),
		srcPatJets_(iConfig.getUntrackedParameter("srcPatJets", (std::string) "selectedPatJets")),
		PatJetScale_(iConfig.getUntrackedParameter("PatJetScale", (std::string) "L3Absolute")),
		BTagAlgorithm_(iConfig.getUntrackedParameter("BTagAlgorithm", (std::string) "trackCountingHighEffBJetTags")),
		jetFlavourTag_(iConfig.getParameter<edm::InputTag>("jetFlavour")),
		PUJetIdDisc_(iConfig.getParameter<edm::InputTag> ("PUJetIdDisc")),
		PUJetIdFlag_(iConfig.getParameter<edm::InputTag> ("PUJetIdFlag"))
 {
	MuonPtCut_ = iConfig.getParameter<double>("MuonPtCut"); //default: 3.0
	MuonEtaCut_ = iConfig.getParameter<double>("MuonEtaCut"); //default: 2.5
	TauPtCut_ = iConfig.getParameter<double>("TauPtCut"); //default: 18.0
	TauEtaCut_ = iConfig.getParameter<double>("TauEtaCut"); //default: 2.4
	ElectronPtCut_ = iConfig.getParameter<double>("ElectronPtCut"); //default: 8.0
	ElectronEtaCut_ = iConfig.getParameter<double>("ElectronEtaCut"); //default: 2.5
	JetPtCut_ = iConfig.getParameter<double>("JetPtCut"); //default: 18.0
	JetEtaCut_ = iConfig.getParameter<double>("JetEtaCut"); //default: 5.2

	primVtxTag_ = iConfig.getParameter<edm::InputTag>("primVtx");
	muonsTag_ = iConfig.getParameter<edm::InputTag>("muons");
	hpsTauProducer_ = iConfig.getParameter<edm::InputTag>("hpsTauProducer");
	PFElectronTag_ = iConfig.getParameter<edm::InputTag>("pfelectrons");

	system("echo 'running system to check directory structure...'");
	system("pwd");
	system("ls");
	system("ls *");
	system("ls ../*/*");

	LumiWeights_ = edm::LumiReWeighting(PUInputFile_, PUInputFile_, PUInputHistoMC_, PUInputHistoData_);
	LumiWeights3D_ = edm::Lumi3DReWeighting(PUInputFile_, PUInputFile_, PUInputHistoMC_, PUInputHistoData_, PUOutputFile_);
	LumiWeights3D_.weight3D_init(1);
	if(PUInputHistoData_p5_ != ""){
		LumiWeights_p5_ = edm::LumiReWeighting(PUInputFile_, PUInputFile_, PUInputHistoMC_, PUInputHistoData_p5_);
		LumiWeights3D_p5_ = edm::Lumi3DReWeighting(PUInputFile_, PUInputFile_, PUInputHistoMC_, PUInputHistoData_p5_, PUOutputFile_);
		LumiWeights3D_p5_.weight3D_init(1);
	}
	if(PUInputHistoData_m5_ != ""){
		LumiWeights_m5_ = edm::LumiReWeighting(PUInputFile_, PUInputFile_, PUInputHistoMC_, PUInputHistoData_m5_);
		LumiWeights3D_m5_ = edm::Lumi3DReWeighting(PUInputFile_, PUInputFile_, PUInputHistoMC_, PUInputHistoData_m5_, PUOutputFile_);
		LumiWeights3D_m5_.weight3D_init(1);
	}
	if((PUInputHistoMCFineBins_ != "") && (PUInputHistoDataFineBins_ != "")){
		LumiWeightsFineBinning_ = edm::LumiReWeighting(PUInputFile_, PUInputFile_, PUInputHistoMCFineBins_, PUInputHistoDataFineBins_);
	}

	// Electron MVA ID

	myManualCatWeightsTrigNoIP2012.clear();
	myManualCatWeightsTrigNoIP2012.push_back(ElectronMVATrigNoIPWeights1_);
	myManualCatWeightsTrigNoIP2012.push_back(ElectronMVATrigNoIPWeights2_);
	myManualCatWeightsTrigNoIP2012.push_back(ElectronMVATrigNoIPWeights3_);
	myManualCatWeightsTrigNoIP2012.push_back(ElectronMVATrigNoIPWeights4_);
	myManualCatWeightsTrigNoIP2012.push_back(ElectronMVATrigNoIPWeights5_);
	myManualCatWeightsTrigNoIP2012.push_back(ElectronMVATrigNoIPWeights6_);

	Bool_t manualCat = true;
	myMVATrigNoIP2012 = new EGammaMvaEleEstimator();
	myMVATrigNoIP2012->initialize("BDT", EGammaMvaEleEstimator::kTrigNoIP, manualCat, myManualCatWeightsTrigNoIP2012);

	myManualCatWeightsTrig2012.clear();
	myManualCatWeightsTrig2012.push_back(ElectronMVATrigWeights1_);
	myManualCatWeightsTrig2012.push_back(ElectronMVATrigWeights2_);
	myManualCatWeightsTrig2012.push_back(ElectronMVATrigWeights3_);
	myManualCatWeightsTrig2012.push_back(ElectronMVATrigWeights4_);
	myManualCatWeightsTrig2012.push_back(ElectronMVATrigWeights5_);
	myManualCatWeightsTrig2012.push_back(ElectronMVATrigWeights6_);

	myMVATrig2012 = new EGammaMvaEleEstimator();
	myMVATrig2012->initialize("BDT", EGammaMvaEleEstimator::kTrig, manualCat, myManualCatWeightsTrig2012);

	myManualCatWeightsNonTrig2012.clear();
	myManualCatWeightsNonTrig2012.push_back(ElectronMVANonTrigWeights1_);
	myManualCatWeightsNonTrig2012.push_back(ElectronMVANonTrigWeights2_);
	myManualCatWeightsNonTrig2012.push_back(ElectronMVANonTrigWeights3_);
	myManualCatWeightsNonTrig2012.push_back(ElectronMVANonTrigWeights4_);
	myManualCatWeightsNonTrig2012.push_back(ElectronMVANonTrigWeights5_);
	myManualCatWeightsNonTrig2012.push_back(ElectronMVANonTrigWeights6_);

	myMVANonTrig2012 = new EGammaMvaEleEstimator();
	myMVANonTrig2012->initialize("BDT", EGammaMvaEleEstimator::kNonTrig, manualCat, myManualCatWeightsNonTrig2012);

}

TauNtuple::~TauNtuple() {
	delete myMVATrigNoIP2012;
	delete myMVATrig2012;
	delete myMVANonTrig2012;
}

bool TauNtuple::isGoodMuon(reco::MuonRef &RefMuon) {
	if (RefMuon.isNonnull()) {
		if (RefMuon->p4().Pt() > MuonPtCut_ && fabs(RefMuon->p4().Eta()) < MuonEtaCut_)
			return true;
	}
	return false;
}

bool TauNtuple::isGoodTau(reco::PFTauRef &RefTau, edm::Handle<reco::PFTauDiscriminator> &Dis) {
	if (RefTau.isNonnull()) {
		if (RefTau->p4().Pt() > TauPtCut_ && fabs(RefTau->p4().Eta()) < TauEtaCut_ && (*Dis)[RefTau])
			return true;
	}
	return false;
}

bool TauNtuple::isGoodElectron(reco::GsfElectronRef &RefElectron) {
	if (RefElectron.isNonnull()) {
		reco::SuperClusterRef refSuperCluster = RefElectron->superCluster();
		if (RefElectron->p4().Et() > ElectronPtCut_ && fabs(refSuperCluster->eta()) < ElectronEtaCut_ && RefElectron->gsfTrack()->trackerExpectedHitsInner().numberOfHits() <= 1)
			return true;
	}
	return false;
}

bool TauNtuple::isGoodVertex(const reco::Vertex &pv) {
	if (!pv.isFake() && pv.ndof() > 4 && fabs(pv.z()) < 24 && pv.position().rho() < 2)
		return true;
	return false;
}

bool TauNtuple::isGoodJet(reco::PFJetRef &RefJet) {
	if (RefJet.isNonnull()) {
		if (RefJet->p4().Pt() > JetPtCut_ && fabs(RefJet->p4().Eta()) < JetEtaCut_)
			return true;
	}
	return false;
}

bool TauNtuple::isGoodJet(pat::JetRef &RefJet) {
	if (RefJet.isNonnull()) {
		if (RefJet->p4().Pt() > JetPtCut_ && fabs(RefJet->p4().Eta()) < JetEtaCut_)
			return true;
	}
	return false;
}

bool TauNtuple::isGoodGenJet(reco::GenJetRef &RefGenJet) {
	if (RefGenJet.isNonnull()) {
		if (RefGenJet->p4().Pt() > std::max(0.,0.8*JetPtCut_) && fabs(RefGenJet->p4().Eta()) < JetEtaCut_)
			return true;
	}
	return false;
}

// member functions
// ------------ method called to produce the data  ------------
void TauNtuple::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
	cnt_++;
	// std::cout<<"Proccess event number ======>>> "<< cnt_ <<std::endl;
	if ((cnt_ % 500) == 0)
		std::cout << "Proccess event number ======>>> " << cnt_ << std::endl;
	ClearEvent();
	using namespace edm;
	DataMCType DMT;
	DataMC_Type_idx = DMT.GetType();
	if (iEvent.isRealData() && !Embedded_) {
		DataMC_Type_idx = DataMCType::Data;
	}
	fillEventInfo(iEvent, iSetup);
	if (doMET_)
		fillMET(iEvent, iSetup);
	edm::Handle<std::vector<reco::Track> > trackCollection;
	iEvent.getByLabel(generalTracks_, trackCollection);
	if (doPrimeVertex_)
		fillPrimeVertex(iEvent, iSetup, trackCollection);
	if (doMuons_)
		fillMuons(iEvent, iSetup, trackCollection);
	if (doElectrons_)
		fillElectrons(iEvent, iSetup, trackCollection);
	if (doPFTaus_)
		fillPFTaus(iEvent, iSetup, trackCollection);
	if (doPFJets_)
		fillPFJets(iEvent, iSetup, trackCollection);
	if (doTracks_)
		fillTracks(trackCollection, iSetup);
	if (doMC_)
		fillMCTruth(iEvent, iSetup);
	if (doTrigger_)
		fillTriggerInfo(iEvent, iSetup);
	output_tree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------

void TauNtuple::fillMCTruth(edm::Event& iEvent, const edm::EventSetup& iSetup) {
	if (!iEvent.isRealData()) {
		TauDecay_CMSSW myTauDecay;
		edm::Handle<reco::GenParticleCollection> genParticles;
		iEvent.getByLabel(gensrc_, genParticles);
		myTauDecay.CheckForSignal(DataMC_Type_idx, genParticles);

		edm::Handle<GenEventInfoProduct> GenEventInfoProduct;
		iEvent.getByLabel("generator", GenEventInfoProduct);
		GenEventInfoProduct_signalProcessID = GenEventInfoProduct->signalProcessID();
		GenEventInfoProduct_weight = GenEventInfoProduct->weight();
		GenEventInfoProduct_weights = GenEventInfoProduct->weights();
		GenEventInfoProduct_qScale = GenEventInfoProduct->qScale();
		GenEventInfoProduct_alphaQCD = GenEventInfoProduct->alphaQCD();
		GenEventInfoProduct_alphaQED = GenEventInfoProduct->alphaQED();

		// info for pdf systematics
		GenEventInfoProduct_id1 = GenEventInfoProduct->pdf()->id.first;
		GenEventInfoProduct_id2 = GenEventInfoProduct->pdf()->id.second;
		GenEventInfoProduct_x1 = GenEventInfoProduct->pdf()->x.first;
		GenEventInfoProduct_x2 = GenEventInfoProduct->pdf()->x.second;
		GenEventInfoProduct_scalePDF = GenEventInfoProduct->pdf()->scalePDF;

		if (do_MCComplete_) {
			std::vector<unsigned int> index;
			for (reco::GenParticleCollection::const_iterator itr = genParticles->begin(); itr != genParticles->end(); ++itr) {
				MC_pdgid.push_back(itr->pdgId());
				MC_charge.push_back(itr->charge());
				std::vector<double> iMC_p4;
				iMC_p4.push_back(itr->p4().E());
				iMC_p4.push_back(itr->p4().Px());
				iMC_p4.push_back(itr->p4().Py());
				iMC_p4.push_back(itr->p4().Pz());

				MC_p4.push_back(iMC_p4);
				MC_midx.push_back(0);
				MC_status.push_back(itr->status());
				MC_childpdgid.push_back(std::vector<int>());
			}
			unsigned int i = 0;
			for (reco::GenParticleCollection::const_iterator itr = genParticles->begin(); itr != genParticles->end(); ++itr, i++) {
				for (unsigned int d = 0; d < itr->numberOfDaughters(); d++) {
					const reco::GenParticle *dau = static_cast<const reco::GenParticle*>(itr->daughter(d));
					MC_childpdgid.at(MC_childpdgid.size() - 1).push_back(dau->pdgId());
					unsigned int j = 0;
					for (reco::GenParticleCollection::const_iterator jtr = genParticles->begin(); jtr != genParticles->end(); ++jtr, j++) {
						if (dau->status() == jtr->status() && dau->p4() == jtr->p4() && dau->pdgId() == jtr->pdgId() && dau->numberOfMothers() == jtr->numberOfMothers()
								&& dau->numberOfDaughters() == jtr->numberOfDaughters()) {
							MC_midx.at(j) = i;
						}
					}
				}
			}
		}
		if (do_MCSummary_) {
			DataMCType DMT;
			for (reco::GenParticleCollection::const_iterator itr = genParticles->begin(); itr != genParticles->end(); ++itr) {
				if (DMT.isSignalParticle(itr->pdgId())) {
					MCSignalParticle_childpdgid.push_back(std::vector<int>());
					MCSignalParticle_pdgid.push_back(itr->pdgId());
					MCSignalParticle_charge.push_back(itr->charge());
					MCSignalParticle_Tauidx.push_back(std::vector<unsigned int>());
					std::vector<double> iSig_Poca;
					iSig_Poca.push_back(itr->vx());
					iSig_Poca.push_back(itr->vy());
					iSig_Poca.push_back(itr->vz());
					MCSignalParticle_Poca.push_back(iSig_Poca);

					std::vector<double> iSig_p4;
					iSig_p4.push_back(itr->p4().E());
					iSig_p4.push_back(itr->p4().Px());
					iSig_p4.push_back(itr->p4().Py());
					iSig_p4.push_back(itr->p4().Pz());
					MCSignalParticle_p4.push_back(iSig_p4);
					// look for daughter tau
					for (unsigned int i = 0; i < itr->numberOfDaughters(); i++) {
						const reco::Candidate *dau = itr->daughter(i);
						MCSignalParticle_childpdgid.at(MCSignalParticle_childpdgid.size() - 1).push_back(dau->pdgId());
						if (abs(dau->pdgId()) == PDGInfo::tau_minus) {
							unsigned int tauidx = MCTauandProd_p4.size();
							MCSignalParticle_Tauidx.at(MCSignalParticle_Tauidx.size() - 1).push_back(tauidx);
							// Analysis the tau decay
							unsigned int JAK_ID, TauBitMask;
							myTauDecay.AnalyzeTau(static_cast<const reco::GenParticle*>(dau), JAK_ID, TauBitMask);
							std::vector<const reco::GenParticle*> TauDecayProducts = myTauDecay.Get_TauDecayProducts();
							MCTauandProd_midx.push_back(myTauDecay.Get_MotherIdx());
							MCTau_JAK.push_back(JAK_ID);
							MCTau_DecayBitMask.push_back(TauBitMask);
							MCTauandProd_pdgid.push_back(std::vector<int>());
							MCTauandProd_charge.push_back(std::vector<int>());
							MCTauandProd_p4.push_back(std::vector<std::vector<double> >());
							MCTauandProd_Vertex.push_back(std::vector<std::vector<double> >());

							for (unsigned int i = 0; i < TauDecayProducts.size(); i++) {
								MCTauandProd_pdgid.at(tauidx).push_back(TauDecayProducts.at(i)->pdgId());
								MCTauandProd_charge.at(tauidx).push_back(TauDecayProducts.at(i)->charge());

								std::vector<double> iTauandProd_p4;
								std::vector<double> iTauandProd_vertex;
								iTauandProd_p4.push_back(TauDecayProducts.at(i)->p4().E());
								iTauandProd_p4.push_back(TauDecayProducts.at(i)->p4().Px());
								iTauandProd_p4.push_back(TauDecayProducts.at(i)->p4().Py());
								iTauandProd_p4.push_back(TauDecayProducts.at(i)->p4().Pz());

								iTauandProd_vertex.push_back(TauDecayProducts.at(i)->vx());
								iTauandProd_vertex.push_back(TauDecayProducts.at(i)->vy());
								iTauandProd_vertex.push_back(TauDecayProducts.at(i)->vz());

								MCTauandProd_p4.at(tauidx).push_back(iTauandProd_p4);
								MCTauandProd_Vertex.at(tauidx).push_back(iTauandProd_vertex);
							}
						}
					}
				}
			}
		}
	}
}

void TauNtuple::fillPrimeVertex(edm::Event& iEvent, const edm::EventSetup& iSetup, edm::Handle<std::vector<reco::Track> > &trackCollection) {
	edm::Handle<reco::VertexCollection> primVtxs;
	iEvent.getByLabel(primVtxTag_, primVtxs);

	float nVtxs = primVtxs->size();
	int ndim = 3;
	nVtxs = primVtxs->size();
	for (int i = 0; i < nVtxs; i++) {
		const reco::Vertex &pv = primVtxs->at(i);
		if (isGoodVertex(pv)) {
			Vtx_isFake.push_back(pv.isFake());
			Vtx_chi2.push_back(pv.chi2());
			Vtx_ndof.push_back(pv.ndof());
			Vtx_x.push_back(pv.x());
			Vtx_y.push_back(pv.y());
			Vtx_z.push_back(pv.z());
			std::vector<std::vector<double> > iVtx_Cov;
			for (int j = 0; j < ndim; j++) {
				iVtx_Cov.push_back(std::vector<double>());
				for (int k = 0; k <= j; k++) {
					iVtx_Cov.at(j).push_back(pv.covariance(j, k));
				}
			}
			Vtx_Cov.push_back(iVtx_Cov);
			std::vector<int> matches;
			std::vector<float> TrackWeights;
			std::vector<std::vector<double> > iVtx_TrackP4;
			for (reco::Vertex::trackRef_iterator iTrack = pv.tracks_begin(); iTrack < pv.tracks_end(); iTrack++) {
				int match(-1);
				reco::TrackRef refTrack = iTrack->castTo<reco::TrackRef>();
				if (refTrack.isNonnull()) {
					std::vector<double> iiVtx_TrackP4;
					double trkEnergy = sqrt(refTrack->px() * refTrack->px() + refTrack->py() * refTrack->py() + refTrack->pz() * refTrack->pz() + 0.13957 * 0.13957);
					iiVtx_TrackP4.push_back(trkEnergy);
					iiVtx_TrackP4.push_back(refTrack->px());
					iiVtx_TrackP4.push_back(refTrack->py());
					iiVtx_TrackP4.push_back(refTrack->pz());
					iVtx_TrackP4.push_back(iiVtx_TrackP4);
					getTrackMatch(trackCollection, refTrack, match);
					matches.push_back(match);
					TrackWeights.push_back(pv.trackWeight(refTrack));
				}
			}

			Vtx_nTrk.push_back(iVtx_TrackP4.size());
			Vtx_TracksP4.push_back(iVtx_TrackP4);
			Vtx_Track_Weights.push_back(TrackWeights);
			Vtx_Track_idx.push_back(matches);
		}
	}
}

void TauNtuple::fillMuons(edm::Event& iEvent, const edm::EventSetup& iSetup, edm::Handle<std::vector<reco::Track> > &trackCollection) {
	edm::Handle<reco::MuonCollection> muonCollection;
	iEvent.getByLabel(muonsTag_, muonCollection);
	int Muon_index = 0;
	for (reco::MuonCollection::const_iterator iMuon = muonCollection->begin(); iMuon != muonCollection->end(); ++iMuon, Muon_index++) {
		reco::MuonRef RefMuon(muonCollection, Muon_index);
		if (isGoodMuon(RefMuon)) {
			std::vector<double> iMuon_Poca;
			iMuon_Poca.push_back(RefMuon->vx());
			iMuon_Poca.push_back(RefMuon->vy());
			iMuon_Poca.push_back(RefMuon->vz());
			Muon_Poca.push_back(iMuon_Poca);
			std::vector<double> iMuon_p4;
			iMuon_p4.push_back(RefMuon->p4().E());
			iMuon_p4.push_back(RefMuon->p4().Px());
			iMuon_p4.push_back(RefMuon->p4().Py());
			iMuon_p4.push_back(RefMuon->p4().Pz());
			Muon_p4.push_back(iMuon_p4);

			const reco::MuonIsolation Iso03 = RefMuon->isolationR03();
			const reco::MuonIsolation Iso05 = RefMuon->isolationR05();

			const reco::MuonPFIsolation PFIso03 = RefMuon->pfIsolationR03();
			const reco::MuonPFIsolation PFIso04 = RefMuon->pfIsolationR04();

			Muon_numberOfChambers.push_back(RefMuon->numberOfChambers());
			Muon_isGlobalMuon.push_back(RefMuon->isGlobalMuon());
			Muon_isPFMuon.push_back(RefMuon->isPFMuon());
			Muon_isStandAloneMuon.push_back(RefMuon->isStandAloneMuon());
			Muon_isTrackerMuon.push_back(RefMuon->isTrackerMuon());
			Muon_isCaloMuon.push_back(RefMuon->isCaloMuon());
			Muon_isQualityValid.push_back(RefMuon->isQualityValid());
			Muon_isTimeValid.push_back(RefMuon->isTimeValid());
			Muon_isIsolationValid.push_back(RefMuon->isIsolationValid());
			Muon_numberOfMatchedStations.push_back(RefMuon->numberOfMatchedStations());
			Muon_numberOfMatches.push_back(RefMuon->numberOfMatches());
			Muon_charge.push_back(RefMuon->charge());

			if (RefMuon->isGlobalMuon()) {
				Muon_normChi2.push_back(RefMuon->globalTrack()->normalizedChi2());
				Muon_hitPattern_numberOfValidMuonHits.push_back(RefMuon->globalTrack()->hitPattern().numberOfValidMuonHits());
				Muon_trackerLayersWithMeasurement.push_back(RefMuon->innerTrack()->hitPattern().trackerLayersWithMeasurement());
				Muon_numberofValidPixelHits.push_back(RefMuon->innerTrack()->hitPattern().numberOfValidPixelHits());

			} else {
				Muon_normChi2.push_back(0);
				Muon_hitPattern_numberOfValidMuonHits.push_back(0);
				Muon_trackerLayersWithMeasurement.push_back(0);
				Muon_numberofValidPixelHits.push_back(0);
			}
			if (RefMuon->isTrackerMuon()) {
				Muon_innerTrack_numberofValidHits.push_back(RefMuon->innerTrack()->numberOfValidHits());
				Muon_hitPattern_pixelLayerwithMeas.push_back(RefMuon->innerTrack()->hitPattern().pixelLayersWithMeasurement());

			} else {
				Muon_innerTrack_numberofValidHits.push_back(0);
				Muon_hitPattern_pixelLayerwithMeas.push_back(0);
			}

			if (RefMuon->isIsolationValid()) {
				Muon_emEt03.push_back(Iso03.emEt);
				Muon_emVetoEt03.push_back(Iso03.emVetoEt);
				Muon_hadEt03.push_back(Iso03.hadEt);
				Muon_hadVetoEt03.push_back(Iso03.hadVetoEt);
				Muon_nJets03.push_back(Iso03.nJets);
				Muon_nTracks03.push_back(Iso03.nTracks);
				Muon_sumPt03.push_back(Iso03.sumPt);
				Muon_trackerVetoPt03.push_back(Iso03.trackerVetoPt);

				Muon_emEt05.push_back(Iso05.emEt);
				Muon_emVetoEt05.push_back(Iso05.emVetoEt);
				Muon_hadEt05.push_back(Iso05.hadEt);
				Muon_hadVetoEt05.push_back(Iso05.hadVetoEt);
				Muon_nJets05.push_back(Iso05.nJets);
				Muon_nTracks05.push_back(Iso05.nTracks);
				Muon_sumPt05.push_back(Iso05.sumPt);
				Muon_trackerVetoPt05.push_back(Iso05.trackerVetoPt);
			} else { // if isolation is not valid use -1 as default
				Muon_emEt03.push_back(-1);
				Muon_emVetoEt03.push_back(-1);
				Muon_hadEt03.push_back(-1);
				Muon_hadVetoEt03.push_back(-1);
				Muon_nJets03.push_back(-1);
				Muon_nTracks03.push_back(-1);
				Muon_sumPt03.push_back(-1);
				Muon_trackerVetoPt03.push_back(-1);

				Muon_emEt05.push_back(-1);
				Muon_emVetoEt05.push_back(-1);
				Muon_hadEt05.push_back(-1);
				Muon_hadVetoEt05.push_back(-1);
				Muon_nJets05.push_back(-1);
				Muon_nTracks05.push_back(-1);
				Muon_sumPt05.push_back(-1);
				Muon_trackerVetoPt05.push_back(-1);
			}

			//--- Fill PFMuonIsolation -----
			if (RefMuon->isPFIsolationValid()) {
				Muon_sumChargedHadronPt03.push_back(PFIso03.sumChargedHadronPt);
				Muon_sumChargedParticlePt03.push_back(PFIso03.sumChargedParticlePt);
				Muon_sumNeutralHadronEt03.push_back(PFIso03.sumNeutralHadronEt);
				Muon_sumNeutralHadronEtHighThreshold03.push_back(PFIso03.sumNeutralHadronEtHighThreshold);
				Muon_sumPhotonEt03.push_back(PFIso03.sumPhotonEt);
				Muon_sumPhotonEtHighThreshold03.push_back(PFIso03.sumPhotonEtHighThreshold);
				Muon_sumPUPt03.push_back(PFIso03.sumPUPt);

				Muon_sumChargedHadronPt04.push_back(PFIso04.sumChargedHadronPt);
				Muon_sumChargedParticlePt04.push_back(PFIso04.sumChargedParticlePt);
				Muon_sumNeutralHadronEt04.push_back(PFIso04.sumNeutralHadronEt);
				Muon_sumNeutralHadronEtHighThreshold04.push_back(PFIso04.sumNeutralHadronEtHighThreshold);
				Muon_sumPhotonEt04.push_back(PFIso04.sumPhotonEt);
				Muon_sumPhotonEtHighThreshold04.push_back(PFIso04.sumPhotonEtHighThreshold);
				Muon_sumPUPt04.push_back(PFIso04.sumPUPt);
			} else { // if isolation is not valid use -1 as default
				Muon_sumChargedHadronPt03.push_back(-1);
				Muon_sumChargedParticlePt03.push_back(-1);
				Muon_sumNeutralHadronEt03.push_back(-1);
				Muon_sumNeutralHadronEtHighThreshold03.push_back(-1);
				Muon_sumPhotonEt03.push_back(-1);
				Muon_sumPhotonEtHighThreshold03.push_back(-1);
				Muon_sumPUPt03.push_back(-1);

				Muon_sumChargedHadronPt04.push_back(-1);
				Muon_sumChargedParticlePt04.push_back(-1);
				Muon_sumNeutralHadronEt04.push_back(-1);
				Muon_sumNeutralHadronEtHighThreshold04.push_back(-1);
				Muon_sumPhotonEt04.push_back(-1);
				Muon_sumPhotonEtHighThreshold04.push_back(-1);
				Muon_sumPUPt04.push_back(-1);
			}

			reco::TrackRef Track = RefMuon->track();
			int ntp = Muon_par.size();
			Muon_par.push_back(std::vector<double>());
			Muon_cov.push_back(std::vector<double>());
			if (Track.isNonnull()) {
				GlobalPoint pvpoint(Track->vx(), Track->vy(), Track->vz());
				edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
				iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transTrackBuilder);
				reco::TransientTrack transTrk = transTrackBuilder->build(Track);
				TrackParticle trackparticle = ParticleBuilder::CreateTrackParticle(transTrk, transTrackBuilder, pvpoint, true, true);
				Muon_trackCharge.push_back(trackparticle.Charge());
				Muon_pdgid.push_back(trackparticle.PDGID());
				Muon_B.push_back(trackparticle.BField());
				Muon_M.push_back(trackparticle.Mass());
				for (int i = 0; i < trackparticle.NParameters(); i++) {
					Muon_par.at(ntp).push_back(trackparticle.Parameter(i));
					for (int j = i; j < trackparticle.NParameters(); j++) {
						Muon_cov.at(ntp).push_back(trackparticle.Covariance(i, j));
					}
				}
			} else {
				Muon_trackCharge.push_back(-999);
				Muon_pdgid.push_back(-999);
				Muon_B.push_back(-999);
				Muon_M.push_back(-999);
			}
			int match;
			getTrackMatch(trackCollection, Track, match);
			Muon_Track_idx.push_back(match);

		}
	}
}

void TauNtuple::fillTracks(edm::Handle<std::vector<reco::Track> > &trackCollection, const edm::EventSetup& iSetup) {
	for (unsigned int iTrack = 0; iTrack < trackCollection->size(); iTrack++) {
		reco::TrackRef Track(trackCollection, iTrack);
		std::vector<double> iTrack_p4;

		//assume pion mass
		double pionmass = PDGInfo::pi_mass();
		iTrack_p4.push_back(sqrt(pow(Track->p(), 2.0) + pow(pionmass, 2.0)));
		iTrack_p4.push_back(Track->px());
		iTrack_p4.push_back(Track->py());
		iTrack_p4.push_back(Track->pz());
		Track_p4.push_back(iTrack_p4);
		std::vector<double> iTrack_Poca;

		iTrack_Poca.push_back(Track->vx());
		iTrack_Poca.push_back(Track->vy());
		iTrack_Poca.push_back(Track->vz());
		Track_Poca.push_back(iTrack_Poca);
		Track_charge.push_back(Track->charge());
		Track_chi2.push_back(Track->chi2());
		Track_ndof.push_back(Track->ndof());
		Track_numberOfLostHits.push_back(Track->numberOfLostHits());
		Track_numberOfValidHits.push_back(Track->numberOfValidHits());
		Track_qualityMask.push_back(Track->qualityMask());

		GlobalPoint pvpoint(Track->vx(), Track->vy(), Track->vz());
		int ntp = Track_par.size();
		Track_par.push_back(std::vector<double>());
		Track_cov.push_back(std::vector<double>());
		edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
		iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transTrackBuilder);
		reco::TransientTrack transTrk = transTrackBuilder->build(Track);
		TrackParticle trackparticle = ParticleBuilder::CreateTrackParticle(transTrk, transTrackBuilder, pvpoint, true, true);
		Track_charge.push_back(trackparticle.Charge());
		Track_pdgid.push_back(trackparticle.PDGID());
		Track_B.push_back(trackparticle.BField());
		Track_M.push_back(trackparticle.Mass());
		for (int i = 0; i < trackparticle.NParameters(); i++) {
			Track_par.at(ntp).push_back(trackparticle.Parameter(i));
			for (int j = i; j < trackparticle.NParameters(); j++) {
				Track_cov.at(ntp).push_back(trackparticle.Covariance(i, j));
			}
		}
	}
}

void TauNtuple::fillPFTaus(edm::Event& iEvent, const edm::EventSetup& iSetup, edm::Handle<std::vector<reco::Track> > &trackCollection) {

	edm::Handle<std::vector<reco::PFTau> > HPStaus;
	iEvent.getByLabel(hpsTauProducer_, HPStaus);

	edm::Handle<std::vector<reco::Photon> > PhotonCollection;
	iEvent.getByLabel("photons", PhotonCollection);

	edm::Handle<reco::PFTauDiscriminator> HPSTightIsoDiscr;
	iEvent.getByLabel(hpsPFTauDiscriminationByTightIsolation_, HPSTightIsoDiscr);
	edm::Handle<reco::PFTauDiscriminator> HPSMediumIsoDiscr;
	iEvent.getByLabel(hpsPFTauDiscriminationByMediumIsolation_, HPSMediumIsoDiscr);
	edm::Handle<reco::PFTauDiscriminator> HPSLooseIsoDiscr;
	iEvent.getByLabel(hpsPFTauDiscriminationByLooseIsolation_, HPSLooseIsoDiscr);
	edm::Handle<reco::PFTauDiscriminator> HPSTightIsoDiscrDBSumPtCorr;
	iEvent.getByLabel(hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr_, HPSTightIsoDiscrDBSumPtCorr);
	edm::Handle<reco::PFTauDiscriminator> HPSMediumIsoDiscrDBSumPtCorr;
	iEvent.getByLabel(hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr_, HPSMediumIsoDiscrDBSumPtCorr);
	edm::Handle<reco::PFTauDiscriminator> HPSLooseIsoDiscrDBSumPtCorr;
	iEvent.getByLabel(hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr_, HPSLooseIsoDiscrDBSumPtCorr);
	edm::Handle<reco::PFTauDiscriminator> HPSVLooseIsoDiscrDBSumPtCorr;
	iEvent.getByLabel(hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr_, HPSVLooseIsoDiscrDBSumPtCorr);
	edm::Handle<reco::PFTauDiscriminator> HPSAgainstElectronsLoose;
	iEvent.getByLabel(hpsPFTauDiscriminationAgainstElectronLoose_, HPSAgainstElectronsLoose);
	edm::Handle<reco::PFTauDiscriminator> HPSAgainstElectronsMedium;
	iEvent.getByLabel(hpsPFTauDiscriminationAgainstElectronMedium_, HPSAgainstElectronsMedium);
	edm::Handle<reco::PFTauDiscriminator> HPSAgainstElectronsTight;
	iEvent.getByLabel(hpsPFTauDiscriminationAgainstElectronTight_, HPSAgainstElectronsTight);
	edm::Handle<reco::PFTauDiscriminator> HPSAgainstMuonLoose;
	iEvent.getByLabel(hpsPFTauDiscriminationAgainstMuonLoose_, HPSAgainstMuonLoose);
	edm::Handle<reco::PFTauDiscriminator> HPSAgainstMuonMedium;
	iEvent.getByLabel(hpsPFTauDiscriminationAgainstMuonMedium_, HPSAgainstMuonMedium);
	edm::Handle<reco::PFTauDiscriminator> HPSAgainstMuonTight;
	iEvent.getByLabel(hpsPFTauDiscriminationAgainstMuonTight_, HPSAgainstMuonTight);
	edm::Handle<reco::PFTauDiscriminator> HPSAgainstMuonLoose2;
	iEvent.getByLabel("hpsPFTauDiscriminationByLooseMuonRejection2", HPSAgainstMuonLoose2);
	edm::Handle<reco::PFTauDiscriminator> HPSAgainstMuonMedium2;
	iEvent.getByLabel("hpsPFTauDiscriminationByMediumMuonRejection2", HPSAgainstMuonMedium2);
	edm::Handle<reco::PFTauDiscriminator> HPSAgainstMuonTight2;
	iEvent.getByLabel("hpsPFTauDiscriminationByTightMuonRejection2", HPSAgainstMuonTight2);
	edm::Handle<reco::PFTauDiscriminator> HPSByDecayModeFinding;
	iEvent.getByLabel(hpsPFTauDiscriminationByDecayModeFinding_, HPSByDecayModeFinding);

	edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByMVA3LooseElectronRejection;
	iEvent.getByLabel("hpsPFTauDiscriminationByMVA3LooseElectronRejection", HPSPFTauDiscriminationByMVA3LooseElectronRejection);
	edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByMVA3MediumElectronRejection;
	iEvent.getByLabel("hpsPFTauDiscriminationByMVA3MediumElectronRejection", HPSPFTauDiscriminationByMVA3MediumElectronRejection);
	edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByMVA3TightElectronRejection;
	iEvent.getByLabel("hpsPFTauDiscriminationByMVA3TightElectronRejection", HPSPFTauDiscriminationByMVA3TightElectronRejection);
	edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByMVA3VTightElectronRejection;
	iEvent.getByLabel("hpsPFTauDiscriminationByMVA3VTightElectronRejection", HPSPFTauDiscriminationByMVA3VTightElectronRejection);

	edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits;
	iEvent.getByLabel("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits", HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits);
	edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits;
	iEvent.getByLabel("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits", HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits);
	edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits;
	iEvent.getByLabel("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits", HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits);
	edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits;
	iEvent.getByLabel("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits", HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits);
	edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByLooseIsolationMVA;
	iEvent.getByLabel("hpsPFTauDiscriminationByLooseIsolationMVA", HPSPFTauDiscriminationByLooseIsolationMVA);
	edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByMediumIsolationMVA;
	iEvent.getByLabel("hpsPFTauDiscriminationByMediumIsolationMVA", HPSPFTauDiscriminationByMediumIsolationMVA);
	edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByTightIsolationMVA;
	iEvent.getByLabel("hpsPFTauDiscriminationByTightIsolationMVA", HPSPFTauDiscriminationByTightIsolationMVA);
	edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByLooseIsolationMVA2;
	iEvent.getByLabel("hpsPFTauDiscriminationByLooseIsolationMVA2", HPSPFTauDiscriminationByLooseIsolationMVA2);
	edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByMediumIsolationMVA2;
	iEvent.getByLabel("hpsPFTauDiscriminationByMediumIsolationMVA2", HPSPFTauDiscriminationByMediumIsolationMVA2);
	edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByTightIsolationMVA2;
	iEvent.getByLabel("hpsPFTauDiscriminationByTightIsolationMVA2", HPSPFTauDiscriminationByTightIsolationMVA2);
	//hpsPFTauDiscriminationByLooseIsolationMVA

	edm::Handle<reco::MuonCollection> muonCollection;
	iEvent.getByLabel(muonsTag_, muonCollection);

	edm::Handle<reco::GsfElectronCollection> ElectronCollection;
	iEvent.getByLabel(PFElectronTag_, ElectronCollection);

	edm::Handle<reco::BeamSpot> beamSpot;
	iEvent.getByLabel(beamSpotTag_, beamSpot);

	for (unsigned iPFTau = 0; iPFTau < HPStaus->size(); ++iPFTau) {
		reco::PFTauRef HPStauCandidate(HPStaus, iPFTau);
		if((*HPSByDecayModeFinding)[HPStauCandidate]){
			PFTau_Photons_p4_inDR05.push_back(std::vector<std::vector<double> >());
			PFTau_photon_hasPixelSeed.push_back(std::vector<int> ());
			PFTau_photon_hadronicOverEm.push_back(std::vector<float> ());
			PFTau_photon_sigmaIetaIeta.push_back(std::vector<float> ());
			PFTau_photon_trkSumPtHollowConeDR04.push_back(std::vector<float> ());
			PFTau_photon_ecalRecHitSumEtConeDR04.push_back(std::vector<float> ());
			PFTau_photon_hcalTowerSumEtConeDR04.push_back(std::vector<float> ());

			PFTau_PionsP4.push_back(std::vector<std::vector<double> >());
			PFTau_PionsCharge.push_back(std::vector<int> ());
			std::vector<double> iPFTau_Poca;
			iPFTau_Poca.push_back(HPStauCandidate->vx());
			iPFTau_Poca.push_back(HPStauCandidate->vy());
			iPFTau_Poca.push_back(HPStauCandidate->vz());
			PFTau_Poca.push_back(iPFTau_Poca);

			std::vector<double> iPFTau_p4;
			iPFTau_p4.push_back(HPStauCandidate->p4().E());
			iPFTau_p4.push_back(HPStauCandidate->p4().Px());
			iPFTau_p4.push_back(HPStauCandidate->p4().Py());
			iPFTau_p4.push_back(HPStauCandidate->p4().Pz());

			PFTau_p4.push_back(iPFTau_p4);

			PFTau_isTightIsolation.push_back((*HPSTightIsoDiscr)[HPStauCandidate]);
			PFTau_isMediumIsolation.push_back((*HPSMediumIsoDiscr)[HPStauCandidate]);
			PFTau_isLooseIsolation.push_back((*HPSLooseIsoDiscr)[HPStauCandidate]);

			PFTau_isTightIsolationDBSumPtCorr.push_back((*HPSTightIsoDiscrDBSumPtCorr)[HPStauCandidate]);
			PFTau_isMediumIsolationDBSumPtCorr.push_back((*HPSMediumIsoDiscrDBSumPtCorr)[HPStauCandidate]);
			PFTau_isLooseIsolationDBSumPtCorr.push_back((*HPSLooseIsoDiscrDBSumPtCorr)[HPStauCandidate]);
			PFTau_isVLooseIsolationDBSumPtCorr.push_back((*HPSVLooseIsoDiscrDBSumPtCorr)[HPStauCandidate]);

			PFTau_isHPSAgainstElectronsLoose.push_back((*HPSAgainstElectronsLoose)[HPStauCandidate]);
			PFTau_isHPSAgainstElectronsMedium.push_back((*HPSAgainstElectronsMedium)[HPStauCandidate]);
			PFTau_isHPSAgainstElectronsTight.push_back((*HPSAgainstElectronsTight)[HPStauCandidate]);
			PFTau_isHPSAgainstMuonLoose.push_back((*HPSAgainstMuonLoose)[HPStauCandidate]);
			PFTau_isHPSAgainstMuonMedium.push_back((*HPSAgainstMuonMedium)[HPStauCandidate]);
			PFTau_isHPSAgainstMuonTight.push_back((*HPSAgainstMuonTight)[HPStauCandidate]);

			PFTau_isHPSAgainstMuonLoose2.push_back((*HPSAgainstMuonLoose2)[HPStauCandidate]);
			PFTau_isHPSAgainstMuonMedium2.push_back((*HPSAgainstMuonMedium2)[HPStauCandidate]);
			PFTau_isHPSAgainstMuonTight2.push_back((*HPSAgainstMuonTight2)[HPStauCandidate]);

			//    PFTau_HPSPFTauDiscriminationByMVA3rawElectronRejection.push_back((*HPSPFTauDiscriminationByMVA3rawElectronRejection)[HPStauCandidate]);
			PFTau_HPSPFTauDiscriminationByMVA3LooseElectronRejection.push_back((*HPSPFTauDiscriminationByMVA3LooseElectronRejection)[HPStauCandidate]);
			PFTau_HPSPFTauDiscriminationByMVA3MediumElectronRejection.push_back((*HPSPFTauDiscriminationByMVA3MediumElectronRejection)[HPStauCandidate]);
			PFTau_HPSPFTauDiscriminationByMVA3TightElectronRejection.push_back((*HPSPFTauDiscriminationByMVA3TightElectronRejection)[HPStauCandidate]);
			PFTau_HPSPFTauDiscriminationByMVA3VTightElectronRejection.push_back((*HPSPFTauDiscriminationByMVA3VTightElectronRejection)[HPStauCandidate]);
			//    PFTau_HPSPFTauDiscriminationByDeadECALElectronRejection.push_back((*HPSPFTauDiscriminationByDeadECALElectronRejection)[HPStauCandidate]);
			PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits.push_back((*HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits)[HPStauCandidate]);
			PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits.push_back((*HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits)[HPStauCandidate]);
			PFTau_HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits.push_back((*HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits)[HPStauCandidate]);
			PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits.push_back((*HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits)[HPStauCandidate]);
			PFTau_HPSPFTauDiscriminationByLooseIsolationMVA.push_back((*HPSPFTauDiscriminationByLooseIsolationMVA)[HPStauCandidate]);
			PFTau_HPSPFTauDiscriminationByMediumIsolationMVA.push_back((*HPSPFTauDiscriminationByMediumIsolationMVA)[HPStauCandidate]);
			PFTau_HPSPFTauDiscriminationByTightIsolationMVA.push_back((*HPSPFTauDiscriminationByTightIsolationMVA)[HPStauCandidate]);
			PFTau_HPSPFTauDiscriminationByLooseIsolationMVA2.push_back((*HPSPFTauDiscriminationByLooseIsolationMVA2)[HPStauCandidate]);
			PFTau_HPSPFTauDiscriminationByMediumIsolationMVA2.push_back((*HPSPFTauDiscriminationByMediumIsolationMVA2)[HPStauCandidate]);
			PFTau_HPSPFTauDiscriminationByTightIsolationMVA2.push_back((*HPSPFTauDiscriminationByTightIsolationMVA2)[HPStauCandidate]);

			PFTau_isHPSByDecayModeFinding.push_back((*HPSByDecayModeFinding)[HPStauCandidate]);
			PFTau_hpsDecayMode.push_back(HPStauCandidate->decayMode());
			PFTau_Charge.push_back(HPStauCandidate->charge());

			////////////////////////////////////////////////////////////////////////////////
			int Ntau = PFTau_daughterTracks.size();
			PFTau_TIP_secondaryVertex_vtxchi2.push_back(std::vector<double>());
			PFTau_TIP_secondaryVertex_vtxndof.push_back(std::vector<double>());
			PFTau_TIP_primaryVertex_vtxchi2.push_back(std::vector<double>());
			PFTau_TIP_primaryVertex_vtxndof.push_back(std::vector<double>());
			PFTau_TIP_primaryVertex_pos.push_back(std::vector<double>());
			PFTau_TIP_primaryVertex_cov.push_back(std::vector<double>());
			PFTau_TIP_secondaryVertex_pos.push_back(std::vector<double>());
			PFTau_TIP_secondaryVertex_cov.push_back(std::vector<double>());
			PFTau_a1_lvp.push_back(std::vector<double>());
			PFTau_a1_cov.push_back(std::vector<double>());

			PFTau_daughterTracks.push_back(std::vector<std::vector<double> >());
			PFTau_daughterTracks_cov.push_back(std::vector<std::vector<double> >());
			PFTau_daughterTracks_charge.push_back(std::vector<int>());
			PFTau_daughterTracks_pdgid.push_back(std::vector<int>());
			PFTau_daughterTracks_B.push_back(std::vector<double>());
			PFTau_daughterTracks_M.push_back(std::vector<double>());
			PFTau_daughterTracks_poca.push_back(std::vector<std::vector<double> >());

			PFTau_3PS_LCchi2.push_back(std::vector<float>());
			PFTau_3PS_has3ProngSolution.push_back(std::vector<int>());
			PFTau_3PS_Tau_LV.push_back(std::vector<std::vector<double> >());

			//
			PFTau_a1_charge.push_back(std::vector<int>());
			PFTau_a1_pdgid.push_back(std::vector<int>());
			PFTau_a1_B.push_back(std::vector<double>());
			PFTau_a1_M.push_back(std::vector<double>());
			// 3PS stuff is not filled at the moment, as it is dependent on HighPt branch of PFTau
			// once HighPt-PFTaus are available, have a look at commit 6dfabd41
			// https://github.com/inugent/TauDataFormat/blob/6dfabd410141b55891f4275145ddb4583e1324a2/TauNtuple/src/TauNtuple.cc
			PFTau_3PS_A1_LV.push_back(std::vector<double>());
			PFTau_3PS_M_A1.push_back(std::vector<double>());
			PFTau_3PS_M_12.push_back(std::vector<double>());
			PFTau_3PS_M_13.push_back(std::vector<double>());
			PFTau_3PS_M_23.push_back(std::vector<double>());
			PFTau_3PS_Tau_Charge.push_back(std::vector<int>());

			PFTau_TIP_flightLength.push_back(std::vector<double>());
			PFTau_TIP_flightLengthSig.push_back(std::vector<double>());

			edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
			iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transTrackBuilder);
			////////////////////////////////////////////////////////
			// Build primary vertex
			std::vector<reco::TrackBaseRef> SignalTracks;
			edm::Handle<reco::VertexCollection> primVtxs;
			iEvent.getByLabel(primVtxTag_, primVtxs);
			reco::Vertex primaryVertex = primVtxs->front();
			// Get tracks form PFTau daugthers
			const reco::PFCandidateRefVector cands = HPStauCandidate->signalPFChargedHadrCands();
			for (reco::PFCandidateRefVector::const_iterator iter = cands.begin(); iter != cands.end(); ++iter) {
				if (iter->get()->trackRef().isNonnull())
					SignalTracks.push_back(reco::TrackBaseRef(iter->get()->trackRef()));
				else if (iter->get()->gsfTrackRef().isNonnull()) {
					SignalTracks.push_back(reco::TrackBaseRef(((iter)->get()->gsfTrackRef())));
				}
			}
			// Get Muon tracks
			if (RemoveMuonTracks_) {
				if (muonCollection.isValid()) {
					for (reco::MuonCollection::size_type iMuon = 0; iMuon < muonCollection->size(); iMuon++) {
						reco::MuonRef RefMuon(muonCollection, iMuon);
						if (isGoodMuon(RefMuon) && RefMuon->track().isNonnull())
							SignalTracks.push_back(reco::TrackBaseRef(RefMuon->track()));
					}
				}
			}
			// Get Electron Tracks
			if (RemoveElectronTracks_) {
				if (ElectronCollection.isValid()) {
					for (reco::GsfElectronCollection::size_type iElectron = 0; iElectron < ElectronCollection->size(); iElectron++) {
						reco::GsfElectronRef RefElectron(ElectronCollection, iElectron);
						if (isGoodElectron(RefElectron) && RefElectron->track().isNonnull())
							SignalTracks.push_back(reco::TrackBaseRef(RefElectron->track()));
					}
				}
			}
			//--------------------  get photons found in the cone of 0.5 around tau candidate ----------------------------------------
			//			std::cout<<" ---->  "<< iPFTau<< "  photon size   "  << PhotonCollection->size() <<std::endl;
			std::vector<std::vector<double> > photonsp4;
			TLorentzVector VisibleTauLV(iPFTau_p4.at(1),iPFTau_p4.at(2),iPFTau_p4.at(3),iPFTau_p4.at(0) );
			for(unsigned int iPhot = 0; iPhot < PhotonCollection->size(); iPhot++){
			  std::vector<double>  iiGammaP4;

			  reco::PhotonRef photon(PhotonCollection,iPhot);
			  TLorentzVector PhotonLV(photon->p4().Px(), photon->p4().Py(), photon->p4().Pz(), photon->p4().E());

			  if(VisibleTauLV.DeltaR(PhotonLV) < 0.5){

				iiGammaP4.push_back(photon->p4().E());
				iiGammaP4.push_back(photon->p4().Px());
				iiGammaP4.push_back(photon->p4().Py());
				iiGammaP4.push_back(photon->p4().Pz());

				PFTau_Photons_p4_inDR05.at(PFTau_Photons_p4_inDR05.size()-1).push_back(iiGammaP4);
				PFTau_photon_hasPixelSeed.at(PFTau_photon_hasPixelSeed.size()-1).push_back(photon->hasPixelSeed());
				PFTau_photon_hadronicOverEm.at(PFTau_photon_hadronicOverEm.size()-1).push_back(photon->hadronicOverEm());
				PFTau_photon_sigmaIetaIeta.at(PFTau_photon_sigmaIetaIeta.size()-1).push_back(photon->sigmaIetaIeta());
				PFTau_photon_trkSumPtHollowConeDR04.at(PFTau_photon_trkSumPtHollowConeDR04.size()-1).push_back(photon->trkSumPtHollowConeDR04());
				PFTau_photon_ecalRecHitSumEtConeDR04.at(PFTau_photon_ecalRecHitSumEtConeDR04.size()-1).push_back(photon->ecalRecHitSumEtConeDR04());
				PFTau_photon_hcalTowerSumEtConeDR04.at(PFTau_photon_hcalTowerSumEtConeDR04.size()-1).push_back(photon->hcalTowerSumEtConeDR04());
			  }

			}

			//PiZero stuff
			std::vector<double> iPFTau_MatchedPFJetP4;
			std::vector<std::vector<double> > iPFGamma_p4;
			std::vector<std::vector<double> > iSCVariables;
			std::vector<std::vector<double> > iPhotVariables;
			std::vector<int> hasSC;
			std::vector<int> hasPhoton;
			double jetphotonEnergyFraction(0);

			edm::Handle<reco::PFJetCollection> JetCollection;
			iEvent.getByLabel(pfjetsTag_, JetCollection);
			unsigned int idx = 0;
			reco::PFJetRef MatchedPFjet = getJetIndexMatchedToGivenHPSTauCandidate(JetCollection, iPFTau_p4, idx);
			if(HPStauCandidate->decayMode()==10 || HPStauCandidate->decayMode()==1){
				if (MatchedPFjet.isNonnull()) {
					iPFTau_MatchedPFJetP4.push_back(MatchedPFjet->p4().E());
					iPFTau_MatchedPFJetP4.push_back(MatchedPFjet->p4().Px());
					iPFTau_MatchedPFJetP4.push_back(MatchedPFjet->p4().Py());
					iPFTau_MatchedPFJetP4.push_back(MatchedPFjet->p4().Pz());

					jetphotonEnergyFraction = MatchedPFjet->photonEnergyFraction();

					std::vector<reco::PFCandidatePtr> PFJetGammas = pfphotons(*MatchedPFjet);

					for (unsigned int iGamma = 0; iGamma < PFJetGammas.size(); iGamma++) {
						std::vector<double> iiPFGamma_p4;
						std::vector<double> iiSCVariables;
						std::vector<double> iiPhotVariables;

						iiPFGamma_p4.push_back(PFJetGammas.at(iGamma)->p4().E());
						iiPFGamma_p4.push_back(PFJetGammas.at(iGamma)->p4().Px());
						iiPFGamma_p4.push_back(PFJetGammas.at(iGamma)->p4().Py());
						iiPFGamma_p4.push_back(PFJetGammas.at(iGamma)->p4().Pz());

						iPFGamma_p4.push_back(iiPFGamma_p4);

						reco::PhotonRef ReferenceToPhoton = PFJetGammas.at(iGamma)->photonRef();
						reco::SuperClusterRef ReferenceToSC = PFJetGammas.at(iGamma)->superClusterRef();

		// 	   hasSC.push_back(ReferenceToSC.isNonnull());
		// 	   hasPhoton.push_back(ReferenceToPhoton.isNonnull());

						if (ReferenceToSC.isNonnull()) {
							iiSCVariables.push_back(ReferenceToSC->energy());
							iiSCVariables.push_back(ReferenceToSC->eta());
							iiSCVariables.push_back(ReferenceToSC->etaWidth());
							iiSCVariables.push_back(ReferenceToSC->phi());
							iiSCVariables.push_back(ReferenceToSC->phiWidth());
						}

						if (ReferenceToPhoton.isNonnull()) {
							iiPhotVariables.push_back(ReferenceToPhoton->caloPosition().X());
							iiPhotVariables.push_back(ReferenceToPhoton->caloPosition().Y());
							iiPhotVariables.push_back(ReferenceToPhoton->caloPosition().Z());
							iiPhotVariables.push_back(ReferenceToPhoton->p4().Px());
							iiPhotVariables.push_back(ReferenceToPhoton->p4().Py());
							iiPhotVariables.push_back(ReferenceToPhoton->p4().Pz());
						}
						iSCVariables.push_back(iiSCVariables);
						iPhotVariables.push_back(iiPhotVariables);
					}
				}
			}

			PFTau_MatchedPFJetP4.push_back(iPFTau_MatchedPFJetP4);
			PFTau_PhotonEnergyFraction.push_back(jetphotonEnergyFraction);

			PFTau_MatchedPFJetGammasP4.push_back(iPFGamma_p4);

	//        PFTau_hasSC.push_back(hasSC);
	//        PFTau_hasPhoton.push_back(hasPhoton);

			PFTau_MatchedPFJetSCVariables.push_back(iSCVariables);
			PFTau_MatchedPFJetPhotonVariables.push_back(iPhotVariables);

			edm::Handle<double> Rhokt6PFJets;
			const edm::InputTag erho("kt6PFJets", "rho");
			iEvent.getByLabel(erho, Rhokt6PFJets);
			PFTau_photon_rho.push_back(*Rhokt6PFJets);

			///////////////////////////////////////////////////////////////////////////////////////////////
			// Get Non-Tau tracks
			reco::TrackCollection nonTauTracks;
			if (trackCollection.isValid()) {
				// remove tau tracks and only tracks associated with the vertex
				unsigned int idx = 0;
				for (reco::TrackCollection::const_iterator iTrk = trackCollection->begin(); iTrk != trackCollection->end(); ++iTrk, idx++) {
					reco::TrackRef tmpRef(trackCollection, idx);
					reco::TrackRef tmpRefForBase = tmpRef;
					bool isSigTrk = false;
					bool fromVertex = false;
					for (unsigned int sigTrk = 0; sigTrk < SignalTracks.size(); sigTrk++) {
						if (reco::TrackBaseRef(tmpRefForBase) == SignalTracks.at(sigTrk)) {
							isSigTrk = true;
							break;
						}
					}
					for (std::vector<reco::TrackBaseRef>::const_iterator vtxTrkRef = primaryVertex.tracks_begin(); vtxTrkRef < primaryVertex.tracks_end(); vtxTrkRef++) {
						if (primaryVertex.trackWeight(*vtxTrkRef) > 0) {
							if ((*vtxTrkRef) == reco::TrackBaseRef(tmpRefForBase)) {
								fromVertex = true;
								break;
							}
						}
					}
					if (!isSigTrk && fromVertex)
						nonTauTracks.push_back(*iTrk);
				}
			}
			///////////////////////////////////////////////////////////////////////////////////////////////
			// Refit the vertex
			TransientVertex transVtx;
			std::vector<reco::TransientTrack> transTracks;
			for (reco::TrackCollection::iterator iter = nonTauTracks.begin(); iter != nonTauTracks.end(); ++iter) {
				transTracks.push_back(transTrackBuilder->build(*iter));
			}
			bool FitOk(true);
			AdaptiveVertexFitter avf;
			avf.setWeightThreshold(0.1); //weight per track. allow almost every fit, else --> exception
			try {
				if (!useBeamSpot_) {
					transVtx = avf.vertex(transTracks);
				} else {
					transVtx = avf.vertex(transTracks, *beamSpot);
				}
			} catch (...) {
				FitOk = false;
			}
			if (FitOk)
				primaryVertex = transVtx;

			//TVector3 pv(primaryVertex.position().x(),primaryVertex.position().y(),primaryVertex.position().z());
			PFTau_TIP_primaryVertex_pos.at(Ntau).push_back(primaryVertex.position().x());
			PFTau_TIP_primaryVertex_pos.at(Ntau).push_back(primaryVertex.position().y());
			PFTau_TIP_primaryVertex_pos.at(Ntau).push_back(primaryVertex.position().z());
			TMatrixTSym<double> pvcov(LorentzVectorParticle::NVertex);
			math::Error<LorentzVectorParticle::NVertex>::type pvCov;
			primaryVertex.fill(pvCov);
			for (int i = 0; i < LorentzVectorParticle::NVertex; i++)
				for (int j = 0; j < LorentzVectorParticle::NVertex; j++) {
					pvcov(i, j) = pvCov(i, j);
					pvcov(j, i) = pvCov(i, j);
				}
			for (int i = 0; i < LorentzVectorParticle::NVertex; i++) {
				for (int j = i; j < LorentzVectorParticle::NVertex; j++) {
					PFTau_TIP_primaryVertex_cov.at(Ntau).push_back(pvcov(i, j));
				}
			}
			double vtxchi2(0), vtxndf(1);
			vtxchi2 = primaryVertex.chi2();
			vtxndf = primaryVertex.ndof();
			PFTau_TIP_primaryVertex_vtxchi2.at(Ntau).push_back(vtxchi2);
			PFTau_TIP_primaryVertex_vtxndof.at(Ntau).push_back(vtxndf);

			///////////////////////////////////
			// if there is a secondary vertex fit it
			if (HPStauCandidate->decayMode() == 10 && (*HPSByDecayModeFinding)[HPStauCandidate]) {
				///////////////////////////////////////////////////////////////////////////////////////////////
				// Get tracks form PFTau daugthers
				std::vector<reco::TransientTrack> transTrk;
				TransientVertex transVtx;
				const reco::PFCandidateRefVector cands = HPStauCandidate->signalPFChargedHadrCands();
				for (reco::PFCandidateRefVector::const_iterator iter = cands.begin(); iter != cands.end(); ++iter) {
					if (iter->get()->trackRef().isNonnull())
						transTrk.push_back(transTrackBuilder->build(iter->get()->trackRef()));
					else if (iter->get()->gsfTrackRef().isNonnull())
						transTrk.push_back(transTrackBuilder->build(iter->get()->gsfTrackRef()));
				}
				///////////////////////////////////////////////////////////////////////////////////////////////
				// Fit the secondary vertex
				bool FitOk(true);
				KalmanVertexFitter kvf(true);
				try {
					transVtx = kvf.vertex(transTrk); //KalmanVertexFitter
				} catch (...) {
					FitOk = false;
				}
				if (!transVtx.hasRefittedTracks())
					FitOk = false;
				if (transVtx.refittedTracks().size() != transTrk.size())
					FitOk = false;
				if (FitOk) {
					reco::Vertex secondaryVertex = transVtx;
					PFTau_TIP_secondaryVertex_pos.at(Ntau).push_back(secondaryVertex.position().x());
					PFTau_TIP_secondaryVertex_pos.at(Ntau).push_back(secondaryVertex.position().y());
					PFTau_TIP_secondaryVertex_pos.at(Ntau).push_back(secondaryVertex.position().z());
					TMatrixTSym<double> svcov(LorentzVectorParticle::NVertex);
					math::Error<LorentzVectorParticle::NVertex>::type svCov;
					secondaryVertex.fill(svCov);
					for (int i = 0; i < LorentzVectorParticle::NVertex; i++)
						for (int j = 0; j < LorentzVectorParticle::NVertex; j++) {
							svcov(i, j) = svCov(i, j);
							svcov(j, i) = svCov(i, j);
						}
					for (int i = 0; i < LorentzVectorParticle::NVertex; i++) {
						for (int j = i; j < LorentzVectorParticle::NVertex; j++) {
							PFTau_TIP_secondaryVertex_cov.at(Ntau).push_back(svcov(i, j));
						}
					}
					GlobalPoint sv(secondaryVertex.position().x(), secondaryVertex.position().y(), secondaryVertex.position().z());
					vtxchi2 = 0;
					vtxndf = 1;
					vtxchi2 = secondaryVertex.chi2();
					vtxndf = secondaryVertex.ndof();
					PFTau_TIP_secondaryVertex_vtxchi2.at(Ntau).push_back(vtxchi2);
					PFTau_TIP_secondaryVertex_vtxndof.at(Ntau).push_back(vtxndf);

					PFTau_TIP_flightLength.at(Ntau).push_back(secondaryVertex.x() - primaryVertex.x());
					PFTau_TIP_flightLength.at(Ntau).push_back(secondaryVertex.y() - primaryVertex.y());
					PFTau_TIP_flightLength.at(Ntau).push_back(secondaryVertex.z() - primaryVertex.z());
					VertexDistance3D vtxdist;
					PFTau_TIP_flightLengthSig.at(Ntau).push_back(vtxdist.distance(primaryVertex, secondaryVertex).significance());

					////////////////////////////////////////////////////////////////////////////////
					LorentzVectorParticle a1;
					std::vector<reco::Track> selectedTracks = secondaryVertex.refittedTracks();
					std::vector<reco::TransientTrack> transTrkVect;
					for (unsigned int i = 0; i != selectedTracks.size(); i++)
						transTrkVect.push_back(transTrackBuilder->build(selectedTracks.at(i)));
					KinematicParticleFactoryFromTransientTrack kinFactory;
					float piMassSigma(sqrt(pow(10., -12.))), piChi(0.0), piNdf(0.0);
					std::vector<RefCountedKinematicParticle> pions;
					for (unsigned int i = 0; i < transTrkVect.size(); i++)
						pions.push_back(kinFactory.particle(transTrkVect.at(i), PDGInfo::pi_mass(), piChi, piNdf, sv, piMassSigma));
					KinematicParticleVertexFitter kpvFitter;
					RefCountedKinematicTree jpTree = kpvFitter.fit(pions);
					jpTree->movePointerToTheTop();
					const KinematicParameters parameters = jpTree->currentParticle()->currentState().kinematicParameters();
					AlgebraicSymMatrix77 cov = jpTree->currentParticle()->currentState().kinematicParametersError().matrix();
					// get pions
					double c(0);
					std::vector<reco::Track> Tracks;
					std::vector<LorentzVectorParticle> ReFitPions;
					for (unsigned int i = 0; i < transTrkVect.size(); i++) {
						std::vector<double> iPionP4;
						std::vector<double> iPionCharge;
						c += transTrkVect.at(i).charge();
						ReFitPions.push_back(ParticleBuilder::CreateLorentzVectorParticle(transTrkVect.at(i), transTrackBuilder, secondaryVertex, true, true));
						iPionP4.push_back(ReFitPions.at(i).LV().E());
						iPionP4.push_back(ReFitPions.at(i).LV().Px());
						iPionP4.push_back(ReFitPions.at(i).LV().Py());
						iPionP4.push_back(ReFitPions.at(i).LV().Pz());

						PFTau_PionsP4.at(PFTau_PionsP4.size()-1).push_back(iPionP4);
						PFTau_PionsCharge.at(PFTau_PionsCharge.size()-1).push_back(transTrkVect.at(i).charge());

					}
					// now covert a1 into LorentzVectorParticle
					TMatrixT<double> a1_par(LorentzVectorParticle::NLorentzandVertexPar, 1);
					TMatrixTSym<double> a1_cov(LorentzVectorParticle::NLorentzandVertexPar);
					for (int i = 0; i < LorentzVectorParticle::NLorentzandVertexPar; i++) {
						a1_par(i, 0) = parameters(i);
						for (int j = 0; j < LorentzVectorParticle::NLorentzandVertexPar; j++) {
							a1_cov(i, j) = cov(i, j);
						}
					}
					a1 = LorentzVectorParticle(a1_par, a1_cov, abs(PDGInfo::a_1_plus) * c, c, transTrackBuilder->field()->inInverseGeV(sv).z());
					PFTau_a1_charge.at(Ntau).push_back(a1.Charge());
					PFTau_a1_pdgid.at(Ntau).push_back(a1.PDGID());
					PFTau_a1_B.at(Ntau).push_back(a1.BField());
					PFTau_a1_M.at(Ntau).push_back(a1.Mass());
					for (int i = 0; i < a1.NParameters(); i++) {
						PFTau_a1_lvp.at(Ntau).push_back(a1.Parameter(i));
						for (int j = i; j < a1.NParameters(); j++) {
							PFTau_a1_cov.at(Ntau).push_back(a1.Covariance(i, j));
						}
					}
				}
			}
			////////////////////////////////////////////////////////////////////////////////
			// Get unfit Tracks
			GlobalPoint pvpoint(primaryVertex.position().x(), primaryVertex.position().y(), primaryVertex.position().z());
			for (reco::PFCandidateRefVector::const_iterator iter = cands.begin(); iter != cands.end(); ++iter) {
				int Npi = PFTau_daughterTracks.at(Ntau).size();
				PFTau_daughterTracks_poca.at(Ntau).push_back(std::vector<double>());
				PFTau_daughterTracks.at(Ntau).push_back(std::vector<double>());
				PFTau_daughterTracks_cov.at(Ntau).push_back(std::vector<double>());
				//
				bool hastrack(false);
				reco::TransientTrack transTrk;
				if (iter->get()->trackRef().isNonnull()) {
					transTrk = transTrackBuilder->build(iter->get()->trackRef());
					hastrack = true;
				}
				//else if(iter->get()->gsfTrackRef().isNonnull()){transTrk=transTrackBuilder->build(iter->get()->gsfTrackRef());hastrack=true;}
				if (hastrack) {
					TrackParticle pion = ParticleBuilder::CreateTrackParticle(transTrk, transTrackBuilder, pvpoint, true, true);
					GlobalPoint pos = transTrk.trajectoryStateClosestToPoint(pvpoint).position();
					PFTau_daughterTracks_poca.at(Ntau).at(Npi).push_back(pos.x());
					PFTau_daughterTracks_poca.at(Ntau).at(Npi).push_back(pos.y());
					PFTau_daughterTracks_poca.at(Ntau).at(Npi).push_back(pos.z());
					PFTau_daughterTracks_charge.at(Ntau).push_back(pion.Charge());
					PFTau_daughterTracks_pdgid.at(Ntau).push_back(pion.PDGID());
					PFTau_daughterTracks_B.at(Ntau).push_back(pion.BField());
					PFTau_daughterTracks_M.at(Ntau).push_back(pion.Mass());
					for (int i = 0; i < pion.NParameters(); i++) {
						PFTau_daughterTracks.at(Ntau).at(Npi).push_back(pion.Parameter(i));
						for (int j = i; j < pion.NParameters(); j++) {
							PFTau_daughterTracks_cov.at(Ntau).at(Npi).push_back(pion.Covariance(i, j));
						}
					}
				}
			}
			////////////////////////////////////////////////////////////////////////////////
			// Get Pi0, Gamma's and Track's
			reco::PFCandidateRefVector GammaCandidate = HPStauCandidate->signalPFGammaCands();
			reco::PFCandidateRefVector ChargedHadrCand = HPStauCandidate->signalPFChargedHadrCands();
			const std::vector<reco::RecoTauPiZero> PiZeroCandiate = HPStauCandidate->signalPiZeroCandidates();

			std::vector<std::vector<double> > iPFTau_PiZeroP4;
			std::vector<int> iPFTau_PiZeroNumOfPhotons;
			std::vector<int> iPFTau_PiZeroNumOfElectrons;

			if (PiZeroCandiate.size() != 0) {
				for (unsigned int Pi0Index = 0; Pi0Index < PiZeroCandiate.size(); Pi0Index++) {
					reco::RecoTauPiZero iPi0 = PiZeroCandiate.at(Pi0Index);
					std::vector<double> iiPFTau_PiZeroP4;

					iiPFTau_PiZeroP4.push_back(iPi0.p4().E());
					iiPFTau_PiZeroP4.push_back(iPi0.p4().Px());
					iiPFTau_PiZeroP4.push_back(iPi0.p4().Py());
					iiPFTau_PiZeroP4.push_back(iPi0.p4().Pz());
					iPFTau_PiZeroP4.push_back(iiPFTau_PiZeroP4);

					iPFTau_PiZeroNumOfPhotons.push_back(iPi0.numberOfGammas());
					iPFTau_PiZeroNumOfElectrons.push_back(iPi0.numberOfElectrons());

				}
			}
			PFTau_PiZeroP4.push_back(iPFTau_PiZeroP4);
			PFTau_PiZeroNumOfPhotons.push_back(iPFTau_PiZeroNumOfPhotons);
			PFTau_PiZeroNumOfElectrons.push_back(iPFTau_PiZeroNumOfElectrons);

			std::vector<std::vector<double> > iPFTau_GammaP4;
			if (GammaCandidate.size() != 0) {
				for (unsigned int iGamma = 0; iGamma < GammaCandidate.size(); iGamma++) {
					reco::PFCandidateRef GammaCand(GammaCandidate, iGamma);
					std::vector<double> iiPFTau_GammaP4;

					iiPFTau_GammaP4.push_back(GammaCand->p4().E());
					iiPFTau_GammaP4.push_back(GammaCand->p4().Px());
					iiPFTau_GammaP4.push_back(GammaCand->p4().Py());
					iiPFTau_GammaP4.push_back(GammaCand->p4().Pz());
					iPFTau_GammaP4.push_back(iiPFTau_GammaP4);
				}
			}
			PFTau_GammaP4.push_back(iPFTau_GammaP4);

			std::vector<std::vector<double> > iPFTau_ChargedHadronP4;
			std::vector<std::vector<int> > iPFTau_ChargedHadronsCharge;
			if (ChargedHadrCand.size() != 0) {
				for (unsigned int iChargedHadron = 0; iChargedHadron < ChargedHadrCand.size(); iChargedHadron++) {
					reco::PFCandidateRef ChargeHadronCand(ChargedHadrCand, iChargedHadron);
					std::vector<double> iiPFTau_ChargedHadronP4;
					std::vector<int> iiPFTau_ChargedHadronsCharge;
					if (ChargedHadrCand.at(iChargedHadron)->trackRef().isNonnull()) {
						iiPFTau_ChargedHadronP4.push_back(sqrt(pow(ChargedHadrCand.at(iChargedHadron)->trackRef()->p(), 2.0) + pow(0.13957018, 2.0)));
						iiPFTau_ChargedHadronP4.push_back(ChargedHadrCand.at(iChargedHadron)->trackRef()->px());
						iiPFTau_ChargedHadronP4.push_back(ChargedHadrCand.at(iChargedHadron)->trackRef()->py());
						iiPFTau_ChargedHadronP4.push_back(ChargedHadrCand.at(iChargedHadron)->trackRef()->pz());
						iiPFTau_ChargedHadronsCharge.push_back(ChargeHadronCand->charge());

						iPFTau_ChargedHadronP4.push_back(iiPFTau_ChargedHadronP4);
						iPFTau_ChargedHadronsCharge.push_back(iiPFTau_ChargedHadronsCharge);
					}
				}
			}
			PFTau_ChargedHadronsP4.push_back(iPFTau_ChargedHadronP4);
			PFTau_ChargedHadronsCharge.push_back(iPFTau_ChargedHadronsCharge);

			////////////////////////////////////////////////////////////////////////////////
			//const reco::PFCandidateRefVector  ChargedHadrCand=HPStauCandidate->signalPFChargedHadrCands();
			std::vector<int> matches;
			for (unsigned int i = 0; i < ChargedHadrCand.size(); i++) {
				reco::TrackRef refTrack = ChargedHadrCand.at(i).get()->trackRef();
				if (refTrack.isNonnull()) {
					int match(-1);
					getTrackMatch(trackCollection, refTrack, match);
					matches.push_back(match);
				}
			}
			PFTau_Track_idx.push_back(matches);
		}
	}
}

void TauNtuple::fillPFJets(edm::Event& iEvent, const edm::EventSetup& iSetup, edm::Handle<std::vector<reco::Track> > &trackCollection) {
	if (!doPatJets_) {
		edm::Handle<reco::PFJetCollection> JetCollection;
		iEvent.getByLabel(pfjetsTag_, JetCollection);

		edm::Handle<edm::ValueMap<float> > puJetIdMva;
		iEvent.getByLabel(PUJetIdDisc_, puJetIdMva);

		edm::Handle<edm::ValueMap<int> > puJetIdFlag;
		iEvent.getByLabel(PUJetIdFlag_, puJetIdFlag);

		JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty();
		if(iEvent.isRealData() || Embedded_){
			jecUnc->setParameters(JECuncData_);
		}else{
			jecUnc->setParameters(JECuncMC_);
		}

		for (reco::PFJetCollection::size_type iPFJet = 0; iPFJet < JetCollection->size(); iPFJet++) {
			reco::PFJetRef PFJet(JetCollection, iPFJet);
			if (isGoodJet(PFJet)) {
				std::vector<double> iPFJet_p4;
				iPFJet_p4.push_back(PFJet->p4().E());
				iPFJet_p4.push_back(PFJet->p4().Px());
				iPFJet_p4.push_back(PFJet->p4().Py());
				iPFJet_p4.push_back(PFJet->p4().Pz());
				PFJet_p4.push_back(iPFJet_p4);
				PFJet_numberOfDaughters.push_back(PFJet->numberOfDaughters());
				PFJet_chargedEmEnergyFraction.push_back(PFJet->chargedEmEnergyFraction());
				PFJet_chargedHadronEnergyFraction.push_back(PFJet->chargedHadronEnergyFraction());
				PFJet_neutralHadronEnergyFraction.push_back(PFJet->neutralHadronEnergyFraction());
				PFJet_neutralEmEnergyFraction.push_back(PFJet->neutralEmEnergyFraction());
				PFJet_chargedEmEnergy.push_back(PFJet->chargedEmEnergy());
				PFJet_chargedHadronEnergy.push_back(PFJet->chargedHadronEnergy());
				PFJet_chargedHadronMultiplicity.push_back(PFJet->chargedHadronMultiplicity());
				PFJet_chargedMuEnergy.push_back(PFJet->chargedMuEnergy());
				PFJet_chargedMultiplicity.push_back(PFJet->chargedMultiplicity());
				PFJet_electronEnergy.push_back(PFJet->electronEnergy());
				PFJet_electronMultiplicity.push_back(PFJet->electronMultiplicity());
				PFJet_HFEMEnergy.push_back(PFJet->HFEMEnergy());
				PFJet_HFEMMultiplicity.push_back(PFJet->HFEMMultiplicity());
				PFJet_HFHadronEnergy.push_back(PFJet->HFHadronEnergy());
				PFJet_HFHadronMultiplicity.push_back(PFJet->HFHadronMultiplicity());
				PFJet_muonEnergy.push_back(PFJet->muonEnergy());
				PFJet_muonMultiplicity.push_back(PFJet->muonMultiplicity());
				PFJet_neutralEmEnergy.push_back(PFJet->neutralEmEnergy());
				PFJet_neutralHadronEnergy.push_back(PFJet->neutralHadronEnergy());
				PFJet_neutralHadronMultiplicity.push_back(PFJet->neutralHadronMultiplicity());
				PFJet_photonEnergy.push_back(PFJet->photonEnergy());
				PFJet_photonMultiplicity.push_back(PFJet->photonMultiplicity());
				PFJet_jetArea.push_back(PFJet->jetArea());
				PFJet_maxDistance.push_back(PFJet->maxDistance());
				PFJet_nConstituents.push_back(PFJet->nConstituents());
				PFJet_pileup.push_back(PFJet->pileup());
				PFJet_etaetaMoment.push_back(PFJet->etaetaMoment());
				PFJet_etaphiMoment.push_back(PFJet->etaphiMoment());
				std::vector<int> matches;
				const edm::ProductID &TrID = trackCollection.id();
				std::vector<std::vector<double> > iPFJet_TrackP4;
				std::vector< std::pair<int, double> > iPFJet_TrackPt;
				unsigned int countTracks = 0;
				for (unsigned i = 0; i < PFJet->numberOfDaughters(); i++) {
					const reco::PFCandidatePtr pfcand = PFJet->getPFConstituent(i);
					reco::TrackRef trackref = pfcand->trackRef();
					if (trackref.isNonnull()) {
						std::vector<double> iiPFJet_TrackP4;
						float trkEnergy = sqrt(trackref->px() * trackref->px() + trackref->py() * trackref->py() + trackref->pz() * trackref->pz() + 0.13957 * 0.13957);
						iiPFJet_TrackP4.push_back(trkEnergy);
						iiPFJet_TrackP4.push_back(trackref->px());
						iiPFJet_TrackP4.push_back(trackref->py());
						iiPFJet_TrackP4.push_back(trackref->pz());
						iPFJet_TrackP4.push_back(iiPFJet_TrackP4);
						iPFJet_TrackPt.push_back( std::make_pair(countTracks,trackref->pt()) );
						if (trackref.id() != TrID)
							continue;
						int match(-1);
						getTrackMatch(trackCollection, trackref, match);
						if (match >= 0) matches.push_back(match);
						countTracks++;
					}
				}
				PFJet_Track_idx.push_back(matches);
				std::sort(iPFJet_TrackPt.begin(), iPFJet_TrackPt.end(), sortIdxByValue());
				std::vector<std::vector<double> > iPFJet_TrackP4_two_leading;
				if(iPFJet_TrackP4.size()>0){
					iPFJet_TrackP4_two_leading.push_back(iPFJet_TrackP4.at(iPFJet_TrackPt.at(0).first));
					if(iPFJet_TrackP4.size()>1){
						iPFJet_TrackP4_two_leading.push_back(iPFJet_TrackP4.at(iPFJet_TrackPt.at(1).first));
					}
				}
				PFJet_TracksP4.push_back(iPFJet_TrackP4_two_leading);
				PFJet_nTrk.push_back(iPFJet_TrackP4.size());

				edm::Handle<std::vector<reco::PFTau> > HPStaus;
				iEvent.getByLabel(hpsTauProducer_, HPStaus);
				int idx = -1;
				reco::PFTauRef MatchedHPSTau = getHPSTauMatchedToJet(HPStaus, iPFJet_p4, idx);
				PFJet_MatchedHPS_idx.push_back(idx);

				/////// PU Jet ID
				float puJetID_discr = (*puJetIdMva)[PFJet];
				int puJetID_idflag = (*puJetIdFlag)[PFJet];

				PFJet_PUJetID_discr.push_back(puJetID_discr);
				PFJet_PUJetID_looseWP.push_back(PileupJetIdentifier::passJetId(puJetID_idflag, PileupJetIdentifier::kLoose));
				PFJet_PUJetID_mediumWP.push_back(PileupJetIdentifier::passJetId(puJetID_idflag, PileupJetIdentifier::kMedium));
				PFJet_PUJetID_tightWP.push_back(PileupJetIdentifier::passJetId(puJetID_idflag, PileupJetIdentifier::kTight));

				if (doBJets_) {
					edm::Handle<reco::JetFloatAssociation::Container> jetDiscriminator;
					edm::InputTag BTagAlgorithmTag = edm::InputTag(BTagAlgorithm_);
					iEvent.getByLabel(BTagAlgorithmTag, jetDiscriminator);
					double bTagValue = reco::JetFloatAssociation::getValue(*jetDiscriminator, *PFJet);
					PFJet_bDiscriminator.push_back(bTagValue);
				}else{
					PFJet_bDiscriminator.push_back(-1);
				}
				//jet flavour (needed for weights)
				if (!iEvent.isRealData()) {
						edm::Handle<reco::JetFlavourMatchingCollection> jetFlavMatch;
						iEvent.getByLabel(jetFlavourTag_, jetFlavMatch);
					PFJet_partonFlavour.push_back((*jetFlavMatch)[reco::JetBaseRef(PFJet)].getFlavour());
				}else{
					PFJet_partonFlavour.push_back(-1);
				}

				jecUnc->setJetEta(PFJet->eta());
				jecUnc->setJetPt(PFJet->pt());
				PFJet_JECuncertainty.push_back(jecUnc->getUncertainty(true));
			}
		}
		if(!iEvent.isRealData() && !Embedded_){

			edm::Handle<reco::GenJetCollection> GenJetCollection;
			iEvent.getByLabel(genjetsTag_, GenJetCollection);

			edm::Handle<reco::GenJetCollection> GenJetNoNuCollection;
			iEvent.getByLabel(genjetsNoNuTag_, GenJetNoNuCollection);

			for (reco::GenJetCollection::size_type iGenJet = 0; iGenJet < GenJetCollection->size(); iGenJet++) {
				reco::GenJetRef RefGenJet(GenJetCollection, iGenJet);
				if(isGoodGenJet(RefGenJet)){
					std::vector<double> iGenJet_p4;
					iGenJet_p4.push_back(RefGenJet->p4().E());
					iGenJet_p4.push_back(RefGenJet->p4().Px());
					iGenJet_p4.push_back(RefGenJet->p4().Py());
					iGenJet_p4.push_back(RefGenJet->p4().Pz());
					PFJet_GenJet_p4.push_back(iGenJet_p4);
					PFJet_GenJet_Constituents_p4.push_back(std::vector<std::vector<double> >());
					for(reco::GenParticleCollection::size_type iGenConst = 0; iGenConst < RefGenJet->getGenConstituents().size(); iGenConst++) {
						std::vector<double> iiGenConst_p4;
						iiGenConst_p4.push_back(RefGenJet->getGenConstituent(iGenConst)->p4().E());
						iiGenConst_p4.push_back(RefGenJet->getGenConstituent(iGenConst)->p4().Px());
						iiGenConst_p4.push_back(RefGenJet->getGenConstituent(iGenConst)->p4().Py());
						iiGenConst_p4.push_back(RefGenJet->getGenConstituent(iGenConst)->p4().Pz());
						PFJet_GenJet_Constituents_p4.at(PFJet_GenJet_Constituents_p4.size()-1).push_back(iiGenConst_p4);
					}
				}
			}
			for (reco::GenJetCollection::size_type iGenJetNoNu = 0; iGenJetNoNu < GenJetNoNuCollection->size(); iGenJetNoNu++) {
				reco::GenJetRef RefGenJetNoNu(GenJetNoNuCollection, iGenJetNoNu);
				if(isGoodGenJet(RefGenJetNoNu)){
					std::vector<double> iGenJetNoNu_p4;
					iGenJetNoNu_p4.push_back(RefGenJetNoNu->p4().E());
					iGenJetNoNu_p4.push_back(RefGenJetNoNu->p4().Px());
					iGenJetNoNu_p4.push_back(RefGenJetNoNu->p4().Py());
					iGenJetNoNu_p4.push_back(RefGenJetNoNu->p4().Pz());
					PFJet_GenJetNoNu_p4.push_back(iGenJetNoNu_p4);
					PFJet_GenJetNoNu_Constituents_p4.push_back(std::vector<std::vector<double> >());
					for(reco::GenParticleCollection::size_type iGenConstNoNu = 0; iGenConstNoNu < RefGenJetNoNu->getGenConstituents().size(); iGenConstNoNu++) {
						std::vector<double> iiGenConstNoNu_p4;
						iiGenConstNoNu_p4.push_back(RefGenJetNoNu->getGenConstituent(iGenConstNoNu)->p4().E());
						iiGenConstNoNu_p4.push_back(RefGenJetNoNu->getGenConstituent(iGenConstNoNu)->p4().Px());
						iiGenConstNoNu_p4.push_back(RefGenJetNoNu->getGenConstituent(iGenConstNoNu)->p4().Py());
						iiGenConstNoNu_p4.push_back(RefGenJetNoNu->getGenConstituent(iGenConstNoNu)->p4().Pz());
						PFJet_GenJetNoNu_Constituents_p4.at(PFJet_GenJetNoNu_Constituents_p4.size()-1).push_back(iiGenConstNoNu_p4);
					}
				}
			}
		}else{
			PFJet_GenJet_p4.push_back(std::vector<double>());
			PFJet_GenJet_Constituents_p4.push_back(std::vector<std::vector<double> >());
			PFJet_GenJetNoNu_p4.push_back(std::vector<double>());
			PFJet_GenJetNoNu_Constituents_p4.push_back(std::vector<std::vector<double> >());
		}
		delete jecUnc;
	} else {
		edm::Handle<pat::JetCollection> jets;
		edm::InputTag labelJets(srcPatJets_);
		iEvent.getByLabel(labelJets, jets);
		//edm::Handle<edm::View<pat::Jet> > jetHandle;
		//iEvent.getByLabel("PatJets", jetHandle);
		//edm::View<pat::Jet> PatJet = *jetHandle;

		edm::Handle<std::vector<reco::PFTau> > HPStaus;
		iEvent.getByLabel(hpsTauProducer_, HPStaus);

		edm::Handle<edm::ValueMap<float> > puJetIdMva;
		iEvent.getByLabel(PUJetIdDisc_, puJetIdMva);

		edm::Handle<edm::ValueMap<int> > puJetIdFlag;
		iEvent.getByLabel(PUJetIdFlag_, puJetIdFlag);

		for (pat::JetCollection::size_type iPatJet = 0; iPatJet < jets->size(); iPatJet++) {
			//for(pat::JetCollection::size_type iPatJet = 0; iPatJet < PatJet.size(); iPatJet++){
			pat::JetRef PatJet(jets, iPatJet);
			if (isGoodJet(PatJet)) {
				std::vector<double> iPatJet_p4;
				iPatJet_p4.push_back(PatJet->correctedP4(PatJetScale_).E());
				iPatJet_p4.push_back(PatJet->correctedP4(PatJetScale_).Px());
				iPatJet_p4.push_back(PatJet->correctedP4(PatJetScale_).Py());
				iPatJet_p4.push_back(PatJet->correctedP4(PatJetScale_).Pz());
				PFJet_p4.push_back(iPatJet_p4);
				PFJet_numberOfDaughters.push_back(PatJet->numberOfDaughters());
				PFJet_chargedEmEnergyFraction.push_back(PatJet->chargedEmEnergyFraction());
				PFJet_chargedHadronEnergyFraction.push_back(PatJet->chargedHadronEnergyFraction());
				PFJet_neutralHadronEnergyFraction.push_back(PatJet->neutralHadronEnergyFraction());
				PFJet_neutralEmEnergyFraction.push_back(PatJet->neutralEmEnergyFraction());
				PFJet_chargedEmEnergy.push_back(PatJet->chargedEmEnergy());
				PFJet_chargedHadronEnergy.push_back(PatJet->chargedHadronEnergy());
				PFJet_chargedHadronMultiplicity.push_back(PatJet->chargedHadronMultiplicity());
				PFJet_chargedMuEnergy.push_back(PatJet->chargedMuEnergy());
				PFJet_chargedMultiplicity.push_back(PatJet->chargedMultiplicity());
				PFJet_electronEnergy.push_back(PatJet->electronEnergy());
				PFJet_electronMultiplicity.push_back(PatJet->electronMultiplicity());
				PFJet_HFEMEnergy.push_back(PatJet->HFEMEnergy());
				PFJet_HFEMMultiplicity.push_back(PatJet->HFEMMultiplicity());
				PFJet_HFHadronEnergy.push_back(PatJet->HFHadronEnergy());
				PFJet_HFHadronMultiplicity.push_back(PatJet->HFHadronMultiplicity());
				PFJet_muonEnergy.push_back(PatJet->muonEnergy());
				PFJet_muonMultiplicity.push_back(PatJet->muonMultiplicity());
				PFJet_neutralEmEnergy.push_back(PatJet->neutralEmEnergy());
				PFJet_neutralHadronEnergy.push_back(PatJet->neutralHadronEnergy());
				PFJet_neutralHadronMultiplicity.push_back(PatJet->neutralHadronMultiplicity());
				PFJet_photonEnergy.push_back(PatJet->photonEnergy());
				PFJet_photonMultiplicity.push_back(PatJet->photonMultiplicity());
				PFJet_jetArea.push_back(PatJet->jetArea());
				PFJet_maxDistance.push_back(PatJet->maxDistance());
				PFJet_nConstituents.push_back(PatJet->nConstituents());
				PFJet_pileup.push_back(PatJet->pileup());
				PFJet_etaetaMoment.push_back(PatJet->etaetaMoment());
				PFJet_etaphiMoment.push_back(PatJet->etaphiMoment());
				std::vector<int> matches;
				const edm::ProductID &TrID = trackCollection.id();
				std::vector<std::vector<double> > iPFJet_TrackP4;
				std::vector< std::pair<int, double> > iPFJet_TrackPt;
				unsigned int countTracks = 0;
				for (unsigned i = 0; i < PatJet->numberOfDaughters(); i++) {
					const reco::PFCandidatePtr pfcand = PatJet->getPFConstituent(i);
					reco::TrackRef trackref = pfcand->trackRef();
					if (trackref.isNonnull()) {
						std::vector<double> iiPFJet_TrackP4;
						float trkEnergy = sqrt(trackref->px() * trackref->px() + trackref->py() * trackref->py() + trackref->pz() * trackref->pz() + 0.13957 * 0.13957);
						iiPFJet_TrackP4.push_back(trkEnergy);
						iiPFJet_TrackP4.push_back(trackref->px());
						iiPFJet_TrackP4.push_back(trackref->py());
						iiPFJet_TrackP4.push_back(trackref->pz());
						iPFJet_TrackP4.push_back(iiPFJet_TrackP4);
						iPFJet_TrackPt.push_back( std::make_pair(countTracks,trackref->pt()) );
						if (trackref.id() != TrID)
							continue;
						int match(-1);
						getTrackMatch(trackCollection, trackref, match);
						if (match >= 0)
							matches.push_back(match);
					}
				}
				PFJet_Track_idx.push_back(matches);
				std::sort(iPFJet_TrackPt.begin(), iPFJet_TrackPt.end(), sortIdxByValue());
				std::vector<std::vector<double> > iPFJet_TrackP4_two_leading;
				if(iPFJet_TrackP4.size()>0){
					iPFJet_TrackP4_two_leading.push_back(iPFJet_TrackP4.at(iPFJet_TrackPt.at(0).first));
					if(iPFJet_TrackP4.size()>1){
						iPFJet_TrackP4_two_leading.push_back(iPFJet_TrackP4.at(iPFJet_TrackPt.at(1).first));
					}
				}
				PFJet_TracksP4.push_back(iPFJet_TrackP4_two_leading);
				PFJet_nTrk.push_back(iPFJet_TrackP4.size());

				int idx = -1;
				reco::PFTauRef MatchedHPSTau = getHPSTauMatchedToJet(HPStaus, iPatJet_p4, idx);
				PFJet_MatchedHPS_idx.push_back(idx);

				/////// PU Jet ID
				float puJetID_discr = (*puJetIdMva)[PatJet];
				int puJetID_idflag = (*puJetIdFlag)[PatJet];

				PFJet_PUJetID_discr.push_back(puJetID_discr);
				PFJet_PUJetID_looseWP.push_back(PileupJetIdentifier::passJetId(puJetID_idflag, PileupJetIdentifier::kLoose));
				PFJet_PUJetID_mediumWP.push_back(PileupJetIdentifier::passJetId(puJetID_idflag, PileupJetIdentifier::kMedium));
				PFJet_PUJetID_tightWP.push_back(PileupJetIdentifier::passJetId(puJetID_idflag, PileupJetIdentifier::kTight));

				///////////////////////////////////////////////
				//
				// B-Tagging
				//
				PFJet_partonFlavour.push_back(PatJet->partonFlavour());
				PFJet_bDiscriminator.push_back(PatJet->bDiscriminator(BTagAlgorithm_));
				//std::vector<float> BTagWeights(0);
				//PFJet_BTagWeight.push_back(BTagWeights);

				//std::cout << "!!!!!!!! b-tagging !!!!!!!!" << std::endl;
				//std::cout << "b-tag algorithm name: " << PatJet->getPairDiscri().first << ", value: " << PatJet->getPairDiscri().second << std::endl;
				//PFJet_bTagAlgorithmName.push_back(PatJet->getPairDiscri().first);
				//PFJet_bTagAlgorithmValue.push_back(PatJet->getPairDiscri().second);
				//std::cout << "name size: " << PFJet_bTagAlgorithmName.size() << ", value size: " << PFJet_bTagAlgorithmValue.size() << std::endl;

			}
		}
	}
}

void TauNtuple::fillElectrons(edm::Event& iEvent, const edm::EventSetup& iSetup, edm::Handle<std::vector<reco::Track> > &trackCollection) {

	edm::Handle<reco::GsfElectronCollection> ElectronCollection;
	iEvent.getByLabel(PFElectronTag_, ElectronCollection);

	edm::Handle<edm::ValueMap<double> > ElectronRegEnergy;
	iEvent.getByLabel("calibratedElectrons","eneRegForGsfEle",ElectronRegEnergy);

	edm::Handle<edm::ValueMap<double> > ElectronRegEnergyError;
	iEvent.getByLabel("calibratedElectrons","eneErrorRegForGsfEle",ElectronRegEnergyError);

	edm::Handle<reco::ConversionCollection> hConversions;
	iEvent.getByLabel("allConversions", hConversions);

	// Isolation variables for cones of dR=0.4 and dR=0.3

	edm::Handle<edm::ValueMap<double> > isoChargedHandle04;
	iEvent.getByLabel(edm::InputTag("elPFIsoValueCharged04PFIdPFIso",""),isoChargedHandle04);
	const edm::ValueMap<double> electronIsoCharged04 = (*isoChargedHandle04.product());
	edm::Handle<edm::ValueMap<double> > isoNeutralHandle04;
	iEvent.getByLabel(edm::InputTag("elPFIsoValueNeutral04PFIdPFIso",""),isoNeutralHandle04);
	const edm::ValueMap<double> electronIsoNeutral04 = (*isoNeutralHandle04.product());
	edm::Handle<edm::ValueMap<double> > isoPhotonHandle04;
	iEvent.getByLabel(edm::InputTag("elPFIsoValueGamma04PFIdPFIso",""),isoPhotonHandle04);
	const edm::ValueMap<double> electronIsoPhoton04 = (*isoPhotonHandle04.product());

	edm::Handle<edm::ValueMap<double> > isoChargedHandle03;
	iEvent.getByLabel(edm::InputTag("elPFIsoValueCharged03PFIdPFIso",""),isoChargedHandle03);
	const edm::ValueMap<double> electronIsoCharged03 = (*isoChargedHandle03.product());
	edm::Handle<edm::ValueMap<double> > isoNeutralHandle03;
	iEvent.getByLabel(edm::InputTag("elPFIsoValueNeutral03PFIdPFIso",""),isoNeutralHandle03);
	const edm::ValueMap<double> electronIsoNeutral03 = (*isoNeutralHandle03.product());
	edm::Handle<edm::ValueMap<double> > isoPhotonHandle03;
	iEvent.getByLabel(edm::InputTag("elPFIsoValueGamma03PFIdPFIso",""),isoPhotonHandle03);
	const edm::ValueMap<double> electronIsoPhoton03 = (*isoPhotonHandle03.product());

	edm::Handle<double> RhoIsolation;
	iEvent.getByLabel(rhoIsolAllInputTag_, RhoIsolation);
	const double *RhoIsolationRef = RhoIsolation.product();

	RhoIsolationAllInputTags = *(RhoIsolationRef);

	/////////////////////////////////////////////////////////////////////////////////////
	//	                                                                              //
	// MVA IDs from https://twiki.cern.ch/twiki/bin/view/CMS/ElectronMVAIDForH2Tau     //
	//                                                                                 //
	// and https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentification //
	//                                                                                 //
	/////////////////////////////////////////////////////////////////////////////////////
	//const reco::GsfElectronCollection theEGamma = *(ElectronCollection.product());

	edm::Handle<reco::VertexCollection> dummyVertexCollection;
	iEvent.getByLabel("offlinePrimaryVertices", dummyVertexCollection);

	reco::Vertex dummy;
	const reco::Vertex *pv = &dummy;
	for (unsigned i = 0; i < dummyVertexCollection->size(); i++) {
		if (isGoodVertex(dummyVertexCollection->at(i))) {
			pv = &dummyVertexCollection->at(i);
			break;
		}
	}

	edm::ESHandle<TransientTrackBuilder> builder;
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
	TransientTrackBuilder theBuilder = *(builder.product());

	EcalClusterLazyTools EcalCluster(iEvent, iSetup, (edm::InputTag) "reducedEcalRecHitsEB", (edm::InputTag) "reducedEcalRecHitsEE");

	Double_t _Rho = 0;
	edm::Handle<double> rhoPtr;
	const edm::InputTag eventrho("kt6PFJets", "rho");
	iEvent.getByLabel(eventrho, rhoPtr);
	_Rho = *rhoPtr;

	double myMVATrig2012Method1 = -1;
	double myMVATrigNoIP2012Method1 = -1;
	double myMVANonTrig2012Method1 = -1;

	///////////////////
	//               //
	// END OF MVA ID //
	//               //
	///////////////////

	for (reco::PFCandidateCollection::size_type iPFElectron = 0; iPFElectron < ElectronCollection->size(); iPFElectron++) {
		reco::GsfElectronRef RefElectron(ElectronCollection, iPFElectron);
		if (isGoodElectron(RefElectron)) {
			std::vector<double> iElectron_Poca;
			iElectron_Poca.push_back(RefElectron->vx());
			iElectron_Poca.push_back(RefElectron->vy());
			iElectron_Poca.push_back(RefElectron->vz());
			Electron_Poca.push_back(iElectron_Poca);
			std::vector<double> iElectron_p4;
			iElectron_p4.push_back(RefElectron->p4().E());
			iElectron_p4.push_back(RefElectron->p4().Px());
			iElectron_p4.push_back(RefElectron->p4().Py());
			iElectron_p4.push_back(RefElectron->p4().Pz());

			Electron_p4.push_back(iElectron_p4);
			Electron_charge.push_back(RefElectron->charge());

			Electron_RegEnergy.push_back(ElectronRegEnergy->get(iPFElectron));
			Electron_RegEnergyError.push_back(ElectronRegEnergyError->get(iPFElectron));

			//////////////////////////
			myMVATrig2012Method1 = -1;
			myMVATrigNoIP2012Method1 = -1;
			myMVANonTrig2012Method1 = -1;
			myMVATrig2012Method1 = myMVATrig2012->mvaValue(*RefElectron, *pv, theBuilder, EcalCluster, false);
			myMVATrigNoIP2012Method1 = myMVATrigNoIP2012->mvaValue(*RefElectron, *pv, _Rho, EcalCluster, false);
			myMVANonTrig2012Method1 = myMVANonTrig2012->mvaValue(*RefElectron, *pv, theBuilder, EcalCluster, false);
			Electron_MVA_Trig_discriminator.push_back(myMVATrig2012Method1);
			Electron_MVA_TrigNoIP_discriminator.push_back(myMVATrigNoIP2012Method1);
			Electron_MVA_NonTrig_discriminator.push_back(myMVANonTrig2012Method1);
			//////////////////////////

			//     reco::GsfElectronRef refGsfElectron = RefElectron->gsfElectronRef();
			Electron_Gsf_deltaEtaEleClusterTrackAtCalo.push_back(RefElectron->deltaEtaEleClusterTrackAtCalo());
			Electron_Gsf_deltaEtaSeedClusterTrackAtCalo.push_back(RefElectron->deltaEtaSeedClusterTrackAtCalo());
			Electron_Gsf_deltaEtaSuperClusterTrackAtVtx.push_back(RefElectron->deltaEtaSuperClusterTrackAtVtx());
			Electron_Gsf_deltaPhiEleClusterTrackAtCalo.push_back(RefElectron->deltaPhiEleClusterTrackAtCalo());
			Electron_Gsf_deltaPhiSeedClusterTrackAtCalo.push_back(RefElectron->deltaPhiSeedClusterTrackAtCalo());
			Electron_Gsf_deltaPhiSuperClusterTrackAtVtx.push_back(RefElectron->deltaPhiSuperClusterTrackAtVtx());
			Electron_Gsf_dr03EcalRecHitSumE.push_back(RefElectron->dr03EcalRecHitSumEt());
			Electron_Gsf_dr03HcalDepth1TowerSumEt.push_back(RefElectron->dr03HcalDepth1TowerSumEt());
			Electron_Gsf_dr03HcalDepth1TowerSumEtBc.push_back(RefElectron->dr03HcalDepth1TowerSumEtBc());
			Electron_Gsf_dr03HcalDepth2TowerSumEt.push_back(RefElectron->dr03HcalDepth2TowerSumEt());
			Electron_Gsf_dr03HcalDepth2TowerSumEtBc.push_back(RefElectron->dr03HcalDepth2TowerSumEtBc());
			Electron_Gsf_dr03HcalTowerSumEt.push_back(RefElectron->dr03HcalTowerSumEt());
			Electron_Gsf_dr03HcalTowerSumEtBc.push_back(RefElectron->dr03HcalTowerSumEtBc());
			Electron_Gsf_dr03TkSumPt.push_back(RefElectron->dr03TkSumPt());
			Electron_Gsf_passingCutBasedPreselection.push_back(RefElectron->passingCutBasedPreselection());
			Electron_Gsf_passingMvaPreselection.push_back(RefElectron->passingMvaPreselection());

			Electron_sigmaIetaIeta.push_back(RefElectron->sigmaIetaIeta());
			Electron_hadronicOverEm.push_back(RefElectron->hadronicOverEm());
			Electron_fbrem.push_back(RefElectron->fbrem());
			Electron_eSuperClusterOverP.push_back(RefElectron->eSuperClusterOverP());
			Electron_ecalEnergy.push_back(RefElectron->ecalEnergy());
			Electron_trackMomentumAtVtx.push_back(RefElectron->ecalEnergy() / RefElectron->eSuperClusterOverP());

			reco::BeamSpot beamSpot;
			edm::Handle<reco::BeamSpot> beamSpotHandle;
			iEvent.getByLabel(beamSpotTag_, beamSpotHandle);
			if (beamSpotHandle.isValid()) {
				beamSpot = *beamSpotHandle;
			} else {
				edm::LogInfo("") << "No beam spot available from EventSetup \n";
			}

			Electron_numberOfMissedHits.push_back(RefElectron->gsfTrack()->trackerExpectedHitsInner().numberOfHits());
			Electron_HasMatchedConversions.push_back(ConversionTools::hasMatchedConversion(ElectronCollection->at(iPFElectron), hConversions, beamSpot.position(), true, 2.0, 1e-6, 0));
			// static bool check = !ConversionTools::hasMatchedConversion(ElectronCollection->at(iPFElectron), hConversions, beamSpot.position());

			Electron_ecalRecHitSumEt03.push_back(RefElectron->isolationVariables03().ecalRecHitSumEt);
			Electron_hcalDepth1TowerSumEt03.push_back(RefElectron->isolationVariables03().hcalDepth1TowerSumEt);
			Electron_hcalDepth1TowerSumEtBc03.push_back(RefElectron->isolationVariables03().hcalDepth1TowerSumEtBc);
			Electron_hcalDepth2TowerSumEt03.push_back(RefElectron->isolationVariables03().hcalDepth2TowerSumEt);
			Electron_hcalDepth2TowerSumEtBc03.push_back(RefElectron->isolationVariables03().hcalDepth2TowerSumEtBc);
			Electron_tkSumPt03.push_back(RefElectron->isolationVariables03().tkSumPt);
			Electron_ecalRecHitSumEt04.push_back(RefElectron->isolationVariables04().ecalRecHitSumEt);
			Electron_hcalDepth1TowerSumEt04.push_back(RefElectron->isolationVariables04().hcalDepth1TowerSumEt);
			Electron_hcalDepth1TowerSumEtBc04.push_back(RefElectron->isolationVariables04().hcalDepth1TowerSumEtBc);
			Electron_hcalDepth2TowerSumEt04.push_back(RefElectron->isolationVariables04().hcalDepth2TowerSumEt);
			Electron_hcalDepth2TowerSumEtBc04.push_back(RefElectron->isolationVariables04().hcalDepth2TowerSumEtBc);
			Electron_tkSumPt04.push_back(RefElectron->isolationVariables04().tkSumPt);

			Electron_chargedHadronIso.push_back(RefElectron->pfIsolationVariables().chargedHadronIso);
			Electron_neutralHadronIso.push_back(RefElectron->pfIsolationVariables().neutralHadronIso);
			Electron_photonIso.push_back(RefElectron->pfIsolationVariables().photonIso);

			Electron_isoDeposits_chargedHadronIso04.push_back(electronIsoCharged04.get(iPFElectron));
			Electron_isoDeposits_neutralHadronIso04.push_back(electronIsoNeutral04.get(iPFElectron));
			Electron_isoDeposits_photonIso04.push_back(electronIsoPhoton04.get(iPFElectron));
			Electron_isoDeposits_chargedHadronIso03.push_back(electronIsoCharged03.get(iPFElectron));
			Electron_isoDeposits_neutralHadronIso03.push_back(electronIsoNeutral03.get(iPFElectron));
			Electron_isoDeposits_photonIso03.push_back(electronIsoPhoton03.get(iPFElectron));

			reco::GsfTrackRef refGsfTrack = RefElectron->gsfTrack();
			Electron_gsftrack_trackerExpectedHitsInner_numberOfLostHits.push_back(refGsfTrack->trackerExpectedHitsInner().numberOfLostHits());

			reco::SuperClusterRef refSuperCluster = RefElectron->superCluster();
			Electron_supercluster_e.push_back(refSuperCluster->energy());
			Electron_supercluster_phi.push_back(refSuperCluster->phi());
			Electron_supercluster_eta.push_back(refSuperCluster->eta());
			Electron_supercluster_centroid_x.push_back(refSuperCluster->x());
			Electron_supercluster_centroid_y.push_back(refSuperCluster->y());
			Electron_supercluster_centroid_z.push_back(refSuperCluster->z());

			int match;
			getTrackMatch(trackCollection, RefElectron, match);
			Electron_Track_idx.push_back(match);

			int ntp = Electron_par.size();
			Electron_par.push_back(std::vector<double>());
			Electron_cov.push_back(std::vector<double>());
			if (refGsfTrack.isNonnull()) {
				GlobalPoint pvpoint(refGsfTrack->vx(), refGsfTrack->vy(), refGsfTrack->vz());
				edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
				iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transTrackBuilder);
				reco::TransientTrack transTrk = transTrackBuilder->build(refGsfTrack);
				TrackParticle trackparticle = ParticleBuilder::CreateTrackParticle(transTrk, transTrackBuilder, pvpoint, false, true);
				Electron_trackCharge.push_back(trackparticle.Charge());
				Electron_pdgid.push_back(trackparticle.PDGID());
				Electron_B.push_back(trackparticle.BField());
				Electron_M.push_back(trackparticle.Mass());
				for (int i = 0; i < trackparticle.NParameters(); i++) {
					Electron_par.at(ntp).push_back(trackparticle.Parameter(i));
					for (int j = i; j < trackparticle.NParameters(); j++) {
						Electron_cov.at(ntp).push_back(trackparticle.Covariance(i, j));
					}
				}
			} else {
				Electron_trackCharge.push_back(-999);
				Electron_pdgid.push_back(-999);
				Electron_B.push_back(-999);
				Electron_M.push_back(-999);
			}
			edm::Handle<double> Rhokt6PFJets;
			const edm::InputTag erho("kt6PFJets", "rho");
			iEvent.getByLabel(erho, Rhokt6PFJets);
			Electron_Rho_kt6PFJets.push_back(*Rhokt6PFJets);
		}
	}
}

void TauNtuple::fillMET(edm::Event& iEvent, const edm::EventSetup& iSetup) {
	if (!doPatMET_) {
		edm::Handle<edm::View<reco::PFMET> > pfMETUncorr;
		iEvent.getByLabel(pfMETUncorr_, pfMETUncorr);

		edm::Handle<std::vector<reco::PFMET> > pfMETCorrT0rt;
		iEvent.getByLabel(pfMETCorrT0rt_, pfMETCorrT0rt);

		edm::Handle<std::vector<reco::PFMET> > pfMETCorrT0rtT1;
		iEvent.getByLabel(pfMETCorrT0rtT1_, pfMETCorrT0rtT1);

		edm::Handle<std::vector<reco::PFMET> > pfMETCorrT0pc;
		iEvent.getByLabel(pfMETCorrT0pc_, pfMETCorrT0pc);

		edm::Handle<std::vector<reco::PFMET> > pfMETCorrT0pcT1;
		iEvent.getByLabel(pfMETCorrT0pcT1_, pfMETCorrT0pcT1);

		edm::Handle<std::vector<reco::PFMET> > pfMETCorrT0rtTxy;
		iEvent.getByLabel(pfMETCorrT0rtTxy_, pfMETCorrT0rtTxy);

		edm::Handle<std::vector<reco::PFMET> > pfMETCorrT0rtT1Txy;
		iEvent.getByLabel(pfMETCorrT0rtT1Txy_, pfMETCorrT0rtT1Txy);

		edm::Handle<std::vector<reco::PFMET> > pfMETCorrT0pcTxy;
		iEvent.getByLabel(pfMETCorrT0pcTxy_, pfMETCorrT0pcTxy);

		edm::Handle<std::vector<reco::PFMET> > pfMETCorrT0pcT1Txy;
		iEvent.getByLabel(pfMETCorrT0pcT1Txy_, pfMETCorrT0pcT1Txy);

		edm::Handle<std::vector<reco::PFMET> > pfMETCorrT1;
		iEvent.getByLabel(pfMETCorrT1_, pfMETCorrT1);

		edm::Handle<std::vector<reco::PFMET> > pfMETCorrT1Txy;
		iEvent.getByLabel(pfMETCorrT1Txy_, pfMETCorrT1Txy);

		edm::Handle<std::vector<reco::CaloMET> > caloMETCorrT1;
		iEvent.getByLabel(caloMETCorrT1_, caloMETCorrT1);

		edm::Handle<std::vector<reco::CaloMET> > caloMETCorrT1T2;
		iEvent.getByLabel(caloMETCorrT1T2_, caloMETCorrT1T2);

//     std::cout<<"-------------------------------------------------------------- >Corrected MET sie"<<CorrectedPFMET->size()<<std::endl;
//     std::cout<<"Corrected  et  pt  and phi "<< CorrectedPFMET->at(0).et() << "  " <<  CorrectedPFMET->at(0).pt()<<"  " <<  CorrectedPFMET->at(0).phi()<<std::endl;
//     std::cout<<"Un Corrected et   pt and phi "<< pfMETUncorr->front().et()<< "  " <<  pfMETUncorr->front().pt()<< "  " << pfMETUncorr->front().phi()<<std::endl;

		TMatrixD sigMat(2, 2);

		MET_Uncorr_et = pfMETUncorr->front().et();
		MET_Uncorr_pt = pfMETUncorr->front().pt();
		MET_Uncorr_phi = pfMETUncorr->front().phi();
		MET_Uncorr_sumET = pfMETUncorr->front().sumEt();
		MET_Uncorr_significance = pfMETUncorr->front().significance();
		sigMat = pfMETUncorr->front().getSignificanceMatrix();
		if (sigMat(0, 1) != sigMat(1, 0))
			std::cout << "WARNING: MET significance matrix not symmetric" << std::endl;
		MET_Uncorr_significance_xx = sigMat(0, 0);
		MET_Uncorr_significance_xy = sigMat(0, 1);
		MET_Uncorr_significance_yy = sigMat(1, 1);
		MET_Uncorr_MuonEtFraction = pfMETUncorr->front().muonEtFraction();
		MET_Uncorr_NeutralEMFraction = pfMETUncorr->front().NeutralEMFraction();
		MET_Uncorr_NeutralHadEtFraction = pfMETUncorr->front().NeutralHadEtFraction();
		MET_Uncorr_Type6EtFraction = pfMETUncorr->front().Type6EtFraction();
		MET_Uncorr_Type7EtFraction = pfMETUncorr->front().Type7EtFraction();

		MET_CorrT0rt_et = pfMETCorrT0rt->front().et();
		MET_CorrT0rt_pt = pfMETCorrT0rt->front().pt();
		MET_CorrT0rt_phi = pfMETCorrT0rt->front().phi();
		MET_CorrT0rt_sumET = pfMETCorrT0rt->front().sumEt();
		MET_CorrT0rt_MuonEtFraction = pfMETCorrT0rt->front().muonEtFraction();
		MET_CorrT0rt_NeutralEMFraction = pfMETCorrT0rt->front().NeutralEMFraction();
		MET_CorrT0rt_NeutralHadEtFraction = pfMETCorrT0rt->front().NeutralHadEtFraction();
		MET_CorrT0rt_Type6EtFraction = pfMETCorrT0rt->front().Type6EtFraction();
		MET_CorrT0rt_Type7EtFraction = pfMETCorrT0rt->front().Type7EtFraction();

		MET_CorrT0rtT1_et = pfMETCorrT0rtT1->front().et();
		MET_CorrT0rtT1_pt = pfMETCorrT0rtT1->front().pt();
		MET_CorrT0rtT1_phi = pfMETCorrT0rtT1->front().phi();
		MET_CorrT0rtT1_sumET = pfMETCorrT0rtT1->front().sumEt();
		MET_CorrT0rtT1_MuonEtFraction = pfMETCorrT0rtT1->front().muonEtFraction();
		MET_CorrT0rtT1_NeutralEMFraction = pfMETCorrT0rtT1->front().NeutralEMFraction();
		MET_CorrT0rtT1_NeutralHadEtFraction = pfMETCorrT0rtT1->front().NeutralHadEtFraction();
		MET_CorrT0rtT1_Type6EtFraction = pfMETCorrT0rtT1->front().Type6EtFraction();
		MET_CorrT0rtT1_Type7EtFraction = pfMETCorrT0rtT1->front().Type7EtFraction();

		MET_CorrT0pc_et = pfMETCorrT0pc->front().et();
		MET_CorrT0pc_pt = pfMETCorrT0pc->front().pt();
		MET_CorrT0pc_phi = pfMETCorrT0pc->front().phi();
		MET_CorrT0pc_sumET = pfMETCorrT0pc->front().sumEt();
		MET_CorrT0pc_MuonEtFraction = pfMETCorrT0pc->front().muonEtFraction();
		MET_CorrT0pc_NeutralEMFraction = pfMETCorrT0pc->front().NeutralEMFraction();
		MET_CorrT0pc_NeutralHadEtFraction = pfMETCorrT0pc->front().NeutralHadEtFraction();
		MET_CorrT0pc_Type6EtFraction = pfMETCorrT0pc->front().Type6EtFraction();
		MET_CorrT0pc_Type7EtFraction = pfMETCorrT0pc->front().Type7EtFraction();

		MET_CorrT0pcT1_et = pfMETCorrT0pcT1->front().et();
		MET_CorrT0pcT1_pt = pfMETCorrT0pcT1->front().pt();
		MET_CorrT0pcT1_phi = pfMETCorrT0pcT1->front().phi();
		MET_CorrT0pcT1_sumET = pfMETCorrT0pcT1->front().sumEt();
		MET_CorrT0pcT1_MuonEtFraction = pfMETCorrT0pcT1->front().muonEtFraction();
		MET_CorrT0pcT1_NeutralEMFraction = pfMETCorrT0pcT1->front().NeutralEMFraction();
		MET_CorrT0pcT1_NeutralHadEtFraction = pfMETCorrT0pcT1->front().NeutralHadEtFraction();
		MET_CorrT0pcT1_Type6EtFraction = pfMETCorrT0pcT1->front().Type6EtFraction();
		MET_CorrT0pcT1_Type7EtFraction = pfMETCorrT0pcT1->front().Type7EtFraction();

		MET_CorrT0rtTxy_et = pfMETCorrT0rtTxy->front().et();
		MET_CorrT0rtTxy_pt = pfMETCorrT0rtTxy->front().pt();
		MET_CorrT0rtTxy_phi = pfMETCorrT0rtTxy->front().phi();
		MET_CorrT0rtTxy_sumET = pfMETCorrT0rtTxy->front().sumEt();
		MET_CorrT0rtTxy_MuonEtFraction = pfMETCorrT0rtTxy->front().muonEtFraction();
		MET_CorrT0rtTxy_NeutralEMFraction = pfMETCorrT0rtTxy->front().NeutralEMFraction();
		MET_CorrT0rtTxy_NeutralHadEtFraction = pfMETCorrT0rtTxy->front().NeutralHadEtFraction();
		MET_CorrT0rtTxy_Type6EtFraction = pfMETCorrT0rtTxy->front().Type6EtFraction();
		MET_CorrT0rtTxy_Type7EtFraction = pfMETCorrT0rtTxy->front().Type7EtFraction();

		MET_CorrT0rtT1Txy_et = pfMETCorrT0rtT1Txy->front().et();
		MET_CorrT0rtT1Txy_pt = pfMETCorrT0rtT1Txy->front().pt();
		MET_CorrT0rtT1Txy_phi = pfMETCorrT0rtT1Txy->front().phi();
		MET_CorrT0rtT1Txy_sumET = pfMETCorrT0rtT1Txy->front().sumEt();
		MET_CorrT0rtT1Txy_MuonEtFraction = pfMETCorrT0rtT1Txy->front().muonEtFraction();
		MET_CorrT0rtT1Txy_NeutralEMFraction = pfMETCorrT0rtT1Txy->front().NeutralEMFraction();
		MET_CorrT0rtT1Txy_NeutralHadEtFraction = pfMETCorrT0rtT1Txy->front().NeutralHadEtFraction();
		MET_CorrT0rtT1Txy_Type6EtFraction = pfMETCorrT0rtT1Txy->front().Type6EtFraction();
		MET_CorrT0rtT1Txy_Type7EtFraction = pfMETCorrT0rtT1Txy->front().Type7EtFraction();

		MET_CorrT0pcTxy_et = pfMETCorrT0pcTxy->front().et();
		MET_CorrT0pcTxy_pt = pfMETCorrT0pcTxy->front().pt();
		MET_CorrT0pcTxy_phi = pfMETCorrT0pcTxy->front().phi();
		MET_CorrT0pcTxy_sumET = pfMETCorrT0pcTxy->front().sumEt();
		MET_CorrT0pcTxy_MuonEtFraction = pfMETCorrT0pcTxy->front().muonEtFraction();
		MET_CorrT0pcTxy_NeutralEMFraction = pfMETCorrT0pcTxy->front().NeutralEMFraction();
		MET_CorrT0pcTxy_NeutralHadEtFraction = pfMETCorrT0pcTxy->front().NeutralHadEtFraction();
		MET_CorrT0pcTxy_Type6EtFraction = pfMETCorrT0pcTxy->front().Type6EtFraction();
		MET_CorrT0pcTxy_Type7EtFraction = pfMETCorrT0pcTxy->front().Type7EtFraction();

		MET_CorrT0pcT1Txy_et = pfMETCorrT0pcT1Txy->front().et();
		MET_CorrT0pcT1Txy_pt = pfMETCorrT0pcT1Txy->front().pt();
		MET_CorrT0pcT1Txy_phi = pfMETCorrT0pcT1Txy->front().phi();
		MET_CorrT0pcT1Txy_sumET = pfMETCorrT0pcT1Txy->front().sumEt();
		MET_CorrT0pcT1Txy_MuonEtFraction = pfMETCorrT0pcT1Txy->front().muonEtFraction();
		MET_CorrT0pcT1Txy_NeutralEMFraction = pfMETCorrT0pcT1Txy->front().NeutralEMFraction();
		MET_CorrT0pcT1Txy_NeutralHadEtFraction = pfMETCorrT0pcT1Txy->front().NeutralHadEtFraction();
		MET_CorrT0pcT1Txy_Type6EtFraction = pfMETCorrT0pcT1Txy->front().Type6EtFraction();
		MET_CorrT0pcT1Txy_Type7EtFraction = pfMETCorrT0pcT1Txy->front().Type7EtFraction();

		MET_CorrT1_et = pfMETCorrT1->front().et();
		MET_CorrT1_pt = pfMETCorrT1->front().pt();
		MET_CorrT1_phi = pfMETCorrT1->front().phi();
		MET_CorrT1_sumET = pfMETCorrT1->front().sumEt();
		MET_CorrT1_MuonEtFraction = pfMETCorrT1->front().muonEtFraction();
		MET_CorrT1_NeutralEMFraction = pfMETCorrT1->front().NeutralEMFraction();
		MET_CorrT1_NeutralHadEtFraction = pfMETCorrT1->front().NeutralHadEtFraction();
		MET_CorrT1_Type6EtFraction = pfMETCorrT1->front().Type6EtFraction();
		MET_CorrT1_Type7EtFraction = pfMETCorrT1->front().Type7EtFraction();

		MET_CorrT1Txy_et = pfMETCorrT1Txy->front().et();
		MET_CorrT1Txy_pt = pfMETCorrT1Txy->front().pt();
		MET_CorrT1Txy_phi = pfMETCorrT1Txy->front().phi();
		MET_CorrT1Txy_sumET = pfMETCorrT1Txy->front().sumEt();
		MET_CorrT1Txy_MuonEtFraction = pfMETCorrT1Txy->front().muonEtFraction();
		MET_CorrT1Txy_NeutralEMFraction = pfMETCorrT1Txy->front().NeutralEMFraction();
		MET_CorrT1Txy_NeutralHadEtFraction = pfMETCorrT1Txy->front().NeutralHadEtFraction();
		MET_CorrT1Txy_Type6EtFraction = pfMETCorrT1Txy->front().Type6EtFraction();
		MET_CorrT1Txy_Type7EtFraction = pfMETCorrT1Txy->front().Type7EtFraction();

		MET_CorrCaloT1_et = caloMETCorrT1->front().et();
		MET_CorrCaloT1_pt = caloMETCorrT1->front().pt();
		MET_CorrCaloT1_phi = caloMETCorrT1->front().phi();
		MET_CorrCaloT1_sumET = caloMETCorrT1->front().sumEt();

		MET_CorrCaloT1T2_et = caloMETCorrT1T2->front().et();
		MET_CorrCaloT1T2_pt = caloMETCorrT1T2->front().pt();
		MET_CorrCaloT1T2_phi = caloMETCorrT1T2->front().phi();
		MET_CorrCaloT1T2_sumET = caloMETCorrT1T2->front().sumEt();

		// MET uncertainties only available in PAT
		MET_Type1Corr_et = -1;
		MET_Type1Corr_pt = -1;
		MET_Type1Corr_phi = -1;
		MET_Type1Corr_sumET = -1;
		MET_Type1Corr_MuonEtFraction = -1;
		MET_Type1Corr_NeutralEMFraction = -1;
		MET_Type1Corr_NeutralHadEtFraction = -1;
		MET_Type1Corr_Type6EtFraction = -1;
		MET_Type1Corr_Type7EtFraction = -1;

		MET_Type1p2Corr_et = -1;
		MET_Type1p2Corr_pt = -1;
		MET_Type1p2Corr_phi = -1;
		MET_Type1p2Corr_sumET = -1;
		MET_Type1p2Corr_MuonEtFraction = -1;
		MET_Type1p2Corr_NeutralEMFraction = -1;
		MET_Type1p2Corr_NeutralHadEtFraction = -1;
		MET_Type1p2Corr_Type6EtFraction = -1;
		MET_Type1p2Corr_Type7EtFraction = -1;

		MET_Type1CorrElectronUp_et = -1;
		MET_Type1CorrElectronDown_et = -1;
		MET_Type1CorrMuonUp_et = -1;
		MET_Type1CorrMuonDown_et = -1;
		MET_Type1CorrTauUp_et = -1;
		MET_Type1CorrTauDown_et = -1;
		MET_Type1CorrJetResUp_et = -1;
		MET_Type1CorrJetResDown_et = -1;
		MET_Type1CorrJetEnUp_et = -1;
		MET_Type1CorrJetEnDown_et = -1;
		MET_Type1CorrUnclusteredUp_et = -1;
		MET_Type1CorrUnclusteredDown_et = -1;

		MET_Type1p2CorrElectronUp_et = -1;
		MET_Type1p2CorrElectronDown_et = -1;
		MET_Type1p2CorrMuonUp_et = -1;
		MET_Type1p2CorrMuonDown_et = -1;
		MET_Type1p2CorrTauUp_et = -1;
		MET_Type1p2CorrTauDown_et = -1;
		MET_Type1p2CorrJetResUp_et = -1;
		MET_Type1p2CorrJetResDown_et = -1;
		MET_Type1p2CorrJetEnUp_et = -1;
		MET_Type1p2CorrJetEnDown_et = -1;
		MET_Type1p2CorrUnclusteredUp_et = -1;
		MET_Type1p2CorrUnclusteredDown_et = -1;

		if (doMVAMET_) {
			edm::Handle<std::vector<reco::PFMET>> pfMETCorrMVA;
			iEvent.getByLabel(pfMETCorrMVA_, pfMETCorrMVA);

			edm::Handle<reco::MuonCollection> mvametMuonCollection;
			iEvent.getByLabel(muonsForPfMetCorrMVA_, mvametMuonCollection);
			edm::Handle<reco::GsfElectronCollection> mvametElectronCollection;
			iEvent.getByLabel(elecsForPfMetCorrMVA_, mvametElectronCollection);
			edm::Handle<std::vector<reco::PFTau> > mvametTauCollection;
			iEvent.getByLabel(tausForPfMetCorrMVA_, mvametTauCollection);

			MET_CorrMVA_et = pfMETCorrMVA->front().et();
			MET_CorrMVA_pt = pfMETCorrMVA->front().pt();
			MET_CorrMVA_phi = pfMETCorrMVA->front().phi();
			MET_CorrMVA_sumET = pfMETCorrMVA->front().sumEt();
			MET_CorrMVA_significance = pfMETCorrMVA->front().significance();
			sigMat = pfMETCorrMVA->front().getSignificanceMatrix();
			if (sigMat(0, 1) != sigMat(1, 0))
				std::cout << "WARNING: MET significance matrix not symmetric" << std::endl;
			MET_CorrMVA_significance_xx = sigMat(0, 0);
			MET_CorrMVA_significance_xy = sigMat(0, 1);
			MET_CorrMVA_significance_yy = sigMat(1, 1);
			MET_CorrMVA_MuonEtFraction = pfMETCorrMVA->front().muonEtFraction();
			MET_CorrMVA_NeutralEMFraction = pfMETCorrMVA->front().NeutralEMFraction();
			MET_CorrMVA_NeutralHadEtFraction = pfMETCorrMVA->front().NeutralHadEtFraction();
			MET_CorrMVA_Type6EtFraction = pfMETCorrMVA->front().Type6EtFraction();
			MET_CorrMVA_Type7EtFraction = pfMETCorrMVA->front().Type7EtFraction();
			// store p4 of leptons used for MVA-MET calculation
			int Muon_index = 0;
			for (reco::MuonCollection::const_iterator iMuon = mvametMuonCollection->begin(); iMuon != mvametMuonCollection->end(); ++iMuon, Muon_index++) {
				reco::MuonRef RefMuon(mvametMuonCollection, Muon_index);
				std::vector<double> iMuon_p4;
				iMuon_p4.push_back(RefMuon->p4().E());
				iMuon_p4.push_back(RefMuon->p4().Px());
				iMuon_p4.push_back(RefMuon->p4().Py());
				iMuon_p4.push_back(RefMuon->p4().Pz());
				MET_CorrMVA_srcMuon_p4.push_back(iMuon_p4);
			}
			for (reco::PFCandidateCollection::size_type iElectron = 0; iElectron < mvametElectronCollection->size(); iElectron++) {
				reco::GsfElectronRef RefElectron(mvametElectronCollection, iElectron);
				std::vector<double> iElectron_p4;
				iElectron_p4.push_back(RefElectron->p4().E());
				iElectron_p4.push_back(RefElectron->p4().Px());
				iElectron_p4.push_back(RefElectron->p4().Py());
				iElectron_p4.push_back(RefElectron->p4().Pz());
				MET_CorrMVA_srcElectron_p4.push_back(iElectron_p4);
			}
			for (unsigned iPFTau = 0; iPFTau < mvametTauCollection->size(); ++iPFTau) {
				reco::PFTauRef HPStauCandidate(mvametTauCollection, iPFTau);
				std::vector<double> iPFTau_p4;
				iPFTau_p4.push_back(HPStauCandidate->p4().E());
				iPFTau_p4.push_back(HPStauCandidate->p4().Px());
				iPFTau_p4.push_back(HPStauCandidate->p4().Py());
				iPFTau_p4.push_back(HPStauCandidate->p4().Pz());
				MET_CorrMVA_srcTau_p4.push_back(iPFTau_p4);
			}

			edm::Handle<std::vector<reco::PFMET>> pfMETCorrMVAMuTau;
			iEvent.getByLabel(pfMETCorrMVAMuTau_, pfMETCorrMVAMuTau);

			edm::Handle<reco::MuonCollection> mvamutaumetMuonCollection;
			iEvent.getByLabel(muonsForPfMetCorrMVAMuTau_, mvamutaumetMuonCollection);
			edm::Handle<std::vector<reco::PFTau> > mvamutaumetTauCollection;
			iEvent.getByLabel(tausForPfMetCorrMVAMuTau_, mvamutaumetTauCollection);

			MET_CorrMVAMuTau_et = pfMETCorrMVAMuTau->front().et();
			MET_CorrMVAMuTau_pt = pfMETCorrMVAMuTau->front().pt();
			MET_CorrMVAMuTau_phi = pfMETCorrMVAMuTau->front().phi();
			MET_CorrMVAMuTau_sumET = pfMETCorrMVAMuTau->front().sumEt();
			MET_CorrMVAMuTau_significance = pfMETCorrMVAMuTau->front().significance();
			sigMat = pfMETCorrMVAMuTau->front().getSignificanceMatrix();
			if (sigMat(0, 1) != sigMat(1, 0))
				std::cout << "WARNING: MET significance matrix not symmetric" << std::endl;
			MET_CorrMVAMuTau_significance_xx = sigMat(0, 0);
			MET_CorrMVAMuTau_significance_xy = sigMat(0, 1);
			MET_CorrMVAMuTau_significance_yy = sigMat(1, 1);
			MET_CorrMVAMuTau_MuonEtFraction = pfMETCorrMVAMuTau->front().muonEtFraction();
			MET_CorrMVAMuTau_NeutralEMFraction = pfMETCorrMVAMuTau->front().NeutralEMFraction();
			MET_CorrMVAMuTau_NeutralHadEtFraction = pfMETCorrMVAMuTau->front().NeutralHadEtFraction();
			MET_CorrMVAMuTau_Type6EtFraction = pfMETCorrMVAMuTau->front().Type6EtFraction();
			MET_CorrMVAMuTau_Type7EtFraction = pfMETCorrMVAMuTau->front().Type7EtFraction();

			// store p4 of leptons used for MVA-MET calculation
			Muon_index = 0;
			for (reco::MuonCollection::const_iterator iMuon = mvamutaumetMuonCollection->begin(); iMuon != mvamutaumetMuonCollection->end(); ++iMuon, Muon_index++) {
				reco::MuonRef RefMuon(mvamutaumetMuonCollection, Muon_index);
				std::vector<double> iMuon_p4;
				iMuon_p4.push_back(RefMuon->p4().E());
				iMuon_p4.push_back(RefMuon->p4().Px());
				iMuon_p4.push_back(RefMuon->p4().Py());
				iMuon_p4.push_back(RefMuon->p4().Pz());
				MET_CorrMVAMuTau_srcMuon_p4.push_back(iMuon_p4);
			}
			for (unsigned iPFTau = 0; iPFTau < mvamutaumetTauCollection->size(); ++iPFTau) {
				reco::PFTauRef HPStauCandidate(mvamutaumetTauCollection, iPFTau);
				std::vector<double> iPFTau_p4;
				iPFTau_p4.push_back(HPStauCandidate->p4().E());
				iPFTau_p4.push_back(HPStauCandidate->p4().Px());
				iPFTau_p4.push_back(HPStauCandidate->p4().Py());
				iPFTau_p4.push_back(HPStauCandidate->p4().Pz());
				MET_CorrMVAMuTau_srcTau_p4.push_back(iPFTau_p4);
			}
		}

	} else {

		edm::Handle<std::vector<pat::MET>> patMETUncorrHandle;
		iEvent.getByLabel(patMETUncorr_, patMETUncorrHandle);
		pat::MET patMETUncorr = patMETUncorrHandle->front();

		edm::Handle<std::vector<pat::MET>> patMETCorrT0rtHandle;
		iEvent.getByLabel(patMETCorrT0rt_, patMETCorrT0rtHandle);
		pat::MET patMETCorrT0rt = patMETCorrT0rtHandle->front();

		edm::Handle<std::vector<pat::MET>> patMETCorrT0rtT1Handle;
		iEvent.getByLabel(patMETCorrT0rtT1_, patMETCorrT0rtT1Handle);
		pat::MET patMETCorrT0rtT1 = patMETCorrT0rtT1Handle->front();

		edm::Handle<std::vector<pat::MET>> patMETCorrT0pcHandle;
		iEvent.getByLabel(patMETCorrT0pc_, patMETCorrT0pcHandle);
		pat::MET patMETCorrT0pc = patMETCorrT0pcHandle->front();

		edm::Handle<std::vector<pat::MET>> patMETCorrT0pcT1Handle;
		iEvent.getByLabel(patMETCorrT0pcT1_, patMETCorrT0pcT1Handle);
		pat::MET patMETCorrT0pcT1 = patMETCorrT0pcT1Handle->front();

		edm::Handle<std::vector<pat::MET>> patMETCorrT0rtTxyHandle;
		iEvent.getByLabel(patMETCorrT0rtTxy_, patMETCorrT0rtTxyHandle);
		pat::MET patMETCorrT0rtTxy = patMETCorrT0rtTxyHandle->front();

		edm::Handle<std::vector<pat::MET>> patMETCorrT0rtT1TxyHandle;
		iEvent.getByLabel(patMETCorrT0rtT1Txy_, patMETCorrT0rtT1TxyHandle);
		pat::MET patMETCorrT0rtT1Txy = patMETCorrT0rtT1TxyHandle->front();

		edm::Handle<std::vector<pat::MET>> patMETCorrT0pcTxyHandle;
		iEvent.getByLabel(patMETCorrT0pcTxy_, patMETCorrT0pcTxyHandle);
		pat::MET patMETCorrT0pcTxy = patMETCorrT0pcTxyHandle->front();

		edm::Handle<std::vector<pat::MET>> patMETCorrT0pcT1TxyHandle;
		iEvent.getByLabel(patMETCorrT0pcT1Txy_, patMETCorrT0pcT1TxyHandle);
		pat::MET patMETCorrT0pcT1Txy = patMETCorrT0pcT1TxyHandle->front();

		edm::Handle<std::vector<pat::MET>> patMETCorrT1Handle;
		iEvent.getByLabel(patMETCorrT1_, patMETCorrT1Handle);
		pat::MET patMETCorrT1 = patMETCorrT1Handle->front();

		edm::Handle<std::vector<pat::MET>> patMETCorrT1TxyHandle;
		iEvent.getByLabel(patMETCorrT1Txy_, patMETCorrT1TxyHandle);
		pat::MET patMETCorrT1Txy = patMETCorrT1TxyHandle->front();

		edm::Handle<std::vector<pat::MET>> patMETCorrCaloT1Handle;
		iEvent.getByLabel(patCaloMETCorrT1_, patMETCorrCaloT1Handle);
		pat::MET patMETCorrCaloT1 = patMETCorrCaloT1Handle->front();

		edm::Handle<std::vector<pat::MET>> patMETCorrCaloT1T2Handle;
		iEvent.getByLabel(patCaloMETCorrT1T2_, patMETCorrCaloT1T2Handle);
		pat::MET patMETCorrCaloT1T2 = patMETCorrCaloT1T2Handle->front();

		TMatrixD sigMat(2, 2);

		MET_Uncorr_et = patMETUncorr.et();
		MET_Uncorr_pt = patMETUncorr.pt();
		MET_Uncorr_phi = patMETUncorr.phi();
		MET_Uncorr_sumET = patMETUncorr.sumEt();
		MET_Uncorr_significance = patMETUncorr.significance();
		sigMat = patMETUncorr.getSignificanceMatrix();
		if (sigMat(0, 1) != sigMat(1, 0))
			std::cout << "WARNING: MET significance matrix not symmetric" << std::endl;
		MET_Uncorr_significance_xx = sigMat(0, 0);
		MET_Uncorr_significance_xy = sigMat(0, 1);
		MET_Uncorr_significance_yy = sigMat(1, 1);
		if (patMETUncorr.isPFMET()) {
			MET_Uncorr_MuonEtFraction = patMETUncorr.MuonEtFraction();
			MET_Uncorr_NeutralEMFraction = patMETUncorr.NeutralEMFraction();
			MET_Uncorr_NeutralHadEtFraction = patMETUncorr.NeutralHadEtFraction();
			MET_Uncorr_Type6EtFraction = patMETUncorr.Type6EtFraction();
			MET_Uncorr_Type7EtFraction = patMETUncorr.Type7EtFraction();
		}

		MET_CorrT0rt_et = patMETCorrT0rt.et();
		MET_CorrT0rt_pt = patMETCorrT0rt.pt();
		MET_CorrT0rt_phi = patMETCorrT0rt.phi();
		MET_CorrT0rt_sumET = patMETCorrT0rt.sumEt();
		if (patMETCorrT0rt.isPFMET()) {
			MET_CorrT0rt_MuonEtFraction = patMETCorrT0rt.MuonEtFraction();
			MET_CorrT0rt_NeutralEMFraction = patMETCorrT0rt.NeutralEMFraction();
			MET_CorrT0rt_NeutralHadEtFraction = patMETCorrT0rt.NeutralHadEtFraction();
			MET_CorrT0rt_Type6EtFraction = patMETCorrT0rt.Type6EtFraction();
			MET_CorrT0rt_Type7EtFraction = patMETCorrT0rt.Type7EtFraction();
		}

		MET_CorrT0rtT1_et = patMETCorrT0rtT1.et();
		MET_CorrT0rtT1_pt = patMETCorrT0rtT1.pt();
		MET_CorrT0rtT1_phi = patMETCorrT0rtT1.phi();
		MET_CorrT0rtT1_sumET = patMETCorrT0rtT1.sumEt();
		if (patMETCorrT0rtT1.isPFMET()) {
			MET_CorrT0rtT1_MuonEtFraction = patMETCorrT0rtT1.MuonEtFraction();
			MET_CorrT0rtT1_NeutralEMFraction = patMETCorrT0rtT1.NeutralEMFraction();
			MET_CorrT0rtT1_NeutralHadEtFraction = patMETCorrT0rtT1.NeutralHadEtFraction();
			MET_CorrT0rtT1_Type6EtFraction = patMETCorrT0rtT1.Type6EtFraction();
			MET_CorrT0rtT1_Type7EtFraction = patMETCorrT0rtT1.Type7EtFraction();
		}

		MET_CorrT0pc_et = patMETCorrT0pc.et();
		MET_CorrT0pc_pt = patMETCorrT0pc.pt();
		MET_CorrT0pc_phi = patMETCorrT0pc.phi();
		MET_CorrT0pc_sumET = patMETCorrT0pc.sumEt();
		if (patMETCorrT0pc.isPFMET()) {
			MET_CorrT0pc_MuonEtFraction = patMETCorrT0pc.MuonEtFraction();
			MET_CorrT0pc_NeutralEMFraction = patMETCorrT0pc.NeutralEMFraction();
			MET_CorrT0pc_NeutralHadEtFraction = patMETCorrT0pc.NeutralHadEtFraction();
			MET_CorrT0pc_Type6EtFraction = patMETCorrT0pc.Type6EtFraction();
			MET_CorrT0pc_Type7EtFraction = patMETCorrT0pc.Type7EtFraction();
		}

		MET_CorrT0pcT1_et = patMETCorrT0pcT1.et();
		MET_CorrT0pcT1_pt = patMETCorrT0pcT1.pt();
		MET_CorrT0pcT1_phi = patMETCorrT0pcT1.phi();
		MET_CorrT0pcT1_sumET = patMETCorrT0pcT1.sumEt();
		if (patMETCorrT0pcT1.isPFMET()) {
			MET_CorrT0pcT1_MuonEtFraction = patMETCorrT0pcT1.MuonEtFraction();
			MET_CorrT0pcT1_NeutralEMFraction = patMETCorrT0pcT1.NeutralEMFraction();
			MET_CorrT0pcT1_NeutralHadEtFraction = patMETCorrT0pcT1.NeutralHadEtFraction();
			MET_CorrT0pcT1_Type6EtFraction = patMETCorrT0pcT1.Type6EtFraction();
			MET_CorrT0pcT1_Type7EtFraction = patMETCorrT0pcT1.Type7EtFraction();
		}

		MET_CorrT0rtTxy_et = patMETCorrT0rtTxy.et();
		MET_CorrT0rtTxy_pt = patMETCorrT0rtTxy.pt();
		MET_CorrT0rtTxy_phi = patMETCorrT0rtTxy.phi();
		MET_CorrT0rtTxy_sumET = patMETCorrT0rtTxy.sumEt();
		if (patMETCorrT0rtTxy.isPFMET()) {
			MET_CorrT0rtTxy_MuonEtFraction = patMETCorrT0rtTxy.MuonEtFraction();
			MET_CorrT0rtTxy_NeutralEMFraction = patMETCorrT0rtTxy.NeutralEMFraction();
			MET_CorrT0rtTxy_NeutralHadEtFraction = patMETCorrT0rtTxy.NeutralHadEtFraction();
			MET_CorrT0rtTxy_Type6EtFraction = patMETCorrT0rtTxy.Type6EtFraction();
			MET_CorrT0rtTxy_Type7EtFraction = patMETCorrT0rtTxy.Type7EtFraction();
		}

		MET_CorrT0rtT1Txy_et = patMETCorrT0rtT1Txy.et();
		MET_CorrT0rtT1Txy_pt = patMETCorrT0rtT1Txy.pt();
		MET_CorrT0rtT1Txy_phi = patMETCorrT0rtT1Txy.phi();
		MET_CorrT0rtT1Txy_sumET = patMETCorrT0rtT1Txy.sumEt();
		if (patMETCorrT0rtT1Txy.isPFMET()) {
			MET_CorrT0rtT1Txy_MuonEtFraction = patMETCorrT0rtT1Txy.MuonEtFraction();
			MET_CorrT0rtT1Txy_NeutralEMFraction = patMETCorrT0rtT1Txy.NeutralEMFraction();
			MET_CorrT0rtT1Txy_NeutralHadEtFraction = patMETCorrT0rtT1Txy.NeutralHadEtFraction();
			MET_CorrT0rtT1Txy_Type6EtFraction = patMETCorrT0rtT1Txy.Type6EtFraction();
			MET_CorrT0rtT1Txy_Type7EtFraction = patMETCorrT0rtT1Txy.Type7EtFraction();
		}

		MET_CorrT0pcTxy_et = patMETCorrT0pcTxy.et();
		MET_CorrT0pcTxy_pt = patMETCorrT0pcTxy.pt();
		MET_CorrT0pcTxy_phi = patMETCorrT0pcTxy.phi();
		MET_CorrT0pcTxy_sumET = patMETCorrT0pcTxy.sumEt();
		if (patMETCorrT0pcTxy.isPFMET()) {
			MET_CorrT0pcTxy_MuonEtFraction = patMETCorrT0pcTxy.MuonEtFraction();
			MET_CorrT0pcTxy_NeutralEMFraction = patMETCorrT0pcTxy.NeutralEMFraction();
			MET_CorrT0pcTxy_NeutralHadEtFraction = patMETCorrT0pcTxy.NeutralHadEtFraction();
			MET_CorrT0pcTxy_Type6EtFraction = patMETCorrT0pcTxy.Type6EtFraction();
			MET_CorrT0pcTxy_Type7EtFraction = patMETCorrT0pcTxy.Type7EtFraction();
		}

		MET_CorrT0pcT1Txy_et = patMETCorrT0pcT1Txy.et();
		MET_CorrT0pcT1Txy_pt = patMETCorrT0pcT1Txy.pt();
		MET_CorrT0pcT1Txy_phi = patMETCorrT0pcT1Txy.phi();
		MET_CorrT0pcT1Txy_sumET = patMETCorrT0pcT1Txy.sumEt();
		if (patMETCorrT0pcT1Txy.isPFMET()) {
			MET_CorrT0pcT1Txy_MuonEtFraction = patMETCorrT0pcT1Txy.MuonEtFraction();
			MET_CorrT0pcT1Txy_NeutralEMFraction = patMETCorrT0pcT1Txy.NeutralEMFraction();
			MET_CorrT0pcT1Txy_NeutralHadEtFraction = patMETCorrT0pcT1Txy.NeutralHadEtFraction();
			MET_CorrT0pcT1Txy_Type6EtFraction = patMETCorrT0pcT1Txy.Type6EtFraction();
			MET_CorrT0pcT1Txy_Type7EtFraction = patMETCorrT0pcT1Txy.Type7EtFraction();
		}

		MET_CorrT1_et = patMETCorrT1.et();
		MET_CorrT1_pt = patMETCorrT1.pt();
		MET_CorrT1_phi = patMETCorrT1.phi();
		MET_CorrT1_sumET = patMETCorrT1.sumEt();
		if (patMETCorrT1.isPFMET()) {
			MET_CorrT1_MuonEtFraction = patMETCorrT1.MuonEtFraction();
			MET_CorrT1_NeutralEMFraction = patMETCorrT1.NeutralEMFraction();
			MET_CorrT1_NeutralHadEtFraction = patMETCorrT1.NeutralHadEtFraction();
			MET_CorrT1_Type6EtFraction = patMETCorrT1.Type6EtFraction();
			MET_CorrT1_Type7EtFraction = patMETCorrT1.Type7EtFraction();
		}

		MET_CorrT1Txy_et = patMETCorrT1Txy.et();
		MET_CorrT1Txy_pt = patMETCorrT1Txy.pt();
		MET_CorrT1Txy_phi = patMETCorrT1Txy.phi();
		MET_CorrT1Txy_sumET = patMETCorrT1Txy.sumEt();
		if (patMETCorrT1Txy.isPFMET()) {
			MET_CorrT1Txy_MuonEtFraction = patMETCorrT1Txy.MuonEtFraction();
			MET_CorrT1Txy_NeutralEMFraction = patMETCorrT1Txy.NeutralEMFraction();
			MET_CorrT1Txy_NeutralHadEtFraction = patMETCorrT1Txy.NeutralHadEtFraction();
			MET_CorrT1Txy_Type6EtFraction = patMETCorrT1Txy.Type6EtFraction();
			MET_CorrT1Txy_Type7EtFraction = patMETCorrT1Txy.Type7EtFraction();
		}

		MET_CorrCaloT1_et = patMETCorrCaloT1.et();
		MET_CorrCaloT1_pt = patMETCorrCaloT1.pt();
		MET_CorrCaloT1_phi = patMETCorrCaloT1.phi();
		MET_CorrCaloT1_sumET = patMETCorrCaloT1.sumEt();

		MET_CorrCaloT1T2_et = patMETCorrCaloT1T2.et();
		MET_CorrCaloT1T2_pt = patMETCorrCaloT1T2.pt();
		MET_CorrCaloT1T2_phi = patMETCorrCaloT1T2.phi();
		MET_CorrCaloT1T2_sumET = patMETCorrCaloT1T2.sumEt();

		if(iEvent.isRealData() || Embedded_){
			MET_Type1Corr_et = -1;
			MET_Type1Corr_pt = -1;
			MET_Type1Corr_phi = -1;
			MET_Type1Corr_sumET = -1;
			MET_Type1Corr_MuonEtFraction = -1;
			MET_Type1Corr_NeutralEMFraction = -1;
			MET_Type1Corr_NeutralHadEtFraction = -1;
			MET_Type1Corr_Type6EtFraction = -1;
			MET_Type1Corr_Type7EtFraction = -1;

			MET_Type1p2Corr_et = -1;
			MET_Type1p2Corr_pt = -1;
			MET_Type1p2Corr_phi = -1;
			MET_Type1p2Corr_sumET = -1;
			MET_Type1p2Corr_MuonEtFraction = -1;
			MET_Type1p2Corr_NeutralEMFraction = -1;
			MET_Type1p2Corr_NeutralHadEtFraction = -1;
			MET_Type1p2Corr_Type6EtFraction = -1;
			MET_Type1p2Corr_Type7EtFraction = -1;

			MET_Type1CorrElectronUp_et = -1;
			MET_Type1CorrElectronDown_et = -1;
			MET_Type1CorrMuonUp_et = -1;
			MET_Type1CorrMuonDown_et = -1;
			MET_Type1CorrTauUp_et = -1;
			MET_Type1CorrTauDown_et = -1;
			MET_Type1CorrJetResUp_et = -1;
			MET_Type1CorrJetResDown_et = -1;
			MET_Type1CorrJetEnUp_et = -1;
			MET_Type1CorrJetEnDown_et = -1;
			MET_Type1CorrUnclusteredUp_et = -1;
			MET_Type1CorrUnclusteredDown_et = -1;

			MET_Type1p2CorrElectronUp_et = -1;
			MET_Type1p2CorrElectronDown_et = -1;
			MET_Type1p2CorrMuonUp_et = -1;
			MET_Type1p2CorrMuonDown_et = -1;
			MET_Type1p2CorrTauUp_et = -1;
			MET_Type1p2CorrTauDown_et = -1;
			MET_Type1p2CorrJetResUp_et = -1;
			MET_Type1p2CorrJetResDown_et = -1;
			MET_Type1p2CorrJetEnUp_et = -1;
			MET_Type1p2CorrJetEnDown_et = -1;
			MET_Type1p2CorrUnclusteredUp_et = -1;
			MET_Type1p2CorrUnclusteredDown_et = -1;
		}else{
			///////////////////////////////
			//
			// MET uncertainties
			//
			// Type1      = T0 + T1 + Txy (recommendation by JetMET POG and Christian Veelken)
			// Type1Type2 = T0 + T1 + Txy + calibration for unclustered energy (better MET response, worse MET resolution)
			//
			edm::Handle<std::vector<pat::MET>> patMETType1CorrHandle;
			iEvent.getByLabel(patMETType1Corr_, patMETType1CorrHandle);
			pat::MET patMETType1Corr = patMETType1CorrHandle->front();

			edm::Handle<std::vector<pat::MET>> patMETType1p2CorrHandle;
			iEvent.getByLabel(patMETType1p2Corr_, patMETType1p2CorrHandle);
			pat::MET patMETType1p2Corr = patMETType1p2CorrHandle->front();

			edm::Handle<std::vector<pat::MET>> patMETType1CorrEleEnUpHandle;
			iEvent.getByLabel(patMETType1CorrEleEnUp_,patMETType1CorrEleEnUpHandle);
			pat::MET patMETType1CorrEleEnUp = patMETType1CorrEleEnUpHandle->front();

			edm::Handle<std::vector<pat::MET>> patMETType1CorrEleEnDownHandle;
			iEvent.getByLabel(patMETType1CorrEleEnDown_,patMETType1CorrEleEnDownHandle);
			pat::MET patMETType1CorrEleEnDown = patMETType1CorrEleEnDownHandle->front();

			edm::Handle<std::vector<pat::MET>> patMETType1CorrMuEnUpHandle;
			iEvent.getByLabel(patMETType1CorrMuEnUp_,patMETType1CorrMuEnUpHandle);
			pat::MET patMETType1CorrMuEnUp = patMETType1CorrMuEnUpHandle->front();

			edm::Handle<std::vector<pat::MET>> patMETType1CorrMuEnDownHandle;
			iEvent.getByLabel(patMETType1CorrMuEnDown_,patMETType1CorrMuEnDownHandle);
			pat::MET patMETType1CorrMuEnDown = patMETType1CorrMuEnDownHandle->front();

			edm::Handle<std::vector<pat::MET>> patMETType1CorrTauEnUpHandle;
			iEvent.getByLabel(patMETType1CorrTauEnUp_,patMETType1CorrTauEnUpHandle);
			pat::MET patMETType1CorrTauEnUp = patMETType1CorrTauEnUpHandle->front();

			edm::Handle<std::vector<pat::MET>> patMETType1CorrTauEnDownHandle;
			iEvent.getByLabel(patMETType1CorrTauEnDown_,patMETType1CorrTauEnDownHandle);
			pat::MET patMETType1CorrTauEnDown = patMETType1CorrTauEnDownHandle->front();

			edm::Handle<std::vector<pat::MET>> patMETType1CorrJetResUpHandle;
			iEvent.getByLabel(patMETType1CorrJetResUp_,patMETType1CorrJetResUpHandle);
			pat::MET patMETType1CorrJetResUp = patMETType1CorrJetResUpHandle->front();

			edm::Handle<std::vector<pat::MET>> patMETType1CorrJetResDownHandle;
			iEvent.getByLabel(patMETType1CorrJetResDown_,patMETType1CorrJetResDownHandle);
			pat::MET patMETType1CorrJetResDown = patMETType1CorrJetResDownHandle->front();

			edm::Handle<std::vector<pat::MET>> patMETType1CorrJetEnUpHandle;
			iEvent.getByLabel(patMETType1CorrJetEnUp_,patMETType1CorrJetEnUpHandle);
			pat::MET patMETType1CorrJetEnUp = patMETType1CorrJetEnUpHandle->front();

			edm::Handle<std::vector<pat::MET>> patMETType1CorrJetEnDownHandle;
			iEvent.getByLabel(patMETType1CorrJetEnDown_,patMETType1CorrJetEnDownHandle);
			pat::MET patMETType1CorrJetEnDown = patMETType1CorrJetEnDownHandle->front();

			edm::Handle<std::vector<pat::MET>> patMETType1CorrUnclusteredUpHandle;
			iEvent.getByLabel(patMETType1CorrUnclusteredUp_,patMETType1CorrUnclusteredUpHandle);
			pat::MET patMETType1CorrUnclusteredUp = patMETType1CorrUnclusteredUpHandle->front();

			edm::Handle<std::vector<pat::MET>> patMETType1CorrUnclusteredDownHandle;
			iEvent.getByLabel(patMETType1CorrUnclusteredDown_,patMETType1CorrUnclusteredDownHandle);
			pat::MET patMETType1CorrUnclusteredDown = patMETType1CorrUnclusteredDownHandle->front();

			edm::Handle<std::vector<pat::MET>> patMETType1p2CorrEleEnUpHandle;
			iEvent.getByLabel(patMETType1p2CorrEleEnUp_,patMETType1p2CorrEleEnUpHandle);
			pat::MET patMETType1p2CorrEleEnUp = patMETType1p2CorrEleEnUpHandle->front();

			edm::Handle<std::vector<pat::MET>> patMETType1p2CorrEleEnDownHandle;
			iEvent.getByLabel(patMETType1p2CorrEleEnDown_,patMETType1p2CorrEleEnDownHandle);
			pat::MET patMETType1p2CorrEleEnDown = patMETType1p2CorrEleEnDownHandle->front();

			edm::Handle<std::vector<pat::MET>> patMETType1p2CorrMuEnUpHandle;
			iEvent.getByLabel(patMETType1p2CorrMuEnUp_,patMETType1p2CorrMuEnUpHandle);
			pat::MET patMETType1p2CorrMuEnUp = patMETType1p2CorrMuEnUpHandle->front();

			edm::Handle<std::vector<pat::MET>> patMETType1p2CorrMuEnDownHandle;
			iEvent.getByLabel(patMETType1p2CorrMuEnDown_,patMETType1p2CorrMuEnDownHandle);
			pat::MET patMETType1p2CorrMuEnDown = patMETType1p2CorrMuEnDownHandle->front();

			edm::Handle<std::vector<pat::MET>> patMETType1p2CorrTauEnUpHandle;
			iEvent.getByLabel(patMETType1p2CorrTauEnUp_,patMETType1p2CorrTauEnUpHandle);
			pat::MET patMETType1p2CorrTauEnUp = patMETType1p2CorrTauEnUpHandle->front();

			edm::Handle<std::vector<pat::MET>> patMETType1p2CorrTauEnDownHandle;
			iEvent.getByLabel(patMETType1p2CorrTauEnDown_,patMETType1p2CorrTauEnDownHandle);
			pat::MET patMETType1p2CorrTauEnDown = patMETType1p2CorrTauEnDownHandle->front();

			edm::Handle<std::vector<pat::MET>> patMETType1p2CorrJetResUpHandle;
			iEvent.getByLabel(patMETType1p2CorrJetResUp_,patMETType1p2CorrJetResUpHandle);
			pat::MET patMETType1p2CorrJetResUp = patMETType1p2CorrJetResUpHandle->front();

			edm::Handle<std::vector<pat::MET>> patMETType1p2CorrJetResDownHandle;
			iEvent.getByLabel(patMETType1p2CorrJetResDown_,patMETType1p2CorrJetResDownHandle);
			pat::MET patMETType1p2CorrJetResDown = patMETType1p2CorrJetResDownHandle->front();

			edm::Handle<std::vector<pat::MET>> patMETType1p2CorrJetEnUpHandle;
			iEvent.getByLabel(patMETType1p2CorrJetEnUp_,patMETType1p2CorrJetEnUpHandle);
			pat::MET patMETType1p2CorrJetEnUp = patMETType1p2CorrJetEnUpHandle->front();

			edm::Handle<std::vector<pat::MET>> patMETType1p2CorrJetEnDownHandle;
			iEvent.getByLabel(patMETType1p2CorrJetEnDown_,patMETType1p2CorrJetEnDownHandle);
			pat::MET patMETType1p2CorrJetEnDown = patMETType1p2CorrJetEnDownHandle->front();

			edm::Handle<std::vector<pat::MET>> patMETType1p2CorrUnclusteredUpHandle;
			iEvent.getByLabel(patMETType1p2CorrUnclusteredUp_,patMETType1p2CorrUnclusteredUpHandle);
			pat::MET patMETType1p2CorrUnclusteredUp = patMETType1p2CorrUnclusteredUpHandle->front();

			edm::Handle<std::vector<pat::MET>> patMETType1p2CorrUnclusteredDownHandle;
			iEvent.getByLabel(patMETType1p2CorrUnclusteredDown_,patMETType1p2CorrUnclusteredDownHandle);
			pat::MET patMETType1p2CorrUnclusteredDown = patMETType1p2CorrUnclusteredDownHandle->front();

			MET_Type1Corr_et = patMETType1Corr.et();
			MET_Type1Corr_pt = patMETType1Corr.pt();
			MET_Type1Corr_phi = patMETType1Corr.phi();
			MET_Type1Corr_sumET = patMETType1Corr.sumEt();
			if(patMETType1Corr.isPFMET()){
				MET_Type1Corr_MuonEtFraction = patMETType1Corr.MuonEtFraction();
				MET_Type1Corr_NeutralEMFraction = patMETType1Corr.NeutralEMFraction();
				MET_Type1Corr_NeutralHadEtFraction = patMETType1Corr.NeutralHadEtFraction();
				MET_Type1Corr_Type6EtFraction = patMETType1Corr.Type6EtFraction();
				MET_Type1Corr_Type7EtFraction = patMETType1Corr.Type7EtFraction();
			}

			MET_Type1p2Corr_et = patMETType1p2Corr.et();
			MET_Type1p2Corr_pt = patMETType1p2Corr.pt();
			MET_Type1p2Corr_phi = patMETType1p2Corr.phi();
			MET_Type1p2Corr_sumET = patMETType1p2Corr.sumEt();
			if(patMETType1p2Corr.isPFMET()){
				MET_Type1p2Corr_MuonEtFraction = patMETType1p2Corr.MuonEtFraction();
				MET_Type1p2Corr_NeutralEMFraction = patMETType1p2Corr.NeutralEMFraction();
				MET_Type1p2Corr_NeutralHadEtFraction = patMETType1p2Corr.NeutralHadEtFraction();
				MET_Type1p2Corr_Type6EtFraction = patMETType1p2Corr.Type6EtFraction();
				MET_Type1p2Corr_Type7EtFraction = patMETType1p2Corr.Type7EtFraction();
			}

			MET_Type1CorrElectronUp_et = patMETType1CorrEleEnUp.et();
			MET_Type1CorrElectronDown_et = patMETType1CorrEleEnDown.et();
			MET_Type1CorrMuonUp_et = patMETType1CorrMuEnUp.et();
			MET_Type1CorrMuonDown_et = patMETType1CorrMuEnDown.et();
			MET_Type1CorrTauUp_et = patMETType1CorrTauEnUp.et();
			MET_Type1CorrTauDown_et = patMETType1CorrTauEnDown.et();
			MET_Type1CorrJetResUp_et = patMETType1CorrJetResUp.et();
			MET_Type1CorrJetResDown_et = patMETType1CorrJetResDown.et();
			MET_Type1CorrJetEnUp_et = patMETType1CorrJetEnUp.et();
			MET_Type1CorrJetEnDown_et = patMETType1CorrJetEnDown.et();
			MET_Type1CorrUnclusteredUp_et = patMETType1CorrUnclusteredUp.et();
			MET_Type1CorrUnclusteredDown_et = patMETType1CorrUnclusteredDown.et();

			MET_Type1p2CorrElectronUp_et = patMETType1p2CorrEleEnUp.et();
			MET_Type1p2CorrElectronDown_et = patMETType1p2CorrEleEnDown.et();
			MET_Type1p2CorrMuonUp_et = patMETType1p2CorrMuEnUp.et();
			MET_Type1p2CorrMuonDown_et = patMETType1p2CorrMuEnDown.et();
			MET_Type1p2CorrTauUp_et = patMETType1p2CorrTauEnUp.et();
			MET_Type1p2CorrTauDown_et = patMETType1p2CorrTauEnDown.et();
			MET_Type1p2CorrJetResUp_et = patMETType1p2CorrJetResUp.et();
			MET_Type1p2CorrJetResDown_et = patMETType1p2CorrJetResDown.et();
			MET_Type1p2CorrJetEnUp_et = patMETType1p2CorrJetEnUp.et();
			MET_Type1p2CorrJetEnDown_et = patMETType1p2CorrJetEnDown.et();
			MET_Type1p2CorrUnclusteredUp_et = patMETType1p2CorrUnclusteredUp.et();
			MET_Type1p2CorrUnclusteredDown_et = patMETType1p2CorrUnclusteredDown.et();
		}

		if (doMVAMET_) {
			edm::Handle<std::vector<pat::MET>> patMETCorrMVAHandle;
			iEvent.getByLabel(patMETCorrMVA_, patMETCorrMVAHandle);
			pat::MET patMETCorrMVA = patMETCorrMVAHandle->front();

			edm::Handle<reco::MuonCollection> mvametMuonCollection;
			iEvent.getByLabel(muonsForPfMetCorrMVA_, mvametMuonCollection);
			edm::Handle<reco::GsfElectronCollection> mvametElectronCollection;
			iEvent.getByLabel(elecsForPfMetCorrMVA_, mvametElectronCollection);
			edm::Handle<std::vector<reco::PFTau> > mvametTauCollection;
			iEvent.getByLabel(tausForPfMetCorrMVA_, mvametTauCollection);

			MET_CorrMVA_et = patMETCorrMVA.et();
			MET_CorrMVA_pt = patMETCorrMVA.pt();
			MET_CorrMVA_phi = patMETCorrMVA.phi();
			MET_CorrMVA_sumET = patMETCorrMVA.sumEt();
			MET_CorrMVA_significance = patMETCorrMVA.significance();
			sigMat = patMETCorrMVA.getSignificanceMatrix();
			if (sigMat(0, 1) != sigMat(1, 0))
				std::cout << "WARNING: MET significance matrix not symmetric" << std::endl;
			MET_CorrMVA_significance_xx = sigMat(0, 0);
			MET_CorrMVA_significance_xy = sigMat(0, 1);
			MET_CorrMVA_significance_yy = sigMat(1, 1);
			MET_CorrMVA_MuonEtFraction = patMETCorrMVA.MuonEtFraction();
			MET_CorrMVA_NeutralEMFraction = patMETCorrMVA.NeutralEMFraction();
			MET_CorrMVA_NeutralHadEtFraction = patMETCorrMVA.NeutralHadEtFraction();
			MET_CorrMVA_Type6EtFraction = patMETCorrMVA.Type6EtFraction();
			MET_CorrMVA_Type7EtFraction = patMETCorrMVA.Type7EtFraction();

			// store p4 of leptons used for MVA-MET calculation
			int Muon_index = 0;
			for (reco::MuonCollection::const_iterator iMuon = mvametMuonCollection->begin(); iMuon != mvametMuonCollection->end(); ++iMuon, Muon_index++) {
				reco::MuonRef RefMuon(mvametMuonCollection, Muon_index);
				std::vector<double> iMuon_p4;
				iMuon_p4.push_back(RefMuon->p4().E());
				iMuon_p4.push_back(RefMuon->p4().Px());
				iMuon_p4.push_back(RefMuon->p4().Py());
				iMuon_p4.push_back(RefMuon->p4().Pz());
				MET_CorrMVA_srcMuon_p4.push_back(iMuon_p4);
			}
			for (reco::PFCandidateCollection::size_type iElectron = 0; iElectron < mvametElectronCollection->size(); iElectron++) {
				reco::GsfElectronRef RefElectron(mvametElectronCollection, iElectron);
				std::vector<double> iElectron_p4;
				iElectron_p4.push_back(RefElectron->p4().E());
				iElectron_p4.push_back(RefElectron->p4().Px());
				iElectron_p4.push_back(RefElectron->p4().Py());
				iElectron_p4.push_back(RefElectron->p4().Pz());
				MET_CorrMVA_srcElectron_p4.push_back(iElectron_p4);
			}
			for (unsigned iPFTau = 0; iPFTau < mvametTauCollection->size(); ++iPFTau) {
				reco::PFTauRef HPStauCandidate(mvametTauCollection, iPFTau);
				std::vector<double> iPFTau_p4;
				iPFTau_p4.push_back(HPStauCandidate->p4().E());
				iPFTau_p4.push_back(HPStauCandidate->p4().Px());
				iPFTau_p4.push_back(HPStauCandidate->p4().Py());
				iPFTau_p4.push_back(HPStauCandidate->p4().Pz());
				MET_CorrMVA_srcTau_p4.push_back(iPFTau_p4);
			}

			edm::Handle<std::vector<pat::MET>> patMETCorrMVAMuTauHandle;
			iEvent.getByLabel(patMETCorrMVAMuTau_, patMETCorrMVAMuTauHandle);
			pat::MET patMETCorrMVAMuTau = patMETCorrMVAMuTauHandle->front();

			edm::Handle<reco::MuonCollection> mvamutaumetMuonCollection;
			iEvent.getByLabel(muonsForPfMetCorrMVAMuTau_, mvamutaumetMuonCollection);
			edm::Handle<std::vector<reco::PFTau> > mvamutaumetTauCollection;
			iEvent.getByLabel(tausForPfMetCorrMVAMuTau_, mvamutaumetTauCollection);

			MET_CorrMVAMuTau_et = patMETCorrMVAMuTau.et();
			MET_CorrMVAMuTau_pt = patMETCorrMVAMuTau.pt();
			MET_CorrMVAMuTau_phi = patMETCorrMVAMuTau.phi();
			MET_CorrMVAMuTau_sumET = patMETCorrMVAMuTau.sumEt();
			MET_CorrMVAMuTau_significance = patMETCorrMVAMuTau.significance();
			sigMat = patMETCorrMVAMuTau.getSignificanceMatrix();
			if (sigMat(0, 1) != sigMat(1, 0))
				std::cout << "WARNING: MET significance matrix not symmetric" << std::endl;
			MET_CorrMVAMuTau_significance_xx = sigMat(0, 0);
			MET_CorrMVAMuTau_significance_xy = sigMat(0, 1);
			MET_CorrMVAMuTau_significance_yy = sigMat(1, 1);
			MET_CorrMVAMuTau_MuonEtFraction = patMETCorrMVAMuTau.MuonEtFraction();
			MET_CorrMVAMuTau_NeutralEMFraction = patMETCorrMVAMuTau.NeutralEMFraction();
			MET_CorrMVAMuTau_NeutralHadEtFraction = patMETCorrMVAMuTau.NeutralHadEtFraction();
			MET_CorrMVAMuTau_Type6EtFraction = patMETCorrMVAMuTau.Type6EtFraction();
			MET_CorrMVAMuTau_Type7EtFraction = patMETCorrMVAMuTau.Type7EtFraction();

			// store p4 of leptons used for MVA-MET calculation
			Muon_index = 0;
			for (reco::MuonCollection::const_iterator iMuon = mvamutaumetMuonCollection->begin(); iMuon != mvamutaumetMuonCollection->end(); ++iMuon, Muon_index++) {
				reco::MuonRef RefMuon(mvamutaumetMuonCollection, Muon_index);
				std::vector<double> iMuon_p4;
				iMuon_p4.push_back(RefMuon->p4().E());
				iMuon_p4.push_back(RefMuon->p4().Px());
				iMuon_p4.push_back(RefMuon->p4().Py());
				iMuon_p4.push_back(RefMuon->p4().Pz());
				MET_CorrMVAMuTau_srcMuon_p4.push_back(iMuon_p4);
			}
			for (unsigned iPFTau = 0; iPFTau < mvamutaumetTauCollection->size(); ++iPFTau) {
				reco::PFTauRef HPStauCandidate(mvamutaumetTauCollection, iPFTau);
				std::vector<double> iPFTau_p4;
				iPFTau_p4.push_back(HPStauCandidate->p4().E());
				iPFTau_p4.push_back(HPStauCandidate->p4().Px());
				iPFTau_p4.push_back(HPStauCandidate->p4().Py());
				iPFTau_p4.push_back(HPStauCandidate->p4().Pz());
				MET_CorrMVAMuTau_srcTau_p4.push_back(iPFTau_p4);
			}
		}
	}
}

void TauNtuple::fillTriggerInfo(edm::Event& iEvent, const edm::EventSetup& iSetup) {
	if (!TriggerOK)
		return;
	edm::Handle<trigger::TriggerEvent> triggerEvent;
	iEvent.getByLabel(TriggerEvent_, triggerEvent);
	edm::Handle<edm::TriggerResults> triggerResults;
	iEvent.getByLabel(TriggerResults_, triggerResults);
	for (unsigned int i = 0; i < MyTriggerInfoNames.size(); i++) {
		HTLTriggerName.push_back(MyTriggerInfoNames.at(i));
		unsigned int triggerIndex = hltConfig_.triggerIndex(HTLTriggerName.at(i));
		TriggerAccept.push_back(triggerResults->accept(triggerIndex));
		TriggerError.push_back(triggerResults->error(triggerIndex));
		TriggerWasRun.push_back(triggerResults->wasrun(triggerIndex));
		////////////////////////////////////////////
		// now get level 1 & HLT prescale
		const std::vector<std::pair<bool, std::string> > level1Seeds = hltConfig_.hltL1GTSeeds(HTLTriggerName.at(i));
		int l1Prescale(-1);
		L1GtUtils l1GtUtils;
		l1GtUtils.retrieveL1EventSetup(iSetup);
		bool isTechbit = false;
		if (level1Seeds.size() == 1) {
			std::vector<std::string> myl1SeedPaths;
			std::stringstream ss(level1Seeds.at(0).second);
			TString buffer = level1Seeds.at(0).second;
			myl1SeedPaths.clear();
			if (!(buffer.Contains("(") || buffer.Contains(")") || buffer.Contains("AND") || buffer.Contains("NOT"))) {
				while (ss.good() && !ss.eof()) {
					ss >> buffer;
					if (!buffer.Contains("OR"))
						myl1SeedPaths.push_back(buffer.Data());
				}
			}
			if (!myl1SeedPaths.empty()) {
				for (unsigned j = 0; j < myl1SeedPaths.size(); j++) {
					int l1TempPrescale(-1);
					int errorCode(0);
					if (level1Seeds.at(0).first) { // technical triggers
						isTechbit = true;
						unsigned int techBit(atoi(myl1SeedPaths.at(j).c_str()));
						const std::string techName(*(triggerMenuLite_->gtTechTrigName(techBit, errorCode)));
						if (errorCode != 0)
							continue;
						if (!l1GtUtils.decision(iEvent, techName, errorCode))
							continue;
						if (errorCode != 0)
							continue;
						l1TempPrescale = l1GtUtils.prescaleFactor(iEvent, techName, errorCode);
						if (errorCode != 0)
							continue;
					} else { // algorithmic triggers
						if (!l1GtUtils.decision(iEvent, myl1SeedPaths.at(j), errorCode))
							continue;
						if (errorCode != 0)
							continue;
						l1TempPrescale = l1GtUtils.prescaleFactor(iEvent, myl1SeedPaths.at(j), errorCode);
						if (errorCode != 0)
							continue;
					}
					if (l1TempPrescale > 0) {
						if (l1Prescale == -1 || l1Prescale > l1TempPrescale)
							l1Prescale = l1TempPrescale;
					}
				}
			}
		}
		HLTPrescale.push_back(hltConfig_.prescaleValue(iEvent, iSetup, HTLTriggerName.at(i)));
		NHLTL1GTSeeds.push_back(level1Seeds.size());
		if (l1Prescale == -1) {
			L1SEEDPrescale.push_back(1);
			L1SEEDInvalidPrescale.push_back(true);
			L1SEEDisTechBit.push_back(isTechbit);
		} else {
			L1SEEDPrescale.push_back((unsigned int) l1Prescale);
			L1SEEDInvalidPrescale.push_back(false);
			L1SEEDisTechBit.push_back(isTechbit);
		}
		////////////////////////////////////
		// Now get Trigger matching
		if (triggerResults->accept(triggerIndex)) {
			std::string filterName_ = "";
			edm::InputTag filterTag;
			std::vector<std::string> filters = hltConfig_.moduleLabels(HTLTriggerName.at(i));
			/*std::cout << "Module names" << std::endl;
			 for(unsigned int ilabel=0; ilabel<filters.size(); ilabel++){
			 std::cout << filters.at(ilabel) << std::endl;
			 }*/
			for (std::vector<std::string>::iterator filter = filters.begin(); filter != filters.end(); ++filter) {
				edm::InputTag testTag(*filter, "", "HLT");
				int testindex = triggerEvent->filterIndex(testTag);
				if (!(testindex >= triggerEvent->sizeFilters())) {
					filterName_ = *filter;
					filterTag = testTag;
				}
			}

			unsigned int index = triggerEvent->filterIndex(filterTag);
			/*std::cout << "TrgPath: " << HTLTriggerName.at(i) << " hltTag_.label(): "
			 << filterTag.label() << "   filter name: "
			 << filterName_ << "  sizeFilters: "
			 << triggerEvent->sizeFilters() << std::endl;*/

			std::vector<float> match;
			// Muons
			edm::Handle<reco::MuonCollection> muonCollection;
			iEvent.getByLabel(muonsTag_, muonCollection);
			TriggerMatch(triggerEvent, index, muonCollection, TriggerMuonMatchingdr_, match);
			TriggerMatchMuon.push_back(match);
			match.clear();
			// Jets
			edm::Handle<reco::PFJetCollection> JetCollection;
			iEvent.getByLabel(pfjetsTag_, JetCollection);
			TriggerMatch(triggerEvent, index, JetCollection, TriggerJetMatchingdr_, match);
			TriggerMatchJet.push_back(match);
			match.clear();
			// Taus
			edm::Handle<std::vector<reco::PFTau> > PFTauCollection;
			iEvent.getByLabel(hpsTauProducer_, PFTauCollection);
			// 	 edm::Handle<reco::PFJetCollection> tauCollection;
			// 	 iEvent.getByLabel(pfjetsTag_, tauCollection);
			TriggerMatch(triggerEvent, index, PFTauCollection, TriggerTauMatchingdr_, match);
			TriggerMatchTau.push_back(match);
			match.clear();

			// Save trigger objects
			std::vector<float> TriggerObj_Pt;
			std::vector<float> TriggerObj_Eta;
			std::vector<float> TriggerObj_Phi;
			std::vector<float> TriggerObj_E;
			std::vector<int> TriggerObj_Id;
			std::vector<trigger::TriggerObject> trgobjs = triggerEvent->getObjects();
			// old version
			/*const trigger::Keys& KEYS(triggerEvent->filterKeys(index));
			 for (unsigned int ipart = 0; ipart < KEYS.size(); ipart++) {
			 TriggerObj_Pt.push_back(trgobjs.at(KEYS.at(ipart)).pt());
			 TriggerObj_Eta.push_back(trgobjs.at(KEYS.at(ipart)).eta());
			 TriggerObj_Phi.push_back(trgobjs.at(KEYS.at(ipart)).phi());
			 TriggerObj_E.push_back(trgobjs.at(KEYS.at(ipart)).energy());
			 std::cout << trgobjs.at(KEYS.at(ipart)).pt() << std::endl;
			 }*/
			// couple trigger objects to their trigger path and save their information *** new version ***
			//std::cout << MyTriggerInfoNames.at(i) << std::endl;
			if (triggerEvent.isValid()) {
				for (unsigned int imodule = 0; imodule < useFilterModules_.size(); imodule++) {
					for (unsigned int ifilter = 0; ifilter < filters.size(); ifilter++) {
						if (filters.at(ifilter) == useFilterModules_.at(imodule)) {
							//std::cout << useFilterModules_.at(imodule) << std::endl;
							edm::InputTag tagEv(useFilterModules_.at(imodule), "", "HLT");
							int ID = triggerEvent->filterIndex(tagEv);
							if (ID != triggerEvent->sizeFilters()) {
								const trigger::Vids& ids(triggerEvent->filterIds(ID));
								const trigger::Keys& keys(triggerEvent->filterKeys(ID));
								//std::cout << "#vids = " << ids.size() << std::endl;
								for (unsigned int ivid = 0; ivid < ids.size(); ivid++) {
									//std::cout << "vid: " << ids.at(ivid) << std::endl;

									// Trigger objects have Ids in range [+81,+96]
									if (ids.at(ivid) >= 81 && ids.at(ivid) <= 96) {
										TriggerObj_Pt.push_back(trgobjs.at(keys.at(ivid)).pt());
										TriggerObj_Eta.push_back(trgobjs.at(keys.at(ivid)).eta());
										TriggerObj_Phi.push_back(trgobjs.at(keys.at(ivid)).phi());
										TriggerObj_E.push_back(trgobjs.at(keys.at(ivid)).energy());
										TriggerObj_Id.push_back(ids.at(ivid));
									}
								}
							}
						}
					}
				}
			} else {
				std::cout << "triggerEvent NOT valid" << std::endl;
			}
			//
			HLTTrigger_objs_Pt.push_back(TriggerObj_Pt);
			HLTTrigger_objs_Eta.push_back(TriggerObj_Eta);
			HLTTrigger_objs_Phi.push_back(TriggerObj_Phi);
			HLTTrigger_objs_E.push_back(TriggerObj_E);
			HLTTrigger_objs_Id.push_back(TriggerObj_Id);
			HLTTrigger_objs_trigger.push_back(MyTriggerInfoNames.at(i));
		}else{
		 TriggerMatchMuon.push_back(std::vector<float>());
		 TriggerMatchJet.push_back(std::vector<float>());
		 TriggerMatchTau.push_back(std::vector<float>());
		 HLTTrigger_objs_Pt.push_back(std::vector<float>());
		 HLTTrigger_objs_Eta.push_back(std::vector<float>());
		 HLTTrigger_objs_Phi.push_back(std::vector<float>());
		 HLTTrigger_objs_E.push_back(std::vector<float>());
		 HLTTrigger_objs_Id.push_back(std::vector<int>());
		 HLTTrigger_objs_trigger.push_back(std::string());
		 }
		////////////////////////////////////
		// Now do L1 TriggerSeeds if requested
		if (doL1Triggers_) {
			for (unsigned j = 0; j < l1TriggerNames_.size(); j++) {
				int errorCode(0);
				L1TriggerName.push_back(l1TriggerNames_.at(j));
				L1TriggerDecision.push_back(l1GtUtils.decision(iEvent, l1TriggerNames_.at(j), errorCode));
				L1ErrorCode.push_back(errorCode);
				L1Prescale.push_back(l1GtUtils.prescaleFactor(iEvent, l1TriggerNames_.at(j), errorCode));
			}
		}else{
			L1TriggerName.push_back(std::string());
			L1TriggerDecision.push_back(false);
			L1ErrorCode.push_back(-1);
			L1Prescale.push_back(1);
		}
	}
}

template<class T>
void TauNtuple::TriggerMatch(edm::Handle<trigger::TriggerEvent> &triggerEvent, unsigned int triggerIndex, T obj, double drmax, std::vector<float> &match) {
	match = std::vector<float>(obj->size(), 999);
	std::vector<trigger::TriggerObject> trgobjs = triggerEvent->getObjects();
	const trigger::Keys& KEYS(triggerEvent->filterKeys(triggerIndex));
	for (unsigned int ipart = 0; ipart < KEYS.size(); ipart++) {
		for (unsigned int i = 0; i < obj->size(); ++i) {
			double dr = reco::deltaR(trgobjs.at(KEYS.at(ipart)).eta(), trgobjs.at(KEYS.at(ipart)).phi(), obj->at(i).eta(), obj->at(i).phi());
			if (dr < drmax) {
				match.at(i) = dr;
				/*	 std::cout << "Found Trigger Match " << i << " " << dr << " " << drmax << " Trigger Obj: " << trgobjs.at(KEYS.at(ipart)).eta()
				 << " " << trgobjs.at(KEYS.at(ipart)).phi() << " " <<  trgobjs.at(KEYS.at(ipart)).energy()
				 << " obj: " <<  obj->at(i).eta() << " " << obj->at(i).phi() << " " << obj->at(i).energy() << std::endl;*/
			}
		}
	}
}

void TauNtuple::fillEventInfo(edm::Event& iEvent, const edm::EventSetup& iSetup) {

	Event_EventNumber = iEvent.id().event();
	Event_RunNumber = iEvent.id().run();
	Event_bunchCrossing = iEvent.bunchCrossing();
	Event_orbitNumber = iEvent.orbitNumber();
	Event_luminosityBlock = iEvent.luminosityBlock();
	Event_isRealData = iEvent.isRealData();

	reco::BeamSpot beamSpot;
	edm::Handle<reco::BeamSpot> beamSpotHandle;
	iEvent.getByLabel(beamSpotTag_, beamSpotHandle);
	if (beamSpotHandle.isValid()) {
		beamSpot = *beamSpotHandle;
		beamspot_par.push_back(beamSpot.x0());
		beamspot_par.push_back(beamSpot.y0());
		beamspot_par.push_back(beamSpot.z0());
		beamspot_par.push_back(beamSpot.sigmaZ());
		beamspot_par.push_back(beamSpot.dxdz());
		beamspot_par.push_back(beamSpot.dydz());
		beamspot_par.push_back(beamSpot.BeamWidthX());
		for (unsigned int i = 0; i < reco::BeamSpot::dimension; i++) {
			for (unsigned int j = i; j < reco::BeamSpot::dimension; j++) {
				beamspot_cov.push_back(beamSpot.covariance(i, j));
			}
		}
		beamspot_emittanceX = beamSpot.emittanceX();
		beamspot_emittanceY = beamSpot.emittanceY();
		beamspot_betaStar = beamSpot.betaStar();
	} else {
		for (unsigned int i = 0; i < reco::BeamSpot::dimension; i++) {
			beamspot_par.push_back(-999);
			for (unsigned int j = i; j < reco::BeamSpot::dimension; j++) {
				beamspot_cov.push_back(-999);
			}
		}
		beamspot_emittanceX = -999;
		beamspot_emittanceY = -999;
		beamspot_betaStar = -999;
	}

	// Add Embedding info
	if (Event_isRealData && Embedded_)
		Event_isRealData = false;
	if (!Event_isRealData && !Embedded_) {
		edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
		iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
		std::vector<PileupSummaryInfo>::const_iterator PVI;
		PileupInfo_TrueNumInteractions_nm1 = -1;
		PileupInfo_TrueNumInteractions_n0 = -1;
		PileupInfo_TrueNumInteractions_np1 = -1;
		for (PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
			int BX = PVI->getBunchCrossing();
			if (BX == -1)
				PileupInfo_TrueNumInteractions_nm1 = PVI->getTrueNumInteractions();
			if (BX == 0)
				PileupInfo_TrueNumInteractions_n0 = PVI->getTrueNumInteractions();
			if (BX == 1)
				PileupInfo_TrueNumInteractions_np1 = PVI->getTrueNumInteractions();
		}
		PUWeight = LumiWeights_.weight( PileupInfo_TrueNumInteractions_n0 );
		PUWeight3D = LumiWeights3D_.weight3D(PileupInfo_TrueNumInteractions_nm1, PileupInfo_TrueNumInteractions_n0, PileupInfo_TrueNumInteractions_np1);
		if(PUInputHistoData_p5_ != ""){
			PUWeight_p5 = LumiWeights_p5_.weight( PileupInfo_TrueNumInteractions_n0 );
			PUWeight3D_p5 = LumiWeights3D_p5_.weight3D(PileupInfo_TrueNumInteractions_nm1, PileupInfo_TrueNumInteractions_n0, PileupInfo_TrueNumInteractions_np1);
		}
		else{
			PUWeight_p5 = -10;
			PUWeight3D_p5 = -10;
		}
		if(PUInputHistoData_m5_ != ""){
			PUWeight_m5 = LumiWeights_m5_.weight( PileupInfo_TrueNumInteractions_n0 );
			PUWeight3D_m5 = LumiWeights3D_m5_.weight3D(PileupInfo_TrueNumInteractions_nm1, PileupInfo_TrueNumInteractions_n0, PileupInfo_TrueNumInteractions_np1);
		}
		else{
			PUWeight_m5 = -10;
			PUWeight3D_m5 = -10;
		}
		if((PUInputHistoMCFineBins_ != "") && (PUInputHistoDataFineBins_ != "")){
			PUWeightFineBins = LumiWeightsFineBinning_.weight( PileupInfo_TrueNumInteractions_n0 );
		}
		else{
			PUWeightFineBins = -10;
		}


	}
	if (!Embedded_) {
		TauSpinnerWeight = 1.;
		SelEffWeight = 1.;
		MinVisPtFilter = 1.;
		KinWeightPt = 1.;
		KinWeightEta = 1.;
		KinWeightMassPt = 1.;
	} else {
		edm::InputTag tauspinner("TauSpinnerReco", "TauSpinnerWT");
		edm::Handle<double> TauSpinnerRecoHandle;
		iEvent.getByLabel(tauspinner, TauSpinnerRecoHandle);
		TauSpinnerWeight = *TauSpinnerRecoHandle;

		edm::InputTag seleffweight("ZmumuEvtSelEffCorrWeightProducer", "weight");
		edm::Handle<double> ZmumuEvtSelEffCorrProducerHandle;
		iEvent.getByLabel(seleffweight, ZmumuEvtSelEffCorrProducerHandle);
		SelEffWeight = *ZmumuEvtSelEffCorrProducerHandle;

		edm::InputTag generator("generator", "minVisPtFilter");
		edm::Handle<GenFilterInfo> generatorHandle;
		iEvent.getByLabel(generator, generatorHandle);
		MinVisPtFilter = generatorHandle->filterEfficiency();

		edm::InputTag kinweightpt("embeddingKineReweightRECembedding", "genTau2PtVsGenTau1Pt");
		edm::Handle<double> embeddingKineReweightRECembeddingPtHandle;
		iEvent.getByLabel(kinweightpt, embeddingKineReweightRECembeddingPtHandle);
		KinWeightPt = *embeddingKineReweightRECembeddingPtHandle;

		edm::InputTag kinweighteta("embeddingKineReweightRECembedding", "genTau2EtaVsGenTau1Eta");
		edm::Handle<double> embeddingKineReweightRECembeddingEtaHandle;
		iEvent.getByLabel(kinweighteta, embeddingKineReweightRECembeddingEtaHandle);
		KinWeightEta = *embeddingKineReweightRECembeddingEtaHandle;

		edm::InputTag kinweightmasspt("embeddingKineReweightRECembedding", "genDiTauMassVsGenDiTauPt");
		edm::Handle<double> embeddingKineReweightRECembeddingMassPtHandle;
		iEvent.getByLabel(kinweightmasspt, embeddingKineReweightRECembeddingEtaHandle);
		KinWeightMassPt = *embeddingKineReweightRECembeddingEtaHandle;
	}
	EmbeddedWeight = TauSpinnerWeight * SelEffWeight * MinVisPtFilter * KinWeightPt * KinWeightEta * KinWeightMassPt;
	if (EmbeddedWeight != TauSpinnerWeight * SelEffWeight * MinVisPtFilter * KinWeightPt * KinWeightEta * KinWeightMassPt) {
		std::cout << "!!! Calculation of embedding weights faulty. Check your code !!!" << std::endl;
	}
}

void TauNtuple::beginJob() {

	std::cout << "----------------------------------- >>>>>>>>>>>>>> TauNtuple begin Job" << std::endl;
	//-------------------------
	//   TString cmd1="pwd";
	//   TString cmd2="ls";
	//   TString cmd3="ls ../";
	//   TString cmd4="ls */";
	//   system(cmd1.Data());
	//   system(cmd2.Data());
	//   system(cmd3.Data());
	//   system(cmd4.Data());
	//-------------------------
	cnt_ = 0;
	output = new TFile("TauNtuple.root", "RECREATE");
	output_tree = new TTree("t", "t");

	output_tree->Branch("DataMC_Type", &DataMC_Type_idx);

	output_tree->Branch("beamspot_par", &beamspot_par);
	output_tree->Branch("beamspot_cov", &beamspot_cov);
	output_tree->Branch("beamspot_emittanceX", &beamspot_emittanceX);
	output_tree->Branch("beamspot_emittanceY", &beamspot_emittanceY);
	output_tree->Branch("beamspot_betaStar", &beamspot_betaStar);

	//=============  Vertex Block ====
	output_tree->Branch("Vtx_chi2", &Vtx_chi2);
	output_tree->Branch("Vtx_nTrk", &Vtx_nTrk);
	output_tree->Branch("Vtx_ndof", &Vtx_ndof);
	output_tree->Branch("Vtx_x", &Vtx_x);
	output_tree->Branch("Vtx_y", &Vtx_y);
	output_tree->Branch("Vtx_z", &Vtx_z);
	output_tree->Branch("Vtx_Cov", &Vtx_Cov);
	output_tree->Branch("Vtx_Track_idx", &Vtx_Track_idx);
	output_tree->Branch("Vtx_Track_Weights", &Vtx_Track_Weights);
	output_tree->Branch("Vtx_isFake", &Vtx_isFake);
	output_tree->Branch("Vtx_TracksP4", &Vtx_TracksP4);

	//=============  Muon Block ====
	output_tree->Branch("isPatMuon", &doPatMuons_);
	output_tree->Branch("Muon_p4", &Muon_p4);
	output_tree->Branch("Muon_Poca", &Muon_Poca);
	output_tree->Branch("Muon_isGlobalMuon", &Muon_isGlobalMuon);
	output_tree->Branch("Muon_isStandAloneMuon", &Muon_isStandAloneMuon);
	output_tree->Branch("Muon_isTrackerMuon", &Muon_isTrackerMuon);
	output_tree->Branch("Muon_isCaloMuon", &Muon_isCaloMuon);
	output_tree->Branch("Muon_isIsolationValid", &Muon_isIsolationValid);
	output_tree->Branch("Muon_isQualityValid", &Muon_isQualityValid);
	output_tree->Branch("Muon_isTimeValid", &Muon_isTimeValid);
	output_tree->Branch("Muon_emEt03", &Muon_emEt03);
	output_tree->Branch("Muon_emVetoEt03", &Muon_emVetoEt03);
	output_tree->Branch("Muon_hadEt03", &Muon_hadEt03);
	output_tree->Branch("Muon_hadVetoEt03", &Muon_hadVetoEt03);
	output_tree->Branch("Muon_nJets03", &Muon_nJets03);
	output_tree->Branch("Muon_nTracks03", &Muon_nTracks03);
	output_tree->Branch("Muon_sumPt03", &Muon_sumPt03);
	output_tree->Branch("Muon_trackerVetoPt03", &Muon_trackerVetoPt03);
	output_tree->Branch("Muon_emEt05", &Muon_emEt05);
	output_tree->Branch("Muon_emVetoEt05", &Muon_emVetoEt05);
	output_tree->Branch("Muon_hadEt05", &Muon_hadEt05);
	output_tree->Branch("Muon_hadVetoEt05", &Muon_hadVetoEt05);
	output_tree->Branch("Muon_nJets05", &Muon_nJets05);
	output_tree->Branch("Muon_nTracks05", &Muon_nTracks05);
	output_tree->Branch("Muon_sumPt05", &Muon_sumPt05);
	output_tree->Branch("Muon_trackerVetoPt05", &Muon_trackerVetoPt05);
	output_tree->Branch("Muon_sumChargedHadronPt03", &Muon_sumChargedHadronPt03);
	output_tree->Branch("Muon_sumChargedParticlePt03", &Muon_sumChargedParticlePt03);
	output_tree->Branch("Muon_sumNeutralHadronEt03", &Muon_sumNeutralHadronEt03);
	output_tree->Branch("Muon_sumNeutralHadronEtHighThreshold03", &Muon_sumNeutralHadronEtHighThreshold03);
	output_tree->Branch("Muon_sumPhotonEt03", &Muon_sumPhotonEt03);
	output_tree->Branch("Muon_sumPhotonEtHighThreshold03", &Muon_sumPhotonEtHighThreshold03);
	output_tree->Branch("Muon_sumPUPt03", &Muon_sumPUPt03);
	output_tree->Branch("Muon_sumChargedHadronPt04", &Muon_sumChargedHadronPt04);
	output_tree->Branch("Muon_sumChargedParticlePt04", &Muon_sumChargedParticlePt04);
	output_tree->Branch("Muon_sumNeutralHadronEt04", &Muon_sumNeutralHadronEt04);
	output_tree->Branch("Muon_sumNeutralHadronEtHighThreshold04", &Muon_sumNeutralHadronEtHighThreshold04);
	output_tree->Branch("Muon_sumPhotonEt04", &Muon_sumPhotonEt04);
	output_tree->Branch("Muon_sumPhotonEtHighThreshold04", &Muon_sumPhotonEtHighThreshold04);
	output_tree->Branch("Muon_sumPUPt04", &Muon_sumPUPt04);
	output_tree->Branch("Muon_Track_idx", &Muon_Track_idx);
	output_tree->Branch("Muon_hitPattern_pixelLayerwithMeas", &Muon_hitPattern_pixelLayerwithMeas);
	output_tree->Branch("Muon_numberOfMatchedStations", &Muon_numberOfMatchedStations);
	output_tree->Branch("Muon_normChi2", &Muon_normChi2);
	output_tree->Branch("Muon_hitPattern_numberOfValidMuonHits", &Muon_hitPattern_numberOfValidMuonHits);
	output_tree->Branch("Muon_innerTrack_numberofValidHits", &Muon_innerTrack_numberofValidHits);
	output_tree->Branch("Muon_numberOfMatches", &Muon_numberOfMatches);
	output_tree->Branch("Muon_numberOfChambers", &Muon_numberOfChambers);
	output_tree->Branch("Muon_isPFMuon", &Muon_isPFMuon);
	output_tree->Branch("Muon_numberofValidPixelHits", &Muon_numberofValidPixelHits);
	output_tree->Branch("Muon_trackerLayersWithMeasurement", &Muon_trackerLayersWithMeasurement);

	output_tree->Branch("Muon_charge", &Muon_charge);
	output_tree->Branch("Muon_trackCharge", &Muon_trackCharge);
	output_tree->Branch("Muon_pdgid", &Muon_pdgid);
	output_tree->Branch("Muon_B", &Muon_B);
	output_tree->Branch("Muon_M", &Muon_M);
	output_tree->Branch("Muon_par", &Muon_par);
	output_tree->Branch("Muon_cov", &Muon_cov);

	//================ Electron block ========
	output_tree->Branch("isPatElectron", &doPatElectrons_);
	output_tree->Branch("Electron_p4", &Electron_p4);
	output_tree->Branch("Electron_Poca", &Electron_Poca);
	output_tree->Branch("Electron_Gsf_deltaEtaEleClusterTrackAtCalo", &Electron_Gsf_deltaEtaEleClusterTrackAtCalo);
	output_tree->Branch("Electron_Gsf_deltaEtaSeedClusterTrackAtCalo", &Electron_Gsf_deltaEtaSeedClusterTrackAtCalo);
	output_tree->Branch("Electron_Gsf_deltaEtaSuperClusterTrackAtVtx", &Electron_Gsf_deltaEtaSuperClusterTrackAtVtx);
	output_tree->Branch("Electron_Gsf_deltaPhiEleClusterTrackAtCalo", &Electron_Gsf_deltaPhiEleClusterTrackAtCalo);
	output_tree->Branch("Electron_Gsf_deltaPhiSeedClusterTrackAtCalo", &Electron_Gsf_deltaPhiSeedClusterTrackAtCalo);
	output_tree->Branch("Electron_Gsf_deltaPhiSuperClusterTrackAtVtx", &Electron_Gsf_deltaPhiSuperClusterTrackAtVtx);
	output_tree->Branch("Electron_Gsf_dr03EcalRecHitSumE", &Electron_Gsf_dr03EcalRecHitSumE);
	output_tree->Branch("Electron_Gsf_dr03HcalDepth1TowerSumEt", &Electron_Gsf_dr03HcalDepth1TowerSumEt);
	output_tree->Branch("Electron_Gsf_dr03HcalDepth1TowerSumEtBc", &Electron_Gsf_dr03HcalDepth1TowerSumEtBc);
	output_tree->Branch("Electron_Gsf_dr03HcalDepth2TowerSumEt", &Electron_Gsf_dr03HcalDepth2TowerSumEt);
	output_tree->Branch("Electron_Gsf_dr03HcalDepth2TowerSumEtBc", &Electron_Gsf_dr03HcalDepth2TowerSumEtBc);
	output_tree->Branch("Electron_Gsf_dr03HcalTowerSumEt", &Electron_Gsf_dr03HcalTowerSumEt);
	output_tree->Branch("Electron_Gsf_dr03HcalTowerSumEtBc", &Electron_Gsf_dr03HcalTowerSumEtBc);
	output_tree->Branch("Electron_Gsf_dr03TkSumPt", &Electron_Gsf_dr03TkSumPt);
	output_tree->Branch("Electron_Gsf_passingCutBasedPreselection", &Electron_Gsf_passingCutBasedPreselection);
	output_tree->Branch("Electron_Gsf_passingMvaPreselection", &Electron_Gsf_passingMvaPreselection);
	output_tree->Branch("Electron_gsftrack_trackerExpectedHitsInner_numberOfLostHits", &Electron_gsftrack_trackerExpectedHitsInner_numberOfLostHits);
	output_tree->Branch("Electron_supercluster_e", &Electron_supercluster_e);
	output_tree->Branch("Electron_supercluster_phi", &Electron_supercluster_phi);
	output_tree->Branch("Electron_supercluster_eta", &Electron_supercluster_eta);
	output_tree->Branch("Electron_supercluster_centroid_x", &Electron_supercluster_centroid_x);
	output_tree->Branch("Electron_supercluster_centroid_y", &Electron_supercluster_centroid_y);
	output_tree->Branch("Electron_supercluster_centroid_z", &Electron_supercluster_centroid_z);
	output_tree->Branch("Electron_Track_idx", &Electron_Track_idx);

	output_tree->Branch("Electron_ecalRecHitSumEt03", &Electron_ecalRecHitSumEt03);
	output_tree->Branch("Electron_hcalDepth1TowerSumEt03", &Electron_hcalDepth1TowerSumEt03);
	output_tree->Branch("Electron_hcalDepth1TowerSumEtBc03", &Electron_hcalDepth1TowerSumEtBc03);
	output_tree->Branch("Electron_hcalDepth2TowerSumEt03", &Electron_hcalDepth2TowerSumEt03);
	output_tree->Branch("Electron_hcalDepth2TowerSumEtBc03", &Electron_hcalDepth2TowerSumEtBc03);
	output_tree->Branch("Electron_tkSumPt03", &Electron_tkSumPt03);
	output_tree->Branch("Electron_ecalRecHitSumEt04", &Electron_ecalRecHitSumEt04);
	output_tree->Branch("Electron_hcalDepth1TowerSumEt04", &Electron_hcalDepth1TowerSumEt04);
	output_tree->Branch("Electron_hcalDepth1TowerSumEtBc04", &Electron_hcalDepth1TowerSumEtBc04);
	output_tree->Branch("Electron_hcalDepth2TowerSumEt04", &Electron_hcalDepth2TowerSumEt04);
	output_tree->Branch("Electron_hcalDepth2TowerSumEtBc04", &Electron_hcalDepth2TowerSumEtBc04);
	output_tree->Branch("Electron_tkSumPt04", &Electron_tkSumPt04);
	output_tree->Branch("Electron_chargedHadronIso", &Electron_chargedHadronIso);
	output_tree->Branch("Electron_neutralHadronIso", &Electron_neutralHadronIso);
	output_tree->Branch("Electron_photonIso", &Electron_photonIso);
	output_tree->Branch("Electron_isoDeposits_chargedHadronIso04", &Electron_isoDeposits_chargedHadronIso04);
	output_tree->Branch("Electron_isoDeposits_neutralHadronIso04", &Electron_isoDeposits_neutralHadronIso04);
	output_tree->Branch("Electron_isoDeposits_photonIso04", &Electron_isoDeposits_photonIso04);
	output_tree->Branch("Electron_isoDeposits_chargedHadronIso03", &Electron_isoDeposits_chargedHadronIso03);
	output_tree->Branch("Electron_isoDeposits_neutralHadronIso03", &Electron_isoDeposits_neutralHadronIso03);
	output_tree->Branch("Electron_isoDeposits_photonIso03", &Electron_isoDeposits_photonIso03);

	output_tree->Branch("Electron_sigmaIetaIeta", &Electron_sigmaIetaIeta);
	output_tree->Branch("Electron_hadronicOverEm", &Electron_hadronicOverEm);
	output_tree->Branch("Electron_fbrem", &Electron_fbrem);
	output_tree->Branch("Electron_eSuperClusterOverP", &Electron_eSuperClusterOverP);
	output_tree->Branch("Electron_ecalEnergy", &Electron_ecalEnergy);
	output_tree->Branch("Electron_trackMomentumAtVtx", &Electron_trackMomentumAtVtx);
	output_tree->Branch("Electron_numberOfMissedHits", &Electron_numberOfMissedHits);
	output_tree->Branch("Electron_HasMatchedConversions", &Electron_HasMatchedConversions);
	output_tree->Branch("RhoIsolationAllInputTags", &RhoIsolationAllInputTags);

	output_tree->Branch("Electron_charge", &Electron_charge);
	output_tree->Branch("Electron_trackCharge", &Electron_trackCharge);
	output_tree->Branch("Electron_pdgid", &Electron_pdgid);
	output_tree->Branch("Electron_B", &Electron_B);
	output_tree->Branch("Electron_M", &Electron_M);
	output_tree->Branch("Electron_par", &Electron_par);
	output_tree->Branch("Electron_cov", &Electron_cov);

	// Electron energy calibration
	output_tree->Branch("Electron_RegEnergy", &Electron_RegEnergy);
	output_tree->Branch("Electron_RegEnergyError", &Electron_RegEnergyError);

	// Electron MVA ID
	output_tree->Branch("Electron_Rho_kt6PFJets", &Electron_Rho_kt6PFJets);
	output_tree->Branch("Electron_MVA_TrigNoIP_discriminator", &Electron_MVA_TrigNoIP_discriminator);
	output_tree->Branch("Electron_MVA_NonTrig_discriminator", &Electron_MVA_NonTrig_discriminator);
	output_tree->Branch("Electron_MVA_Trig_discriminator", &Electron_MVA_Trig_discriminator);

	//================  PFTau block ==========
	output_tree->Branch("PFTau_p4", &PFTau_p4);
	output_tree->Branch("PFTau_Poca", &PFTau_Poca);
	output_tree->Branch("PFTau_isTightIsolation", &PFTau_isTightIsolation);
	output_tree->Branch("PFTau_isMediumIsolation", &PFTau_isMediumIsolation);
	output_tree->Branch("PFTau_isLooseIsolation", &PFTau_isLooseIsolation);
	output_tree->Branch("PFTau_isTightIsolationDBSumPtCorr", &PFTau_isTightIsolationDBSumPtCorr);
	output_tree->Branch("PFTau_isMediumIsolationDBSumPtCorr", &PFTau_isMediumIsolationDBSumPtCorr);
	output_tree->Branch("PFTau_isLooseIsolationDBSumPtCorr", &PFTau_isLooseIsolationDBSumPtCorr);
	output_tree->Branch("PFTau_isVLooseIsolationDBSumPtCorr", &PFTau_isVLooseIsolationDBSumPtCorr);
	output_tree->Branch("PFTau_isHPSAgainstElectronsLoose", &PFTau_isHPSAgainstElectronsLoose);
	output_tree->Branch("PFTau_isHPSAgainstElectronsMedium", &PFTau_isHPSAgainstElectronsMedium);
	output_tree->Branch("PFTau_isHPSAgainstElectronsTight", &PFTau_isHPSAgainstElectronsTight);
	output_tree->Branch("PFTau_isHPSAgainstMuonLoose", &PFTau_isHPSAgainstMuonLoose);
	output_tree->Branch("PFTau_isHPSAgainstMuonMedium", &PFTau_isHPSAgainstMuonMedium);
	output_tree->Branch("PFTau_isHPSAgainstMuonTight", &PFTau_isHPSAgainstMuonTight);
	output_tree->Branch("PFTau_isHPSAgainstMuonLoose2", &PFTau_isHPSAgainstMuonLoose2);
	output_tree->Branch("PFTau_isHPSAgainstMuonMedium2", &PFTau_isHPSAgainstMuonMedium2);
	output_tree->Branch("PFTau_isHPSAgainstMuonTight2", &PFTau_isHPSAgainstMuonTight2);
	output_tree->Branch("PFTau_isHPSByDecayModeFinding", &PFTau_isHPSByDecayModeFinding);

	//  output_tree->Branch("PFTau_HPSPFTauDiscriminationByMVA3rawElectronRejection",&PFTau_HPSPFTauDiscriminationByMVA3rawElectronRejection);
	output_tree->Branch("PFTau_HPSPFTauDiscriminationByMVA3LooseElectronRejection", &PFTau_HPSPFTauDiscriminationByMVA3LooseElectronRejection);
	output_tree->Branch("PFTau_HPSPFTauDiscriminationByMVA3MediumElectronRejection", &PFTau_HPSPFTauDiscriminationByMVA3MediumElectronRejection);
	output_tree->Branch("PFTau_HPSPFTauDiscriminationByMVA3TightElectronRejection", &PFTau_HPSPFTauDiscriminationByMVA3TightElectronRejection);
	output_tree->Branch("PFTau_HPSPFTauDiscriminationByMVA3VTightElectronRejection", &PFTau_HPSPFTauDiscriminationByMVA3VTightElectronRejection);
	//   output_tree->Branch("PFTau_HPSPFTauDiscriminationByDeadECALElectronRejection",&PFTau_HPSPFTauDiscriminationByDeadECALElectronRejection);
	output_tree->Branch("PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits", &PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits);
	output_tree->Branch("PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits", &PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits);
	output_tree->Branch("PFTau_HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits", &PFTau_HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits);
	output_tree->Branch("PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits", &PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits);
	output_tree->Branch("PFTau_HPSPFTauDiscriminationByLooseIsolationMVA", &PFTau_HPSPFTauDiscriminationByLooseIsolationMVA);
	output_tree->Branch("PFTau_HPSPFTauDiscriminationByMediumIsolationMVA", &PFTau_HPSPFTauDiscriminationByMediumIsolationMVA);
	output_tree->Branch("PFTau_HPSPFTauDiscriminationByTightIsolationMVA", &PFTau_HPSPFTauDiscriminationByTightIsolationMVA);

	output_tree->Branch("PFTau_HPSPFTauDiscriminationByLooseIsolationMVA2", &PFTau_HPSPFTauDiscriminationByLooseIsolationMVA2);
	output_tree->Branch("PFTau_HPSPFTauDiscriminationByMediumIsolationMVA2", &PFTau_HPSPFTauDiscriminationByMediumIsolationMVA2);
	output_tree->Branch("PFTau_HPSPFTauDiscriminationByTightIsolationMVA2", &PFTau_HPSPFTauDiscriminationByTightIsolationMVA2);

	output_tree->Branch("PFTau_hpsDecayMode", &PFTau_hpsDecayMode);
	output_tree->Branch("PFTau_Charge", &PFTau_Charge);
	output_tree->Branch("PFTau_Track_idx", &PFTau_Track_idx);

	output_tree->Branch("PFTau_TIP_primaryVertex_pos", &PFTau_TIP_primaryVertex_pos);
	output_tree->Branch("PFTau_TIP_primaryVertex_cov", &PFTau_TIP_primaryVertex_cov);
	output_tree->Branch("PFTau_TIP_secondaryVertex_pos", &PFTau_TIP_secondaryVertex_pos);
	output_tree->Branch("PFTau_TIP_secondaryVertex_cov", &PFTau_TIP_secondaryVertex_cov);
	output_tree->Branch("PFTau_TIP_secondaryVertex_vtxchi2", &PFTau_TIP_secondaryVertex_vtxchi2);
	output_tree->Branch("PFTau_TIP_secondaryVertex_vtxndof", &PFTau_TIP_secondaryVertex_vtxndof);
	output_tree->Branch("PFTau_TIP_primaryVertex_vtxchi2", &PFTau_TIP_primaryVertex_vtxchi2);
	output_tree->Branch("PFTau_TIP_primaryVertex_vtxndof", &PFTau_TIP_primaryVertex_vtxndof);
	output_tree->Branch("PFTau_TIP_flightLength", &PFTau_TIP_flightLength);
	output_tree->Branch("PFTau_TIP_flightLengthSig", &PFTau_TIP_flightLengthSig);

	output_tree->Branch("PFTau_a1_lvp", &PFTau_a1_lvp);
	output_tree->Branch("PFTau_a1_cov", &PFTau_a1_cov);
	output_tree->Branch("PFTau_a1_charge", &PFTau_a1_charge);
	output_tree->Branch("PFTau_a1_pdgid", &PFTau_a1_pdgid);
	output_tree->Branch("PFTau_a1_B", &PFTau_a1_B);
	output_tree->Branch("PFTau_a1_M", &PFTau_a1_M);

	output_tree->Branch("PFTau_daughterTracks", &PFTau_daughterTracks);
	output_tree->Branch("PFTau_daughterTracks_cov", &PFTau_daughterTracks_cov);
	output_tree->Branch("PFTau_daughterTracks_charge", &PFTau_daughterTracks_charge);
	output_tree->Branch("PFTau_daughterTracks_pdgid", &PFTau_daughterTracks_pdgid);
	output_tree->Branch("PFTau_daughterTracks_B", &PFTau_daughterTracks_B);
	output_tree->Branch("PFTau_daughterTracks_M", &PFTau_daughterTracks_M);
	output_tree->Branch("PFTau_daughterTracks_poca", &PFTau_daughterTracks_poca);

	output_tree->Branch("PFTau_3PS_A1_LV", &PFTau_3PS_A1_LV);
	output_tree->Branch("PFTau_3PS_M_A1", &PFTau_3PS_M_A1);
	output_tree->Branch("PFTau_3PS_M_12", &PFTau_3PS_M_12);
	output_tree->Branch("PFTau_3PS_M_13", &PFTau_3PS_M_13);
	output_tree->Branch("PFTau_3PS_M_23", &PFTau_3PS_M_23);
	output_tree->Branch("PFTau_3PS_Tau_Charge", &PFTau_3PS_Tau_Charge);
	output_tree->Branch("PFTau_3PS_LCchi2", &PFTau_3PS_LCchi2);
	output_tree->Branch("PFTau_3PS_has3ProngSolution", &PFTau_3PS_has3ProngSolution);
	output_tree->Branch("PFTau_3PS_Tau_LV", &PFTau_3PS_Tau_LV);
	output_tree->Branch("PFTau_PionsP4", &PFTau_PionsP4);
	output_tree->Branch("PFTau_PionsCharge", &PFTau_PionsCharge);



	output_tree->Branch("PFTau_PiZeroP4", &PFTau_PiZeroP4);
	output_tree->Branch("PFTau_PiZeroNumOfPhotons", &PFTau_PiZeroNumOfPhotons);
	output_tree->Branch("PFTau_PiZeroNumOfElectrons", &PFTau_PiZeroNumOfElectrons);
	output_tree->Branch("PFTau_ChargedHadronsP4", &PFTau_ChargedHadronsP4);
	output_tree->Branch("PFTau_ChargedHadronsCharge", &PFTau_ChargedHadronsCharge);
	output_tree->Branch("PFTau_GammaP4", &PFTau_GammaP4);
	output_tree->Branch("PFTau_Photons_p4_inDR05", &PFTau_Photons_p4_inDR05);

	output_tree->Branch("PFTau_MatchedPFJetP4", &PFTau_MatchedPFJetP4);
	output_tree->Branch("PFTau_MatchedPFJetGammasP4", &PFTau_MatchedPFJetGammasP4);
	output_tree->Branch("PFTau_MatchedPFJetSCVariables", &PFTau_MatchedPFJetSCVariables);
	output_tree->Branch("PFTau_MatchedPFJetPhotonVariables", &PFTau_MatchedPFJetPhotonVariables);

	output_tree->Branch("PFTau_PhotonEnergyFraction", &PFTau_PhotonEnergyFraction);

	output_tree->Branch("PFTau_photon_hasPixelSeed", &PFTau_photon_hasPixelSeed);
	output_tree->Branch("PFTau_photon_hadronicOverEm", &PFTau_photon_hadronicOverEm);
	output_tree->Branch("PFTau_photon_sigmaIetaIeta", &PFTau_photon_sigmaIetaIeta);
	output_tree->Branch("PFTau_photon_trkSumPtHollowConeDR04", &PFTau_photon_trkSumPtHollowConeDR04);
	output_tree->Branch("PFTau_photon_ecalRecHitSumEtConeDR04", &PFTau_photon_ecalRecHitSumEtConeDR04);
	output_tree->Branch("PFTau_photon_hcalTowerSumEtConeDR04", &PFTau_photon_hcalTowerSumEtConeDR04);
	output_tree->Branch("PFTau_photon_rho", &PFTau_photon_rho);


//    output_tree->Branch("PFTau_hasSC",&PFTau_hasSC);    
//    output_tree->Branch("PFTau_hasPhoton",&PFTau_hasPhoton);    
	//=======  PFJets ===
	output_tree->Branch("isPatJet", &doPatJets_);
	output_tree->Branch("PFJet_p4", &PFJet_p4);
	output_tree->Branch("PFJet_chargedEmEnergy", &PFJet_chargedEmEnergy);
	output_tree->Branch("PFJet_chargedHadronEnergy", &PFJet_chargedHadronEnergy);
	output_tree->Branch("PFJet_chargedHadronMultiplicity", &PFJet_chargedHadronMultiplicity);
	output_tree->Branch("PFJet_chargedMuEnergy", &PFJet_chargedMuEnergy);
	output_tree->Branch("PFJet_chargedMultiplicity", &PFJet_chargedMultiplicity);
	output_tree->Branch("PFJet_electronEnergy", &PFJet_electronEnergy);
	output_tree->Branch("PFJet_electronMultiplicity", &PFJet_electronMultiplicity);
	output_tree->Branch("PFJet_HFEMEnergy", &PFJet_HFEMEnergy);
	output_tree->Branch("PFJet_HFEMMultiplicity", &PFJet_HFEMMultiplicity);
	output_tree->Branch("PFJet_HFHadronEnergy", &PFJet_HFHadronEnergy);
	output_tree->Branch("PFJet_HFHadronMultiplicity", &PFJet_HFHadronMultiplicity);
	output_tree->Branch("PFJet_muonEnergy", &PFJet_muonEnergy);
	output_tree->Branch("PFJet_muonMultiplicity", &PFJet_muonMultiplicity);
	output_tree->Branch("PFJet_neutralEmEnergy", &PFJet_neutralEmEnergy);
	output_tree->Branch("PFJet_neutralHadronEnergy", &PFJet_neutralHadronEnergy);
	output_tree->Branch("PFJet_neutralHadronMultiplicity", &PFJet_neutralHadronMultiplicity);
	output_tree->Branch("PFJet_photonEnergy", &PFJet_photonEnergy);
	output_tree->Branch("PFJet_photonMultiplicity", &PFJet_photonMultiplicity);
	output_tree->Branch("PFJet_jetArea", &PFJet_jetArea);
	output_tree->Branch("PFJet_maxDistance", &PFJet_maxDistance);
	output_tree->Branch("PFJet_nConstituents", &PFJet_nConstituents);
	output_tree->Branch("PFJet_pileup", &PFJet_pileup);
	output_tree->Branch("PFJet_etaetaMoment", &PFJet_etaetaMoment);
	output_tree->Branch("PFJet_etaphiMoment", &PFJet_etaphiMoment);
	output_tree->Branch("PFJet_Track_idx", &PFJet_Track_idx);
	output_tree->Branch("PFJet_MatchedHPS_idx", &PFJet_MatchedHPS_idx);
	output_tree->Branch("PFJet_numberOfDaughters", &PFJet_numberOfDaughters);
	output_tree->Branch("PFJet_chargedEmEnergyFraction", &PFJet_chargedEmEnergyFraction);
	output_tree->Branch("PFJet_chargedHadronEnergyFraction", &PFJet_chargedHadronEnergyFraction);
	output_tree->Branch("PFJet_neutralHadronEnergyFraction", &PFJet_neutralHadronEnergyFraction);
	output_tree->Branch("PFJet_neutralEmEnergyFraction", &PFJet_neutralEmEnergyFraction);

	output_tree->Branch("PFJet_PUJetID_discr", &PFJet_PUJetID_discr);
	output_tree->Branch("PFJet_PUJetID_looseWP", &PFJet_PUJetID_looseWP);
	output_tree->Branch("PFJet_PUJetID_mediumWP", &PFJet_PUJetID_mediumWP);
	output_tree->Branch("PFJet_PUJetID_tightWP", &PFJet_PUJetID_tightWP);

	output_tree->Branch("PFJet_partonFlavour", &PFJet_partonFlavour);
	output_tree->Branch("PFJet_bDiscriminator", &PFJet_bDiscriminator);
	//output_tree->Branch("PFJet_BTagWeight", &PFJet_BTagWeight);
	//output_tree->Branch("PFJet_bTagAlgorithmName",&PFJet_bTagAlgorithmName);
	//output_tree->Branch("PFJet_bTagAlgorithmValue",&PFJet_bTagAlgorithmValue);

	output_tree->Branch("PFJet_TracksP4", &PFJet_TracksP4);
	output_tree->Branch("PFJet_nTrk", &PFJet_nTrk);
	output_tree->Branch("PFJet_JECuncertainty", &PFJet_JECuncertainty);
	output_tree->Branch("PFJet_GenJet_p4", &PFJet_GenJet_p4);
	output_tree->Branch("PFJet_GenJet_Constituents_p4", &PFJet_GenJet_Constituents_p4);
	output_tree->Branch("PFJet_GenJetNoNu_p4", &PFJet_GenJetNoNu_p4);
	output_tree->Branch("PFJet_GenJetNoNu_Constituents_p4", &PFJet_GenJetNoNu_Constituents_p4);

	//================  MET block ==========
	output_tree->Branch("isPatMET", &doPatMET_);

	output_tree->Branch("MET_Uncorr_et", &MET_Uncorr_et);
	output_tree->Branch("MET_Uncorr_pt", &MET_Uncorr_pt);
	output_tree->Branch("MET_Uncorr_phi", &MET_Uncorr_phi);
	output_tree->Branch("MET_Uncorr_sumET", &MET_Uncorr_sumET);
	output_tree->Branch("MET_Uncorr_significance", &MET_Uncorr_significance);
	output_tree->Branch("MET_Uncorr_significance_xx", &MET_Uncorr_significance_xx);
	output_tree->Branch("MET_Uncorr_significance_xy", &MET_Uncorr_significance_xy);
	output_tree->Branch("MET_Uncorr_significance_yy", &MET_Uncorr_significance_yy);
	output_tree->Branch("MET_Uncorr_MuonEtFraction", &MET_Uncorr_MuonEtFraction);
	output_tree->Branch("MET_Uncorr_NeutralEMFraction", &MET_Uncorr_NeutralEMFraction);
	output_tree->Branch("MET_Uncorr_NeutralHadEtFraction", &MET_Uncorr_NeutralHadEtFraction);
	output_tree->Branch("MET_Uncorr_Type6EtFraction", &MET_Uncorr_Type6EtFraction);
	output_tree->Branch("MET_Uncorr_Type7EtFraction", &MET_Uncorr_Type7EtFraction);

	output_tree->Branch("MET_CorrT0rt_et", &MET_CorrT0rt_et);
	output_tree->Branch("MET_CorrT0rt_pt", &MET_CorrT0rt_pt);
	output_tree->Branch("MET_CorrT0rt_phi", &MET_CorrT0rt_phi);
	output_tree->Branch("MET_CorrT0rt_sumET", &MET_CorrT0rt_sumET);
	output_tree->Branch("MET_CorrT0rt_MuonEtFraction", &MET_CorrT0rt_MuonEtFraction);
	output_tree->Branch("MET_CorrT0rt_NeutralEMFraction", &MET_CorrT0rt_NeutralEMFraction);
	output_tree->Branch("MET_CorrT0rt_NeutralHadEtFraction", &MET_CorrT0rt_NeutralHadEtFraction);
	output_tree->Branch("MET_CorrT0rt_Type6EtFraction", &MET_CorrT0rt_Type6EtFraction);
	output_tree->Branch("MET_CorrT0rt_Type7EtFraction", &MET_CorrT0rt_Type7EtFraction);

	output_tree->Branch("MET_CorrT0rtT1_et", &MET_CorrT0rtT1_et);
	output_tree->Branch("MET_CorrT0rtT1_pt", &MET_CorrT0rtT1_pt);
	output_tree->Branch("MET_CorrT0rtT1_phi", &MET_CorrT0rtT1_phi);
	output_tree->Branch("MET_CorrT0rtT1_sumET", &MET_CorrT0rtT1_sumET);
	output_tree->Branch("MET_CorrT0rtT1_MuonEtFraction", &MET_CorrT0rtT1_MuonEtFraction);
	output_tree->Branch("MET_CorrT0rtT1_NeutralEMFraction", &MET_CorrT0rtT1_NeutralEMFraction);
	output_tree->Branch("MET_CorrT0rtT1_NeutralHadEtFraction", &MET_CorrT0rtT1_NeutralHadEtFraction);
	output_tree->Branch("MET_CorrT0rtT1_Type6EtFraction", &MET_CorrT0rtT1_Type6EtFraction);
	output_tree->Branch("MET_CorrT0rtT1_Type7EtFraction", &MET_CorrT0rtT1_Type7EtFraction);

	output_tree->Branch("MET_CorrT0pc_et", &MET_CorrT0pc_et);
	output_tree->Branch("MET_CorrT0pc_pt", &MET_CorrT0pc_pt);
	output_tree->Branch("MET_CorrT0pc_phi", &MET_CorrT0pc_phi);
	output_tree->Branch("MET_CorrT0pc_sumET", &MET_CorrT0pc_sumET);
	output_tree->Branch("MET_CorrT0pc_MuonEtFraction", &MET_CorrT0pc_MuonEtFraction);
	output_tree->Branch("MET_CorrT0pc_NeutralEMFraction", &MET_CorrT0pc_NeutralEMFraction);
	output_tree->Branch("MET_CorrT0pc_NeutralHadEtFraction", &MET_CorrT0pc_NeutralHadEtFraction);
	output_tree->Branch("MET_CorrT0pc_Type6EtFraction", &MET_CorrT0pc_Type6EtFraction);
	output_tree->Branch("MET_CorrT0pc_Type7EtFraction", &MET_CorrT0pc_Type7EtFraction);

	output_tree->Branch("MET_CorrT0pcT1_et", &MET_CorrT0pcT1_et);
	output_tree->Branch("MET_CorrT0pcT1_pt", &MET_CorrT0pcT1_pt);
	output_tree->Branch("MET_CorrT0pcT1_phi", &MET_CorrT0pcT1_phi);
	output_tree->Branch("MET_CorrT0pcT1_sumET", &MET_CorrT0pcT1_sumET);
	output_tree->Branch("MET_CorrT0pcT1_MuonEtFraction", &MET_CorrT0pcT1_MuonEtFraction);
	output_tree->Branch("MET_CorrT0pcT1_NeutralEMFraction", &MET_CorrT0pcT1_NeutralEMFraction);
	output_tree->Branch("MET_CorrT0pcT1_NeutralHadEtFraction", &MET_CorrT0pcT1_NeutralHadEtFraction);
	output_tree->Branch("MET_CorrT0pcT1_Type6EtFraction", &MET_CorrT0pcT1_Type6EtFraction);
	output_tree->Branch("MET_CorrT0pcT1_Type7EtFraction", &MET_CorrT0pcT1_Type7EtFraction);

	output_tree->Branch("MET_CorrT0rtTxy_et", &MET_CorrT0rtTxy_et);
	output_tree->Branch("MET_CorrT0rtTxy_pt", &MET_CorrT0rtTxy_pt);
	output_tree->Branch("MET_CorrT0rtTxy_phi", &MET_CorrT0rtTxy_phi);
	output_tree->Branch("MET_CorrT0rtTxy_sumET", &MET_CorrT0rtTxy_sumET);
	output_tree->Branch("MET_CorrT0rtTxy_MuonEtFraction", &MET_CorrT0rtTxy_MuonEtFraction);
	output_tree->Branch("MET_CorrT0rtTxy_NeutralEMFraction", &MET_CorrT0rtTxy_NeutralEMFraction);
	output_tree->Branch("MET_CorrT0rtTxy_NeutralHadEtFraction", &MET_CorrT0rtTxy_NeutralHadEtFraction);
	output_tree->Branch("MET_CorrT0rtTxy_Type6EtFraction", &MET_CorrT0rtTxy_Type6EtFraction);
	output_tree->Branch("MET_CorrT0rtTxy_Type7EtFraction", &MET_CorrT0rtTxy_Type7EtFraction);

	output_tree->Branch("MET_CorrT0rtT1Txy_et", &MET_CorrT0rtT1Txy_et);
	output_tree->Branch("MET_CorrT0rtT1Txy_pt", &MET_CorrT0rtT1Txy_pt);
	output_tree->Branch("MET_CorrT0rtT1Txy_phi", &MET_CorrT0rtT1Txy_phi);
	output_tree->Branch("MET_CorrT0rtT1Txy_sumET", &MET_CorrT0rtT1Txy_sumET);
	output_tree->Branch("MET_CorrT0rtT1Txy_MuonEtFraction", &MET_CorrT0rtT1Txy_MuonEtFraction);
	output_tree->Branch("MET_CorrT0rtT1Txy_NeutralEMFraction", &MET_CorrT0rtT1Txy_NeutralEMFraction);
	output_tree->Branch("MET_CorrT0rtT1Txy_NeutralHadEtFraction", &MET_CorrT0rtT1Txy_NeutralHadEtFraction);
	output_tree->Branch("MET_CorrT0rtT1Txy_Type6EtFraction", &MET_CorrT0rtT1Txy_Type6EtFraction);
	output_tree->Branch("MET_CorrT0rtT1Txy_Type7EtFraction", &MET_CorrT0rtT1Txy_Type7EtFraction);

	output_tree->Branch("MET_CorrT0pcTxy_et", &MET_CorrT0pcTxy_et);
	output_tree->Branch("MET_CorrT0pcTxy_pt", &MET_CorrT0pcTxy_pt);
	output_tree->Branch("MET_CorrT0pcTxy_phi", &MET_CorrT0pcTxy_phi);
	output_tree->Branch("MET_CorrT0pcTxy_sumET", &MET_CorrT0pcTxy_sumET);
	output_tree->Branch("MET_CorrT0pcTxy_MuonEtFraction", &MET_CorrT0pcTxy_MuonEtFraction);
	output_tree->Branch("MET_CorrT0pcTxy_NeutralEMFraction", &MET_CorrT0pcTxy_NeutralEMFraction);
	output_tree->Branch("MET_CorrT0pcTxy_NeutralHadEtFraction", &MET_CorrT0pcTxy_NeutralHadEtFraction);
	output_tree->Branch("MET_CorrT0pcTxy_Type6EtFraction", &MET_CorrT0pcTxy_Type6EtFraction);
	output_tree->Branch("MET_CorrT0pcTxy_Type7EtFraction", &MET_CorrT0pcTxy_Type7EtFraction);

	output_tree->Branch("MET_CorrT0pcT1Txy_et", &MET_CorrT0pcT1Txy_et);
	output_tree->Branch("MET_CorrT0pcT1Txy_pt", &MET_CorrT0pcT1Txy_pt);
	output_tree->Branch("MET_CorrT0pcT1Txy_phi", &MET_CorrT0pcT1Txy_phi);
	output_tree->Branch("MET_CorrT0pcT1Txy_sumET", &MET_CorrT0pcT1Txy_sumET);
	output_tree->Branch("MET_CorrT0pcT1Txy_MuonEtFraction", &MET_CorrT0pcT1Txy_MuonEtFraction);
	output_tree->Branch("MET_CorrT0pcT1Txy_NeutralEMFraction", &MET_CorrT0pcT1Txy_NeutralEMFraction);
	output_tree->Branch("MET_CorrT0pcT1Txy_NeutralHadEtFraction", &MET_CorrT0pcT1Txy_NeutralHadEtFraction);
	output_tree->Branch("MET_CorrT0pcT1Txy_Type6EtFraction", &MET_CorrT0pcT1Txy_Type6EtFraction);
	output_tree->Branch("MET_CorrT0pcT1Txy_Type7EtFraction", &MET_CorrT0pcT1Txy_Type7EtFraction);

	output_tree->Branch("MET_CorrT1_et", &MET_CorrT1_et);
	output_tree->Branch("MET_CorrT1_pt", &MET_CorrT1_pt);
	output_tree->Branch("MET_CorrT1_phi", &MET_CorrT1_phi);
	output_tree->Branch("MET_CorrT1_sumET", &MET_CorrT1_sumET);
	output_tree->Branch("MET_CorrT1_MuonEtFraction", &MET_CorrT1_MuonEtFraction);
	output_tree->Branch("MET_CorrT1_NeutralEMFraction", &MET_CorrT1_NeutralEMFraction);
	output_tree->Branch("MET_CorrT1_NeutralHadEtFraction", &MET_CorrT1_NeutralHadEtFraction);
	output_tree->Branch("MET_CorrT1_Type6EtFraction", &MET_CorrT1_Type6EtFraction);
	output_tree->Branch("MET_CorrT1_Type7EtFraction", &MET_CorrT1_Type7EtFraction);

	output_tree->Branch("MET_CorrT1Txy_et", &MET_CorrT1Txy_et);
	output_tree->Branch("MET_CorrT1Txy_pt", &MET_CorrT1Txy_pt);
	output_tree->Branch("MET_CorrT1Txy_phi", &MET_CorrT1Txy_phi);
	output_tree->Branch("MET_CorrT1Txy_sumET", &MET_CorrT1Txy_sumET);
	output_tree->Branch("MET_CorrT1Txy_MuonEtFraction", &MET_CorrT1Txy_MuonEtFraction);
	output_tree->Branch("MET_CorrT1Txy_NeutralEMFraction", &MET_CorrT1Txy_NeutralEMFraction);
	output_tree->Branch("MET_CorrT1Txy_NeutralHadEtFraction", &MET_CorrT1Txy_NeutralHadEtFraction);
	output_tree->Branch("MET_CorrT1Txy_Type6EtFraction", &MET_CorrT1Txy_Type6EtFraction);
	output_tree->Branch("MET_CorrT1Txy_Type7EtFraction", &MET_CorrT1Txy_Type7EtFraction);

	output_tree->Branch("MET_CorrCaloT1_et", &MET_CorrCaloT1_et);
	output_tree->Branch("MET_CorrCaloT1_pt", &MET_CorrCaloT1_pt);
	output_tree->Branch("MET_CorrCaloT1_phi", &MET_CorrCaloT1_phi);
	output_tree->Branch("MET_CorrCaloT1_sumET", &MET_CorrCaloT1_sumET);

	output_tree->Branch("MET_CorrCaloT1T2_et", &MET_CorrCaloT1T2_et);
	output_tree->Branch("MET_CorrCaloT1T2_pt", &MET_CorrCaloT1T2_pt);
	output_tree->Branch("MET_CorrCaloT1T2_phi", &MET_CorrCaloT1T2_phi);
	output_tree->Branch("MET_CorrCaloT1T2_sumET", &MET_CorrCaloT1T2_sumET);

	output_tree->Branch("MET_CorrMVA_et", &MET_CorrMVA_et);
	output_tree->Branch("MET_CorrMVA_pt", &MET_CorrMVA_pt);
	output_tree->Branch("MET_CorrMVA_phi", &MET_CorrMVA_phi);
	output_tree->Branch("MET_CorrMVA_sumET", &MET_CorrMVA_sumET);
	output_tree->Branch("MET_CorrMVA_significance", &MET_CorrMVA_significance);
	output_tree->Branch("MET_CorrMVA_significance_xx", &MET_CorrMVA_significance_xx);
	output_tree->Branch("MET_CorrMVA_significance_xy", &MET_CorrMVA_significance_xy);
	output_tree->Branch("MET_CorrMVA_significance_yy", &MET_CorrMVA_significance_yy);
	output_tree->Branch("MET_CorrMVA_MuonEtFraction", &MET_CorrMVA_MuonEtFraction);
	output_tree->Branch("MET_CorrMVA_NeutralEMFraction", &MET_CorrMVA_NeutralEMFraction);
	output_tree->Branch("MET_CorrMVA_NeutralHadEtFraction", &MET_CorrMVA_NeutralHadEtFraction);
	output_tree->Branch("MET_CorrMVA_Type6EtFraction", &MET_CorrMVA_Type6EtFraction);
	output_tree->Branch("MET_CorrMVA_Type7EtFraction", &MET_CorrMVA_Type7EtFraction);
	output_tree->Branch("MET_CorrMVA_srcMuon_p4", &MET_CorrMVA_srcMuon_p4);
	output_tree->Branch("MET_CorrMVA_srcElectron_p4", &MET_CorrMVA_srcElectron_p4);
	output_tree->Branch("MET_CorrMVA_srcTau_p4", &MET_CorrMVA_srcTau_p4);

	output_tree->Branch("MET_CorrMVAMuTau_et", &MET_CorrMVAMuTau_et);
	output_tree->Branch("MET_CorrMVAMuTau_pt", &MET_CorrMVAMuTau_pt);
	output_tree->Branch("MET_CorrMVAMuTau_phi", &MET_CorrMVAMuTau_phi);
	output_tree->Branch("MET_CorrMVAMuTau_sumET", &MET_CorrMVAMuTau_sumET);
	output_tree->Branch("MET_CorrMVAMuTau_significance", &MET_CorrMVAMuTau_significance);
	output_tree->Branch("MET_CorrMVAMuTau_significance_xx", &MET_CorrMVAMuTau_significance_xx);
	output_tree->Branch("MET_CorrMVAMuTau_significance_xy", &MET_CorrMVAMuTau_significance_xy);
	output_tree->Branch("MET_CorrMVAMuTau_significance_yy", &MET_CorrMVAMuTau_significance_yy);
	output_tree->Branch("MET_CorrMVAMuTau_MuonEtFraction", &MET_CorrMVAMuTau_MuonEtFraction);
	output_tree->Branch("MET_CorrMVAMuTau_NeutralEMFraction", &MET_CorrMVAMuTau_NeutralEMFraction);
	output_tree->Branch("MET_CorrMVAMuTau_NeutralHadEtFraction", &MET_CorrMVAMuTau_NeutralHadEtFraction);
	output_tree->Branch("MET_CorrMVAMuTau_Type6EtFraction", &MET_CorrMVAMuTau_Type6EtFraction);
	output_tree->Branch("MET_CorrMVAMuTau_Type7EtFraction", &MET_CorrMVAMuTau_Type7EtFraction);
	output_tree->Branch("MET_CorrMVAMuTau_srcMuon_p4", &MET_CorrMVAMuTau_srcMuon_p4);
	output_tree->Branch("MET_CorrMVAMuTau_srcTau_p4", &MET_CorrMVAMuTau_srcTau_p4);

	output_tree->Branch("MET_Type1Corr_et", &MET_Type1Corr_et);
	output_tree->Branch("MET_Type1Corr_pt", &MET_Type1Corr_pt);
	output_tree->Branch("MET_Type1Corr_phi", &MET_Type1Corr_phi);
	output_tree->Branch("MET_Type1Corr_sumET", &MET_Type1Corr_sumET);
	output_tree->Branch("MET_Type1Corr_MuonEtFraction", &MET_Type1Corr_MuonEtFraction);
	output_tree->Branch("MET_Type1Corr_NeutralEMFraction", &MET_Type1Corr_NeutralEMFraction);
	output_tree->Branch("MET_Type1Corr_NeutralHadEtFraction", &MET_Type1Corr_NeutralHadEtFraction);
	output_tree->Branch("MET_Type1Corr_Type6EtFraction", &MET_Type1Corr_Type6EtFraction);
	output_tree->Branch("MET_Type1Corr_Type7EtFraction", &MET_Type1Corr_Type7EtFraction);

	output_tree->Branch("MET_Type1p2Corr_et", &MET_Type1p2Corr_et);
	output_tree->Branch("MET_Type1p2Corr_pt", &MET_Type1p2Corr_pt);
	output_tree->Branch("MET_Type1p2Corr_phi", &MET_Type1p2Corr_phi);
	output_tree->Branch("MET_Type1p2Corr_sumET", &MET_Type1p2Corr_sumET);
	output_tree->Branch("MET_Type1p2Corr_MuonEtFraction", &MET_Type1p2Corr_MuonEtFraction);
	output_tree->Branch("MET_Type1p2Corr_NeutralEMFraction", &MET_Type1p2Corr_NeutralEMFraction);
	output_tree->Branch("MET_Type1p2Corr_NeutralHadEtFraction", &MET_Type1p2Corr_NeutralHadEtFraction);
	output_tree->Branch("MET_Type1p2Corr_Type6EtFraction", &MET_Type1p2Corr_Type6EtFraction);
	output_tree->Branch("MET_Type1p2Corr_Type7EtFraction", &MET_Type1p2Corr_Type7EtFraction);

	output_tree->Branch("MET_Type1CorrElectronUp_et", &MET_Type1CorrElectronUp_et);
	output_tree->Branch("MET_Type1CorrElectronDown_et", &MET_Type1CorrElectronDown_et);
	output_tree->Branch("MET_Type1CorrMuonUp_et", &MET_Type1CorrMuonUp_et);
	output_tree->Branch("MET_Type1CorrMuonDown_et", &MET_Type1CorrMuonDown_et);
	output_tree->Branch("MET_Type1CorrTauUp_et", &MET_Type1CorrTauUp_et);
	output_tree->Branch("MET_Type1CorrTauDown_et", &MET_Type1CorrTauDown_et);
	output_tree->Branch("MET_Type1CorrJetResUp_et", &MET_Type1CorrJetResUp_et);
	output_tree->Branch("MET_Type1CorrJetResDown_et", &MET_Type1CorrJetResDown_et);
	output_tree->Branch("MET_Type1CorrJetEnUp_et", &MET_Type1CorrJetEnUp_et);
	output_tree->Branch("MET_Type1CorrJetEnDown_et", &MET_Type1CorrJetEnDown_et);
	output_tree->Branch("MET_Type1CorrUnclusteredUp_et", &MET_Type1CorrUnclusteredUp_et);
	output_tree->Branch("MET_Type1CorrUnclusteredDown_et", &MET_Type1CorrUnclusteredDown_et);
	output_tree->Branch("MET_Type1p2CorrElectronUp_et", &MET_Type1p2CorrElectronUp_et);
	output_tree->Branch("MET_Type1p2CorrElectronDown_et", &MET_Type1p2CorrElectronDown_et);
	output_tree->Branch("MET_Type1p2CorrMuonUp_et", &MET_Type1p2CorrMuonUp_et);
	output_tree->Branch("MET_Type1p2CorrMuonDown_et", &MET_Type1p2CorrMuonDown_et);
	output_tree->Branch("MET_Type1p2CorrTauUp_et", &MET_Type1p2CorrTauUp_et);
	output_tree->Branch("MET_Type1p2CorrTauDown_et", &MET_Type1p2CorrTauDown_et);
	output_tree->Branch("MET_Type1p2CorrJetResUp_et", &MET_Type1p2CorrJetResUp_et);
	output_tree->Branch("MET_Type1p2CorrJetResDown_et", &MET_Type1p2CorrJetResDown_et);
	output_tree->Branch("MET_Type1p2CorrJetEnUp_et", &MET_Type1p2CorrJetEnUp_et);
	output_tree->Branch("MET_Type1p2CorrJetEnDown_et", &MET_Type1p2CorrJetEnDown_et);
	output_tree->Branch("MET_Type1p2CorrUnclusteredUp_et", &MET_Type1p2CorrUnclusteredUp_et);
	output_tree->Branch("MET_Type1p2CorrUnclusteredDown_et", &MET_Type1p2CorrUnclusteredDown_et);

	//=============== Event Block ==============
	output_tree->Branch("Event_EventNumber", &Event_EventNumber);
	output_tree->Branch("Event_RunNumber", &Event_RunNumber);
	output_tree->Branch("Event_bunchCrossing", &Event_bunchCrossing);
	output_tree->Branch("Event_orbitNumber", &Event_orbitNumber);
	output_tree->Branch("Event_luminosityBlock", &Event_luminosityBlock);
	output_tree->Branch("Event_isRealData", &Event_isRealData);

	output_tree->Branch("PileupInfo_TrueNumInteractions_nm1", &PileupInfo_TrueNumInteractions_nm1);
	output_tree->Branch("PileupInfo_TrueNumInteractions_n0", &PileupInfo_TrueNumInteractions_n0);
	output_tree->Branch("PileupInfo_TrueNumInteractions_np1", &PileupInfo_TrueNumInteractions_np1);
	output_tree->Branch("PUWeight",&PUWeight);
	output_tree->Branch("PUWeight_p5",&PUWeight_p5);
	output_tree->Branch("PUWeight_m5",&PUWeight_m5);
	output_tree->Branch("PUWeight3D", &PUWeight3D);
	output_tree->Branch("PUWeight3D_p5", &PUWeight3D_p5);
	output_tree->Branch("PUWeight3D_m5", &PUWeight3D_m5);
	output_tree->Branch("PUWeightFineBins",&PUWeightFineBins);

	// for embbeded samples
	output_tree->Branch("TauSpinnerWeight", &TauSpinnerWeight);
	output_tree->Branch("SelEffWeight", &SelEffWeight);
	output_tree->Branch("MinVisPtFilter", &MinVisPtFilter);
	output_tree->Branch("KinWeightPt", &KinWeightPt);
	output_tree->Branch("KinWeightEta", &KinWeightEta);
	output_tree->Branch("KinWeightMassPt", &KinWeightMassPt);
	output_tree->Branch("EmbeddedWeight", &EmbeddedWeight);

	//=============== Track Block ==============
	output_tree->Branch("Track_p4", &Track_p4);
	output_tree->Branch("Track_Poca", &Track_Poca);
	output_tree->Branch("Track_chi2", &Track_chi2);
	output_tree->Branch("Track_ndof", &Track_ndof);
	output_tree->Branch("Track_numberOfLostHits", &Track_numberOfLostHits);
	output_tree->Branch("Track_numberOfValidHits", &Track_numberOfValidHits);
	output_tree->Branch("Track_qualityMask", &Track_qualityMask);

	output_tree->Branch("Track_charge", &Track_charge);
	output_tree->Branch("Track_pdgid", &Track_pdgid);
	output_tree->Branch("Track_B", &Track_B);
	output_tree->Branch("Track_M", &Track_M);
	output_tree->Branch("Track_par", &Track_par);
	output_tree->Branch("Track_cov", &Track_cov);

	//=============== MC Block ==============

	output_tree->Branch("GenEventInfoProduct_signalProcessID", &GenEventInfoProduct_signalProcessID);
	output_tree->Branch("GenEventInfoProduct_weight", &GenEventInfoProduct_weight);
	output_tree->Branch("GenEventInfoProduct_weights", &GenEventInfoProduct_weights);
	output_tree->Branch("GenEventInfoProduct_qScale", &GenEventInfoProduct_qScale);
	output_tree->Branch("GenEventInfoProduct_alphaQED", &GenEventInfoProduct_alphaQED);
	output_tree->Branch("GenEventInfoProduct_alphaQCD", &GenEventInfoProduct_alphaQCD);

	output_tree->Branch("GenEventInfoProduct_id1", &GenEventInfoProduct_id1);
	output_tree->Branch("GenEventInfoProduct_id2", &GenEventInfoProduct_id2);
	output_tree->Branch("GenEventInfoProduct_x1", &GenEventInfoProduct_x1);
	output_tree->Branch("GenEventInfoProduct_x2", &GenEventInfoProduct_x2);
	output_tree->Branch("GenEventInfoProduct_scalePDF", &GenEventInfoProduct_scalePDF);

	if (do_MCComplete_) {
		output_tree->Branch("MC_p4", &MC_p4);
		output_tree->Branch("MC_pdgid", &MC_pdgid);
		output_tree->Branch("MC_charge", &MC_charge);
		output_tree->Branch("MC_midx", &MC_midx);
		output_tree->Branch("MC_childpdgid", &MC_childpdgid);
		output_tree->Branch("MC_status", &MC_status);
	}
	if (do_MCSummary_) {
		output_tree->Branch("MCSignalParticle_p4", &MCSignalParticle_p4);
		output_tree->Branch("MCSignalParticle_pdgid", &MCSignalParticle_pdgid);
		output_tree->Branch("MCSignalParticle_childpdgid", &MCSignalParticle_childpdgid);
		output_tree->Branch("MCSignalParticle_charge", &MCSignalParticle_charge);
		output_tree->Branch("MCSignalParticle_Poca", &MCSignalParticle_Poca);
		output_tree->Branch("MCSignalParticle_Tauidx", &MCSignalParticle_Tauidx);
		output_tree->Branch("MCTauandProd_p4", &MCTauandProd_p4);
		output_tree->Branch("MCTauandProd_Vertex", &MCTauandProd_Vertex);

		output_tree->Branch("MCTauandProd_pdgid", &MCTauandProd_pdgid);
		output_tree->Branch("MCTauandProd_midx", &MCTauandProd_midx);
		output_tree->Branch("MCTauandProd_charge", &MCTauandProd_charge);
		output_tree->Branch("MCTau_JAK", &MCTau_JAK);
		output_tree->Branch("MCTau_DecayBitMask", &MCTau_DecayBitMask);
	}

	//================= Trigger Block ===============
	output_tree->Branch("HTLTriggerName", &HTLTriggerName);
	output_tree->Branch("TriggerAccept", &TriggerAccept);
	output_tree->Branch("TriggerError", &TriggerError);
	output_tree->Branch("TriggerWasRun", &TriggerWasRun);
	output_tree->Branch("HLTPrescale", &HLTPrescale);
	output_tree->Branch("NHLTL1GTSeeds", &NHLTL1GTSeeds);
	output_tree->Branch("L1SEEDPrescale", &L1SEEDPrescale);
	output_tree->Branch("L1SEEDInvalidPrescale", &L1SEEDInvalidPrescale);
	output_tree->Branch("L1SEEDisTechBit", &L1SEEDisTechBit);
	output_tree->Branch("TriggerMatchMuon", &TriggerMatchMuon);
	//output_tree->Branch("TriggerMatchElectron", &TriggerMatchElectron);
	output_tree->Branch("TriggerMatchJet", &TriggerMatchJet);
	output_tree->Branch("TriggerMatchTau", &TriggerMatchTau);
	output_tree->Branch("HLTTrigger_objs_Pt", &HLTTrigger_objs_Pt);
	output_tree->Branch("HLTTrigger_objs_Eta", &HLTTrigger_objs_Eta);
	output_tree->Branch("HLTTrigger_objs_Phi", &HLTTrigger_objs_Phi);
	output_tree->Branch("HLTTrigger_objs_E", &HLTTrigger_objs_E);
	output_tree->Branch("HLTTrigger_objs_Id", &HLTTrigger_objs_Id);
	output_tree->Branch("HLTTrigger_objs_trigger", &HLTTrigger_objs_trigger);

	output_tree->Branch("L1TriggerName", &L1TriggerName);
	output_tree->Branch("L1TriggerDecision", &L1TriggerDecision);
	output_tree->Branch("L1ErrorCode", &L1ErrorCode);
	output_tree->Branch("L1Prescale", &L1Prescale);

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// std::vector<bool> TauNtuple::CheckTauDiscriminators(std::vector<edm::Handle<reco::PFTauDiscriminator> > tauDiscriminators, reco::PFTauRef tauRef)
//
// checks that tau candidate pass kinematic fit and KinFit quality criteria
// returns vector of bool variables, first variable true/false if tau candidate pass/fail kinematic fit
// second variable is true/false if  refitted tau candidate pass/fail quality requirements
std::vector<bool> TauNtuple::CheckTauDiscriminators(std::vector<edm::Handle<reco::PFTauDiscriminator> > tauDiscriminators, const reco::PFTauRef tauRef) {
	std::vector<bool> output_pair;
	bool discriminateByKinFit = false;
	bool discriminateByKinQC = false;
	int iDiscr = 0;
	for (std::vector<edm::Handle<reco::PFTauDiscriminator> >::const_iterator discr = tauDiscriminators.begin(); discr != tauDiscriminators.end(); ++discr) {
		iDiscr = iDiscr + (**discr)[tauRef];
	}
	if (iDiscr == 1)
		discriminateByKinFit = true;
	if (iDiscr == 2) {
		discriminateByKinFit = true;
		discriminateByKinQC = true;
	}
	output_pair.push_back(discriminateByKinFit);
	output_pair.push_back(discriminateByKinQC);

	return output_pair;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// reco::PFTauRef TauNtuple::getHPSTauMatchedToJet(edm::Handle<std::vector<reco::PFTau> > & HPStaus,   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  &Jet, unsigned int &match)
//
// finds HPS tau candidate for a given KinFit tau candidate
// the closest by deltaR HPS candidate is accepted
reco::PFTauRef TauNtuple::getMatchedHPSTau(edm::Handle<std::vector<reco::PFTau> > & HPStaus, std::vector<double> &UnmodifiedTau, int &match) {
	TLorentzVector TauVisible;
	TauVisible.SetE(UnmodifiedTau.at(0));
	TauVisible.SetPx(UnmodifiedTau.at(1));
	TauVisible.SetPy(UnmodifiedTau.at(2));
	TauVisible.SetPz(UnmodifiedTau.at(3));
	reco::PFTauRef MatchedHPSTau;
	double deltaR = 0.5; // only match if distance is less than 0.5
	match = -1;
	for (unsigned int iTau = 0; iTau < HPStaus->size(); ++iTau) {
		reco::PFTauRef HPStauCandidate(HPStaus, iTau);
		double dr = sqrt(pow(DeltaPhi(HPStauCandidate->p4().Phi(), TauVisible.Phi()), 2) + pow(HPStauCandidate->p4().Eta() - TauVisible.Eta(), 2));
		if (dr < deltaR) {
			deltaR = dr;
			match = iTau;
			MatchedHPSTau = HPStauCandidate;
		}

	}
	return MatchedHPSTau;
}

std::vector<reco::PFCandidatePtr> TauNtuple::pfCandidates(const reco::PFJet& jet, int particleId, bool sort) {
	std::vector<reco::PFCandidatePtr> pfCands = jet.getPFConstituents();
	std::vector<reco::PFCandidatePtr> selectedPFCands;

	for (std::vector<reco::PFCandidatePtr>::const_iterator iter = pfCands.begin(); iter != pfCands.end(); ++iter) {
		reco::PFCandidatePtr ptr(*iter);
		if (ptr->particleId() == particleId)
			selectedPFCands.push_back(ptr);
	}
	return selectedPFCands;
}

std::vector<reco::PFCandidatePtr> TauNtuple::pfCandidates(const reco::PFJet& jet, const std::vector<int>& particleIds, bool sort) {
	std::vector<reco::PFCandidatePtr> output;
	// Get each desired candidate type, unsorted for now
	for (std::vector<int>::const_iterator particleId = particleIds.begin(); particleId != particleIds.end(); ++particleId) {
		std::vector<reco::PFCandidatePtr> selectedPFCands = pfCandidates(jet, *particleId, false);
		output.insert(output.end(), selectedPFCands.begin(), selectedPFCands.end());
	}

	return output;
}

std::vector<reco::PFCandidatePtr> TauNtuple::pfphotons(const reco::PFJet& jet, bool sort) {
	return pfCandidates(jet, reco::PFCandidate::gamma, sort);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// reco::PFTauRef TauNtuple::getHPSTauMatchedToJet(edm::Handle<std::vector<reco::PFTau> > & HPStaus,   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  &Jet, unsigned int &match)
//
// finds HPS tau candidate for a given KinFit tau candidate
// the closest by deltaR HPS candidate is accepted
reco::PFTauRef TauNtuple::getHPSTauMatchedToJet(edm::Handle<std::vector<reco::PFTau> > & HPStaus, std::vector<double> &Jet, int &match) {
	TLorentzVector Jetp4;
	Jetp4.SetE(Jet.at(0));
	Jetp4.SetPx(Jet.at(1));
	Jetp4.SetPy(Jet.at(2));
	Jetp4.SetPz(Jet.at(3));

	reco::PFTauRef MatchedHPSTau;
	double deltaR = 0.5; // only match if distance is less than 0.5
	match = -1;
	for (unsigned int iTau = 0; iTau < HPStaus->size(); ++iTau) {
		reco::PFTauRef HPStauCandidate(HPStaus, iTau);
		double dr = sqrt(pow(DeltaPhi(HPStauCandidate->p4().Phi(), Jetp4.Phi()), 2) + pow(HPStauCandidate->p4().Eta() - Jetp4.Eta(), 2));
		if (dr < deltaR) {
			deltaR = dr;
			match = iTau;
			MatchedHPSTau = HPStauCandidate;
		}

	}
	return MatchedHPSTau;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// reco::PFTauRef TauNtuple::getHPSTauMatchedToJet(edm::Handle<std::vector<reco::PFTau> > & HPStaus,   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  &Jet, unsigned int &match)
//
// finds HPS tau candidate for a given KinFit tau candidate
// the closest by deltaR HPS candidate is accepted
reco::PFJetRef TauNtuple::getJetIndexMatchedToGivenHPSTauCandidate(edm::Handle<std::vector<reco::PFJet> > & PFJets, std::vector<double> &Tau, unsigned int &match) {
	reco::PFJetRef MatchedPFJet;
	double deltaR = 999;

	TLorentzVector TauLV;
	TauLV.SetE(Tau.at(0));
	TauLV.SetPx(Tau.at(1));
	TauLV.SetPy(Tau.at(2));
	TauLV.SetPz(Tau.at(3));

	match = 0;

	for (unsigned int iJet = 0; iJet < PFJets->size(); ++iJet) {
		reco::PFJetRef JetRefCandidate(PFJets, iJet);
		double dr = sqrt(pow(DeltaPhi(JetRefCandidate->p4().Phi(), TauLV.Phi()), 2) + pow(JetRefCandidate->p4().Eta() - TauLV.Eta(), 2));
		if (dr < deltaR) {
			deltaR = dr;
			match = iJet;
			MatchedPFJet = JetRefCandidate;
		}

	}
	return MatchedPFJet;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// bool TauNtuple::getTrackMatch(edm::Handle< std::vector<reco::Track>  > &trackCollection, reco::TrackRef &refTrack, int &match)
//
// finds track match for a given TrackRef
// returns true  if the matching track is found in the collection and sets match to the index of the found track
// retruns false if on match is found in the collection and match is set to -1.
bool TauNtuple::getTrackMatch(edm::Handle<std::vector<reco::Track> > &trackCollection, reco::TrackRef &refTrack, int &match) {
	match = -1;
	for (unsigned int iTrack = 0; iTrack < trackCollection->size(); iTrack++) {
		reco::TrackRef Track(trackCollection, iTrack);
		if (refTrack == Track) {
			match = iTrack;
			return true;
		}
	}
	return false;
}

bool TauNtuple::getTrackMatch(edm::Handle<std::vector<reco::Track> > &trackCollection, reco::GsfElectronRef &refElectron, int &match) {
	match = -1;
	for (unsigned int iTrack = 0; iTrack < trackCollection->size(); iTrack++) {
		reco::TrackRef Track(trackCollection, iTrack);
		if(Track == refElectron->closestCtfTrackRef() && refElectron->shFracInnerHits()>0.45){
			match = iTrack;
			return true;
		}
	}
	return false;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// double TauNtuple::DeltaPhi(double phi1, double phi2)
//
// Calculates the the difference between two phi angles
// 
double TauNtuple::DeltaPhi(double phi1, double phi2) {
	double dphi = fabs(phi1 - phi2);
	if (dphi > TMath::Pi())
		dphi = 2 * TMath::Pi() - dphi;
	double sign = 1;
	if (phi1 - phi2 < 0)
		sign = -1;
	return dphi * sign;
}

// ------------ method called once each job just after ending the event loop  ------------
void TauNtuple::endJob() {
	std::cout << " No Of event processed: " << cnt_ << std::endl;
	output->Write();
	output->Close();
}

// ------------ method called when starting to processes a run  ------------
void TauNtuple::beginRun(edm::Run& Run, edm::EventSetup const& Setup) {
	bool changed(true);
	TriggerOK = true;
	if (hltConfig_.init(Run, Setup, processName_, changed)) {
		// if init returns TRUE, initialisation has succeeded!
		if (changed) {
			// The HLT config has actually changed wrt the previous Run, hence rebook your
			// histograms or do anything else dependent on the revised HLT config
			//     std::cout << "Initalizing HLTConfigProvider"  << std::endl;
		}
	} else {
		// if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
		// with the file and/or code and needs to be investigated!
		//  std::cout << " HLT config extraction failure with process name " << processName_ << std::endl;
		// In this case, all access methods will return empty values!
		TriggerOK = false;
	}
	if (!Run.getByLabel(l1GtTriggerMenuLite_.label(), triggerMenuLite_)) {
		std::cout << " l1GtTrigger config extraction failure " << std::endl;
		TriggerOK = false;
	}
}

// ------------ method called when ending the processing of a run  ------------
void TauNtuple::endRun(edm::Run&, edm::EventSetup const&) {
}

// ------------ method called when starting to processes a luminosity block  ------------
void TauNtuple::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a luminosity block  ------------
void TauNtuple::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TauNtuple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

void TauNtuple::ClearEvent() {
	beamspot_par.clear();
	beamspot_cov.clear();

	Vtx_chi2.clear();
	Vtx_nTrk.clear();
	Vtx_ndof.clear();
	Vtx_y.clear();
	Vtx_x.clear();
	Vtx_z.clear();
	Vtx_Cov.clear();
	Vtx_Track_idx.clear();
	Vtx_Track_Weights.clear();
	Vtx_isFake.clear();
	Vtx_TracksP4.clear();
	//=======  Muons ===
	Muon_p4.clear();
	Muon_Poca.clear();
	Muon_isGlobalMuon.clear();
	Muon_isPFMuon.clear();
	Muon_isStandAloneMuon.clear();
	Muon_isTrackerMuon.clear();
	Muon_isCaloMuon.clear();
	Muon_isIsolationValid.clear();
	Muon_isQualityValid.clear();
	Muon_isTimeValid.clear();

	Muon_emEt03.clear();
	Muon_emVetoEt03.clear();
	Muon_hadEt03.clear();
	Muon_hadVetoEt03.clear();
	Muon_nJets03.clear();
	Muon_nTracks03.clear();
	Muon_sumPt03.clear();
	Muon_trackerVetoPt03.clear();

	Muon_emEt05.clear();
	Muon_emVetoEt05.clear();
	Muon_hadEt05.clear();
	Muon_hadVetoEt05.clear();
	Muon_nJets05.clear();
	Muon_nTracks05.clear();
	Muon_sumPt05.clear();
	Muon_trackerVetoPt05.clear();

	Muon_sumChargedHadronPt03.clear();
	Muon_sumChargedParticlePt03.clear();
	Muon_sumNeutralHadronEt03.clear();
	Muon_sumNeutralHadronEtHighThreshold03.clear();
	Muon_sumPhotonEt03.clear();
	Muon_sumPhotonEtHighThreshold03.clear();
	Muon_sumPUPt03.clear();

	Muon_sumChargedHadronPt04.clear();
	Muon_sumChargedParticlePt04.clear();
	Muon_sumNeutralHadronEt04.clear();
	Muon_sumNeutralHadronEtHighThreshold04.clear();
	Muon_sumPhotonEt04.clear();
	Muon_sumPhotonEtHighThreshold04.clear();
	Muon_sumPUPt04.clear();

	Muon_numberOfChambers.clear();
	Muon_Track_idx.clear();

	Muon_charge.clear();
	Muon_trackCharge.clear();
	Muon_pdgid.clear();
	Muon_B.clear();
	Muon_M.clear();
	Muon_par.clear();
	Muon_cov.clear();

	Muon_hitPattern_pixelLayerwithMeas.clear();
	Muon_numberOfMatchedStations.clear();
	Muon_normChi2.clear();
	Muon_hitPattern_numberOfValidMuonHits.clear();
	Muon_innerTrack_numberofValidHits.clear();
	Muon_numberOfMatches.clear();
	Muon_numberofValidPixelHits.clear();
	Muon_trackerLayersWithMeasurement.clear();

	//======= PFTaus ===
	PFTau_p4.clear();
	PFTau_Poca.clear();
	PFTau_isTightIsolation.clear();
	PFTau_isMediumIsolation.clear();
	PFTau_isLooseIsolation.clear();
	PFTau_isTightIsolationDBSumPtCorr.clear();
	PFTau_isMediumIsolationDBSumPtCorr.clear();
	PFTau_isLooseIsolationDBSumPtCorr.clear();
	PFTau_isVLooseIsolationDBSumPtCorr.clear();
	PFTau_isHPSAgainstElectronsLoose.clear();
	PFTau_isHPSAgainstElectronsMedium.clear();
	PFTau_isHPSAgainstElectronsTight.clear();
	PFTau_isHPSAgainstMuonLoose.clear();
	PFTau_isHPSAgainstMuonMedium.clear();
	PFTau_isHPSAgainstMuonTight.clear();
	PFTau_isHPSAgainstMuonLoose2.clear();
	PFTau_isHPSAgainstMuonMedium2.clear();
	PFTau_isHPSAgainstMuonTight2.clear();
	PFTau_isHPSByDecayModeFinding.clear();

	//  PFTau_HPSPFTauDiscriminationByMVA3rawElectronRejection.clear();
	PFTau_HPSPFTauDiscriminationByMVA3LooseElectronRejection.clear();
	PFTau_HPSPFTauDiscriminationByMVA3MediumElectronRejection.clear();
	PFTau_HPSPFTauDiscriminationByMVA3TightElectronRejection.clear();
	PFTau_HPSPFTauDiscriminationByMVA3VTightElectronRejection.clear();
	//  PFTau_HPSPFTauDiscriminationByDeadECALElectronRejection.clear();
	PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits.clear();
	PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits.clear();
	PFTau_HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits.clear();
	PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits.clear();
	PFTau_HPSPFTauDiscriminationByLooseIsolationMVA.clear();
	PFTau_HPSPFTauDiscriminationByMediumIsolationMVA.clear();
	PFTau_HPSPFTauDiscriminationByTightIsolationMVA.clear();

	PFTau_HPSPFTauDiscriminationByLooseIsolationMVA2.clear();
	PFTau_HPSPFTauDiscriminationByMediumIsolationMVA2.clear();
	PFTau_HPSPFTauDiscriminationByTightIsolationMVA2.clear();

	PFTau_hpsDecayMode.clear();
	PFTau_Charge.clear();
	PFTau_Track_idx.clear();

	PFTau_TIP_primaryVertex_pos.clear();
	PFTau_TIP_primaryVertex_cov.clear();
	PFTau_TIP_secondaryVertex_pos.clear();
	PFTau_TIP_secondaryVertex_cov.clear();
	PFTau_TIP_secondaryVertex_vtxchi2.clear();
	PFTau_TIP_secondaryVertex_vtxndof.clear();
	PFTau_TIP_primaryVertex_vtxchi2.clear();
	PFTau_TIP_primaryVertex_vtxndof.clear();
	PFTau_TIP_flightLength.clear();
	PFTau_TIP_flightLengthSig.clear();

	PFTau_a1_lvp.clear();
	PFTau_a1_cov.clear();
	PFTau_a1_charge.clear();
	PFTau_a1_pdgid.clear();
	PFTau_a1_B.clear();
	PFTau_a1_M.clear();

	PFTau_daughterTracks.clear();
	PFTau_daughterTracks_cov.clear();
	PFTau_daughterTracks_charge.clear();
	PFTau_daughterTracks_pdgid.clear();
	PFTau_daughterTracks_B.clear();
	PFTau_daughterTracks_M.clear();
	PFTau_daughterTracks_poca.clear();

	PFTau_3PS_A1_LV.clear();
	PFTau_3PS_M_A1.clear();
	PFTau_3PS_M_12.clear();
	PFTau_3PS_M_13.clear();
	PFTau_3PS_M_23.clear();
	PFTau_3PS_Tau_Charge.clear();
	PFTau_3PS_LCchi2.clear();
	PFTau_3PS_has3ProngSolution.clear();
	PFTau_3PS_Tau_LV.clear();

	PFTau_PionsP4.clear();
	PFTau_PionsCharge.clear();


	PFTau_PiZeroP4.clear();
	PFTau_PiZeroNumOfPhotons.clear();
	PFTau_PiZeroNumOfElectrons.clear();
	PFTau_ChargedHadronsP4.clear();
	PFTau_ChargedHadronsCharge.clear();
	PFTau_GammaP4.clear();
	PFTau_Photons_p4_inDR05.clear();
	PFTau_MatchedPFJetP4.clear();
	PFTau_MatchedPFJetGammasP4.clear();
	PFTau_MatchedPFJetSCVariables.clear();
	PFTau_MatchedPFJetPhotonVariables.clear();
	PFTau_PhotonEnergyFraction.clear();

	PFTau_photon_hasPixelSeed.clear();
	PFTau_photon_hadronicOverEm.clear();
	PFTau_photon_sigmaIetaIeta.clear();
	PFTau_photon_trkSumPtHollowConeDR04.clear();
	PFTau_photon_ecalRecHitSumEtConeDR04.clear();
	PFTau_photon_hcalTowerSumEtConeDR04.clear();
	PFTau_photon_rho.clear();


//   PFTau_hasSC.clear();
//   PFTau_hasPhoton.clear();

	//=======  Electrons ===
	Electron_p4.clear();
	Electron_Poca.clear();

	Electron_charge.clear();
	Electron_trackCharge.clear();
	Electron_pdgid.clear();
	Electron_B.clear();
	Electron_M.clear();
	Electron_par.clear();
	Electron_cov.clear();

	Electron_Gsf_deltaEtaEleClusterTrackAtCalo.clear();
	Electron_Gsf_deltaEtaSeedClusterTrackAtCalo.clear();
	Electron_Gsf_deltaEtaSuperClusterTrackAtVtx.clear();
	Electron_Gsf_deltaPhiEleClusterTrackAtCalo.clear();
	Electron_Gsf_deltaPhiSeedClusterTrackAtCalo.clear();
	Electron_Gsf_deltaPhiSuperClusterTrackAtVtx.clear();
	Electron_Gsf_dr03EcalRecHitSumE.clear();
	Electron_Gsf_dr03HcalDepth1TowerSumEt.clear();
	Electron_Gsf_dr03HcalDepth1TowerSumEtBc.clear();
	Electron_Gsf_dr03HcalDepth2TowerSumEt.clear();
	Electron_Gsf_dr03HcalDepth2TowerSumEtBc.clear();
	Electron_Gsf_dr03HcalTowerSumEt.clear();
	Electron_Gsf_dr03HcalTowerSumEtBc.clear();
	Electron_Gsf_dr03TkSumPt.clear();
	Electron_Gsf_passingCutBasedPreselection.clear();
	Electron_Gsf_passingMvaPreselection.clear();
	Electron_gsftrack_trackerExpectedHitsInner_numberOfLostHits.clear();
	Electron_supercluster_e.clear();
	Electron_supercluster_phi.clear();
	Electron_supercluster_eta.clear();
	Electron_supercluster_centroid_x.clear();
	Electron_supercluster_centroid_y.clear();
	Electron_supercluster_centroid_z.clear();
	Electron_Track_idx.clear();

	Electron_ecalRecHitSumEt03.clear();
	Electron_hcalDepth1TowerSumEt03.clear();
	Electron_hcalDepth1TowerSumEtBc03.clear();
	Electron_hcalDepth2TowerSumEt03.clear();
	Electron_hcalDepth2TowerSumEtBc03.clear();
	Electron_tkSumPt03.clear();
	Electron_ecalRecHitSumEt04.clear();
	Electron_hcalDepth1TowerSumEt04.clear();
	Electron_hcalDepth1TowerSumEtBc04.clear();
	Electron_hcalDepth2TowerSumEt04.clear();
	Electron_hcalDepth2TowerSumEtBc04.clear();
	Electron_tkSumPt04.clear();

	Electron_chargedHadronIso.clear();
	Electron_neutralHadronIso.clear();
	Electron_photonIso.clear();
	Electron_isoDeposits_chargedHadronIso04.clear();
	Electron_isoDeposits_neutralHadronIso04.clear();
	Electron_isoDeposits_photonIso04.clear();
	Electron_isoDeposits_chargedHadronIso03.clear();
	Electron_isoDeposits_neutralHadronIso03.clear();
	Electron_isoDeposits_photonIso03.clear();

	Electron_sigmaIetaIeta.clear();
	Electron_hadronicOverEm.clear();
	Electron_fbrem.clear();
	Electron_eSuperClusterOverP.clear();
	Electron_ecalEnergy.clear();
	Electron_trackMomentumAtVtx.clear();
	Electron_numberOfMissedHits.clear();
	Electron_HasMatchedConversions.clear();

	Electron_RegEnergy.clear();
	Electron_RegEnergyError.clear();

	Electron_Rho_kt6PFJets.clear();
	Electron_MVA_TrigNoIP_discriminator.clear();
	Electron_MVA_NonTrig_discriminator.clear();
	Electron_MVA_Trig_discriminator.clear();

	//=======  PFJets ===
	PFJet_p4.clear();
	PFJet_chargedEmEnergy.clear();
	PFJet_chargedHadronEnergy.clear();
	PFJet_chargedHadronMultiplicity.clear();
	PFJet_chargedMuEnergy.clear();
	PFJet_chargedMultiplicity.clear();
	PFJet_electronEnergy.clear();
	PFJet_electronMultiplicity.clear();
	PFJet_HFEMEnergy.clear();
	PFJet_HFEMMultiplicity.clear();
	PFJet_HFHadronEnergy.clear();
	PFJet_HFHadronMultiplicity.clear();
	PFJet_muonEnergy.clear();
	PFJet_muonMultiplicity.clear();
	PFJet_neutralEmEnergy.clear();
	PFJet_neutralHadronEnergy.clear();
	PFJet_neutralHadronMultiplicity.clear();
	PFJet_photonEnergy.clear();
	PFJet_photonMultiplicity.clear();
	PFJet_jetArea.clear();
	PFJet_maxDistance.clear();
	PFJet_nConstituents.clear();
	PFJet_pileup.clear();
	PFJet_etaetaMoment.clear();
	PFJet_etaphiMoment.clear();
	PFJet_Track_idx.clear();
	PFJet_MatchedHPS_idx.clear();

	PFJet_numberOfDaughters.clear();
	PFJet_chargedEmEnergyFraction.clear();
	PFJet_chargedHadronEnergyFraction.clear();
	PFJet_neutralHadronEnergyFraction.clear();
	PFJet_neutralEmEnergyFraction.clear();

	PFJet_PUJetID_discr.clear();
	PFJet_PUJetID_looseWP.clear();
	PFJet_PUJetID_mediumWP.clear();
	PFJet_PUJetID_tightWP.clear();

	PFJet_partonFlavour.clear();
	PFJet_bDiscriminator.clear();
	//PFJet_BTagWeight.clear();

	//PFJet_bTagAlgorithmName.clear();
	//PFJet_bTagAlgorithmValue.clear();

	PFJet_TracksP4.clear();
	PFJet_nTrk.clear();
	PFJet_JECuncertainty.clear();
	PFJet_GenJet_p4.clear();
	PFJet_GenJetNoNu_p4.clear();
	PFJet_GenJet_Constituents_p4.clear();
	PFJet_GenJetNoNu_Constituents_p4.clear();

	//=======  MET =======
	MET_CorrMVA_srcMuon_p4.clear();
	MET_CorrMVA_srcElectron_p4.clear();
	MET_CorrMVA_srcTau_p4.clear();
	MET_CorrMVAMuTau_srcMuon_p4.clear();
	MET_CorrMVAMuTau_srcTau_p4.clear();

	//=============== Track Block ==============
	Track_p4.clear();
	Track_Poca.clear();
	Track_chi2.clear();
	Track_ndof.clear();
	Track_numberOfLostHits.clear();
	Track_numberOfValidHits.clear();
	Track_qualityMask.clear();

	Track_charge.clear();
	Track_pdgid.clear();
	Track_B.clear();
	Track_M.clear();
	Track_par.clear();
	Track_cov.clear();

	//=============== Event Block ==============
	PUWeight = 0;
	PUWeight_p5 = 0;
	PUWeight_m5 = 0;
	PUWeight3D = 0;
	PUWeight3D_p5 = 0;
	PUWeight3D_m5 = 0;
	PUWeightFineBins = 0;

	//=============== MC Block ==============

	GenEventInfoProduct_weights.clear();
	GenEventInfoProduct_signalProcessID = 0;
	GenEventInfoProduct_weight = 0;
	GenEventInfoProduct_qScale = 0;
	GenEventInfoProduct_alphaQED = 0;
	GenEventInfoProduct_alphaQCD = 0;

	GenEventInfoProduct_id1 = 0;
	GenEventInfoProduct_id2 = 0;
	GenEventInfoProduct_x1 = 0;
	GenEventInfoProduct_x2 = 0;
	GenEventInfoProduct_scalePDF = 0;

	if (do_MCComplete_) {
		MC_p4.clear();
		MC_pdgid.clear();
		MC_charge.clear();
		MC_midx.clear();
		MC_status.clear();
		MC_childpdgid.clear();
	}
	if (do_MCSummary_) {
		MCSignalParticle_p4.clear();
		MCSignalParticle_pdgid.clear();
		MCSignalParticle_charge.clear();
		MCSignalParticle_Poca.clear();
		MCSignalParticle_Tauidx.clear();
		MCTauandProd_p4.clear();
		MCTauandProd_Vertex.clear();
		MCTauandProd_pdgid.clear();
		MCTauandProd_midx.clear();
		MCTauandProd_charge.clear();
		MCTau_JAK.clear();
		MCTau_DecayBitMask.clear();
		MCSignalParticle_childpdgid.clear();
	}

	//======================Trigger Block ======================
	HTLTriggerName.clear();
	HLTPrescale.clear();
	NHLTL1GTSeeds.clear();
	L1SEEDPrescale.clear();
	L1SEEDInvalidPrescale.clear();
	L1SEEDisTechBit.clear();
	TriggerAccept.clear();
	TriggerMatchMuon.clear();
	//TriggerMatchElectron.clear();
	TriggerMatchJet.clear();
	TriggerMatchTau.clear();
	TriggerError.clear();
	TriggerWasRun.clear();
	HLTTrigger_objs_Pt.clear();
	HLTTrigger_objs_Eta.clear();
	HLTTrigger_objs_Phi.clear();
	HLTTrigger_objs_E.clear();
	HLTTrigger_objs_Id.clear();
	HLTTrigger_objs_trigger.clear();

	L1TriggerName.clear();
	L1TriggerDecision.clear();
	L1ErrorCode.clear();
	L1Prescale.clear();

}

//define this as a plug-in
DEFINE_FWK_MODULE(TauNtuple);

