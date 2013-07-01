#include "TauDataFormat/TauNtuple/interface/TauNtuple.h"
#include "TauDataFormat/TauNtuple/interface/TauDecay_CMSSW.h"
#include "Validation/EventGenerator/interface/PdtPdgMini.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"
#include <vector>
#include <map>
#include "TMatrixT.h"

//#include "CMGTools/External/plugins/PileupJetIdProducer.cc"
// #include "CMGTools/External/interface/PileupJetIdentifier.h"
// #include "CMGTools/External/interface/PileupJetIdAlgo.h"

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
#include "RecoTauTag/ImpactParameter/interface/ParticleBuilder.h"
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

#include "DataFormats/TauReco/interface/PFTauTransverseImpactParameter.h"
#include "DataFormats/TauReco/interface/PFTauTransverseImpactParameterFwd.h"
#include "DataFormats/TauReco/interface/PFTau3ProngSummary.h"
#include "DataFormats/TauReco/interface/PFTau3ProngSummaryFwd.h"

TauNtuple::TauNtuple(const edm::ParameterSet& iConfig):
  primVtxTag_( iConfig.getParameter<edm::InputTag>( "primVtx" ) ),
  muonsTag_(iConfig.getParameter<edm::InputTag>( "muons" )),
  hpsTauProducer_( iConfig.getParameter<edm::InputTag>( "hpsTauProducer" ) ),
  PFTauTIPTag_(iConfig.getParameter<edm::InputTag>("PFTauTIPTag")),
  PFTau3PSTag_(iConfig.getParameter<edm::InputTag>("PFTau3PSTag")),
  hpsPFTauDiscriminationByTightIsolation_( iConfig.getParameter<edm::InputTag>( "hpsPFTauDiscriminationByTightIsolation" ) ),
  hpsPFTauDiscriminationByMediumIsolation_( iConfig.getParameter<edm::InputTag>( "hpsPFTauDiscriminationByMediumIsolation" ) ),
  hpsPFTauDiscriminationByLooseIsolation_( iConfig.getParameter<edm::InputTag>( "hpsPFTauDiscriminationByLooseIsolation" ) ),
  hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr_( iConfig.getParameter<edm::InputTag>( "hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr" ) ),
  hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr_( iConfig.getParameter<edm::InputTag>( "hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr" ) ),
  hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr_( iConfig.getParameter<edm::InputTag>( "hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr" ) ),
  hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr_( iConfig.getParameter<edm::InputTag>( "hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr" ) ),
  hpsPFTauDiscriminationAgainstElectronLoose_( iConfig.getParameter<edm::InputTag>( "hpsPFTauDiscriminationAgainstElectronLoose" ) ),
  hpsPFTauDiscriminationAgainstElectronMedium_( iConfig.getParameter<edm::InputTag>( "hpsPFTauDiscriminationAgainstElectronMedium" ) ),
  hpsPFTauDiscriminationAgainstElectronTight_( iConfig.getParameter<edm::InputTag>( "hpsPFTauDiscriminationAgainstElectronTight" ) ),
  hpsPFTauDiscriminationAgainstMuonLoose_( iConfig.getParameter<edm::InputTag>( "hpsPFTauDiscriminationAgainstMuonLoose" ) ),
  hpsPFTauDiscriminationAgainstMuonMedium_( iConfig.getParameter<edm::InputTag>( "hpsPFTauDiscriminationAgainstMuonMedium" ) ),
  hpsPFTauDiscriminationAgainstMuonTight_( iConfig.getParameter<edm::InputTag>( "hpsPFTauDiscriminationAgainstMuonTight" ) ),
  hpsPFTauDiscriminationByDecayModeFinding_( iConfig.getParameter<edm::InputTag>( "hpsPFTauDiscriminationByDecayModeFinding" ) ),
  pfMETTag_( iConfig.getParameter<edm::InputTag>( "pfMet" ) ),
  pfjetsTag_( iConfig.getParameter<edm::InputTag>( "pfjets" ) ),
  PFElectronTag_( iConfig.getParameter<edm::InputTag>( "pfelectrons" ) ),
  rhoIsolAllInputTag_( iConfig.getParameter<edm::InputTag>( "RhoIsolAllInputTag" ) ),
  generalTracks_(iConfig.getParameter<edm::InputTag>( "generalTracks" )),
  gensrc_(iConfig.getParameter<edm::InputTag>( "gensrc" )),
  GenEventInfo_(iConfig.getParameter<edm::InputTag>("GenEventInfo")),
  ElectronMVAWeights1_(iConfig.getUntrackedParameter<std::string>("EleMVAWeights1")), // Electron MVA ID 
  ElectronMVAWeights2_(iConfig.getUntrackedParameter<std::string>("EleMVAWeights2")), //  |              
  ElectronMVAWeights3_(iConfig.getUntrackedParameter<std::string>("EleMVAWeights3")), //  |              
  ElectronMVAWeights4_(iConfig.getUntrackedParameter<std::string>("EleMVAWeights4")), // \ /             
  ElectronMVAWeights5_(iConfig.getUntrackedParameter<std::string>("EleMVAWeights5")),
  ElectronMVAWeights6_(iConfig.getUntrackedParameter<std::string>("EleMVAWeights6")),
  ElectronMVAPtCut_(iConfig.getParameter<double>("ElectronMVAPtCut")),
  discriminators_( iConfig.getParameter< std::vector<std::string> >("discriminators") ),
  ScaleFactor_(iConfig.getUntrackedParameter<std::string>("ScaleFactor")),
  PUInputFile_(iConfig.getUntrackedParameter<std::string>("PUInputFile")),
  PUInputHistoMC_(iConfig.getUntrackedParameter<std::string>("PUInputHistoMC")),
  PUInputHistoData_(iConfig.getUntrackedParameter<std::string>("PUInputHistoData")),
  PUOutputFile_(iConfig.getUntrackedParameter("PUOutputFile",(std::string)("Weight3D.root"))),
  do_MCSummary_(iConfig.getUntrackedParameter("do_MCSummary",(bool)(true))),
  do_MCComplete_(iConfig.getUntrackedParameter("do_MCComplete",(bool)(true))),
  processName_(iConfig.getUntrackedParameter("TriggerProcessName",(std::string)"HLT")),
  TriggerInfoName_( iConfig.getParameter<edm::InputTag>("TriggerInfoName")),
  TriggerEvent_( iConfig.getParameter<edm::InputTag>("TriggerEvent")),
  TriggerResults_( iConfig.getParameter<edm::InputTag>("TriggerResults")),
  l1GtTriggerMenuLite_(iConfig.getParameter< edm::InputTag >("L1GtTriggerMenuLite")),
  doL1Triggers_(iConfig.getUntrackedParameter("doL1Triggers_",(bool)(false))),
  l1TriggerNames_(iConfig.getParameter< std::vector<std::string> >("l1TriggerNames")),
  TriggerJetMatchingdr_(iConfig.getUntrackedParameter("TriggerJetMatchingdr",(double)0.3)),
  TriggerMuonMatchingdr_(iConfig.getUntrackedParameter("TriggerMuonMatchingdr",(double)0.3)),
  TriggerElectronMatchingdr_(iConfig.getUntrackedParameter("TriggerElectronMatchingdr",(double)0.3)),
  TriggerTauMatchingdr_(iConfig.getUntrackedParameter("TriggerTauMatchingdr",(double)0.3)),
  doBJets_(iConfig.getUntrackedParameter("doBJets",(bool)(false))),
  doPFJets_(iConfig.getUntrackedParameter("doPFJets",(bool)(true))),
  doMuons_(iConfig.getUntrackedParameter("doMuons",(bool)(true))),
  doElectrons_(iConfig.getUntrackedParameter("doElectrons",(bool)(true))),
  doPFTaus_(iConfig.getUntrackedParameter("doPFTaus",(bool)(true))),
  doTracks_(iConfig.getUntrackedParameter("doTrack",(bool)(true))),
  doTrigger_(iConfig.getUntrackedParameter("doTrigger",(bool)(true))),
  doPrimeVertex_(iConfig.getUntrackedParameter("doPrimeVertex",(bool)(true))),
  doMET_(iConfig.getUntrackedParameter("doMET",(bool)(true))),
  doMC_(iConfig.getUntrackedParameter("doMC",(bool)(true))),
  doPatJets_(iConfig.getUntrackedParameter("doPatJets",(bool)(false))),
  doPatElectrons_(iConfig.getUntrackedParameter("doPatElectrons",(bool)(false))),
  doPatMuons_(iConfig.getUntrackedParameter("doPatMuons",(bool)(false))),
  doPatMET_(iConfig.getUntrackedParameter("doPatMET",(bool)(false))),
  srcPatJets_(iConfig.getUntrackedParameter("srcPatJets",(std::string)"selectedPatJets")),
  PatJetScale_(iConfig.getUntrackedParameter("PatJetScale",(std::string)"L3Absolute")),
  BTagAlgorithim_(iConfig.getUntrackedParameter("BTagAlgorithim",(std::string)"trackCountingHighEffBJetTags")),
  srcPatMET_(iConfig.getUntrackedParameter("srcPatMET",(std::string)"patMETsPF"))
{   
  system("echo 'running system to check directory structure...'");
  system("pwd");
  system("ls");
  system("ls *");
  system("ls ../*/*");

  LumiWeights_ = edm::Lumi3DReWeighting(PUInputFile_,PUInputFile_, PUInputHistoMC_, PUInputHistoData_,PUOutputFile_);
  LumiWeights_.weight3D_init(1);
  
  // Electron MVA ID
  
  myManualCatWeightsTrigNoIP2012.clear();
  myManualCatWeightsTrigNoIP2012.push_back(ElectronMVAWeights1_);
  myManualCatWeightsTrigNoIP2012.push_back(ElectronMVAWeights2_);
  myManualCatWeightsTrigNoIP2012.push_back(ElectronMVAWeights3_);
  myManualCatWeightsTrigNoIP2012.push_back(ElectronMVAWeights4_);
  myManualCatWeightsTrigNoIP2012.push_back(ElectronMVAWeights5_);
  myManualCatWeightsTrigNoIP2012.push_back(ElectronMVAWeights6_);
  
  Bool_t manualCat = true;
  myMVATrigNoIP2012 = new EGammaMvaEleEstimator();
  myMVATrigNoIP2012->initialize("BDT",
  			EGammaMvaEleEstimator::kTrigNoIP,//kNonTrig,
  			manualCat,
  			myManualCatWeightsTrigNoIP2012);
  
} 



TauNtuple::~TauNtuple()
{  
}   





// member functions
// ------------ method called to produce the data  ------------
void TauNtuple::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{ 
  cnt_++;
  // std::cout<<"Proccess event number ======>>> "<< cnt_ <<std::endl;
  if( (cnt_%500)==0) std::cout<<"Proccess event number ======>>> "<< cnt_ <<std::endl;
  ClearEvent();
  using namespace edm;
  DataMCType DMT;
  DataMC_Type_idx=DMT.GetType();
  if(iEvent.isRealData()){
    DataMC_Type_idx=DataMCType::Data;
  }
  fillEventInfo(iEvent, iSetup);  
  if(doMET_)fillMET(iEvent, iSetup);
  edm::Handle< std::vector<reco::Track>  > trackCollection;
  iEvent.getByLabel(generalTracks_,trackCollection);
  if(doPrimeVertex_)fillPrimeVertex(iEvent,iSetup,trackCollection);
  if(doMuons_)fillMuons(iEvent,iSetup,trackCollection);
  if(doElectrons_)fillElectrons(iEvent,iSetup,trackCollection);
  if(doPFTaus_)fillPFTaus(iEvent,iSetup,trackCollection);
  if(doPFJets_)fillPFJets(iEvent,iSetup,trackCollection);
  if(doTracks_)fillTracks(trackCollection);
  if(doMC_)fillMCTruth(iEvent,iSetup);
  if(doTrigger_)fillTriggerInfo(iEvent,iSetup);
  output_tree->Fill();
}      


// ------------ method called once each job just before starting event loop  ------------
                          
void TauNtuple::fillMCTruth(edm::Event& iEvent, const edm::EventSetup& iSetup){
  if(!iEvent.isRealData()){
    TauDecay_CMSSW myTauDecay;
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel(gensrc_, genParticles);
    myTauDecay.CheckForSignal(DataMC_Type_idx,genParticles);


    edm::Handle<GenEventInfoProduct> GenEventInfoProduct;   
    iEvent.getByLabel("generator",GenEventInfoProduct);            
    GenEventInfoProduct_signalProcessID=GenEventInfoProduct->signalProcessID();   
    GenEventInfoProduct_weight=GenEventInfoProduct->weight();
    GenEventInfoProduct_weights=GenEventInfoProduct->weights();        
    GenEventInfoProduct_qScale=GenEventInfoProduct->qScale();    
    GenEventInfoProduct_alphaQCD=GenEventInfoProduct->alphaQCD();      
    GenEventInfoProduct_alphaQED=GenEventInfoProduct->alphaQED();  

    if(do_MCComplete_){
      std::vector<unsigned int> index;
      for(reco::GenParticleCollection::const_iterator itr = genParticles->begin(); itr!= genParticles->end(); ++itr){
	MC_pdgid.push_back(itr->pdgId());
	MC_charge.push_back(itr->charge());
	std::vector<float> iMC_p4;
	iMC_p4.push_back(itr->p4().E());
	iMC_p4.push_back(itr->p4().Px());
	iMC_p4.push_back(itr->p4().Py());
	iMC_p4.push_back(itr->p4().Pz());

	MC_p4.push_back(iMC_p4);
	MC_midx.push_back(0);
      }   
      unsigned int i=0; 
      for(reco::GenParticleCollection::const_iterator itr = genParticles->begin(); itr!= genParticles->end(); ++itr,i++){
	for(unsigned int d = 0; d <itr->numberOfDaughters(); d++){
	  const reco::GenParticle *dau=static_cast<const reco::GenParticle*>(itr->daughter(d));
	  unsigned int j=0;
	  for(reco::GenParticleCollection::const_iterator jtr = genParticles->begin(); jtr!= genParticles->end(); ++jtr,j++){
	    if(dau->status()==jtr->status() && dau->p4()==jtr->p4() && dau->pdgId()==jtr->pdgId() && dau->numberOfMothers()==jtr->numberOfMothers() && dau->numberOfDaughters()==jtr->numberOfDaughters()){
	      MC_midx.at(j)=i;
	    }
	  }
	}
      }
    }
    if(do_MCSummary_){
      DataMCType DMT;
      for(reco::GenParticleCollection::const_iterator itr = genParticles->begin(); itr!= genParticles->end(); ++itr){
	if(DMT.isSignalParticle(itr->pdgId())){
	  MC_childpdgid.push_back(std::vector<int>());
	  MCSignalParticle_pdgid.push_back(itr->pdgId());
	  MCSignalParticle_charge.push_back(itr->charge());
	  MCSignalParticle_Tauidx.push_back(std::vector<unsigned int>());
	  std::vector<float> iSig_Poca;
	  iSig_Poca.push_back(itr->vx());
	  iSig_Poca.push_back(itr->vy());
	  iSig_Poca.push_back(itr->vz());
	  MCSignalParticle_Poca.push_back(iSig_Poca);

	  std::vector<float> iSig_p4;
	  iSig_p4.push_back(itr->p4().E());
	  iSig_p4.push_back(itr->p4().Px());
	  iSig_p4.push_back(itr->p4().Py());
	  iSig_p4.push_back(itr->p4().Pz());
	  MCSignalParticle_p4.push_back(iSig_p4);
	  // look for daughter tau
	  for(unsigned int i = 0; i <itr->numberOfDaughters(); i++){
	    const reco::Candidate *dau=itr->daughter(i);
	    MC_childpdgid.at(MC_childpdgid.size()-1).push_back(dau->pdgId());
	    if(abs(dau->pdgId())==PDGInfo::tau_minus){
	      unsigned int tauidx=MCTauandProd_p4.size();
	      MCSignalParticle_Tauidx.at(MCSignalParticle_Tauidx.size()-1).push_back(tauidx);
	      // Analysis the tau decay
	      unsigned int JAK_ID,TauBitMask;
	      myTauDecay.AnalyzeTau(static_cast<const reco::GenParticle*>(dau),JAK_ID,TauBitMask);
	      std::vector<const reco::GenParticle* > TauDecayProducts=myTauDecay.Get_TauDecayProducts();
	      MCTauandProd_midx.push_back(myTauDecay.Get_MotherIdx());
	      MCTau_JAK.push_back(JAK_ID);
	      MCTau_DecayBitMask.push_back(TauBitMask);
	      MCTauandProd_pdgid.push_back(std::vector<int>());
	      MCTauandProd_charge.push_back(std::vector<int>());
	      MCTauandProd_p4.push_back(std::vector<std::vector<float> >());
	      MCTauandProd_Vertex.push_back(std::vector<std::vector<float> >());
	     
	      for(unsigned int i=0;i<TauDecayProducts.size();i++){
		MCTauandProd_pdgid.at(tauidx).push_back(TauDecayProducts.at(i)->pdgId());
		MCTauandProd_charge.at(tauidx).push_back(TauDecayProducts.at(i)->charge());
		
		std::vector<float > iTauandProd_p4;
		std::vector<float > iTauandProd_vertex;
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

void
TauNtuple::fillPrimeVertex(edm::Event& iEvent, const edm::EventSetup& iSetup,edm::Handle< std::vector<reco::Track>  > &trackCollection){
  edm::Handle<reco::VertexCollection> primVtxs;
  iEvent.getByLabel( primVtxTag_, primVtxs);

  float nVtxs=primVtxs->size();
  int ndim=3;
  nVtxs = primVtxs->size();
  for(int i=0;i<nVtxs;i++){
    const reco::Vertex &pv = primVtxs->at(i);
    Vtx_isFake.push_back(pv.isFake());
    Vtx_chi2.push_back(pv.chi2());
    Vtx_ndof.push_back(pv.ndof());
    Vtx_x.push_back(pv.x());
    Vtx_y.push_back(pv.y());
    Vtx_z.push_back(pv.z());
    std::vector<std::vector<float> > iVtx_Cov;
    for(int j=0;j<ndim;j++){
      iVtx_Cov.push_back(std::vector<float>());
      for(int k=0;k<=j;k++){
        iVtx_Cov.at(j).push_back(pv.covariance(j,k));
      }
    }
    Vtx_Cov.push_back(iVtx_Cov);
    std::vector<int> matches;
    std::vector<float> TrackWeights;
    std::vector<std::vector<float> > iVtx_TrackP4;
    // Vtx_TracksP4.push_back(std::vector<std::vector<float>  >());
    for(reco::Vertex::trackRef_iterator iTrack=pv.tracks_begin(); iTrack<pv.tracks_end();iTrack++){
      int match(-1);
      reco::TrackRef refTrack=iTrack->castTo<reco::TrackRef>();
      if( refTrack.isNonnull() ) {
	std::vector<float> iiVtx_TrackP4;
	float trkEnergy = sqrt(refTrack->px()*refTrack->px() + refTrack->py()*refTrack->py() + refTrack->pz()*refTrack->pz() + 0.13957*0.13957);
	iiVtx_TrackP4.push_back(trkEnergy);
	iiVtx_TrackP4.push_back(refTrack->px());
	iiVtx_TrackP4.push_back(refTrack->py());
	iiVtx_TrackP4.push_back(refTrack->pz());
	iVtx_TrackP4.push_back(iiVtx_TrackP4);
	getTrackMatch(trackCollection,refTrack,match);
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

void 
TauNtuple::fillMuons(edm::Event& iEvent, const edm::EventSetup& iSetup,edm::Handle< std::vector<reco::Track>  > &trackCollection){
  edm::Handle< reco::MuonCollection > muonCollection;
  iEvent.getByLabel(muonsTag_,  muonCollection);
  int Muon_index =0;
  for(reco::MuonCollection::const_iterator iMuon = muonCollection->begin(); iMuon!= muonCollection->end(); ++iMuon, Muon_index++){
    reco::MuonRef RefMuon(muonCollection, Muon_index);
    std::vector<float> iMuon_Poca;
    iMuon_Poca.push_back(RefMuon->vx());
    iMuon_Poca.push_back(RefMuon->vy());
    iMuon_Poca.push_back(RefMuon->vz());
    Muon_Poca.push_back(iMuon_Poca);
    std::vector<float> iMuon_p4;
    iMuon_p4.push_back(RefMuon->p4().E());
    iMuon_p4.push_back(RefMuon->p4().Px());
    iMuon_p4.push_back(RefMuon->p4().Py());
    iMuon_p4.push_back(RefMuon->p4().Pz());
    Muon_p4.push_back(iMuon_p4);

    const reco::MuonIsolation  Iso03 = RefMuon->isolationR03(); 
    const reco::MuonIsolation  Iso05 = RefMuon->isolationR05(); 

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
    Muon_Charge.push_back(RefMuon->charge());
    Muon_numberOfMatchedStations.push_back(RefMuon->numberOfMatchedStations());
    Muon_numberOfMatches.push_back(RefMuon->numberOfMatches());



    if(RefMuon->isGlobalMuon()){
      Muon_normChi2.push_back(RefMuon->globalTrack()->normalizedChi2());
      Muon_hitPattern_numberOfValidMuonHits.push_back(RefMuon->globalTrack()->hitPattern().numberOfValidMuonHits());
      Muon_trackerLayersWithMeasurement.push_back(RefMuon->innerTrack()->hitPattern().trackerLayersWithMeasurement());
      Muon_numberofValidPixelHits.push_back(RefMuon->innerTrack()->hitPattern().numberOfValidPixelHits());
//       Muon_dz.push_back(RefMuon->innerTrack()->dz(vertex->position()));
//       Muon_dxy.push_back(RefMuon->innerTrack()->dxy(vertex->position()));


    }
    else{
      Muon_normChi2.push_back(0);
      Muon_hitPattern_numberOfValidMuonHits.push_back(0);
    }
    if(RefMuon->isTrackerMuon()){
      Muon_innerTrack_numberofValidHits.push_back(RefMuon->innerTrack()->numberOfValidHits());
      Muon_hitPattern_pixelLayerwithMeas.push_back(RefMuon->innerTrack()->hitPattern().pixelLayersWithMeasurement());

    }
    else{
      Muon_innerTrack_numberofValidHits.push_back(0);
      Muon_hitPattern_pixelLayerwithMeas.push_back(0);
    }

    if(RefMuon->isIsolationValid()){
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
    }
    else{// if isolation is not valid use -1 as default
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
    if(RefMuon->isPFIsolationValid() ){
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
    }
    else{// if isolation is not valid use -1 as default
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



    reco::TrackRef refTrack=RefMuon->track();
    int match;
    getTrackMatch(trackCollection,refTrack,match);
    Muon_Track_idx.push_back(match);

  }


}


void 
 TauNtuple::fillTracks(edm::Handle< std::vector<reco::Track>  > &trackCollection){
   for(unsigned int iTrack = 0; iTrack < trackCollection->size(); iTrack++) {
     reco::TrackRef Track(trackCollection, iTrack);
     std::vector<float> iTrack_p4;
    
     //assume pion mass
     float pionmass=PDGInfo::pi_mass();
     iTrack_p4.push_back(sqrt(pow(Track->p(),2.0)+pow(pionmass,2.0)));
     iTrack_p4.push_back(Track->px());
     iTrack_p4.push_back(Track->py());
     iTrack_p4.push_back(Track->pz());
     Track_p4.push_back(iTrack_p4);
     std::vector<float> iTrack_Poca;

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

     //Track par.
     std::vector<float> iTrack_par;
     std::vector<std::vector<float> > iTrack_parCov;
     for(int j=0;j<reco::Track::dimension;j++){
       iTrack_par.push_back(Track->parameters()(j));
       iTrack_parCov.push_back(std::vector<float>());
       for(int k=0;k<=j;k++){
	 iTrack_parCov.at(j).push_back(Track->covariance()(j,k));
       }
     }
     Track_par.push_back(iTrack_par);
     Track_parCov.push_back(iTrack_parCov);
   }
 }



 void  TauNtuple::fillPFTaus(edm::Event& iEvent, const edm::EventSetup& iSetup,edm::Handle< std::vector<reco::Track>  > &trackCollection){

   edm::Handle<std::vector<reco::PFTau> > HPStaus;
   iEvent.getByLabel(hpsTauProducer_, HPStaus);

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
   iEvent.getByLabel("hpsPFTauDiscriminationByLooseMuonRejection2", HPSAgainstMuonMedium2);
   edm::Handle<reco::PFTauDiscriminator> HPSAgainstMuonTight2;
   iEvent.getByLabel("hpsPFTauDiscriminationByLooseMuonRejection2", HPSAgainstMuonTight2);
   edm::Handle<reco::PFTauDiscriminator> HPSByDecayModeFinding;
   iEvent.getByLabel(hpsPFTauDiscriminationByDecayModeFinding_, HPSByDecayModeFinding);

//    edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByMVA3rawElectronRejection;
//    iEvent.getByLabel("hpsPFTauDiscriminationByMVA3rawElectronRejection", HPSPFTauDiscriminationByMVA3rawElectronRejection);
//    edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByMVA3LooseElectronRejection;
//    iEvent.getByLabel("hpsPFTauDiscriminationByMVA3LooseElectronRejection", HPSPFTauDiscriminationByMVA3LooseElectronRejection);
//    edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByMVA3MediumElectronRejection;
//    iEvent.getByLabel("hpsPFTauDiscriminationByMVA3MediumElectronRejection", HPSPFTauDiscriminationByMVA3MediumElectronRejection);
//    edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByMVA3TightElectronRejection;
//    iEvent.getByLabel("hpsPFTauDiscriminationByMVA3TightElectronRejection", HPSPFTauDiscriminationByMVA3TightElectronRejection);
//    edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByMVA3VTightElectronRejection;
//    iEvent.getByLabel("hpsPFTauDiscriminationByMVA3VTightElectronRejection", HPSPFTauDiscriminationByMVA3VTightElectronRejection);
//    edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByMVA3rawElectronRejection;
//    iEvent.getByLabel("hpsPFTauDiscriminationByMVA3rawElectronRejection", HPSPFTauDiscriminationByMVA3rawElectronRejection);
 
   edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByMVA3LooseElectronRejection;
   iEvent.getByLabel("hpsPFTauDiscriminationByMVA3LooseElectronRejection", HPSPFTauDiscriminationByMVA3LooseElectronRejection);
   edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByMVA3MediumElectronRejection;
   iEvent.getByLabel("hpsPFTauDiscriminationByMVA3MediumElectronRejection", HPSPFTauDiscriminationByMVA3MediumElectronRejection);
   edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByMVA3TightElectronRejection;
   iEvent.getByLabel("hpsPFTauDiscriminationByMVA3TightElectronRejection", HPSPFTauDiscriminationByMVA3TightElectronRejection);
   edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByMVA3VTightElectronRejection;
   iEvent.getByLabel("hpsPFTauDiscriminationByMVA3VTightElectronRejection", HPSPFTauDiscriminationByMVA3VTightElectronRejection);

//    edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByDeadECALElectronRejection;
//    iEvent.getByLabel("hpsPFTauDiscriminationByDeadECALElectronRejection", HPSPFTauDiscriminationByDeadECALElectronRejection);

   edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits;
   iEvent.getByLabel("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits", HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits);
   edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits;
   iEvent.getByLabel("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits", HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits);
   edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits;
   iEvent.getByLabel("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits", HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits);
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


   edm::Handle<edm::AssociationVector<reco::PFTauRefProd, std::vector<reco::PFTauTransverseImpactParameterRef> > > TIPAV;
   iEvent.getByLabel(PFTauTIPTag_,TIPAV);

   edm::Handle<edm::AssociationVector<reco::PFTauRefProd,  std::vector<reco::PFTau3ProngSummaryRef> > > ThreePSV;
   iEvent.getByLabel(PFTau3PSTag_,ThreePSV);

   for ( unsigned iPFTau = 0; iPFTau < HPStaus->size(); ++iPFTau ){
     reco::PFTauRef HPStauCandidate(HPStaus, iPFTau);
     if(HPStauCandidate->p4().Pt()>18 && fabs(HPStauCandidate->p4().Eta())<2.2){
     std::vector<float> iPFTau_Poca;
     iPFTau_Poca.push_back(HPStauCandidate->vx());
     iPFTau_Poca.push_back(HPStauCandidate->vy());
     iPFTau_Poca.push_back(HPStauCandidate->vz());
     PFTau_Poca.push_back(iPFTau_Poca);
     
     std::vector<float> iPFTau_p4;
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
     int Ntau=PFTau_daughterTracks.size();
     PFTau_TIP_secondaryVertex_vtxchi2.push_back(std::vector<float>());
     PFTau_TIP_secondaryVertex_vtxndof.push_back(std::vector<float>());
     PFTau_TIP_primaryVertex_vtxchi2.push_back(std::vector<float>());
     PFTau_TIP_primaryVertex_vtxndof.push_back(std::vector<float>());
     PFTau_TIP_primaryVertex_pos.push_back(std::vector<float>());
     PFTau_TIP_primaryVertex_cov.push_back(std::vector<float>());
     PFTau_TIP_secondaryVertex_pos.push_back(std::vector<float>());
     PFTau_TIP_secondaryVertex_cov.push_back(std::vector<float>());
     PFTau_a1_lvp.push_back(std::vector<float>());
     PFTau_a1_cov.push_back(std::vector<float>());

     PFTau_daughterTracks.push_back(std::vector<std::vector<float> >());
     PFTau_daughterTracks_cov.push_back(std::vector<std::vector<float> >());
     PFTau_daughterTracks_charge.push_back(std::vector<int>());
     PFTau_daughterTracks_pdgid.push_back(std::vector<int>());
     PFTau_daughterTracks_B.push_back(std::vector<float>());
     PFTau_daughterTracks_M.push_back(std::vector<float>());
     PFTau_daughterTracks_poca.push_back(std::vector<std::vector<float> >());

     PFTau_3PS_LCchi2.push_back(std::vector<float>());
     PFTau_3PS_has3ProngSolution.push_back(std::vector<int>());
     PFTau_3PS_Tau_LV.push_back(std::vector<std::vector<float> >());

     //
     PFTau_a1_charge.push_back(std::vector<int>());
     PFTau_a1_pdgid.push_back(std::vector<int>());
     PFTau_a1_B.push_back(std::vector<float>());
     PFTau_a1_M.push_back(std::vector<float>());

     PFTau_3PS_A1_LV.push_back(std::vector<float>());
     PFTau_3PS_M_A1.push_back(std::vector<float>());
     PFTau_3PS_M_12.push_back(std::vector<float>());
     PFTau_3PS_M_13.push_back(std::vector<float>());
     PFTau_3PS_M_23.push_back(std::vector<float>());
     PFTau_3PS_Tau_Charge.push_back(std::vector<int>());

     PFTau_TIP_flightLength.push_back(std::vector<float>());
     PFTau_TIP_flightLengthSig.push_back(std::vector<float>());

     const reco::PFTauTransverseImpactParameterRef theTIP=TIPAV->value(HPStauCandidate.key());
     //fill primary and Secondary vertex
     const reco::VertexRef primaryVertex=theTIP->primaryVertex();
     reco::Vertex pvtx=(*primaryVertex);
     edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
     iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder);
     //TVector3 pv(primaryVertex->position().x(),primaryVertex->position().y(),primaryVertex->position().z());
     PFTau_TIP_primaryVertex_pos.at(Ntau).push_back(primaryVertex->position().x());
     PFTau_TIP_primaryVertex_pos.at(Ntau).push_back(primaryVertex->position().y());
     PFTau_TIP_primaryVertex_pos.at(Ntau).push_back(primaryVertex->position().z());
     TMatrixTSym<double> pvcov(LorentzVectorParticle::NVertex);
     math::Error<LorentzVectorParticle::NVertex>::type pvCov;
     primaryVertex->fill(pvCov);
     for(int i = 0; i<LorentzVectorParticle::NVertex; i++)for(int j = 0; j<LorentzVectorParticle::NVertex; j++){pvcov(i,j)=pvCov(i,j);pvcov(j,i)=pvCov(i,j);}
     for(int i = 0; i<LorentzVectorParticle::NVertex; i++){
       for(int j = i; j<LorentzVectorParticle::NVertex; j++){
         PFTau_TIP_primaryVertex_pos.at(Ntau).push_back(pvcov(i,j));
       }
     }
     double vtxchi2(0), vtxndf(1);
     vtxchi2=pvtx.chi2();
     vtxndf=pvtx.ndof();
     PFTau_TIP_primaryVertex_vtxchi2.at(Ntau).push_back(vtxchi2);
     PFTau_TIP_primaryVertex_vtxndof.at(Ntau).push_back(vtxndf);

     if(theTIP->hasSecondaryVertex()){
       const reco::VertexRef secondaryVertex=theTIP->secondaryVertex();
       reco::Vertex svtx=(*secondaryVertex);
       //TVector3 sv(secondaryVertex->position().x(),secondaryVertex->position().y(),secondaryVertex->position().z());
       PFTau_TIP_secondaryVertex_pos.at(Ntau).push_back(secondaryVertex->position().x());
       PFTau_TIP_secondaryVertex_pos.at(Ntau).push_back(secondaryVertex->position().y());
       PFTau_TIP_secondaryVertex_pos.at(Ntau).push_back(secondaryVertex->position().z());
       TMatrixTSym<double> svcov(LorentzVectorParticle::NVertex);
       math::Error<LorentzVectorParticle::NVertex>::type svCov;
       secondaryVertex->fill(svCov);
       for(int i = 0; i<LorentzVectorParticle::NVertex; i++)for(int j = 0; j<LorentzVectorParticle::NVertex; j++){svcov(i,j)=svCov(i,j);svcov(j,i)=svCov(i,j);}
       for(int i = 0; i<LorentzVectorParticle::NVertex; i++){
	 for(int j = i; j<LorentzVectorParticle::NVertex; j++){
	   PFTau_TIP_secondaryVertex_cov.at(Ntau).push_back(svcov(i,j));
	 }
       }
       const reco::VertexRef secVtx=theTIP->secondaryVertex();
       GlobalPoint sv(secVtx->position().x(), secVtx->position().y(), secVtx->position().z());
       vtxchi2=0; vtxndf=1;
       vtxchi2=secVtx->chi2();
       vtxndf=secVtx->ndof();
       PFTau_TIP_secondaryVertex_vtxchi2.at(Ntau).push_back(vtxchi2);
       PFTau_TIP_secondaryVertex_vtxndof.at(Ntau).push_back(vtxndf);

       PFTau_TIP_flightLength.at(Ntau).push_back(theTIP->flightLength().X());
       PFTau_TIP_flightLength.at(Ntau).push_back(theTIP->flightLength().Y());
       PFTau_TIP_flightLength.at(Ntau).push_back(theTIP->flightLength().Z());
       PFTau_TIP_flightLengthSig.at(Ntau).push_back(theTIP->flightLengthSig());

       ////////////////////////////////////////////////////////////////////////////////
       LorentzVectorParticle a1;
       std::vector<reco::Track> selectedTracks=secVtx->refittedTracks();
       std::vector<reco::TransientTrack> transTrkVect;
       for(unsigned int i = 0; i!=selectedTracks.size();i++) transTrkVect.push_back(transTrackBuilder->build(selectedTracks.at(i)));
       KinematicParticleFactoryFromTransientTrack kinFactory;
       float piMassSigma(sqrt(pow(10.,-12.))), piChi(0.0), piNdf(0.0);
       std::vector<RefCountedKinematicParticle> pions;
       for(unsigned int i = 0; i<transTrkVect.size();i++) pions.push_back(kinFactory.particle(transTrkVect.at(i),PDGInfo::pi_mass(),piChi,piNdf,sv,piMassSigma));
       KinematicParticleVertexFitter kpvFitter;
       RefCountedKinematicTree jpTree = kpvFitter.fit(pions);
       jpTree->movePointerToTheTop();
       const KinematicParameters parameters = jpTree->currentParticle()->currentState().kinematicParameters();
       AlgebraicSymMatrix77 cov=jpTree->currentParticle()->currentState().kinematicParametersError().matrix();
       // get pions                                
       double c(0);
       std::vector<reco::Track> Tracks;
       std::vector<LorentzVectorParticle> ReFitPions;
       for(unsigned int i=0;i<transTrkVect.size();i++){
	 c+=transTrkVect.at(i).charge();
	 ReFitPions.push_back(ParticleBuilder::CreateLorentzVectorParticle(transTrkVect.at(i),transTrackBuilder,svtx,true,true));
       }
       // now covert a1 into LorentzVectorParticle
       TMatrixT<double>    a1_par(LorentzVectorParticle::NLorentzandVertexPar,1);
       TMatrixTSym<double> a1_cov(LorentzVectorParticle::NLorentzandVertexPar);
       for(int i = 0; i<LorentzVectorParticle::NLorentzandVertexPar; i++){a1_par(i,0)=parameters(i);for(int j = 0; j<LorentzVectorParticle::NLorentzandVertexPar; j++){a1_cov(i,j)=cov(i,j);} }
       a1=LorentzVectorParticle(a1_par,a1_cov,abs(PDGInfo::a_1_plus)*c,c,transTrackBuilder->field()->inInverseGeV(sv).z());
       PFTau_a1_charge.at(Ntau).push_back(a1.Charge());
       PFTau_a1_pdgid.at(Ntau).push_back(a1.PDGID());
       PFTau_a1_B.at(Ntau).push_back(a1.BField());
       PFTau_a1_M.at(Ntau).push_back(a1.Mass());
       for(int i=0;i<a1.NParameters();i++){
         PFTau_a1_lvp.at(Ntau).push_back(a1.Parameter(i));
         for(int j=i;j<a1.NParameters();j++){
           PFTau_a1_cov.at(Ntau).push_back(a1.Covariance(i,j));
         }
       }
     }
     ////////////////////////////////////////////////////////////////////////////////
     // Get unfit Tracks
     GlobalPoint pvpoint(primaryVertex->position().x(),primaryVertex->position().y(),primaryVertex->position().z());
     const std::vector<edm::Ptr<reco::PFCandidate> > cands = HPStauCandidate->signalPFChargedHadrCands();
     for (std::vector<edm::Ptr<reco::PFCandidate> >::const_iterator iter = cands.begin(); iter!=cands.end(); ++iter) {
       int Npi=PFTau_daughterTracks.at(Ntau).size();
       PFTau_daughterTracks_poca.at(Ntau).push_back(std::vector<float>());
       PFTau_daughterTracks.at(Ntau).push_back(std::vector<float>());
       PFTau_daughterTracks_cov.at(Ntau).push_back(std::vector<float>());
       //
       bool hastrack(false);
       reco::TransientTrack transTrk;
       if(iter->get()->trackRef().isNonnull()){transTrk=transTrackBuilder->build(iter->get()->trackRef());hastrack=true;}
       //else if(iter->get()->gsfTrackRef().isNonnull()){transTrk=transTrackBuilder->build(iter->get()->gsfTrackRef());hastrack=true;}
       if(hastrack){
         TrackParticle pion=ParticleBuilder::CreateTrackParticle(transTrk,transTrackBuilder,pvpoint,true,true);
         GlobalPoint pos=transTrk.trajectoryStateClosestToPoint(pvpoint).position();
         PFTau_daughterTracks_poca.at(Ntau).at(Npi).push_back(pos.x());
         PFTau_daughterTracks_poca.at(Ntau).at(Npi).push_back(pos.y());
         PFTau_daughterTracks_poca.at(Ntau).at(Npi).push_back(pos.z());
	 PFTau_daughterTracks_charge.at(Ntau).push_back(pion.Charge());
	 PFTau_daughterTracks_pdgid.at(Ntau).push_back(pion.PDGID());
	 PFTau_daughterTracks_B.at(Ntau).push_back(pion.BField());
	 PFTau_daughterTracks_M.at(Ntau).push_back(pion.Mass());
	 for(int i=0;i<pion.NParameters();i++){
	   PFTau_daughterTracks.at(Ntau).at(Npi).push_back(pion.Parameter(i));
	   for(int j=i;j<pion.NParameters();j++){
	     PFTau_daughterTracks_cov.at(Ntau).at(Npi).push_back(pion.Covariance(i,j));
	   }
	 }
       }
     }
     ////////////////////////////////////////////////////////////////////////////////
     const reco::PFTau3ProngSummaryRef the3PS=ThreePSV->value(HPStauCandidate.key());
     std::vector<float> A1_LV;
     if(theTIP->hasSecondaryVertex()){
       A1_LV.push_back(the3PS->A1_LV().E());
       A1_LV.push_back(the3PS->A1_LV().Px());
       A1_LV.push_back(the3PS->A1_LV().Py());
       A1_LV.push_back(the3PS->A1_LV().Pz());
       PFTau_3PS_A1_LV.at(Ntau)=A1_LV;
       PFTau_3PS_M_A1.at(Ntau).push_back(the3PS->M_A1());
       PFTau_3PS_M_12.at(Ntau).push_back(the3PS->M_12());
       PFTau_3PS_M_13.at(Ntau).push_back(the3PS->M_13());
       PFTau_3PS_M_23.at(Ntau).push_back(the3PS->M_23());
       PFTau_3PS_Tau_Charge.at(Ntau).push_back(the3PS->Tau_Charge());
       for(unsigned int i=0;i<reco::PFTau3ProngSummary::nsolutions;i++){
	 PFTau_3PS_LCchi2.at(Ntau).push_back(the3PS->Solution_Chi2(i));
	 PFTau_3PS_has3ProngSolution.at(Ntau).push_back((int)the3PS->has3ProngSolution(i));
	 std::vector<float> Tau_LV;
	 Tau_LV.push_back(the3PS->Tau(i).E());
	 Tau_LV.push_back(the3PS->Tau(i).Px());
	 Tau_LV.push_back(the3PS->Tau(i).Py());
	 Tau_LV.push_back(the3PS->Tau(i).Pz());
	 PFTau_3PS_Tau_LV.at(Ntau).push_back(Tau_LV);
       }
     }
     else{
       PFTau_3PS_M_A1.at(Ntau).push_back(0);
       PFTau_3PS_M_12.at(Ntau).push_back(0);
       PFTau_3PS_M_13.at(Ntau).push_back(0);
       PFTau_3PS_M_23.at(Ntau).push_back(0);
       PFTau_3PS_Tau_Charge.at(Ntau).push_back(0);
     }
     ////////////////////////////////////////////////////////////////////////////////
     const std::vector<edm::Ptr<reco::PFCandidate> >  ChargedHadrCand=HPStauCandidate->signalPFChargedHadrCands();
     std::vector<int> matches;
     for(unsigned int i=0; i<ChargedHadrCand.size();i++){
       reco::TrackRef refTrack=ChargedHadrCand.at(i).get()->trackRef();
       if( refTrack.isNonnull() ) {
	 int match(-1);
	 getTrackMatch(trackCollection,refTrack,match);
	 matches.push_back(match);
       }
     }
     PFTau_Track_idx.push_back(matches);
   }
   }
 }


void TauNtuple::fillPFJets(edm::Event& iEvent, const edm::EventSetup& iSetup,edm::Handle< std::vector<reco::Track>  > &trackCollection){
  if(!doPatJets_){
    edm::Handle<reco::PFJetCollection> JetCollection;
    iEvent.getByLabel(pfjetsTag_,  JetCollection);
    for(reco::PFJetCollection::size_type iPFJet = 0; iPFJet < JetCollection->size(); iPFJet++) {
      reco::PFJetRef PFJet(JetCollection, iPFJet);
      std::vector<float> iPFJet_Poca;
      iPFJet_Poca.push_back(PFJet->vx());
      iPFJet_Poca.push_back(PFJet->vy());
      iPFJet_Poca.push_back(PFJet->vz());
      PFJet_Poca.push_back(iPFJet_Poca);
      
      std::vector<float> iPFJet_p4;
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
      std::vector<std::vector<float> > iPFJet_TrackP4;
      for (unsigned i = 0;  i <  PFJet->numberOfDaughters (); i++) {
	const reco::PFCandidatePtr pfcand = PFJet->getPFConstituent(i);
	reco::TrackRef trackref = pfcand->trackRef();
	if( trackref.isNonnull() ) {
	  std::vector<float> iiPFJet_TrackP4;
	  float trkEnergy = sqrt(trackref->px()*trackref->px() + trackref->py()*trackref->py() + trackref->pz()*trackref->pz() + 0.13957*0.13957);
	  iiPFJet_TrackP4.push_back(trkEnergy);
	  iiPFJet_TrackP4.push_back(trackref->px());
	  iiPFJet_TrackP4.push_back(trackref->py());
	  iiPFJet_TrackP4.push_back(trackref->pz());
	  iPFJet_TrackP4.push_back(iiPFJet_TrackP4);
	  if(trackref.id() != TrID) continue;
	  int match(-1);
	  getTrackMatch(trackCollection,trackref,match);
	  if(match>=0)matches.push_back(match);
	}
      }
      PFJet_Track_idx.push_back(matches);

      PFJet_TracksP4.push_back(iPFJet_TrackP4);
      PFJet_nTrk.push_back(iPFJet_TrackP4.size());

      edm::Handle<std::vector<reco::PFTau> > HPStaus;
      iEvent.getByLabel(hpsTauProducer_, HPStaus);
      unsigned int idx =0; 
      reco::PFTauRef MatchedHPSTau = getHPSTauMatchedToJet(HPStaus,iPFJet_p4,idx);
      PFJet_MatchedHPS_idx.push_back(idx);
    }
  }
  else{
    edm::Handle<pat::JetCollection> jets;
    edm::InputTag labelJets(srcPatJets_);
    iEvent.getByLabel(labelJets, jets);   
    //edm::Handle<edm::View<pat::Jet> > jetHandle;
    //iEvent.getByLabel("PatJets", jetHandle);
    //edm::View<pat::Jet> PatJet = *jetHandle;
    
    for(pat::JetCollection::size_type iPatJet = 0; iPatJet < jets->size(); iPatJet++) {
    //for(pat::JetCollection::size_type iPatJet = 0; iPatJet < PatJet.size(); iPatJet++){
      pat::JetRef PatJet(jets, iPatJet);
      std::vector<float> iPatJet_Poca;
      iPatJet_Poca.push_back(PatJet->vx());
      iPatJet_Poca.push_back(PatJet->vy());
      iPatJet_Poca.push_back(PatJet->vz());
      PFJet_Poca.push_back(iPatJet_Poca);
 
      std::vector<float> iPatJet_p4;
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
      for (unsigned i = 0;  i <  PatJet->numberOfDaughters (); i++) {
	const reco::PFCandidatePtr pfcand = PatJet->getPFConstituent(i);
	reco::TrackRef trackref = pfcand->trackRef();
	if( trackref.isNonnull() ) {
	  if(trackref.id() != TrID) continue;
	  int match(-1);
	  getTrackMatch(trackCollection,trackref,match);
	  if(match>=0)matches.push_back(match);
	}
      }
      PFJet_Track_idx.push_back(matches);
      edm::Handle<std::vector<reco::PFTau> > HPStaus;
      iEvent.getByLabel(hpsTauProducer_, HPStaus);
      unsigned int idx =0;
      reco::PFTauRef MatchedHPSTau = getHPSTauMatchedToJet(HPStaus,iPatJet_p4,idx);
      PFJet_MatchedHPS_idx.push_back(idx);
    
      ///////////////////////////////////////////////
      //
      // B-Tagging
      //
      PFJet_partonFlavour.push_back(PatJet->partonFlavour());
      PFJet_bDiscriminator.push_back(PatJet->bDiscriminator(BTagAlgorithim_));
      std::vector<float> BTagWeights(0);
      PFJet_BTagWeight.push_back(BTagWeights);
      
      //std::cout << "!!!!!!!! b-tagging !!!!!!!!" << std::endl;
      //std::cout << "b-tag algorithm name: " << PatJet->getPairDiscri().first << ", value: " << PatJet->getPairDiscri().second << std::endl;
      //PFJet_bTagAlgorithmName.push_back(PatJet->getPairDiscri().first);
      //PFJet_bTagAlgorithmValue.push_back(PatJet->getPairDiscri().second);
      //std::cout << "name size: " << PFJet_bTagAlgorithmName.size() << ", value size: " << PFJet_bTagAlgorithmValue.size() << std::endl;

    }
  }
}



 void TauNtuple::fillElectrons(edm::Event& iEvent, const edm::EventSetup& iSetup,edm::Handle< std::vector<reco::Track>  > &trackCollection){

   edm::Handle<reco::GsfElectronCollection> ElectronCollection;
   iEvent.getByLabel(PFElectronTag_, ElectronCollection);

   edm::Handle<reco::ConversionCollection> hConversions;
   iEvent.getByLabel("allConversions", hConversions);

   edm::Handle<double> RhoIsolation;
   iEvent.getByLabel(rhoIsolAllInputTag_, RhoIsolation);
   const double *RhoIsolationRef = RhoIsolation.product();


   RhoIsolationAllInputTags  = *(RhoIsolationRef);


   
 for(reco::PFCandidateCollection::size_type iPFElectron = 0; iPFElectron < ElectronCollection->size(); iPFElectron++) {

     reco::GsfElectronRef RefElectron(ElectronCollection, iPFElectron);
     std::vector<float> iElectron_Poca;
     iElectron_Poca.push_back(RefElectron->vx());
     iElectron_Poca.push_back(RefElectron->vy());
     iElectron_Poca.push_back(RefElectron->vz());
     Electron_Poca.push_back(iElectron_Poca);
     std::vector<float> iElectron_p4;
     iElectron_p4.push_back(RefElectron->p4().E());
     iElectron_p4.push_back(RefElectron->p4().Px());
     iElectron_p4.push_back(RefElectron->p4().Py());
     iElectron_p4.push_back(RefElectron->p4().Pz());

     Electron_p4.push_back(iElectron_p4);
     Electron_Charge.push_back(RefElectron->charge());

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
     Electron_trackMomentumAtVtx.push_back(RefElectron->ecalEnergy()/RefElectron->eSuperClusterOverP());


     reco::BeamSpot beamSpot;
     edm::Handle<reco::BeamSpot> beamSpotHandle;
     iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
    
     if ( beamSpotHandle.isValid() )
       {
	 beamSpot = *beamSpotHandle;
	 
       } else
       {
	 edm::LogInfo("")
	   << "No beam spot available from EventSetup \n";
       }



     Electron_numberOfMissedHits.push_back(RefElectron->gsfTrack()->trackerExpectedHitsInner().numberOfHits());
     Electron_HasMatchedConversions.push_back(ConversionTools::hasMatchedConversion(ElectronCollection->at(iPFElectron), hConversions, beamSpot.position(),true,2.0,1e-6,0));
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


     reco::GsfTrackRef    refGsfTrack = RefElectron->gsfTrack();
     Electron_gsftrack_trackerExpectedHitsInner_numberOfLostHits.push_back(refGsfTrack->trackerExpectedHitsInner().numberOfLostHits());
     
     reco::SuperClusterRef refSuperCluster = RefElectron->superCluster();
     Electron_supercluster_e.push_back(refSuperCluster->energy());
     Electron_supercluster_phi.push_back(refSuperCluster->phi());
     Electron_supercluster_eta.push_back(refSuperCluster->eta()); 
     Electron_supercluster_centroid_x.push_back(refSuperCluster->x());
     Electron_supercluster_centroid_y.push_back(refSuperCluster->y());
     Electron_supercluster_centroid_z.push_back(refSuperCluster->z());
       
     int match;
     getTrackMatch(trackCollection,refGsfTrack,match);
     Electron_Track_idx.push_back(match);
     
     edm::Handle<double> Rhokt6PFJets;
     const edm::InputTag erho("kt6PFJets","rho");
     iEvent.getByLabel(erho, Rhokt6PFJets);
     Electron_Rho_kt6PFJets.push_back(*Rhokt6PFJets);
     
     ////////////////////////////////////////////////////////////////////////////////
     //	                                                                           //
     // MVA ID from https://twiki.cern.ch/twiki/bin/view/CMS/ElectronMVAIDForH2Tau //
     //                                                                            //
     ////////////////////////////////////////////////////////////////////////////////
     
     double myMVATrigNoIP2012Method1 = -1.;
     
     if(RefElectron->pt()>=ElectronMVAPtCut_){
    	 Bool_t validKF = true;
         reco::TrackRef myTrackRef = RefElectron->closestCtfTrackRef();
         validKF = (myTrackRef.isAvailable());
         validKF = (myTrackRef.isNonnull());
         
         EcalClusterLazyTools myEcalCluster(iEvent, iSetup, (edm::InputTag)"reducedEcalRecHitsEB", (edm::InputTag)"reducedEcalRecHitsEE");
         
         Double_t Var_spp;
         std::vector<float> vCov = myEcalCluster.localCovariances(*(RefElectron->superCluster()->seed()));
         if(!isnan(vCov[2]))Var_spp = sqrt(vCov[2]);
         else Var_spp = 0.;
         
         Double_t Var_R9 = myEcalCluster.e3x3(*(RefElectron->superCluster()->seed()))/RefElectron->superCluster()->rawEnergy();
         
         edm::Handle<double> Rho;
         const edm::InputTag eventrho("kt6PFJets","rho");
         iEvent.getByLabel(eventrho, Rho);
         Double_t Var_rho = *Rho;
         
         Bool_t printDebug = false;
         
         myMVATrigNoIP2012Method1 =  myMVATrigNoIP2012->mvaValue( (Double_t)RefElectron->fbrem(),
        		 (Double_t)((validKF) ? myTrackRef->normalizedChi2() : 0),
        		 (Int_t)((validKF) ? myTrackRef->hitPattern().trackerLayersWithMeasurement() : -1.),
        		 (Double_t)RefElectron->gsfTrack()->normalizedChi2(),
        		 (Double_t)RefElectron->deltaEtaSuperClusterTrackAtVtx(),
        		 (Double_t)RefElectron->deltaPhiSuperClusterTrackAtVtx(),
        		 (Double_t)RefElectron->deltaEtaSeedClusterTrackAtCalo(),
        		 (Double_t)RefElectron->sigmaIetaIeta(),
        		 (Double_t)Var_spp,
        		 (Double_t)RefElectron->superCluster()->etaWidth(),
        		 (Double_t)RefElectron->superCluster()->phiWidth(),
        		 (Double_t)((RefElectron->e5x5()) !=0. ? 1.-(RefElectron->e1x5()/RefElectron->e5x5()) : -1.),
        		 (Double_t)Var_R9,
        		 (Double_t)RefElectron->hadronicOverEm(),
        		 (Double_t)RefElectron->eSuperClusterOverP(),
        		 (Double_t)((1.0/RefElectron->superCluster()->energy())-(1.0/RefElectron->gsfTrack()->p())),
        		 (Double_t)RefElectron->eEleClusterOverPout(),
        		 (Double_t)Var_rho,
        		 (Double_t)(RefElectron->superCluster()->preshowerEnergy()/RefElectron->superCluster()->rawEnergy()),
        		 (Double_t)RefElectron->superCluster()->eta(),
        		 (Double_t)RefElectron->pt(),
        		 (Bool_t)printDebug);
     }
     
     Electron_MVA_discriminator.push_back(myMVATrigNoIP2012Method1);
     
   }
 }

void TauNtuple::fillMET(edm::Event& iEvent, const edm::EventSetup& iSetup){
  if(!doPatMET_){
    edm::Handle<edm::View<reco::PFMET> > pfMEThandle;
    iEvent.getByLabel(pfMETTag_, pfMEThandle);


    edm::Handle<std::vector<reco::PFMET>  > CorrectedPFMET;
    iEvent.getByLabel("pfType1CorrectedMet", CorrectedPFMET);
//     std::cout<<"-------------------------------------------------------------- >Corrected MET sie"<<CorrectedPFMET->size()<<std::endl;
//     std::cout<<"Corrected  et  pt  and phi "<< CorrectedPFMET->at(0).et() << "  " <<  CorrectedPFMET->at(0).pt()<<"  " <<  CorrectedPFMET->at(0).phi()<<std::endl;
//     std::cout<<"Un Corrected et   pt and phi "<< pfMEThandle->front().et()<< "  " <<  pfMEThandle->front().pt()<< "  " << pfMEThandle->front().phi()<<std::endl;



    MET_et=pfMEThandle->front().et();
    MET_pt=pfMEThandle->front().pt();
    MET_phi=pfMEThandle->front().phi(); 
    MET_sumET=pfMEThandle->front().sumEt();
    MET_metSignificance=-1;
    MET_MuonEtFraction=-1;
    MET_NeutralEMFraction=-1;
    MET_NeutralHadEtFraction=-1;
    MET_Type6EtFraction=-1;
    MET_Type7EtFraction=-1;

    MET_Corr_et= CorrectedPFMET->at(0).et();
    MET_Corr_pt= CorrectedPFMET->at(0).pt();
    MET_Corr_phi=CorrectedPFMET->at(0).phi();
    MET_Corr_SumET=CorrectedPFMET->at(0).sumEt();



  }
  else{
    edm::Handle<pat::MET> PatMET;
    edm::InputTag labelMET(srcPatMET_);
    iEvent.getByLabel(labelMET, PatMET);
    MET_et=PatMET->et();
    MET_phi=PatMET->phi();
    MET_sumET=PatMET->sumEt();
    MET_metSignificance=PatMET->metSignificance();
    MET_MuonEtFraction=PatMET->MuonEtFraction();
    MET_NeutralEMFraction=PatMET->NeutralEMFraction();
    MET_NeutralHadEtFraction=PatMET->NeutralHadEtFraction();
    MET_Type6EtFraction=PatMET->Type6EtFraction();
    MET_Type7EtFraction=PatMET->Type7EtFraction();
  }
}

void TauNtuple::fillTriggerInfo(edm::Event& iEvent, const edm::EventSetup& iSetup){
   if(!TriggerOK) return;
   edm::Handle<trigger::TriggerEvent> triggerEvent;
   iEvent.getByLabel(TriggerEvent_,triggerEvent);
   edm::Handle<edm::TriggerResults> triggerResults;
   iEvent.getByLabel(TriggerResults_, triggerResults);
   edm::Handle<std::vector<std::string> > MyTriggerInfoNames;
   iEvent.getByLabel(TriggerInfoName_, MyTriggerInfoNames);
   for(unsigned int i=0; i<MyTriggerInfoNames->size();i++){
     HTLTriggerName.push_back(MyTriggerInfoNames->at(i));
     unsigned int triggerIndex = hltConfig_.triggerIndex(HTLTriggerName.at(i));
     TriggerAccept.push_back(triggerResults->accept(triggerIndex));
     TriggerError.push_back(triggerResults->error(triggerIndex));
     TriggerWasRun.push_back(triggerResults->wasrun(triggerIndex));
     ////////////////////////////////////////////
     // now get level 1 & HLT prescale
     const std::vector< std::pair < bool, std::string > > level1Seeds = hltConfig_.hltL1GTSeeds(HTLTriggerName.at(i));
     int l1Prescale(-1);
     L1GtUtils l1GtUtils;
     l1GtUtils.retrieveL1EventSetup(iSetup);
     bool isTechbit=false;
     if(level1Seeds.size() == 1){
       std::vector<std::string> myl1SeedPaths;
       std::stringstream ss(level1Seeds.at(0).second);
       TString       buffer=level1Seeds.at(0).second;
       myl1SeedPaths.clear();
       if(!(buffer.Contains("(") || buffer.Contains(")") || buffer.Contains("AND") || buffer.Contains("NOT") )){
	 while(ss.good() && ! ss.eof()){
	   ss >> buffer;
	   if(!buffer.Contains("OR")) myl1SeedPaths.push_back(buffer.Data());
	 }
       }
       if(!myl1SeedPaths.empty()){
	 for(unsigned j=0; j<myl1SeedPaths.size(); j++){
	   int l1TempPrescale(-1);
	   int errorCode(0);
	   if(level1Seeds.at(0).first) { // technical triggers
	     isTechbit=true;
	     unsigned int techBit(atoi(myl1SeedPaths.at(j).c_str()));
	     const std::string techName(*(triggerMenuLite_->gtTechTrigName(techBit, errorCode)));
	     if(errorCode != 0) continue;
	     if(!l1GtUtils.decision(iEvent,techName,errorCode)) continue;
	     if(errorCode != 0) continue;
	     l1TempPrescale = l1GtUtils.prescaleFactor(iEvent,techName,errorCode);
	     if (errorCode != 0) continue; 
	   }
	   else{ // algorithmic triggers
	     if(!l1GtUtils.decision(iEvent,myl1SeedPaths.at(j),errorCode)) continue;
	     if(errorCode != 0) continue;
	     l1TempPrescale = l1GtUtils.prescaleFactor(iEvent,myl1SeedPaths.at(j),errorCode);
	     if (errorCode != 0) continue;
	   }
	   if(l1TempPrescale > 0){
	     if( l1Prescale == -1 || l1Prescale > l1TempPrescale) l1Prescale = l1TempPrescale;
	   }
	 }
       }
     }
     HLTPrescale.push_back(hltConfig_.prescaleValue(iEvent, iSetup, HTLTriggerName.at(i)));
     NHLTL1GTSeeds.push_back(level1Seeds.size());
     if(l1Prescale==-1){
       L1SEEDPrescale.push_back(1);
       L1SEEDInvalidPrescale.push_back(true);
       L1SEEDisTechBit.push_back(isTechbit);
     }
     else{
       L1SEEDPrescale.push_back((unsigned int)l1Prescale);
       L1SEEDInvalidPrescale.push_back(false);
       L1SEEDisTechBit.push_back(isTechbit);
     }
     ////////////////////////////////////
     // Now get Trigger matching
     if(triggerResults->accept(triggerIndex)){
       std::string filterName_="";
       edm::InputTag filterTag;
       std::vector<std::string> filters = hltConfig_.moduleLabels(HTLTriggerName.at(i));
       for(std::vector<std::string>::iterator filter =
	     filters.begin(); filter!= filters.end(); ++filter ) {
	 edm::InputTag testTag(*filter,"","HLT");
	 int testindex = triggerEvent->filterIndex(testTag);
	 if ( !(testindex >= triggerEvent->sizeFilters()) ) {
	   filterName_ = *filter;
	   filterTag=testTag;
	 }
       }       
       
       unsigned int index=triggerEvent->filterIndex(filterTag);
       /*std::cout << "TrgPath: " << HTLTriggerName.at(i) << " hltTag_.label(): "  
		 << filterTag.label() << "   filter name: "   
		 << filterName_ << "  sizeFilters: "   
		 << triggerEvent->sizeFilters() << std::endl;*/
       
       
       
       std::vector<float> match;
       // Muons
       edm::Handle< reco::MuonCollection > muonCollection;
       iEvent.getByLabel(muonsTag_,muonCollection);
       TriggerMatch(triggerEvent,index,muonCollection,TriggerMuonMatchingdr_,match);
       MuonTriggerMatch.push_back(match);
       match.clear();
       // Jets 
       edm::Handle<reco::PFJetCollection> JetCollection;
       iEvent.getByLabel(pfjetsTag_,  JetCollection);
       TriggerMatch(triggerEvent,index,JetCollection,TriggerJetMatchingdr_,match);
       JetTriggerMatch.push_back(match);
       match.clear();
       // Taus
       edm::Handle<std::vector<reco::PFTau> > PFTauCollection;
       iEvent.getByLabel(hpsTauProducer_, PFTauCollection);
       // 	 edm::Handle<reco::PFJetCollection> tauCollection;
       // 	 iEvent.getByLabel(pfjetsTag_, tauCollection);
       TriggerMatch(triggerEvent,index,PFTauCollection,TriggerTauMatchingdr_,match);
       TauTriggerMatch.push_back(match);
       match.clear();
       
       // Save trigger objects
       std::vector<float> TriggerObj_Pt;
       std::vector<float> TriggerObj_Eta;
       std::vector<float>  TriggerObj_Phi;
       std::vector<trigger::TriggerObject> trgobjs=triggerEvent->getObjects();
       const trigger::Keys& KEYS(triggerEvent->filterKeys(index));
       for(unsigned int ipart=0; ipart<KEYS.size();ipart++){
	 TriggerObj_Pt.push_back(trgobjs.at(KEYS.at(ipart)).pt());
	 TriggerObj_Eta.push_back(trgobjs.at(KEYS.at(ipart)).eta());
	 TriggerObj_Phi.push_back(trgobjs.at(KEYS.at(ipart)).phi());
       }       
       HLTTrigger_objs_Pt.push_back(TriggerObj_Pt);
       HLTTrigger_objs_Eta.push_back(TriggerObj_Eta);
       HLTTrigger_objs_Phi.push_back(TriggerObj_Phi);
     }
     else{
       MuonTriggerMatch.push_back(std::vector<float>());
       JetTriggerMatch.push_back(std::vector<float>());
       TauTriggerMatch.push_back(std::vector<float>());
       HLTTrigger_objs_Pt.push_back(std::vector<float>());
       HLTTrigger_objs_Eta.push_back(std::vector<float>());
       HLTTrigger_objs_Phi.push_back(std::vector<float>());
     }
     //////////////////////////////////// 
     // Now do L1 TriggerSeeds if requested
     if(doL1Triggers_){
       for(unsigned j=0; j<l1TriggerNames_.size(); j++){
	 int errorCode(0);
         L1TriggerName.push_back(l1TriggerNames_.at(j));
         L1TriggerDecision.push_back(l1GtUtils.decision(iEvent,l1TriggerNames_.at(j),errorCode));
         L1ErrorCode.push_back(errorCode);
         L1Prescale.push_back(l1GtUtils.prescaleFactor(iEvent,l1TriggerNames_.at(j),errorCode));
       }
     }
   }
}


 template <class T>
 void TauNtuple::TriggerMatch(edm::Handle<trigger::TriggerEvent> &triggerEvent,unsigned int triggerIndex,T obj,
			      double drmax,std::vector<float> &match){
   match=std::vector<float>(obj->size(),999);
   std::vector<trigger::TriggerObject> trgobjs=triggerEvent->getObjects();
   const trigger::Keys& KEYS(triggerEvent->filterKeys(triggerIndex));
   for(unsigned int ipart=0; ipart<KEYS.size();ipart++){
     for(unsigned int i=0; i< obj->size(); ++i ){
       double dr = reco::deltaR(trgobjs.at(KEYS.at(ipart)).eta(),trgobjs.at(KEYS.at(ipart)).phi(),obj->at(i).eta(),obj->at(i).phi());
       if(dr<drmax){
	 match.at(i)=dr;
	 /*	 std::cout << "Found Trigger Match " << i << " " << dr << " " << drmax << " Trigger Obj: " << trgobjs.at(KEYS.at(ipart)).eta()
		   << " " << trgobjs.at(KEYS.at(ipart)).phi() << " " <<  trgobjs.at(KEYS.at(ipart)).energy() 
		   << " obj: " <<  obj->at(i).eta() << " " << obj->at(i).phi() << " " << obj->at(i).energy() << std::endl;*/
       }
     }
   }
 }



 void 
 TauNtuple::fillEventInfo(edm::Event& iEvent, const edm::EventSetup& iSetup){

   Event_EventNumber=iEvent.id().event();
   Event_RunNumber=iEvent.id().run();
   Event_bunchCrossing=iEvent.bunchCrossing(); 
   Event_orbitNumber=iEvent.orbitNumber();
   Event_luminosityBlock=iEvent.luminosityBlock(); 
   Event_isRealData=iEvent.isRealData();
   if(!Event_isRealData){
     edm::Handle<std::vector< PileupSummaryInfo > > PupInfo; 
     iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo); 
     std::vector<PileupSummaryInfo>::const_iterator PVI; 
     PileupInfo_NumInteractions_nm1 = -1;
     PileupInfo_NumInteractions_n0  = -1; 
     PileupInfo_NumInteractions_np1 = -1; 
     for(PVI = PupInfo->begin(); PVI != PupInfo->end();++PVI) { 
       int BX = PVI->getBunchCrossing(); 
       if(BX == -1) PileupInfo_NumInteractions_nm1 =  PVI->getPU_NumInteractions(); 
       if(BX == 0)  PileupInfo_NumInteractions_n0  =  PVI->getPU_NumInteractions();  
       if(BX == 1)  PileupInfo_NumInteractions_np1 =  PVI->getPU_NumInteractions(); 
     } 
     EvtWeight3D = LumiWeights_.weight3D( PileupInfo_NumInteractions_nm1,PileupInfo_NumInteractions_n0,PileupInfo_NumInteractions_np1);}
 }





 void 
 TauNtuple::beginJob()
 {

   std::cout<<"----------------------------------- >>>>>>>>>>>>>> TauNtuple begin Job"<<std::endl;
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
   cnt_=0;
   output = new TFile("TauNtuple.root","RECREATE");
   output_tree = new TTree("t","t");

   output_tree->Branch("DataMC_Type",&DataMC_Type_idx);

   //=============  Vertex Block ====
   output_tree->Branch("Vtx_chi2",&Vtx_chi2);
   output_tree->Branch("Vtx_nTrk",&Vtx_nTrk);
   output_tree->Branch("Vtx_ndof",&Vtx_ndof);
   output_tree->Branch("Vtx_x",&Vtx_x);
   output_tree->Branch("Vtx_y",&Vtx_y);
   output_tree->Branch("Vtx_z",&Vtx_z);
   output_tree->Branch("Vtx_Cov",&Vtx_Cov);
   output_tree->Branch("Vtx_Track_idx",&Vtx_Track_idx);
   output_tree->Branch("Vtx_Track_Weights",&Vtx_Track_Weights);
   output_tree->Branch("Vtx_isFake",&Vtx_isFake);
   output_tree->Branch("Vtx_TracksP4",&Vtx_TracksP4);


   //=============  Muon Block ====
   output_tree->Branch("isPatMuon",&doPatMuons_);
   output_tree->Branch("Muon_p4",&Muon_p4);		       
   output_tree->Branch("Muon_Poca",&Muon_Poca);	       
   output_tree->Branch("Muon_isGlobalMuon",&Muon_isGlobalMuon);      
   output_tree->Branch("Muon_isStandAloneMuon",&Muon_isStandAloneMuon);  
   output_tree->Branch("Muon_isTrackerMuon",&Muon_isTrackerMuon);     
   output_tree->Branch("Muon_isCaloMuon",&Muon_isCaloMuon);	       
   output_tree->Branch("Muon_isIsolationValid",&Muon_isIsolationValid);  
   output_tree->Branch("Muon_isQualityValid",&Muon_isQualityValid);    
   output_tree->Branch("Muon_isTimeValid",&Muon_isTimeValid);       
   output_tree->Branch("Muon_emEt03",&Muon_emEt03);          
   output_tree->Branch("Muon_emVetoEt03",&Muon_emVetoEt03);      
   output_tree->Branch("Muon_hadEt03",&Muon_hadEt03);         
   output_tree->Branch("Muon_hadVetoEt03",&Muon_hadVetoEt03);     
   output_tree->Branch("Muon_nJets03",&Muon_nJets03);         
   output_tree->Branch("Muon_nTracks03",&Muon_nTracks03);       
   output_tree->Branch("Muon_sumPt03",&Muon_sumPt03);         
   output_tree->Branch("Muon_trackerVetoPt03",&Muon_trackerVetoPt03); 
   output_tree->Branch("Muon_emEt05",&Muon_emEt05);          
   output_tree->Branch("Muon_emVetoEt05",&Muon_emVetoEt05);      
   output_tree->Branch("Muon_hadEt05",&Muon_hadEt05);         
   output_tree->Branch("Muon_hadVetoEt05",&Muon_hadVetoEt05);     
   output_tree->Branch("Muon_nJets05",&Muon_nJets05);         
   output_tree->Branch("Muon_nTracks05",&Muon_nTracks05);       
   output_tree->Branch("Muon_sumPt05",&Muon_sumPt05);         
   output_tree->Branch("Muon_trackerVetoPt05",&Muon_trackerVetoPt05); 
   output_tree->Branch("Muon_sumChargedHadronPt03",&Muon_sumChargedHadronPt03);           
   output_tree->Branch("Muon_sumChargedParticlePt03",&Muon_sumChargedParticlePt03);        
   output_tree->Branch("Muon_sumNeutralHadronEt03",&Muon_sumNeutralHadronEt03);      
   output_tree->Branch("Muon_sumNeutralHadronEtHighThreshold03",&Muon_sumNeutralHadronEtHighThreshold03);
   output_tree->Branch("Muon_sumPhotonEt03",&Muon_sumPhotonEt03);                
   output_tree->Branch("Muon_sumPhotonEtHighThreshold03",&Muon_sumPhotonEtHighThreshold03);   
   output_tree->Branch("Muon_sumPUPt03",&Muon_sumPUPt03);                        
   output_tree->Branch("Muon_sumChargedHadronPt04",&Muon_sumChargedHadronPt04);           
   output_tree->Branch("Muon_sumChargedParticlePt04",&Muon_sumChargedParticlePt04);      
   output_tree->Branch("Muon_sumNeutralHadronEt04",&Muon_sumNeutralHadronEt04);      
   output_tree->Branch("Muon_sumNeutralHadronEtHighThreshold04",&Muon_sumNeutralHadronEtHighThreshold04);
   output_tree->Branch("Muon_sumPhotonEt04",&Muon_sumPhotonEt04);                    
   output_tree->Branch("Muon_sumPhotonEtHighThreshold04",&Muon_sumPhotonEtHighThreshold04);      
   output_tree->Branch("Muon_sumPUPt04",&Muon_sumPUPt04);                      
   output_tree->Branch("Muon_Track_idx",&Muon_Track_idx);
   output_tree->Branch("Muon_hitPattern_pixelLayerwithMeas",&Muon_hitPattern_pixelLayerwithMeas);
   output_tree->Branch("Muon_numberOfMatchedStations",&Muon_numberOfMatchedStations);
   output_tree->Branch("Muon_normChi2",&Muon_normChi2);
   output_tree->Branch("Muon_hitPattern_numberOfValidMuonHits",&Muon_hitPattern_numberOfValidMuonHits);
   output_tree->Branch("Muon_innerTrack_numberofValidHits",&Muon_innerTrack_numberofValidHits);
   output_tree->Branch("Muon_numberOfMatches",&Muon_numberOfMatches);
   output_tree->Branch("Muon_Charge",&Muon_Charge);
   output_tree->Branch("Muon_numberOfChambers",&Muon_numberOfChambers);
   output_tree->Branch("Muon_isPFMuon",&Muon_isPFMuon);
   output_tree->Branch("Muon_numberofValidPixelHits",&Muon_numberofValidPixelHits);
   output_tree->Branch("Muon_trackerLayersWithMeasurement",&Muon_trackerLayersWithMeasurement);


   //================ Electron block ========
   output_tree->Branch("isPatElectron",&doPatElectrons_);
   output_tree->Branch("Electron_p4",&Electron_p4);
   output_tree->Branch("Electron_Poca",&Electron_Poca);
   output_tree->Branch("Electron_Charge",&Electron_Charge);
   output_tree->Branch("Electron_Gsf_deltaEtaEleClusterTrackAtCalo",&Electron_Gsf_deltaEtaEleClusterTrackAtCalo);
   output_tree->Branch("Electron_Gsf_deltaEtaSeedClusterTrackAtCalo",&Electron_Gsf_deltaEtaSeedClusterTrackAtCalo);
   output_tree->Branch("Electron_Gsf_deltaEtaSuperClusterTrackAtVtx",&Electron_Gsf_deltaEtaSuperClusterTrackAtVtx);
   output_tree->Branch("Electron_Gsf_deltaPhiEleClusterTrackAtCalo",&Electron_Gsf_deltaPhiEleClusterTrackAtCalo);
   output_tree->Branch("Electron_Gsf_deltaPhiSeedClusterTrackAtCalo",&Electron_Gsf_deltaPhiSeedClusterTrackAtCalo);
   output_tree->Branch("Electron_Gsf_deltaPhiSuperClusterTrackAtVtx",&Electron_Gsf_deltaPhiSuperClusterTrackAtVtx);
   output_tree->Branch("Electron_Gsf_dr03EcalRecHitSumE",&Electron_Gsf_dr03EcalRecHitSumE);
   output_tree->Branch("Electron_Gsf_dr03HcalDepth1TowerSumEt",&Electron_Gsf_dr03HcalDepth1TowerSumEt);
   output_tree->Branch("Electron_Gsf_dr03HcalDepth1TowerSumEtBc",&Electron_Gsf_dr03HcalDepth1TowerSumEtBc);
   output_tree->Branch("Electron_Gsf_dr03HcalDepth2TowerSumEt",&Electron_Gsf_dr03HcalDepth2TowerSumEt);
   output_tree->Branch("Electron_Gsf_dr03HcalDepth2TowerSumEtBc",&Electron_Gsf_dr03HcalDepth2TowerSumEtBc);
   output_tree->Branch("Electron_Gsf_dr03HcalTowerSumEt",&Electron_Gsf_dr03HcalTowerSumEt);
   output_tree->Branch("Electron_Gsf_dr03HcalTowerSumEtBc",&Electron_Gsf_dr03HcalTowerSumEtBc);
   output_tree->Branch("Electron_Gsf_dr03TkSumPt",&Electron_Gsf_dr03TkSumPt);
   output_tree->Branch("Electron_Gsf_passingCutBasedPreselection",&Electron_Gsf_passingCutBasedPreselection);
   output_tree->Branch("Electron_Gsf_passingMvaPreselection",&Electron_Gsf_passingMvaPreselection);
   output_tree->Branch("Electron_gsftrack_trackerExpectedHitsInner_numberOfLostHits",&Electron_gsftrack_trackerExpectedHitsInner_numberOfLostHits);
   output_tree->Branch("Electron_supercluster_e",&Electron_supercluster_e);
   output_tree->Branch("Electron_supercluster_phi",&Electron_supercluster_phi);
   output_tree->Branch("Electron_supercluster_eta",&Electron_supercluster_eta);
   output_tree->Branch("Electron_supercluster_centroid_x",&Electron_supercluster_centroid_x);
   output_tree->Branch("Electron_supercluster_centroid_y",&Electron_supercluster_centroid_y);
   output_tree->Branch("Electron_supercluster_centroid_z",&Electron_supercluster_centroid_z);
   output_tree->Branch("Electron_Track_idx",&Electron_Track_idx);

   output_tree->Branch("Electron_ecalRecHitSumEt03",&Electron_ecalRecHitSumEt03);	   
   output_tree->Branch("Electron_hcalDepth1TowerSumEt03",&Electron_hcalDepth1TowerSumEt03);   
   output_tree->Branch("Electron_hcalDepth1TowerSumEtBc03",&Electron_hcalDepth1TowerSumEtBc03); 
   output_tree->Branch("Electron_hcalDepth2TowerSumEt03",&Electron_hcalDepth2TowerSumEt03);   
   output_tree->Branch("Electron_hcalDepth2TowerSumEtBc03",&Electron_hcalDepth2TowerSumEtBc03); 
   output_tree->Branch("Electron_tkSumPt03",&Electron_tkSumPt03);		   
   output_tree->Branch("Electron_ecalRecHitSumEt04",&Electron_ecalRecHitSumEt04);	   
   output_tree->Branch("Electron_hcalDepth1TowerSumEt04",&Electron_hcalDepth1TowerSumEt04);   
   output_tree->Branch("Electron_hcalDepth1TowerSumEtBc04",&Electron_hcalDepth1TowerSumEtBc04); 
   output_tree->Branch("Electron_hcalDepth2TowerSumEt04",&Electron_hcalDepth2TowerSumEt04);   
   output_tree->Branch("Electron_hcalDepth2TowerSumEtBc04",&Electron_hcalDepth2TowerSumEtBc04); 
   output_tree->Branch("Electron_tkSumPt04",&Electron_tkSumPt04);		   
   output_tree->Branch("Electron_chargedHadronIso",&Electron_chargedHadronIso);	   
   output_tree->Branch("Electron_neutralHadronIso",&Electron_neutralHadronIso);	   
   output_tree->Branch("Electron_photonIso",&Electron_photonIso);                

   output_tree->Branch("Electron_sigmaIetaIeta",&Electron_sigmaIetaIeta); 	   
   output_tree->Branch("Electron_hadronicOverEm",&Electron_hadronicOverEm);   
   output_tree->Branch("Electron_fbrem",&Electron_fbrem); 		   
   output_tree->Branch("Electron_eSuperClusterOverP",&Electron_eSuperClusterOverP);    
   output_tree->Branch("Electron_ecalEnergy",&Electron_ecalEnergy); 	   
   output_tree->Branch("Electron_trackMomentumAtVtx",&Electron_trackMomentumAtVtx);    
   output_tree->Branch("Electron_numberOfMissedHits",&Electron_numberOfMissedHits);    
   output_tree->Branch("Electron_HasMatchedConversions",&Electron_HasMatchedConversions); 
   output_tree->Branch("RhoIsolationAllInputTags",&RhoIsolationAllInputTags); 
   
   output_tree->Branch("Electron_Track_dR",&Electron_Track_dR);
   // Electron MVA ID
   output_tree->Branch("Electron_Rho_kt6PFJets",&Electron_Rho_kt6PFJets);
   output_tree->Branch("Electron_MVA_discriminator",&Electron_MVA_discriminator);

   //================  PFTau block ==========
   output_tree->Branch("PFTau_p4",&PFTau_p4);
   output_tree->Branch("PFTau_Poca",&PFTau_Poca);
   output_tree->Branch("PFTau_isTightIsolation",&PFTau_isTightIsolation);
   output_tree->Branch("PFTau_isMediumIsolation",&PFTau_isMediumIsolation);
   output_tree->Branch("PFTau_isLooseIsolation",&PFTau_isLooseIsolation);
   output_tree->Branch("PFTau_isTightIsolationDBSumPtCorr",&PFTau_isTightIsolationDBSumPtCorr); 
   output_tree->Branch("PFTau_isMediumIsolationDBSumPtCorr",&PFTau_isMediumIsolationDBSumPtCorr);
   output_tree->Branch("PFTau_isLooseIsolationDBSumPtCorr",&PFTau_isLooseIsolationDBSumPtCorr); 
   output_tree->Branch("PFTau_isVLooseIsolationDBSumPtCorr",&PFTau_isVLooseIsolationDBSumPtCorr);
   output_tree->Branch("PFTau_isHPSAgainstElectronsLoose",&PFTau_isHPSAgainstElectronsLoose); 
   output_tree->Branch("PFTau_isHPSAgainstElectronsMedium",&PFTau_isHPSAgainstElectronsMedium);
   output_tree->Branch("PFTau_isHPSAgainstElectronsTight",&PFTau_isHPSAgainstElectronsTight); 
   output_tree->Branch("PFTau_isHPSAgainstMuonLoose",&PFTau_isHPSAgainstMuonLoose);      
   output_tree->Branch("PFTau_isHPSAgainstMuonMedium",&PFTau_isHPSAgainstMuonMedium);    
   output_tree->Branch("PFTau_isHPSAgainstMuonTight",&PFTau_isHPSAgainstMuonTight);      
   output_tree->Branch("PFTau_isHPSAgainstMuonLoose2",&PFTau_isHPSAgainstMuonLoose2);      
   output_tree->Branch("PFTau_isHPSAgainstMuonMedium2",&PFTau_isHPSAgainstMuonMedium2);    
   output_tree->Branch("PFTau_isHPSAgainstMuonTight2",&PFTau_isHPSAgainstMuonTight2);      
   output_tree->Branch("PFTau_isHPSByDecayModeFinding",&PFTau_isHPSByDecayModeFinding);    

   //  output_tree->Branch("PFTau_HPSPFTauDiscriminationByMVA3rawElectronRejection",&PFTau_HPSPFTauDiscriminationByMVA3rawElectronRejection);		   
   output_tree->Branch("PFTau_HPSPFTauDiscriminationByMVA3LooseElectronRejection",&PFTau_HPSPFTauDiscriminationByMVA3LooseElectronRejection);		   
   output_tree->Branch("PFTau_HPSPFTauDiscriminationByMVA3MediumElectronRejection",&PFTau_HPSPFTauDiscriminationByMVA3MediumElectronRejection);		   
   output_tree->Branch("PFTau_HPSPFTauDiscriminationByMVA3TightElectronRejection",&PFTau_HPSPFTauDiscriminationByMVA3TightElectronRejection);		   
   output_tree->Branch("PFTau_HPSPFTauDiscriminationByMVA3VTightElectronRejection",&PFTau_HPSPFTauDiscriminationByMVA3VTightElectronRejection);		   
   //   output_tree->Branch("PFTau_HPSPFTauDiscriminationByDeadECALElectronRejection",&PFTau_HPSPFTauDiscriminationByDeadECALElectronRejection);		   
   output_tree->Branch("PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits",&PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits);  
   output_tree->Branch("PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits",&PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits); 
   output_tree->Branch("PFTau_HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits",&PFTau_HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits);  
   output_tree->Branch("PFTau_HPSPFTauDiscriminationByLooseIsolationMVA",&PFTau_HPSPFTauDiscriminationByLooseIsolationMVA);			   
   output_tree->Branch("PFTau_HPSPFTauDiscriminationByMediumIsolationMVA",&PFTau_HPSPFTauDiscriminationByMediumIsolationMVA);                      
   output_tree->Branch("PFTau_HPSPFTauDiscriminationByTightIsolationMVA",&PFTau_HPSPFTauDiscriminationByTightIsolationMVA);                      

   output_tree->Branch("PFTau_HPSPFTauDiscriminationByLooseIsolationMVA2",&PFTau_HPSPFTauDiscriminationByLooseIsolationMVA2);			   
   output_tree->Branch("PFTau_HPSPFTauDiscriminationByMediumIsolationMVA2",&PFTau_HPSPFTauDiscriminationByMediumIsolationMVA2);                      
   output_tree->Branch("PFTau_HPSPFTauDiscriminationByTightIsolationMVA2",&PFTau_HPSPFTauDiscriminationByTightIsolationMVA2);                      

   output_tree->Branch("PFTau_hpsDecayMode",&PFTau_hpsDecayMode);
   output_tree->Branch("PFTau_Charge",&PFTau_Charge);
   output_tree->Branch("PFTau_Track_idx",&PFTau_Track_idx);

   output_tree->Branch("PFTau_TIP_primaryVertex_pos",&PFTau_TIP_primaryVertex_pos);
   output_tree->Branch("PFTau_TIP_primaryVertex_cov",&PFTau_TIP_primaryVertex_cov);
   output_tree->Branch("PFTau_TIP_secondaryVertex_pos",&PFTau_TIP_secondaryVertex_pos);
   output_tree->Branch("PFTau_TIP_secondaryVertex_cov",&PFTau_TIP_secondaryVertex_cov);
   output_tree->Branch("PFTau_TIP_secondaryVertex_vtxchi2",&PFTau_TIP_secondaryVertex_vtxchi2);
   output_tree->Branch("PFTau_TIP_secondaryVertex_vtxndof",&PFTau_TIP_secondaryVertex_vtxndof);
   output_tree->Branch("PFTau_TIP_primaryVertex_vtxchi2",&PFTau_TIP_secondaryVertex_vtxchi2);
   output_tree->Branch("PFTau_TIP_primaryVertex_vtxndof",&PFTau_TIP_secondaryVertex_vtxndof);

   output_tree->Branch("PFTau_a1_lvp",&PFTau_a1_lvp);
   output_tree->Branch("PFTau_a1_cov",&PFTau_a1_cov);
   output_tree->Branch("PFTau_a1_charge",&PFTau_a1_charge);
   output_tree->Branch("PFTau_a1_pdgid",&PFTau_a1_pdgid);
   output_tree->Branch("PFTau_a1_B",&PFTau_a1_B);
   output_tree->Branch("PFTau_a1_M",&PFTau_a1_M);

   output_tree->Branch("PFTau_daughterTracks",&PFTau_daughterTracks);
   output_tree->Branch("PFTau_daughterTracks_cov",&PFTau_daughterTracks_cov);
   output_tree->Branch("PFTau_daughterTracks_charge",&PFTau_daughterTracks_charge);
   output_tree->Branch("PFTau_daughterTracks_pdgid",&PFTau_daughterTracks_pdgid);
   output_tree->Branch("PFTau_daughterTracks_B",&PFTau_daughterTracks_B);
   output_tree->Branch("PFTau_daughterTracks_M",&PFTau_daughterTracks_M);
   output_tree->Branch("PFTau_daughterTracks_poca",&PFTau_daughterTracks_poca);

   output_tree->Branch("PFTau_3PS_A1_LV",&PFTau_3PS_A1_LV);
   output_tree->Branch("PFTau_3PS_M_A1",&PFTau_3PS_M_A1);
   output_tree->Branch("PFTau_3PS_M_12",&PFTau_3PS_M_12);
   output_tree->Branch("PFTau_3PS_M_13",&PFTau_3PS_M_13);
   output_tree->Branch("PFTau_3PS_M_23",&PFTau_3PS_M_23);
   output_tree->Branch("PFTau_3PS_Tau_Charge",&PFTau_3PS_Tau_Charge);
   output_tree->Branch("PFTau_3PS_LCchi2",&PFTau_3PS_LCchi2);
   output_tree->Branch("PFTau_3PS_has3ProngSolution",&PFTau_3PS_has3ProngSolution);
   output_tree->Branch("PFTau_3PS_Tau_LV",&PFTau_3PS_Tau_LV);

  //=======  PFJets ===
   output_tree->Branch("isPatJet",&doPatJets_);
   output_tree->Branch("PFJet_p4",&PFJet_p4);
   output_tree->Branch("PFJet_Poca",&PFJet_Poca);
   output_tree->Branch("PFJet_chargedEmEnergy",&PFJet_chargedEmEnergy);
   output_tree->Branch("PFJet_chargedHadronEnergy",&PFJet_chargedHadronEnergy);
   output_tree->Branch("PFJet_chargedHadronMultiplicity",&PFJet_chargedHadronMultiplicity);
   output_tree->Branch("PFJet_chargedMuEnergy",&PFJet_chargedMuEnergy);
   output_tree->Branch("PFJet_chargedMultiplicity",&PFJet_chargedMultiplicity);
   output_tree->Branch("PFJet_electronEnergy",&PFJet_electronEnergy);
   output_tree->Branch("PFJet_electronMultiplicity",&PFJet_electronMultiplicity);
   output_tree->Branch("PFJet_HFEMEnergy",&PFJet_HFEMEnergy);
   output_tree->Branch("PFJet_HFEMMultiplicity",&PFJet_HFEMMultiplicity);
   output_tree->Branch("PFJet_HFHadronEnergy",&PFJet_HFHadronEnergy);
   output_tree->Branch("PFJet_HFHadronMultiplicity",&PFJet_HFHadronMultiplicity);
   output_tree->Branch("PFJet_muonEnergy",&PFJet_muonEnergy);
   output_tree->Branch("PFJet_muonMultiplicity",&PFJet_muonMultiplicity);
   output_tree->Branch("PFJet_neutralEmEnergy",&PFJet_neutralEmEnergy);
   output_tree->Branch("PFJet_neutralHadronEnergy",&PFJet_neutralHadronEnergy);
   output_tree->Branch("PFJet_neutralHadronMultiplicity",&PFJet_neutralHadronMultiplicity);
   output_tree->Branch("PFJet_photonEnergy",&PFJet_photonEnergy);
   output_tree->Branch("PFJet_photonMultiplicity",&PFJet_photonMultiplicity);
   output_tree->Branch("PFJet_jetArea",&PFJet_jetArea); 
   output_tree->Branch("PFJet_maxDistance",&PFJet_maxDistance);
   output_tree->Branch("PFJet_nConstituents",&PFJet_nConstituents);
   output_tree->Branch("PFJet_pileup",&PFJet_pileup);  
   output_tree->Branch("PFJet_etaetaMoment",&PFJet_etaetaMoment);
   output_tree->Branch("PFJet_etaphiMoment",&PFJet_etaphiMoment);
   output_tree->Branch("PFJet_Track_idx",&PFJet_Track_idx);
   output_tree->Branch("PFJet_MatchedHPS_idx",&PFJet_MatchedHPS_idx);
   output_tree->Branch("PFJet_numberOfDaughters",&PFJet_numberOfDaughters);
   output_tree->Branch("PFJet_chargedEmEnergyFraction",&PFJet_chargedEmEnergyFraction);
   output_tree->Branch("PFJet_chargedHadronEnergyFraction",&PFJet_chargedHadronEnergyFraction);
   output_tree->Branch("PFJet_neutralHadronEnergyFraction",&PFJet_neutralHadronEnergyFraction);
   output_tree->Branch("PFJet_neutralEmEnergyFraction",&PFJet_neutralEmEnergyFraction);

   output_tree->Branch("PFJet_partonFlavour",&PFJet_partonFlavour);
   output_tree->Branch("PFJet_bDiscriminator",&PFJet_bDiscriminator);
   output_tree->Branch("PFJet_BTagWeight",&PFJet_BTagWeight);
   //output_tree->Branch("PFJet_bTagAlgorithmName",&PFJet_bTagAlgorithmName);
   //output_tree->Branch("PFJet_bTagAlgorithmValue",&PFJet_bTagAlgorithmValue);

   output_tree->Branch("PFJet_TracksP4",&PFJet_TracksP4);
   output_tree->Branch("PFJet_nTrk",&PFJet_nTrk);



   //================  MET block ==========
   output_tree->Branch("isPatMET",&doPatMET_);
   output_tree->Branch("MET_et",&MET_et);
   output_tree->Branch("MET_pt",&MET_pt);
   output_tree->Branch("MET_phi",&MET_phi);
   output_tree->Branch("MET_sumET",&MET_sumET);
   output_tree->Branch("MET_metSignificance",&MET_metSignificance);
   output_tree->Branch("MET_MuonEtFraction",&MET_MuonEtFraction);
   output_tree->Branch("MET_NeutralEMFraction",&MET_NeutralEMFraction);
   output_tree->Branch("MET_NeutralHadEtFraction",&MET_NeutralHadEtFraction);
   output_tree->Branch("MET_Type6EtFraction",&MET_Type6EtFraction);
   output_tree->Branch("MET_Type7EtFraction",&MET_Type7EtFraction);

   output_tree->Branch("MET_Corr_et",&MET_Corr_et);
   output_tree->Branch("MET_Corr_pt",&MET_Corr_pt);
   output_tree->Branch("MET_Corr_phi",&MET_Corr_phi);
   output_tree->Branch("MET_Corr_SumET",&MET_Corr_SumET);

   //=============== Event Block ==============
   output_tree->Branch("Event_EventNumber",&Event_EventNumber);
   output_tree->Branch("Event_RunNumber",&Event_RunNumber);
   output_tree->Branch("Event_bunchCrossing",&Event_bunchCrossing);
   output_tree->Branch("Event_orbitNumber",&Event_orbitNumber);
   output_tree->Branch("Event_luminosityBlock",&Event_luminosityBlock);	 
   output_tree->Branch("Event_isRealData",&Event_isRealData);

   output_tree->Branch("PileupInfo_NumInteractions_nm1",&PileupInfo_NumInteractions_nm1);
   output_tree->Branch("PileupInfo_NumInteractions_n0",&PileupInfo_NumInteractions_n0);
   output_tree->Branch("PileupInfo_NumInteractions_np1",&PileupInfo_NumInteractions_np1);
   output_tree->Branch("EvtWeight3D",&EvtWeight3D);

   //=============== Track Block ==============
   output_tree->Branch("Track_p4",&Track_p4);
   output_tree->Branch("Track_Poca",&Track_Poca);
   output_tree->Branch("Track_charge",&Track_charge);
   output_tree->Branch("Track_chi2",&Track_chi2);
   output_tree->Branch("Track_ndof",&Track_ndof);
   output_tree->Branch("Track_numberOfLostHits",&Track_numberOfLostHits);
   output_tree->Branch("Track_numberOfValidHits",&Track_numberOfValidHits);
   output_tree->Branch("Track_qualityMask",&Track_qualityMask);
   output_tree->Branch("Track_par",&Track_par);
   output_tree->Branch("Track_parCov",&Track_parCov);

   //=============== MC Block ==============

   output_tree->Branch("GenEventInfoProduct_signalProcessID",&GenEventInfoProduct_signalProcessID);
   output_tree->Branch("GenEventInfoProduct_weight",&GenEventInfoProduct_weight);
   output_tree->Branch("GenEventInfoProduct_weights",&GenEventInfoProduct_weights);
   output_tree->Branch("GenEventInfoProduct_qScale",&GenEventInfoProduct_qScale);
   output_tree->Branch("GenEventInfoProduct_alphaQED",&GenEventInfoProduct_alphaQED);
   output_tree->Branch("GenEventInfoProduct_alphaQCD",&GenEventInfoProduct_alphaQCD);

   if(do_MCComplete_){
     output_tree->Branch("MC_p4",&MC_p4);
     output_tree->Branch("MC_pdgid",&MC_pdgid);
     output_tree->Branch("MC_charge",&MC_charge);
     output_tree->Branch("MC_midx",&MC_midx);
     output_tree->Branch("MC_childpdgid",&MC_childpdgid);
   }
   if(do_MCSummary_){
     output_tree->Branch("MCSignalParticle_p4",&MCSignalParticle_p4);
     output_tree->Branch("MCSignalParticle_pdgid",&MCSignalParticle_pdgid);
     output_tree->Branch("MCSignalParticle_charge",&MCSignalParticle_charge);
     output_tree->Branch("MCSignalParticle_Poca",&MCSignalParticle_Poca);
     output_tree->Branch("MCSignalParticle_Tauidx",&MCSignalParticle_Tauidx);
     output_tree->Branch("MCTauandProd_p4",&MCTauandProd_p4);
     output_tree->Branch("MCTauandProd_Vertex",&MCTauandProd_Vertex);

     output_tree->Branch("MCTauandProd_pdgid",&MCTauandProd_pdgid);
     output_tree->Branch("MCTauandProd_midx",&MCTauandProd_midx);
     output_tree->Branch("MCTauandProd_charge",&MCTauandProd_charge);
     output_tree->Branch("MCTau_JAK",&MCTau_JAK);   
     output_tree->Branch("MCTau_DecayBitMask",&MCTau_DecayBitMask); 
   }

   //================= Trigger Block ===============
   output_tree->Branch("HTLTriggerName",&HTLTriggerName);
   output_tree->Branch("TriggerAccept",&TriggerAccept);
   output_tree->Branch("TriggerError",&TriggerError);
   output_tree->Branch("TriggerWasRun",&TriggerWasRun);
   output_tree->Branch("HLTPrescale",&HLTPrescale);
   output_tree->Branch("NHLTL1GTSeeds",&NHLTL1GTSeeds);
  output_tree->Branch("L1SEEDPrescale",&L1SEEDPrescale);
   output_tree->Branch("L1SEEDInvalidPrescale",&L1SEEDInvalidPrescale);
   output_tree->Branch("L1SEEDisTechBit",&L1SEEDisTechBit);
   output_tree->Branch("MuonTriggerMatch",&MuonTriggerMatch);
   output_tree->Branch("ElectronTriggerMatch",&ElectronTriggerMatch);
   output_tree->Branch("JetTriggerMatch",&JetTriggerMatch);
   output_tree->Branch("TauTriggerMatch",&TauTriggerMatch);
   output_tree->Branch("HLTTrigger_objs_Pt",&HLTTrigger_objs_Pt);
   output_tree->Branch("HLTTrigger_objs_Eta",&HLTTrigger_objs_Eta);
   output_tree->Branch("HLTTrigger_objs_Phi",&HLTTrigger_objs_Phi);

   output_tree->Branch("L1TriggerName",&L1TriggerName);
   output_tree->Branch("L1TriggerDecision",&L1TriggerDecision);
   output_tree->Branch("L1ErrorCode",&L1ErrorCode);
   output_tree->Branch("L1Prescale",&L1Prescale);

 } 


 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 //
 // std::vector<bool> TauNtuple::CheckTauDiscriminators(std::vector<edm::Handle<reco::PFTauDiscriminator> > tauDiscriminators, reco::PFTauRef tauRef)
 //
 // checks that tau candidate pass kinematic fit and KinFit quality criteria 
 // returns vector of bool variables, first variable true/false if tau candidate pass/fail kinematic fit
 // second variable is true/false if  refitted tau candidate pass/fail quality requirements
 std::vector<bool>  TauNtuple::CheckTauDiscriminators(std::vector<edm::Handle<reco::PFTauDiscriminator> > tauDiscriminators, const reco::PFTauRef tauRef){
   std::vector<bool> output_pair;
   bool discriminateByKinFit = false;
   bool discriminateByKinQC  = false;
   int iDiscr =0;
   for (std::vector<edm::Handle<reco::PFTauDiscriminator> >::const_iterator discr = tauDiscriminators.begin(); discr!=tauDiscriminators.end(); ++discr) {
     iDiscr = iDiscr + (**discr)[tauRef];
     }
   if(iDiscr == 1 )discriminateByKinFit=true;
   if(iDiscr == 2 ){discriminateByKinFit=true;discriminateByKinQC=true;}
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
 reco::PFTauRef TauNtuple::getMatchedHPSTau(edm::Handle<std::vector<reco::PFTau> > & HPStaus,   std::vector<float>   &UnmodifiedTau, unsigned int &match){
   TLorentzVector TauVisible;
   TauVisible.SetE(UnmodifiedTau.at(0));
   TauVisible.SetPx(UnmodifiedTau.at(1));
   TauVisible.SetPy(UnmodifiedTau.at(2));
   TauVisible.SetPz(UnmodifiedTau.at(3));
   reco::PFTauRef MatchedHPSTau;
   double deltaR = 999;
   match=0;
    for ( unsigned int iTau = 0; iTau < HPStaus->size(); ++iTau ) {
     reco::PFTauRef HPStauCandidate(HPStaus, iTau);
     double dr=sqrt( pow(DeltaPhi(HPStauCandidate->p4().Phi(),TauVisible.Phi()),2) + pow(HPStauCandidate->p4().Eta() - TauVisible.Eta(),2));
     if(dr < deltaR){
       deltaR = dr;
       match=iTau;
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
 reco::PFTauRef TauNtuple::getHPSTauMatchedToJet(edm::Handle<std::vector<reco::PFTau> > & HPStaus,   std::vector<float>   &Jet, unsigned int &match){
   TLorentzVector Jetp4;;
   Jetp4.SetE(Jet.at(0));
   Jetp4.SetPx(Jet.at(1));
   Jetp4.SetPy(Jet.at(2));
   Jetp4.SetPz(Jet.at(3));

   reco::PFTauRef MatchedHPSTau;
   double deltaR = 999;
   match=0;
    for ( unsigned int iTau = 0; iTau < HPStaus->size(); ++iTau ) {
     reco::PFTauRef HPStauCandidate(HPStaus, iTau);
     double dr=sqrt( pow(DeltaPhi(HPStauCandidate->p4().Phi(),Jetp4.Phi()),2) + pow(HPStauCandidate->p4().Eta() - Jetp4.Eta(),2));
     if(dr < deltaR){
       deltaR = dr;
       match=iTau;
       MatchedHPSTau = HPStauCandidate;
     }

   }
   return MatchedHPSTau;
 }




 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 //
 // bool TauNtuple::getTrackMatch(edm::Handle< std::vector<reco::Track>  > &trackCollection, reco::TrackRef &refTrack, int &match)
 //
 // finds track match for a given TrackRef
 // returns true  if the matching track is found in the collection and sets match to the index of the found track
 // retruns false if on match is found in the collection and match is set to -1.
 bool TauNtuple::getTrackMatch(edm::Handle< std::vector<reco::Track>  > &trackCollection, reco::TrackRef &refTrack, int &match){
   match=-1;
   for(unsigned int iTrack = 0; iTrack < trackCollection->size(); iTrack++) {
     reco::TrackRef Track(trackCollection, iTrack);
     if(refTrack==Track){
       match=iTrack;
       return true;  
     }
   }
   return false;
 }
 
 bool TauNtuple::getTrackMatch(edm::Handle< std::vector<reco::Track> > &trackCollection, reco::GsfTrackRef &refTrack, int &match){
	 match=-1;
	 for(unsigned int iTrack = 0; iTrack < trackCollection->size(); iTrack++){
		 reco::TrackRef Track(trackCollection, iTrack);
		 double dr = TMath::Sqrt(TMath::Power(refTrack->eta()-Track->eta(),2)+TMath::Power(refTrack->phi()-Track->phi(),2));
		 Electron_Track_dR.push_back(dr);
		 if(dr<0.1){
			 match=iTrack;
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
double TauNtuple::DeltaPhi(double phi1, double phi2){
  double dphi=fabs(phi1-phi2);
  if (dphi>TMath::Pi())   dphi=2*TMath::Pi()-dphi;
  double sign=1;
  if(phi1-phi2<0) sign=-1;
  return dphi*sign;
}


// ------------ method called once each job just after ending the event loop  ------------
void TauNtuple::endJob() {
  std::cout<<" No Of event processed: "<< cnt_ << std::endl;
  output->Write();
  output->Close();
}  

// ------------ method called when starting to processes a run  ------------
void TauNtuple::beginRun(edm::Run& Run, edm::EventSetup const& Setup){ 
  bool changed(true);
  TriggerOK=true;
  if (hltConfig_.init(Run,Setup,processName_,changed)) {
    // if init returns TRUE, initialisation has succeeded!
    if (changed) {
      // The HLT config has actually changed wrt the previous Run, hence rebook your
      // histograms or do anything else dependent on the revised HLT config
      //     std::cout << "Initalizing HLTConfigProvider"  << std::endl;
    }
  } 
  else{
    // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
    // with the file and/or code and needs to be investigated!
    //  std::cout << " HLT config extraction failure with process name " << processName_ << std::endl;
    // In this case, all access methods will return empty values!
    TriggerOK=false;
  }
  if(!Run.getByLabel(l1GtTriggerMenuLite_.label(), triggerMenuLite_)){
    std::cout << " l1GtTrigger config extraction failure "  << std::endl;
    TriggerOK=false;
  }
} 

// ------------ method called when ending the processing of a run  ------------
void 
TauNtuple::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TauNtuple::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TauNtuple::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TauNtuple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
} 

void 
TauNtuple::ClearEvent(){
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
  Muon_Charge.clear();
  Muon_Track_idx.clear();

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

  //=======  Electrons ===
  Electron_p4.clear();
  Electron_Poca.clear();
  Electron_Charge.clear();
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


  Electron_sigmaIetaIeta.clear();          	   
  Electron_hadronicOverEm.clear();          	   
  Electron_fbrem.clear();          		   
  Electron_eSuperClusterOverP.clear();             
  Electron_ecalEnergy.clear();          	   
  Electron_trackMomentumAtVtx.clear();           
  Electron_numberOfMissedHits.clear();           
  Electron_HasMatchedConversions.clear();    

  Electron_Track_dR.clear();
  Electron_MVA_discriminator.clear();


  //=======  PFJets ===
   PFJet_p4.clear();
   PFJet_Poca.clear();
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

   PFJet_partonFlavour.clear();
   PFJet_bDiscriminator.clear();
   PFJet_BTagWeight.clear();

   //PFJet_bTagAlgorithmName.clear();
   //PFJet_bTagAlgorithmValue.clear();

   PFJet_TracksP4.clear();
   PFJet_nTrk.clear();

   //=============== Track Block ==============
   Track_p4.clear();
   Track_Poca.clear();
   Track_charge.clear();
   Track_chi2.clear();
   Track_ndof.clear();
   Track_numberOfLostHits.clear();
   Track_numberOfValidHits.clear();
   Track_qualityMask.clear();
   Track_par.clear();
   Track_parCov.clear();

   // Event Block
   EvtWeight3D=0;


   //=============== MC Block ==============
   
   GenEventInfoProduct_weights.clear();
   GenEventInfoProduct_signalProcessID=0;
   GenEventInfoProduct_weight=0;
   GenEventInfoProduct_qScale=0;
   GenEventInfoProduct_alphaQED=0;
   GenEventInfoProduct_alphaQCD=0;


  if(do_MCComplete_){
    MC_p4.clear();
    MC_pdgid.clear();
    MC_charge.clear();
    MC_midx.clear();
  }
  if(do_MCSummary_){
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
  MC_childpdgid.clear();
  }

  //======================Trigger Block ======================
  HTLTriggerName.clear();
  HLTPrescale.clear();
  NHLTL1GTSeeds.clear();
  L1SEEDPrescale.clear();
  L1SEEDInvalidPrescale.clear();
  L1SEEDisTechBit.clear();
  TriggerAccept.clear();
  MuonTriggerMatch.clear();
  ElectronTriggerMatch.clear();
  JetTriggerMatch.clear();
  TauTriggerMatch.clear();
  TriggerError.clear();
  TriggerWasRun.clear();
  HLTTrigger_objs_Pt.clear();
  HLTTrigger_objs_Eta.clear();
  HLTTrigger_objs_Phi.clear();

  L1TriggerName.clear();
  L1TriggerDecision.clear();
  L1ErrorCode.clear();
  L1Prescale.clear();


}








