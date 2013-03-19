#include "TauDataFormat/TauNtuple/interface/TauNtuple.h"
#include "TauDataFormat/TauNtuple/interface/TauDecay_CMSSW.h"
#include "TauDataFormat/TauNtuple/interface/PdtPdgMini.h"
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

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/ParticleFlowReco/interface/GsfPFRecTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include <DataFormats/EgammaReco/interface/SuperCluster.h>

#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/PatCandidates/interface/MET.h>

TauNtuple::TauNtuple(const edm::ParameterSet& iConfig):
  primVtxTag_( iConfig.getParameter<edm::InputTag>( "primVtx" ) ),
  muonsTag_(iConfig.getParameter<edm::InputTag>( "muons" )),
  hpsTauProducer_( iConfig.getParameter<edm::InputTag>( "hpsTauProducer" ) ),
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
  //  kinTausTag_( iConfig.getParameter<edm::InputTag>( "kinematicTaus" ) ),
  KinFitAdvanced_( iConfig.getParameter<edm::InputTag>( "kinematicTausAdvanced" ) ),
  tauPrimaryVtx_( iConfig.getParameter<edm::InputTag>( "tauPrimaryVtx" ) ),
  pfjetsTag_( iConfig.getParameter<edm::InputTag>( "pfjets" ) ),
  PFElectronTag_( iConfig.getParameter<edm::InputTag>( "pfelectrons" ) ),
  rhoIsolAllInputTag_( iConfig.getParameter<edm::InputTag>( "RhoIsolAllInputTag" ) ),
  generalTracks_(iConfig.getParameter<edm::InputTag>( "generalTracks" )),
  gensrc_(iConfig.getParameter<edm::InputTag>( "gensrc" )),
  GenEventInfo_(iConfig.getParameter<edm::InputTag>("GenEventInfo")),
  discriminators_( iConfig.getParameter< std::vector<std::string> >("discriminators") ),
  ScaleFactor_(iConfig.getUntrackedParameter<std::string>("ScaleFactor")),
  PUInputFile_(iConfig.getUntrackedParameter<std::string>("PUInputFile")),
  PUInputHistoMC_(iConfig.getUntrackedParameter<std::string>("PUInputHistoMC")),
  PUInputHistoData_(iConfig.getUntrackedParameter<std::string>("PUInputHistoData")),
  PUOutputFile_(iConfig.getUntrackedParameter("PUOutputFile",(std::string)("Weight3D.root"))),
  do_MCSummary_(iConfig.getUntrackedParameter("do_MCSummary",(bool)(true))),
  do_MCComplete_(iConfig.getUntrackedParameter("do_MCComplete",(bool)(false))),
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
  doKinFitTaus_(iConfig.getUntrackedParameter("doKinFitTaus",(bool)(true))),
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
  if(doKinFitTaus_)fillKinFitTaus(iEvent,iSetup,trackCollection);
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
	    if(abs(dau->pdgId())==PdtPdgMini::tau_minus){
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
    for(reco::Vertex::trackRef_iterator iTrack=pv.tracks_begin(); iTrack<pv.tracks_end();iTrack++){
      int match(-1);
      reco::TrackRef refTrack=iTrack->castTo<reco::TrackRef>();
      if( refTrack.isNonnull() ) {
	getTrackMatch(trackCollection,refTrack,match);
	matches.push_back(match);
	TrackWeights.push_back(pv.trackWeight(refTrack));
      }
    }
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
     float pionmass=0.13957018;
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



 void 
 TauNtuple::fillPFTaus(edm::Event& iEvent, const edm::EventSetup& iSetup,edm::Handle< std::vector<reco::Track>  > &trackCollection){
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

   for ( unsigned iPFTau = 0; iPFTau < HPStaus->size(); ++iPFTau ) {


 
     reco::PFTauRef HPStauCandidate(HPStaus, iPFTau);
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


     PFTau_isHPSByDecayModeFinding.push_back((*HPSByDecayModeFinding)[HPStauCandidate]);


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




     PFTau_hpsDecayMode.push_back(HPStauCandidate->decayMode());
     PFTau_Charge.push_back(HPStauCandidate->charge());


     reco::PFCandidateRefVector ChargedHadrCand=HPStauCandidate->signalPFChargedHadrCands();
     std::vector<int> matches;
     for(unsigned int i=0; i<ChargedHadrCand.size();i++){
       reco::PFCandidateRef Cand(ChargedHadrCand,i);
       reco::TrackRef refTrack=Cand.get()->trackRef();
       if( refTrack.isNonnull() ) {
	 int match(-1);
	 getTrackMatch(trackCollection,refTrack,match);
	 matches.push_back(match);
       }
     }
     PFTau_Track_idx.push_back(matches);
   }
 }


void  TauNtuple::fillKinFitTaus(edm::Event& iEvent, const edm::EventSetup& iSetup,edm::Handle< std::vector<reco::Track>  > &trackCollection){

  //======== Get Reduced Vertex ================
   edm::Handle<reco::VertexCollection> RedprimVtxs;
   iEvent.getByLabel( tauPrimaryVtx_, RedprimVtxs);

   unsigned int nReducedVtxs=RedprimVtxs->size();
   int ndim=3;
   if(nReducedVtxs>1) nReducedVtxs=1;
   for(unsigned int i=0;i<nReducedVtxs;i++){
     const reco::Vertex &pv = RedprimVtxs->at(i);
     ReducedVtx_isFake.push_back(pv.isFake());
     ReducedVtx_chi2.push_back(pv.chi2());
     ReducedVtx_ndof.push_back(pv.ndof());
     ReducedVtx_x.push_back(pv.x());
     ReducedVtx_y.push_back(pv.y());
     ReducedVtx_z.push_back(pv.z());
     std::vector<std::vector<float> > iReducedVtx_Cov;
     for(int j=0;j<ndim;j++){
       iReducedVtx_Cov.push_back(std::vector<float>());
       for(int k=0;k<=j;k++){
 	iReducedVtx_Cov.at(j).push_back(pv.covariance(j,k));
       }
     }
     ReducedVtx_Cov.push_back(iReducedVtx_Cov);
     std::vector<int> matches;
     for(reco::Vertex::trackRef_iterator iTrack=pv.tracks_begin(); iTrack<pv.tracks_end();iTrack++){
       int match(-1);
       reco::TrackRef refTrack=iTrack->castTo<reco::TrackRef>();
       if( refTrack.isNonnull() ) {
 	getTrackMatch(trackCollection,refTrack,match);
 	matches.push_back(match);
       }
     }
     ReducedVtx_Track_idx.push_back(matches);
   }
  

  //========= HPS taus for matching issues                                                                                                                                                                                                  
  edm::Handle<std::vector<reco::PFTau> > HPStaus;
  iEvent.getByLabel(hpsTauProducer_, HPStaus);
  //========= HPS taus for matching issues       
  

  //========= Get Tau Discriminators =========//
  std::vector<edm::Handle<reco::PFTauDiscriminator> > tauDiscriminators;
  for(std::vector<std::string>::const_iterator discr=discriminators_.begin(); discr!=discriminators_.end(); ++discr) {
    edm::Handle<reco::PFTauDiscriminator> tmpHandle;
    iEvent.getByLabel("KinematicTauBasicProducer", *discr, tmpHandle);
    tauDiscriminators.push_back(tmpHandle);
  }


  //======== Get Tau Collection =====//
  edm::Handle<reco::PFTauCollection> tauCollection;
  iEvent.getByLabel(pfjetsTag_, tauCollection);


  //================== KinematicFit Info ===================
  edm::Handle<SelectedKinematicDecayCollection> selected;
  iEvent.getByLabel(KinFitAdvanced_, selected);

  //start loop over KinFit decays 
  unsigned int tauindex=0;
  unsigned int ntaus=0;
  for(SelectedKinematicDecayCollection::const_iterator decay = selected->begin(); decay != selected->end(); ++decay, tauindex++){
    ntaus = tauindex;

    SelectedKinematicDecay KFTau=(*decay);
    
      std::vector<std::vector<float> > iKFTau_TauVis_p4;
      std::vector<float> iKFTau_a1Initial_p4;
      std::vector<std::vector<float> > iKFTau_TauFit_p4;
      std::vector<std::vector<float> > iKFTau_Neutrino_p4;

      std::vector<std::vector<float> > iKFTau_TauFitInitial_p4;
      std::vector<std::vector<float> > iKFTau_NeutrinoInitial_p4;

      std::vector<float> iKFTau_Fit_PrimaryVertex;
      std::vector<float> iKFTau_Fit_InitialPrimaryVertex;
      std::vector<float> iKFTau_Fit_InitialPrimaryVertexReFit;
      std::vector<float> iKFTau_Fit_InitialSecondaryVertex;
      std::vector<std::vector<float> > iKFTau_Fit_SecondaryVertex;
      std::vector<float> iKFTau_Fit_TauPrimVtx;

      iKFTau_a1Initial_p4.push_back(KFTau.Initial_a1_p4().E());
      iKFTau_a1Initial_p4.push_back(KFTau.Initial_a1_p4().Px());
      iKFTau_a1Initial_p4.push_back(KFTau.Initial_a1_p4().Py());
      iKFTau_a1Initial_p4.push_back(KFTau.Initial_a1_p4().Pz());


      iKFTau_Fit_PrimaryVertex.push_back(KFTau.PrimaryVertexReFitAndRotated().position().x());
      iKFTau_Fit_PrimaryVertex.push_back(KFTau.PrimaryVertexReFitAndRotated().position().y());
      iKFTau_Fit_PrimaryVertex.push_back(KFTau.PrimaryVertexReFitAndRotated().position().z());
      
      iKFTau_Fit_InitialPrimaryVertex.push_back(KFTau.InitialPrimaryVertex().position().x());
      iKFTau_Fit_InitialPrimaryVertex.push_back(KFTau.InitialPrimaryVertex().position().y());
      iKFTau_Fit_InitialPrimaryVertex.push_back(KFTau.InitialPrimaryVertex().position().z());

      iKFTau_Fit_InitialPrimaryVertexReFit.push_back(KFTau.InitialPrimaryVertexReFit().position().x());
      iKFTau_Fit_InitialPrimaryVertexReFit.push_back(KFTau.InitialPrimaryVertexReFit().position().y());
      iKFTau_Fit_InitialPrimaryVertexReFit.push_back(KFTau.InitialPrimaryVertexReFit().position().z());

      
      iKFTau_Fit_InitialSecondaryVertex.push_back(KFTau.InitialSecondaryVertex().position().x());
      iKFTau_Fit_InitialSecondaryVertex.push_back(KFTau.InitialSecondaryVertex().position().y());
      iKFTau_Fit_InitialSecondaryVertex.push_back(KFTau.InitialSecondaryVertex().position().z());
      
      // KFTau_Fit_TauCharge.push_back();


      
      KFTau_Daughter_pdgid.push_back(std::vector<int>());
      KFTau_Daughter_charge.push_back(std::vector<int>());

      KFTau_Daughter_par.push_back(std::vector<std::vector<float>  >());
      KFTau_Daughter_parCov.push_back(std::vector<std::vector<float>  >());
      KFTau_Daughter_inputpar.push_back(std::vector<std::vector<float> > ());
      KFTau_Daughter_inputparCov.push_back(std::vector<std::vector<float> > ());
    

      std::vector<float>   iKFTau_Daughter_par;
      std::vector<float>   iKFTau_Daughter_parCov;
      std::vector<float>   iKFTau_Daughter_inputpar;
      std::vector<float>   iKFTau_Daughter_inputparCov;

        

   
      int taucharge;
      //    Start loop over ambiguity points.

      //   std::vector<float> 
	KFTau_discriminatorByKFit.push_back(std::vector<int>());
	KFTau_discriminatorByQC.push_back(std::vector<int>());

	KFTau_Fit_Chi2Prob.push_back(std::vector<float>());
	KFTau_Fit_BDTVal.push_back(std::vector<float>());
	KFTau_Fit_ndf.push_back(std::vector<float>());
	KFTau_Fit_csum.push_back(std::vector<int>());
	KFTau_Fit_iterations.push_back(std::vector<int>());

	KFTau_Fit_TauEnergyFraction.push_back(std::vector<float>());
	KFTau_Fit_RefitVisibleMass.push_back(std::vector<float>());
	KFTau_Fit_chi2.push_back(std::vector<float>());
	KFTau_Fit_PV_PV_significance.push_back(std::vector<float>());
	KFTau_Fit_SV_PV_significance.push_back(std::vector<float>());


// 	KFTau_pions.push_back(std::vector<std::vector<float> >());
// 	KFTau_Initial_pions.push_back(std::vector<std::vector<float> > ());


      for(unsigned int ambiguity=0; ambiguity<SelectedKinematicDecay::NAmbiguity; ambiguity++){
	




	iKFTau_TauVis_p4.push_back(std::vector<float>());
	iKFTau_TauFit_p4.push_back(std::vector<float>());
	iKFTau_Neutrino_p4.push_back(std::vector<float>());
	iKFTau_TauFitInitial_p4.push_back(std::vector<float>());
	iKFTau_NeutrinoInitial_p4.push_back(std::vector<float>());
	iKFTau_Fit_SecondaryVertex.push_back(std::vector<float>());
	
	
	
	



	KFTau_Fit_chi2.at(tauindex).push_back(KFTau.chi2(ambiguity));
	KFTau_Fit_ndf.at(tauindex).push_back(KFTau.ndf(ambiguity));
	KFTau_Fit_csum.at(tauindex).push_back(KFTau.csum(ambiguity));
	KFTau_Fit_iterations.at(tauindex).push_back(KFTau.iterations(ambiguity));
	KFTau_Fit_BDTVal.at(tauindex).push_back(KFTau.BDTVal(ambiguity));
	
	KFTau_Fit_TauEnergyFraction.at(tauindex).push_back(KFTau.energyTFraction(ambiguity));
	KFTau_Fit_RefitVisibleMass.at(tauindex).push_back(KFTau.a1Mass(ambiguity));
	KFTau_Fit_Chi2Prob.at(tauindex).push_back(KFTau.chi2prob(ambiguity));
	KFTau_Fit_PV_PV_significance.at(tauindex).push_back(KFTau.vtxSignPVRotPVRed(ambiguity));
	KFTau_Fit_SV_PV_significance.at(tauindex).push_back(KFTau.vtxSignPVRotSV(ambiguity));
	std::cout<<"KFTau.vtxSignPVRotSV(ambiguity) "<<KFTau.vtxSignPVRotSV(ambiguity)<<std::endl;
	
	iKFTau_TauVis_p4.at(ambiguity).push_back(KFTau.a1_p4(ambiguity).E());
	iKFTau_TauVis_p4.at(ambiguity).push_back(KFTau.a1_p4(ambiguity).Px());
	iKFTau_TauVis_p4.at(ambiguity).push_back(KFTau.a1_p4(ambiguity).Py());
	iKFTau_TauVis_p4.at(ambiguity).push_back(KFTau.a1_p4(ambiguity).Pz());
	
	//	std::cout<<"ambiguity: "<<ambiguity<<"  tauPt  "<<KFTau.Tau(ambiguity).Pt()<<std::endl;
	iKFTau_TauFit_p4.at(ambiguity).push_back(KFTau.Tau(ambiguity).E());
	iKFTau_TauFit_p4.at(ambiguity).push_back(KFTau.Tau(ambiguity).Px());
	iKFTau_TauFit_p4.at(ambiguity).push_back(KFTau.Tau(ambiguity).Py());
	iKFTau_TauFit_p4.at(ambiguity).push_back(KFTau.Tau(ambiguity).Pz());

	iKFTau_Neutrino_p4.at(ambiguity).push_back(KFTau.Neutrino(ambiguity).E());
	iKFTau_Neutrino_p4.at(ambiguity).push_back(KFTau.Neutrino(ambiguity).Px());
	iKFTau_Neutrino_p4.at(ambiguity).push_back(KFTau.Neutrino(ambiguity).Py());
	iKFTau_Neutrino_p4.at(ambiguity).push_back(KFTau.Neutrino(ambiguity).Pz());
	
	
	iKFTau_TauFitInitial_p4.at(ambiguity).push_back(KFTau.InitialTauGuess(ambiguity).E());
	iKFTau_TauFitInitial_p4.at(ambiguity).push_back(KFTau.InitialTauGuess(ambiguity).Px());
	iKFTau_TauFitInitial_p4.at(ambiguity).push_back(KFTau.InitialTauGuess(ambiguity).Py());
	iKFTau_TauFitInitial_p4.at(ambiguity).push_back(KFTau.InitialTauGuess(ambiguity).Pz());
	
	iKFTau_NeutrinoInitial_p4.at(ambiguity).push_back(KFTau.InitialNeutrinoGuess(ambiguity).E());
	iKFTau_NeutrinoInitial_p4.at(ambiguity).push_back(KFTau.InitialNeutrinoGuess(ambiguity).Px());
	iKFTau_NeutrinoInitial_p4.at(ambiguity).push_back(KFTau.InitialNeutrinoGuess(ambiguity).Py());
	iKFTau_NeutrinoInitial_p4.at(ambiguity).push_back(KFTau.InitialNeutrinoGuess(ambiguity).Pz());
	
	
	iKFTau_Fit_SecondaryVertex.at(ambiguity).push_back(KFTau.SecondaryVertex(ambiguity).position().x());
	iKFTau_Fit_SecondaryVertex.at(ambiguity).push_back(KFTau.SecondaryVertex(ambiguity).position().y());
	iKFTau_Fit_SecondaryVertex.at(ambiguity).push_back(KFTau.SecondaryVertex(ambiguity).position().z());
	


	if(KFTau.discriminators(ambiguity).find("PFRecoTauDiscriminationByKinematicFit")->second) KFTau_discriminatorByKFit.at(tauindex).push_back(1);
	else KFTau_discriminatorByKFit.at(tauindex).push_back(0);
	
	if(KFTau.discriminators(ambiguity).find("PFRecoTauDiscriminationByKinematicFitQuality")->second){ KFTau_discriminatorByQC.at(tauindex).push_back(1); }
	else{ KFTau_discriminatorByQC.at(tauindex).push_back(0); }
	

	//	std::cout<<"ambiguity "<< ambiguity << "  ByKfit  "<<KFTau.discriminators(ambiguity).find("PFRecoTauDiscriminationByKinematicFit")->second<<	"  byQC  "<< KFTau.discriminators(ambiguity).find("PFRecoTauDiscriminationByKinematicFitQuality")->second<<" vector size   " <<KFTau_discriminatorByQC.at(ambiguity).size() <<"check ambig size  " << KFTau_Fit_chi2.at(ambiguity).size() <<std::endl;

	const SelectedKinematicParticleCollection& Particles =decay->particles(ambiguity);
	
	for(std::vector<SelectedKinematicParticle>::const_iterator iParticle = Particles.begin(); iParticle != Particles.end(); ++iParticle){
	    if(fabs(iParticle->pdgid())==15){

	      taucharge = iParticle->charge();

	    }

	  if(ambiguity == 0){ // supposed to  be the same for all ambuguity points

	    //First Store Tau's

	    int d_pdgid=0;
// 	    if(iParticle->name()=="neutrino")                              d_pdgid=PdtPdgMini::nu_tau;
// 	    else if(iParticle->name()=="tau"  && iParticle->charge()==1)   d_pdgid=PdtPdgMini::tau_plus;
// 	    else if(iParticle->name()=="tau"  && iParticle->charge()==-1)  d_pdgid=PdtPdgMini::tau_minus;
// 	    else if(iParticle->name()=="a1"   && iParticle->charge()==1)   d_pdgid=PdtPdgMini::a_1_plus;
// 	    else if(iParticle->name()=="a1"   && iParticle->charge()==-1)  d_pdgid=PdtPdgMini::a_1_minus;
// 	    else if(iParticle->name()=="pion" && iParticle->charge()==1)   d_pdgid=PdtPdgMini::pi_plus;
// 	    else if(iParticle->name()=="pion" && iParticle->charge()==-1)  d_pdgid=PdtPdgMini::pi_minus;

	    d_pdgid = iParticle->pdgid();
// 	    if(iParticle->name()=="neutrino")                              d_pdgid=PdtPdgMini::nu_tau;
// 	    else if(iParticle->name()=="tau"  && iParticle->charge()==1)   d_pdgid=PdtPdgMini::tau_plus;
// 	    else if(iParticle->name()=="tau"  && iParticle->charge()==-1)  d_pdgid=PdtPdgMini::tau_minus;
// 	    else if(iParticle->name()=="a1"   && iParticle->charge()==1)   d_pdgid=PdtPdgMini::a_1_plus;
// 	    else if(iParticle->name()=="a1"   && iParticle->charge()==-1)  d_pdgid=PdtPdgMini::a_1_minus;
// 	    else if(iParticle->name()=="pion" && iParticle->charge()==1)   d_pdgid=PdtPdgMini::pi_plus;
// 	    else if(iParticle->name()=="pion" && iParticle->charge()==-1)  d_pdgid=PdtPdgMini::pi_minus;
	    

	    
	    KFTau_Daughter_pdgid.at(ntaus).push_back(d_pdgid);
	    KFTau_Daughter_charge.at(ntaus).push_back(iParticle->charge());

	    
	  }
	  

 	  for(int j=0;j<iParticle->matrix().GetNrows();j++){
 	    iKFTau_Daughter_par.push_back(iParticle->parameters()(j));

 	    iKFTau_Daughter_inputpar.push_back(iParticle->input_parameters()(j));
 	    for(int k=0;k<=j;k++){

 	      iKFTau_Daughter_parCov.push_back(iParticle->matrix()(j,k));
 	      iKFTau_Daughter_inputparCov.push_back(iParticle->input_matrix()(j,k));
 	    }
 	  }
	  


// 	  KFTau_Daughter_par.at(ntaus).push_back(iKFTau_Daughter_par);
// 	  KFTau_Daughter_parCov.at(ntaus).push_back(iKFTau_Daughter_parCov);
// 	  KFTau_Daughter_inputpar.at(ntaus).push_back(iKFTau_Daughter_inputpar);
// 	  KFTau_Daughter_inputparCov.at(ntaus).push_back(iKFTau_Daughter_inputparCov);


	}

      }
  
      // std::cout<<"Size of QC disc ================= "<<KFTau_discriminatorByQC.at(0).size()<<std::end;

	KFTau_pions.push_back(std::vector<std::vector<float> >());
	KFTau_Initial_pions.push_back(std::vector<std::vector<float> > ());
	unsigned int TauAmbiguityForPions;
	if(KFTau.Pions(0).size()!=0){TauAmbiguityForPions =0;}
	else if(KFTau.Pions(1).size()!=0){TauAmbiguityForPions=1;}
	else if(KFTau.Pions(2).size()!=0){TauAmbiguityForPions=2;}
	for(unsigned int ipion = 0; ipion < KFTau.Pions(TauAmbiguityForPions).size(); ipion++){
	  
	  std::vector<float > iiKFTau_pions;
	  std::vector<float > iiKFTau_Initial_pions;
	  
	  iiKFTau_pions.push_back(KFTau.Pions(TauAmbiguityForPions).at(ipion).E());
	  iiKFTau_pions.push_back(KFTau.Pions(TauAmbiguityForPions).at(ipion).Px());
	  iiKFTau_pions.push_back(KFTau.Pions(TauAmbiguityForPions).at(ipion).Py());
	  iiKFTau_pions.push_back(KFTau.Pions(TauAmbiguityForPions).at(ipion).Pz());
	  //    std::cout<<"particle parameters: "<< KFTau.Pions(TauAmbiguityForPions).at(ipion).E()<<std::endl;
	  
	  iiKFTau_Initial_pions.push_back(KFTau.InitialPions().at(ipion).E());
	  iiKFTau_Initial_pions.push_back(KFTau.InitialPions().at(ipion).Px());
	  iiKFTau_Initial_pions.push_back(KFTau.InitialPions().at(ipion).Py());
	  iiKFTau_Initial_pions.push_back(KFTau.InitialPions().at(ipion).Pz());
	  
	  KFTau_pions.at(ntaus).push_back(iiKFTau_pions);
	  KFTau_Initial_pions.at(ntaus).push_back(iiKFTau_Initial_pions);

	  //	  std::cout<<"npions  "<<iiKFTau_pions.size()<<"  "<<KFTau_Initial_pions.at(ntaus).size()<<" ipion" << ipion<< "  "<< ntaus <<std::endl;
	}
	//	std::cout<<"ntau "<< ntaus  <<" " <<selected->size()<<" KFTau.Pions(1).size() " <<KFTau.Pions(TauAmbiguityForPions).size() <<std::endl;
	//	std::cout<<" 0--->"<<KFTau_pions.at(ntaus).at(TauAmbiguityForPions).at(3) <<std::endl;
      KFTau_TauVis_p4.push_back(iKFTau_TauVis_p4);
      KFTau_TauFit_p4.push_back(iKFTau_TauFit_p4);
      KFTau_Neutrino_p4.push_back(iKFTau_Neutrino_p4);
      KFTau_TauFitInitial_p4.push_back(iKFTau_TauFitInitial_p4);
      KFTau_NeutrinoInitial_p4.push_back(iKFTau_NeutrinoInitial_p4);

      KFTau_Fit_SecondaryVertex.push_back(iKFTau_Fit_SecondaryVertex);
      KFTau_a1Initial_p4.push_back(iKFTau_a1Initial_p4);
      KFTau_Fit_PrimaryVertex.push_back(iKFTau_Fit_PrimaryVertex);
      KFTau_Fit_InitialPrimaryVertex.push_back(iKFTau_Fit_InitialPrimaryVertex);
      KFTau_Fit_InitialPrimaryVertexReFit.push_back(iKFTau_Fit_InitialPrimaryVertexReFit);
      KFTau_Fit_InitialSecondaryVertex.push_back(iKFTau_Fit_InitialSecondaryVertex);

      KFTau_Fit_charge.push_back(taucharge);

      // std::cout<<" discr  "<<KFTau_discriminatorByQC.size() <<"check ambig size  " << iKFTau_TauVis_p4.at(1).size() <<" taucharge  "<< "taucharge<<taucharge size  " <<KFTau_Fit_charge.size() <<std::endl;
      //-------------------------------------
  
      
      //Match to tau collection to find discriminants
      unsigned int index = 0;
      //  bool discriminatorByKFit(false),discriminatorByQC(false);

      TLorentzVector FitTau(iKFTau_TauFit_p4.at(0).at(1),iKFTau_TauFit_p4.at(0).at(2),iKFTau_TauFit_p4.at(0).at(3),iKFTau_TauFit_p4.at(0).at(0));

      double dP=0.01;
      //     for(reco::PFTauCollection::const_iterator tau = tauCollection->begin(); tau != tauCollection->end(); ++tau, index++) {
      //       reco::PFTauRef tauRef(tauCollection, index);
      //       TLorentzVector CollTau(tauRef->alternatLorentzVect().Px(),tauRef->alternatLorentzVect().Py(),tauRef->alternatLorentzVect().Pz(),tauRef->alternatLorentzVect().E());
      //       double deltaP=sqrt(pow(CollTau.Px()-FitTau.Px(),2.0)+pow(CollTau.Py()-FitTau.Py(),2.0)+pow(CollTau.Pz()-FitTau.Pz(),2.0));
      //       if(deltaP<dP){
      // 	dP=deltaP;
      // 	std::vector<bool> discriminatorPair = CheckTauDiscriminators(tauDiscriminators,tauRef);
      // 	discriminatorByKFit=discriminatorPair.at(0);
      // 	discriminatorByQC=discriminatorPair.at(1);
      //       }
      //     }
      /*std::cout << "drmatch " << dP << " " << (int)discriminatorByKFit << " " << (int)discriminatorByQC 
	<< "iKFTau_TauFit_p4 E: " << iKFTau_TauFit_p4.at(0) 
	<< "iKFTau_Neutrino_p4 E: " << iKFTau_Neutrino_p4.at(0) 
	<< "iKFTau_TauVis_p4 E: " << iKFTau_TauVis_p4.at(0) 
	<< "tauRef->alternatLorentzVect().E()" << E << std::endl;*/
      //     if(dP<0.001){
      //       KFTau_discriminatorByKFit.push_back(discriminatorByKFit);
      //       KFTau_discriminatorByQC.push_back(discriminatorByQC);
      //     }
      //     else{
      //       KFTau_discriminatorByKFit.push_back(false);
      //       KFTau_discriminatorByQC.push_back(false);
      //     }

      unsigned int idx =0;
      unsigned int convergedTau;
      if(iKFTau_TauVis_p4.at(0).at(0)!=0)convergedTau =0;
      if(iKFTau_TauVis_p4.at(1).at(0)!=0)convergedTau =1;
      if(iKFTau_TauVis_p4.at(2).at(0)!=0)convergedTau =2;

      reco::PFTauRef MatchedHPSTau = getMatchedHPSTau(HPStaus,iKFTau_TauVis_p4.at(convergedTau),idx);
      KFTau_MatchedHPS_idx.push_back(idx);
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
      for (unsigned i = 0;  i <  PFJet->numberOfDaughters (); i++) {
	const reco::PFCandidatePtr pfcand = PFJet->getPFConstituent(i);
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
      reco::PFTauRef MatchedHPSTau = getHPSTauMatchedToJet(HPStaus,iPFJet_p4,idx);
      PFJet_MatchedHPS_idx.push_back(idx);
      
    }
  }
  else{
    edm::Handle<pat::JetCollection> jets;
    edm::InputTag labelJets(srcPatJets_);
    iEvent.getByLabel(labelJets, jets);   
    
    for(pat::JetCollection::size_type iPatJet = 0; iPatJet < jets->size(); iPatJet++) {
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
     Electron_HasMatchedConversions.push_back(ConversionTools::hasMatchedConversion(ElectronCollection->at(iPFElectron), hConversions, beamSpot.position()));
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
       
     //reco::TrackRef refTrack=static_cast<reco::TrackRef>(refGsfTrack);//RefElectron->trackRef();
     //int match;
     //getTrackMatch(trackCollection,refTrack,match);
     Electron_Track_idx.push_back(-1);
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
       if(doKinFitTaus_){
	 edm::Handle<std::vector<reco::PFTau> > PFTauCollection;
	 iEvent.getByLabel(hpsTauProducer_, PFTauCollection);
// 	 edm::Handle<reco::PFJetCollection> tauCollection;
// 	 iEvent.getByLabel(pfjetsTag_, tauCollection);
	 TriggerMatch(triggerEvent,index,PFTauCollection,TriggerTauMatchingdr_,match);
	 TauTriggerMatch.push_back(match);
	 match.clear();
       }
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
//       std::cout<<"Size of QC disc  1 ================= "<<KFTau_discriminatorByQC.at(0).size()<<std::end;
//       std::cout<<"Size of Fit disc 2 ================= "<<KFTau_discriminatorByKFit.at(0).size()<<std::end;
   //================  KinFitTaus block ==========
   output_tree->Branch("KFTau_discriminatorByKFit",&KFTau_discriminatorByKFit);
   output_tree->Branch("KFTau_discriminatorByQC",&KFTau_discriminatorByQC);
   output_tree->Branch("KFTau_nKinTaus",&KFTau_nKinTaus);	 

   output_tree->Branch("KFTau_TauFitInitial_p4",&KFTau_TauFitInitial_p4);
   output_tree->Branch("KFTau_NeutrinoInitial_p4",&KFTau_NeutrinoInitial_p4);
   output_tree->Branch("KFTau_a1Initial_p4",&KFTau_a1Initial_p4);
   output_tree->Branch("KFTau_pions",&KFTau_pions);
   output_tree->Branch("KFTau_Initial_pions",&KFTau_Initial_pions);


   output_tree->Branch("KFTau_Fit_PrimaryVertex",&KFTau_Fit_PrimaryVertex);
   output_tree->Branch("KFTau_Fit_InitialPrimaryVertex",&KFTau_Fit_InitialPrimaryVertex);
   output_tree->Branch("KFTau_Fit_InitialPrimaryVertexReFit",&KFTau_Fit_InitialPrimaryVertexReFit);

   output_tree->Branch("KFTau_Fit_SecondaryVertex",&KFTau_Fit_SecondaryVertex);
   output_tree->Branch("KFTau_Fit_InitialSecondaryVertex",&KFTau_Fit_InitialSecondaryVertex);
   output_tree->Branch("KFTau_Fit_Chi2Prob",&KFTau_Fit_Chi2Prob);
   	    
   output_tree->Branch("KFTau_TauVis_p4",&KFTau_TauVis_p4);    
   output_tree->Branch("KFTau_TauFit_p4",&KFTau_TauFit_p4);    
   output_tree->Branch("KFTau_Neutrino_p4",&KFTau_Neutrino_p4);  
   output_tree->Branch("KFTau_MatchedHPS_idx",&KFTau_MatchedHPS_idx);
   output_tree->Branch("KFTau_Track_idx",&PFTau_Track_idx);
   output_tree->Branch("KFTau_indexOfFitInfo",&KFTau_indexOfFitInfo);

   output_tree->Branch("KFTau_Fit_TauPrimVtx",&KFTau_Fit_TauPrimVtx);
   output_tree->Branch("KFTau_Fit_chi2",&KFTau_Fit_chi2);
   output_tree->Branch("KFTau_Fit_ndf",&KFTau_Fit_ndf);
   //  output_tree->Branch("KFTau_Fit_ambiguity",&KFTau_Fit_ambiguity);
   output_tree->Branch("KFTau_Fit_charge",&KFTau_Fit_charge);
   output_tree->Branch("KFTau_Fit_csum",&KFTau_Fit_csum);
   output_tree->Branch("KFTau_Fit_iterations",&KFTau_Fit_iterations);
   output_tree->Branch("KFTau_Fit_BDTVal",&KFTau_Fit_BDTVal);
   output_tree->Branch("KFTau_Fit_TauEnergyFraction",&KFTau_Fit_TauEnergyFraction);
   output_tree->Branch("KFTau_Fit_RefitVisibleMass",&KFTau_Fit_RefitVisibleMass);
   output_tree->Branch("KFTau_Fit_Chi2Prob",&KFTau_Fit_Chi2Prob);
   output_tree->Branch("KFTau_Fit_PV_PV_significance",&KFTau_Fit_PV_PV_significance);
   output_tree->Branch("KFTau_Fit_SV_PV_significance",&KFTau_Fit_SV_PV_significance);

   output_tree->Branch("KFTau_Daughter_pdgid",&KFTau_Daughter_pdgid);
   output_tree->Branch("KFTau_Daughter_charge",&KFTau_Daughter_charge);
   output_tree->Branch("KFTau_Daughter_ambiguity",&KFTau_Daughter_ambiguity);

   output_tree->Branch("KFTau_Daughter_par",&KFTau_Daughter_par);
   output_tree->Branch("KFTau_Daughter_parCov",&KFTau_Daughter_parCov);
   output_tree->Branch("KFTau_Daughter_inputpar",&KFTau_Daughter_inputpar);
   output_tree->Branch("KFTau_Daughter_inputparCov",&KFTau_Daughter_inputparCov);

   output_tree->Branch("ReducedVtx_chi2",&ReducedVtx_chi2);
   output_tree->Branch("ReducedVtx_nTrk",&ReducedVtx_nTrk);
   output_tree->Branch("ReducedVtx_ndof",&ReducedVtx_ndof);
   output_tree->Branch("ReducedVtx_x",&ReducedVtx_x);
   output_tree->Branch("ReducedVtx_y",&ReducedVtx_y);
   output_tree->Branch("ReducedVtx_z",&ReducedVtx_z);
   output_tree->Branch("ReducedVtx_Cov",&ReducedVtx_Cov);
   output_tree->Branch("ReducedVtx_Track_idx",&ReducedVtx_Track_idx);
   output_tree->Branch("ReducedVtx_isFake",&ReducedVtx_isFake);


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

  //======= KinFitTaus ===
  KFTau_discriminatorByKFit.clear();
  KFTau_discriminatorByQC.clear();
  KFTau_TauVis_p4.clear();
  KFTau_TauFit_p4.clear();
  KFTau_TauFitInitial_p4.clear();
  KFTau_Neutrino_p4.clear();
  KFTau_NeutrinoInitial_p4.clear();
  KFTau_MatchedHPS_idx.clear();
  KFTau_Track_idx.clear();
  KFTau_indexOfFitInfo.clear();
  KFTau_a1Initial_p4.clear();


  KFTau_pions.clear();
  KFTau_Initial_pions.clear();

  KFTau_Fit_PrimaryVertex.clear();
  KFTau_Fit_InitialPrimaryVertex.clear();
  KFTau_Fit_InitialPrimaryVertexReFit.clear();
  KFTau_Fit_SecondaryVertex.clear();
  KFTau_Fit_InitialSecondaryVertex.clear();
  KFTau_Fit_BDTVal.clear();

  KFTau_Fit_chi2.clear();	    
  KFTau_Fit_ndf.clear();	    
  // KFTau_Fit_ambiguity.clear();
  KFTau_Fit_charge.clear();   
  KFTau_Fit_csum.clear();     
  KFTau_Fit_iterations.clear();
  KFTau_Fit_TauPrimVtx.clear();

  KFTau_Fit_TauEnergyFraction.clear();
  KFTau_Fit_RefitVisibleMass.clear();
  KFTau_Fit_Chi2Prob.clear();
  KFTau_Fit_PV_PV_significance.clear();
  KFTau_Fit_SV_PV_significance.clear();
  
  KFTau_Daughter_pdgid.clear();
  KFTau_Daughter_charge.clear();
  KFTau_Daughter_ambiguity.clear();

  KFTau_Daughter_par.clear();
  KFTau_Daughter_parCov.clear();
  KFTau_Daughter_inputpar.clear();
  KFTau_Daughter_inputparCov.clear();

  ReducedVtx_chi2.clear();
  ReducedVtx_nTrk.clear();
  ReducedVtx_ndof.clear();
  ReducedVtx_y.clear();
  ReducedVtx_x.clear();
  ReducedVtx_z.clear();
  ReducedVtx_Cov.clear();
  ReducedVtx_Track_idx.clear();
  ReducedVtx_isFake.clear();



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
  Electron_trackMomentumAtVtx  .clear();           
  Electron_numberOfMissedHits  .clear();           
  Electron_HasMatchedConversions.clear();          



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








