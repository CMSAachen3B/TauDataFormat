#include "TauDataFormat/TauNtuple/interface/TauNtuple.h"
#include "TauDataFormat/TauNtuple/interface/TauDecay_CMSSW.h"
#include "TauDataFormat/TauNtuple/interface/PdtPdgMini.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"
#include <vector>
#include "TMatrixT.h"




#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include <sys/types.h>
#include <dirent.h>



TauNtuple::TauNtuple(const edm::ParameterSet& iConfig):
  primVtxTag_( iConfig.getParameter<edm::InputTag>( "primVtx" ) ),
  muonsTag_(iConfig.getParameter<edm::InputTag>( "muons" )),
  hpsTauProducer_( iConfig.getParameter<edm::InputTag>( "hpsTauProducer" ) ),
  hpsPFTauDiscriminationByTightIsolation_( iConfig.getParameter<edm::InputTag>( "hpsPFTauDiscriminationByTightIsolation" ) ),
  hpsPFTauDiscriminationByMediumIsolation_( iConfig.getParameter<edm::InputTag>( "hpsPFTauDiscriminationByMediumIsolation" ) ),
  hpsPFTauDiscriminationByLooseIsolation_( iConfig.getParameter<edm::InputTag>( "hpsPFTauDiscriminationByLooseIsolation" ) ),
  pfMETTag_( iConfig.getParameter<edm::InputTag>( "pfMet" ) ),
  kinTausTag_( iConfig.getParameter<edm::InputTag>( "kinematicTaus" ) ),
  KinFitAdvanced_( iConfig.getParameter<edm::InputTag>( "kinematicTausAdvanced" ) ),
  pfjetsTag_( iConfig.getParameter<edm::InputTag>( "pfjets" ) ),
  generalTracks_(iConfig.getParameter<edm::InputTag>( "generalTracks" )),
  gensrc_(iConfig.getParameter<edm::InputTag>( "gensrc" )),
  GenEventInfo_(iConfig.getParameter<edm::InputTag>("GenEventInfo")),
  discriminators_( iConfig.getParameter< std::vector<std::string> >("discriminators") ),
  ScaleFactor_(iConfig.getUntrackedParameter<std::string>("ScaleFactor")),
  PUInputFile_(iConfig.getUntrackedParameter<std::string>("PUInputFile")),
  PUInputHistoMC_(iConfig.getUntrackedParameter<std::string>("PUInputHistoMC")),
  PUInputHistoData_(iConfig.getUntrackedParameter<std::string>("PUInputHistoData")),
  do_MCComplete_(iConfig.getUntrackedParameter("do_MCComplete",(bool)(false))),
  do_MCSummary_(iConfig.getUntrackedParameter("do_MCSummary",(bool)(true)))

{   
  LumiWeights_ = edm::Lumi3DReWeighting(PUInputFile_,PUInputFile_, PUInputHistoMC_, PUInputHistoData_);
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
  iEvent_=&iEvent;
  DataMCType DMT;
  DataMC_Type_idx=DMT.GetType();
  if(iEvent.isRealData()){
    DataMC_Type_idx=DataMCType::Data;
  }
  fillEventInfo(iEvent, iSetup);
  fillMET(iEvent, iSetup);
  //  std::cout<<" filKinTaus 1"<<std::endl;
  // Get trackcollection for matching to objects
  edm::Handle< std::vector<reco::Track>  > trackCollection;
  iEvent_->getByLabel(generalTracks_,  trackCollection);
  fillPrimeVertex(iEvent, iSetup,trackCollection);
  fillMuons(iEvent, iSetup,trackCollection);
  fillPFTaus(iEvent, iSetup,trackCollection);
  //  fillPFJets(iEvent, iSetup,trackCollection);
  fillKinFitTaus(iEvent, iSetup,trackCollection); 
  fillTracks(trackCollection);
  fillMCTruth(iEvent, iSetup);
  output_tree->Fill();
}      

// ------------ method called once each job just before starting event loop  ------------
                          
void TauNtuple::fillMCTruth(edm::Event& iEvent, const edm::EventSetup& iSetup){
  if(!iEvent.isRealData()){
    TauDecay_CMSSW myTauDecay;
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent_->getByLabel(gensrc_, genParticles);
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
	  // flag to only select particles that has a daughter tau
	  bool hastaudaughter=false;
	  // look for daughter tau
	  for(unsigned int i = 0; i <itr->numberOfDaughters(); i++){
	    const reco::Candidate *dau=itr->daughter(i);
	    if(abs(dau->pdgId())==PdtPdgMini::tau_minus){
	      if(!hastaudaughter){
		hastaudaughter=true;
		//Fill information as signal particle
		MCSignalParticle_pdgid.push_back(itr->pdgId());
		MCSignalParticle_charge.push_back(itr->charge());
		MCSignalParticle_Tauidx.push_back(std::vector<unsigned int>());
		std::vector<double> iSig_Poca;
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
	      }
	      unsigned int tauidx=MCTauandProd_p4.size();
	      MCSignalParticle_Tauidx.at(MCSignalParticle_Tauidx.size()-1).push_back(tauidx);
	      // Analysis the tau decay
	      unsigned int JAK_ID,TauBitMask;
	      //MCTauandProd_p4.push_back(std::vector<std::vector<float> >());
	      //MCTauandProd_pdgid.push_back(std::vector<int>());
	      //MCTauandProd_charge.push_back(std::vector<int>());
	      myTauDecay.AnalyzeTau(static_cast<const reco::GenParticle*>(dau),JAK_ID,TauBitMask);
	      std::vector<const reco::GenParticle* > TauDecayProducts=myTauDecay.Get_TauDecayProducts();
	      MCTauandProd_midx.push_back(myTauDecay.Get_MotherIdx());
	      MCTau_JAK.push_back(JAK_ID);
	      MCTau_DecayBitMask.push_back(TauBitMask);
	      MCTauandProd_pdgid.push_back(std::vector<int>());
	      MCTauandProd_charge.push_back(std::vector<int>());
	      MCTauandProd_p4.push_back(std::vector<std::vector<float> >());
	     
	      for(unsigned int i=0;i<TauDecayProducts.size();i++){
		MCTauandProd_pdgid.at(tauidx).push_back(TauDecayProducts.at(i)->pdgId());
		MCTauandProd_charge.at(tauidx).push_back(TauDecayProducts.at(i)->charge());
		std::vector<float > iTauandProd_p4;

		iTauandProd_p4.push_back(TauDecayProducts.at(i)->p4().E());
		iTauandProd_p4.push_back(TauDecayProducts.at(i)->p4().Px());
		iTauandProd_p4.push_back(TauDecayProducts.at(i)->p4().Py());
		iTauandProd_p4.push_back(TauDecayProducts.at(i)->p4().Pz());
	
		MCTauandProd_p4.at(tauidx).push_back(iTauandProd_p4);
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
  iEvent_->getByLabel( primVtxTag_, primVtxs);

  float nVtxs=primVtxs->size();
  int ndim=3;
  nVtxs = primVtxs->size();
  for(int i=0;i<nVtxs;i++){
    const reco::Vertex &pv = (*primVtxs)[i];
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
    std::vector<bool> found;
    for(reco::Vertex::trackRef_iterator iTrack=pv.tracks_begin(); iTrack<pv.tracks_end();iTrack++){
      int match;
      reco::TrackRef refTrack=iTrack->castTo<reco::TrackRef>();
      getTrackMatch(trackCollection,refTrack,match);
      matches.push_back(match);
    }
    Vtx_Track_idx.push_back(matches);
  }
}

void 
TauNtuple::fillMuons(edm::Event& iEvent, const edm::EventSetup& iSetup,edm::Handle< std::vector<reco::Track>  > &trackCollection){
  edm::Handle< reco::MuonCollection > muonCollection;
  iEvent_->getByLabel(muonsTag_,  muonCollection);
  int Muon_index =0;
  for(reco::MuonCollection::const_iterator iMuon = muonCollection->begin(); iMuon!= muonCollection->end(); ++iMuon, Muon_index++){
    reco::MuonRef RefMuon(muonCollection, Muon_index);
    std::vector<double> iMuon_Poca;
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

    Muon_numberOfChambers.push_back(RefMuon->numberOfChambers());
    Muon_isGlobalMuon.push_back(RefMuon->isGlobalMuon());
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
  iEvent_->getByLabel(hpsTauProducer_, HPStaus);
  
  edm::Handle<reco::PFTauDiscriminator> HPSTightIsoDiscr;
  iEvent_->getByLabel(hpsPFTauDiscriminationByTightIsolation_, HPSTightIsoDiscr);
  edm::Handle<reco::PFTauDiscriminator> HPSMediumIsoDiscr;
  iEvent_->getByLabel(hpsPFTauDiscriminationByMediumIsolation_, HPSMediumIsoDiscr);
  edm::Handle<reco::PFTauDiscriminator> HPSLooseIsoDiscr;
  iEvent_->getByLabel(hpsPFTauDiscriminationByLooseIsolation_, HPSLooseIsoDiscr);
  
  for ( unsigned iPFTau = 0; iPFTau < HPStaus->size(); ++iPFTau ) {
    
    reco::PFTauRef HPStauCandidate(HPStaus, iPFTau);
    std::vector<float> iPFTau_p4;
    iPFTau_p4.push_back(HPStauCandidate->p4().E());
    iPFTau_p4.push_back(HPStauCandidate->p4().Px());
    iPFTau_p4.push_back(HPStauCandidate->p4().Py());
    iPFTau_p4.push_back(HPStauCandidate->p4().Pz());

    PFTau_p4.push_back(iPFTau_p4);
    
    PFTau_isTightIsolation.push_back((*HPSTightIsoDiscr)[HPStauCandidate]);
    PFTau_isMediumIsolation.push_back((*HPSMediumIsoDiscr)[HPStauCandidate]);
    PFTau_isLooseIsolation.push_back((*HPSLooseIsoDiscr)[HPStauCandidate]);
    PFTau_hpsDecayMode.push_back(HPStauCandidate->decayMode());
    PFTau_Charge.push_back(HPStauCandidate->charge());
    

    reco::PFCandidateRefVector ChargedHadrCand=HPStauCandidate->signalPFChargedHadrCands();
    std::vector<int> matches;
    std::vector<bool> found;
    for(unsigned int i=0; i<ChargedHadrCand.size();i++){
      reco::PFCandidateRef Cand(ChargedHadrCand,i);
      reco::TrackRef refTrack=Cand.get()->trackRef();
      int match;
      found.push_back(getTrackMatch(trackCollection,refTrack,match));
      matches.push_back(match);
    }
    PFTau_Track_idx.push_back(matches);
  }
}


void 
TauNtuple::fillKinFitTaus(edm::Event& iEvent, const edm::EventSetup& iSetup,edm::Handle< std::vector<reco::Track>  > &trackCollection){
  // printf("fillKinFittaus 2  \n");
  edm::Handle<reco::PFTauCollection> tauCollection;
  iEvent.getByLabel(kinTausTag_, tauCollection);
  
  //========= HPS taus for matching issues
  edm::Handle<std::vector<reco::PFTau> > HPStaus;
  iEvent_->getByLabel(hpsTauProducer_, HPStaus);
  //========= HPS taus for matching issues

  std::vector<edm::Handle<reco::PFTauDiscriminator> > tauDiscriminators;
  for(std::vector<std::string>::const_iterator discr=discriminators_.begin(); discr!=discriminators_.end(); ++discr) {
    edm::Handle<reco::PFTauDiscriminator> tmpHandle;
    iEvent.getByLabel("KinematicTauBasicProducer", *discr, tmpHandle);
    tauDiscriminators.push_back(tmpHandle);
  }
  unsigned int index = 0;
  int NoRefittedTaus =0;
  for(reco::PFTauCollection::const_iterator tau = tauCollection->begin(); tau != tauCollection->end(); ++tau, index++) {
    reco::PFTauRef tauRef(tauCollection, index);
    
    std::vector<bool> discriminatorPair = CheckTauDiscriminators(tauDiscriminators,tauRef);
    KFTau_discriminatorByKFit.push_back(discriminatorPair.at(0));
    KFTau_discriminatorByQC.push_back(discriminatorPair.at(1));
    //---------------------- count number of taus passed KF ----------------
    if(discriminatorPair.at(0)){ 
      NoRefittedTaus++;
      KFTau_indexOfFitInfo.push_back(NoRefittedTaus - 1);
    }else{
      KFTau_indexOfFitInfo.push_back(-1);
    }
    //---------------------- count number of taus passed KF ----------------
    std::vector<float> iKFTau_TauVis_p4;
    std::vector<float> iKFTau_TauFit_p4;
    std::vector<float> iKFTau_Neutrino_p4;
    
    iKFTau_TauVis_p4.push_back(tauRef->p4().E());
    iKFTau_TauVis_p4.push_back(tauRef->p4().Px());
    iKFTau_TauVis_p4.push_back(tauRef->p4().Py());
    iKFTau_TauVis_p4.push_back(tauRef->p4().Pz()); 
    
    iKFTau_TauFit_p4.push_back(tauRef->alternatLorentzVect().E());
    iKFTau_TauFit_p4.push_back(tauRef->alternatLorentzVect().Px());
    iKFTau_TauFit_p4.push_back(tauRef->alternatLorentzVect().Py());
    iKFTau_TauFit_p4.push_back(tauRef->alternatLorentzVect().Pz());
    
    iKFTau_Neutrino_p4.push_back(tauRef->alternatLorentzVect().E() - tauRef->p4().E());
    iKFTau_Neutrino_p4.push_back(tauRef->alternatLorentzVect().Px()- tauRef->p4().Px());
    iKFTau_Neutrino_p4.push_back(tauRef->alternatLorentzVect().Py()- tauRef->p4().Py());
    iKFTau_Neutrino_p4.push_back(tauRef->alternatLorentzVect().Pz()- tauRef->p4().Pz());
    
    KFTau_TauFit_p4.push_back(iKFTau_TauFit_p4);
    KFTau_TauVis_p4.push_back(iKFTau_TauVis_p4);
    KFTau_Neutrino_p4.push_back(iKFTau_Neutrino_p4);
    
    TLorentzVector vecKf;
    double px  = iKFTau_TauFit_p4.at(1);
    double py  = iKFTau_TauFit_p4.at(2);
    double pz  = iKFTau_TauFit_p4.at(3);
    double e   = iKFTau_TauFit_p4.at(0);
    vecKf.SetPxPyPzE(px,py,pz,e);
    
    
    unsigned int idx =0;
    reco::PFTauRef MatchedHPSTau = getMatchedHPSTau(HPStaus,iKFTau_TauVis_p4,idx);
    KFTau_MatchedHPS_idx.push_back(idx);
    
    
  }
  KFTau_nKinTaus = NoRefittedTaus;
  //   std::cout<<"fillKinFittaus  " <<std::endl;
  //    printf("fillKinFittaus 3  \n");

  //================== KinematicFit Info ===================
  edm::Handle<SelectedKinematicDecayCollection> selected;
  iEvent.getByLabel(KinFitAdvanced_, selected);
  for(SelectedKinematicDecayCollection::const_iterator decay = selected->begin(); decay != selected->end(); ++decay){
    //  std::vector< SelectedKinematicParticle const* > Particles;
    const SelectedKinematicParticleCollection& Particles =decay->particles();
    //    Particles = decay->particles();
    //----------------- Store Quality values

    KFTau_Fit_TauEnergyFraction.push_back(decay->energyTFraction());
    KFTau_Fit_RefitVisibleMass.push_back(decay->a1Mass());
    KFTau_Fit_Chi2.push_back(decay->chi2prob());
    KFTau_Fit_PV_PV_significance.push_back(decay->vtxSignPVRotPVRed());
    KFTau_Fit_SV_PV_significance.push_back(decay->vtxSignPVRotSV());
    //    printf("decay->chi2prob()   %f  \n",decay->chi2prob());
    //   std::cout<<"decay chi2  " << decay->chi2prob()<<std::endl;
    //----------------- Store Quality values
    for(std::vector<SelectedKinematicParticle>::const_iterator iParticle = Particles.begin(); iParticle != Particles.end(); ++iParticle){
      //      std::cout<<iParticle->name()<<std::endl;
      if(iParticle->name()=="tau"){
	TVectorT<double> intauParam ;
	intauParam.ResizeTo(7);
	intauParam=iParticle->SelectedKinematicParticle::input_parameters();
	
	unsigned int PrimaryVertexIndex =999;
	float diff=999.;
	for(unsigned int iPrVtx=0; iPrVtx <Vtx_x.size(); iPrVtx++ ){
	  
	  float delta = sqrt(pow(Vtx_x.at(iPrVtx)-intauParam[0],2)
 		             + pow(Vtx_y.at(iPrVtx)-intauParam[1],2)
                             + pow(Vtx_z.at(iPrVtx)-intauParam[2],2));
	  
	  if(delta < diff){
	    diff = delta;
	    PrimaryVertexIndex = iPrVtx;
	  }
	}
	
	std::vector<float> iKFTau_Fit_TauPrimVtx;
	iKFTau_Fit_TauPrimVtx.push_back(iParticle->vertex().X());
	iKFTau_Fit_TauPrimVtx.push_back(iParticle->vertex().Y());
	iKFTau_Fit_TauPrimVtx.push_back(iParticle->vertex().Z());
	KFTau_Fit_TauPrimVtx.push_back(iKFTau_Fit_TauPrimVtx);
	KFTau_Fit_IndexToPrimVertexVector.push_back(PrimaryVertexIndex);
	KFTau_Fit_chi2.push_back(decay->chi2());
	KFTau_Fit_ndf.push_back(decay->ndf());
	KFTau_Fit_ambiguity.push_back(iParticle->ambiguity());
	KFTau_Fit_charge.push_back(iParticle->charge());
	KFTau_Fit_csum.push_back(decay->csum());
	KFTau_Fit_iterations.push_back(decay->iterations());

	//	std::cout<<"tau   charge " << iParticle->charge()<<std::endl;
      }
    }
  }
}

void 
TauNtuple::fillPFJets(edm::Event& iEvent, const edm::EventSetup& iSetup,edm::Handle< std::vector<reco::Track>  > &trackCollection){
  edm::Handle<reco::PFJetCollection> JetCollection;
  iEvent_->getByLabel(pfjetsTag_,  JetCollection);

  for(reco::PFJetCollection::size_type iPFJet = 0; iPFJet < JetCollection->size(); iPFJet++) {
    reco::PFJetRef PFJet(JetCollection, iPFJet);
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
    PFJet_PFJet_neutralEmEnergyFraction.push_back(PFJet->neutralEmEnergyFraction());
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
    std::vector<bool> found;
                                             
    edm::Handle< std::vector<reco::Track>  > trackCollection;
    iEvent_->getByLabel(generalTracks_,  trackCollection);
    const edm::ProductID &TrID = trackCollection.id();
    
    reco::TrackRefVector refTracks;
    refTracks.reserve( PFJet->chargedMultiplicity() );
    for (unsigned i = 0;  i <  PFJet->numberOfDaughters (); i++) {
      const reco::PFCandidatePtr pfcand = PFJet->getPFConstituent(i);
      reco::TrackRef trackref = pfcand->trackRef();
      if( trackref.isNonnull() ) {
	if(trackref.id() != TrID) continue;
	refTracks.push_back( trackref );
      }
    }

    //  reco::TrackRefVector refTracks = PFJet->getTrackRefs();
    //std::vector<reco::TrackRef> refTracks = PFJet->getTrackRefs();   // <-------  fixing bug with ID uncoincedens, use std::vector<reco::TrackRef> instead reco::TrackRefVector refTracks

    getTrackMatch(trackCollection,refTracks,matches,found);
    PFJet_Track_idx.push_back(matches);
    edm::Handle<std::vector<reco::PFTau> > HPStaus;
    iEvent_->getByLabel(hpsTauProducer_, HPStaus);
    unsigned int idx =0; 
    reco::PFTauRef MatchedHPSTau = getHPSTauMatchedToJet(HPStaus,iPFJet_p4,idx);
    PFJet_MatchedHPS_idx.push_back(idx);

  }

}


void 
TauNtuple::fillElectrons(edm::Event& iEvent, const edm::EventSetup& iSetup,edm::Handle< std::vector<reco::Track>  > &trackCollection){}
void 

TauNtuple::fillMET(edm::Event& iEvent, const edm::EventSetup& iSetup){

  edm::Handle< edm::View<reco::PFMET> > pfMEThandle;
  iEvent.getByLabel(pfMETTag_, pfMEThandle);
  MET_et=pfMEThandle->front().et();
  MET_phi=pfMEThandle->front().phi();
  MET_sumET=pfMEThandle->front().sumEt();
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
  output_tree->Branch("Vtx_isFake",&Vtx_isFake);
  
  //=============  Muon Block ====
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
  output_tree->Branch("Muon_Track_idx",&Muon_Track_idx);
  output_tree->Branch("Muon_hitPattern_pixelLayerwithMeas",&Muon_hitPattern_pixelLayerwithMeas);
  output_tree->Branch("Muon_numberOfMatchedStations",&Muon_numberOfMatchedStations);
  output_tree->Branch("Muon_normChi2",&Muon_normChi2);
  output_tree->Branch("Muon_hitPattern_numberOfValidMuonHits",&Muon_hitPattern_numberOfValidMuonHits);
  output_tree->Branch("Muon_innerTrack_numberofValidHits",&Muon_innerTrack_numberofValidHits);
  output_tree->Branch("Muon_numberOfMatches",&Muon_numberOfMatches);
  output_tree->Branch("Muon_Charge",&Muon_Charge);
  output_tree->Branch("Muon_numberOfChambers",&Muon_numberOfChambers);
  
  //================  PFTau block ==========
  output_tree->Branch("PFTau_p4",&PFTau_p4);
  output_tree->Branch("PFTau_isTightIsolation",&PFTau_isTightIsolation);
  output_tree->Branch("PFTau_isMediumIsolation",&PFTau_isMediumIsolation);
  output_tree->Branch("PFTau_isLooseIsolation",&PFTau_isLooseIsolation);
  output_tree->Branch("PFTau_hpsDecayMode",&PFTau_hpsDecayMode);
  output_tree->Branch("PFTau_Charge",&PFTau_Charge);
  output_tree->Branch("PFTau_Track_idx",&PFTau_Track_idx);

  //================  KinFitTaus block ==========
  output_tree->Branch("KFTau_discriminatorByKFit",&KFTau_discriminatorByKFit);
  output_tree->Branch("KFTau_discriminatorByQC",&KFTau_discriminatorByQC);
  output_tree->Branch("KFTau_nKinTaus",&KFTau_nKinTaus);	    	    
  output_tree->Branch("KFTau_TauVis_p4",&KFTau_TauVis_p4);    
  output_tree->Branch("KFTau_TauFit_p4",&KFTau_TauFit_p4);    
  output_tree->Branch("KFTau_Neutrino_p4",&KFTau_Neutrino_p4);  
  output_tree->Branch("KFTau_MatchedHPS_idx",&KFTau_MatchedHPS_idx);
  output_tree->Branch("KFTau_Track_idx",&PFTau_Track_idx);
  output_tree->Branch("KFTau_indexOfFitInfo",&KFTau_indexOfFitInfo);

  output_tree->Branch("KFTau_Fit_TauPrimVtx",&KFTau_Fit_TauPrimVtx);
  output_tree->Branch("KFTau_Fit_IndexToPrimVertexVector",&KFTau_Fit_IndexToPrimVertexVector);
  output_tree->Branch("KFTau_Fit_chi2",&KFTau_Fit_chi2);
  output_tree->Branch("KFTau_Fit_ndf",&KFTau_Fit_ndf);
  output_tree->Branch("KFTau_Fit_ambiguity",&KFTau_Fit_ambiguity);
  output_tree->Branch("KFTau_Fit_charge",&KFTau_Fit_charge);
  output_tree->Branch("KFTau_Fit_csum",&KFTau_Fit_csum);
  output_tree->Branch("KFTau_Fit_iterations",&KFTau_Fit_iterations);

  output_tree->Branch("KFTau_Fit_TauEnergyFraction",&KFTau_Fit_TauEnergyFraction);
  output_tree->Branch("KFTau_Fit_RefitVisibleMass",&KFTau_Fit_RefitVisibleMass);
  output_tree->Branch("KFTau_Fit_Chi2",&KFTau_Fit_Chi2);
  output_tree->Branch("KFTau_Fit_PV_PV_significance",&KFTau_Fit_PV_PV_significance);
  output_tree->Branch("KFTau_Fit_SV_PV_significance",&KFTau_Fit_SV_PV_significance);
  


 //=======  PFJets ===
  output_tree->Branch("PFJet_p4",&PFJet_p4);
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
  output_tree->Branch("PFJet_PFJet_neutralEmEnergyFraction",&PFJet_PFJet_neutralEmEnergyFraction);

  //================  MET block ==========
  output_tree->Branch("MET_et",&MET_et);
  output_tree->Branch("MET_phi",&MET_phi);
  output_tree->Branch("MET_sumET",&MET_sumET);

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
  }
  if(do_MCSummary_){
    output_tree->Branch("MCSignalParticle_p4",&MCSignalParticle_p4);
    output_tree->Branch("MCSignalParticle_pdgid",&MCSignalParticle_pdgid);
    output_tree->Branch("MCSignalParticle_charge",&MCSignalParticle_charge);
    output_tree->Branch("MCSignalParticle_Poca",&MCSignalParticle_Poca);
    output_tree->Branch("MCSignalParticle_Tauidx",&MCSignalParticle_Tauidx);
    output_tree->Branch("MCTauandProd_p4",&MCTauandProd_p4);
    output_tree->Branch("MCTauandProd_pdgid",&MCTauandProd_pdgid);
    output_tree->Branch("MCTauandProd_midx",&MCTauandProd_midx);
    output_tree->Branch("MCTauandProd_charge",&MCTauandProd_charge);
    output_tree->Branch("MCTau_JAK",&MCTau_JAK);   
    output_tree->Branch("MCTau_DecayBitMask",&MCTau_DecayBitMask); 
  }
} 


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// std::vector<bool> TauNtuple::CheckTauDiscriminators(std::vector<edm::Handle<reco::PFTauDiscriminator> > tauDiscriminators, reco::PFTauRef tauRef)
//
// checks that tau candidate pass kinematic fit and KinFit quality criteria 
// returns vector of bool variables, first variable true/false if tau candidate pass/fail kinematic fit
// second variable is true/false if  refitted tau candidate pass/fail quality requirements
std::vector<bool>
TauNtuple::CheckTauDiscriminators(std::vector<edm::Handle<reco::PFTauDiscriminator> > tauDiscriminators, reco::PFTauRef tauRef){
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
reco::PFTauRef 
TauNtuple::getMatchedHPSTau(edm::Handle<std::vector<reco::PFTau> > & HPStaus,   std::vector<float>   &UnmodifiedTau, unsigned int &match){
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
reco::PFTauRef 
TauNtuple::getHPSTauMatchedToJet(edm::Handle<std::vector<reco::PFTau> > & HPStaus,   std::vector<float>   &Jet, unsigned int &match){
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// bool TauNtuple::getTrackMatch(edm::Handle< std::vector<reco::Track>  > &trackCollection, reco::TrackRefVector &refTracks, std::vector<int> &matches,std::vector<bool> &found)
//
// Finds vector of track matches
// calls  getTrackMatch(edm::Handle< std::vector<reco::Track>  > &trackCollection, reco::TrackRef &refTrack, int &match) to match individual tracks
// returns is true if all tracks are found, false otherwise.
bool TauNtuple::getTrackMatch(edm::Handle< std::vector<reco::Track>  > &trackCollection, reco::TrackRefVector &refTracks, std::vector<int> &matches,std::vector<bool> &found){
  int match(-1);
  bool flag=true;
  matches.clear();
  found.clear();
  for(unsigned int i=0; i<refTracks.size();i++){
    reco::TrackRef refTrack(trackCollection,i);
    found.push_back(getTrackMatch(trackCollection,refTrack,match));
    matches.push_back(match);
    if(found.at(i)==false)flag=false;
  }
  return flag;
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
void 
TauNtuple::endJob() {
  std::cout<<" No Of event processed: "<< cnt_ << std::endl;
  output->Write();
  output->Close();
}  

// ------------ method called when starting to processes a run  ------------
void   
TauNtuple::beginRun(edm::Run& Run, edm::EventSetup const& Setup){ 
//------------------- printf out triggeer menu
//   std::string processName = "HLT";
//   bool changed(true);
//   if (hltConfig_.init(Run,Setup,processName,changed)) {
//     if (changed) {
//       const std::string &  DSName = hltConfig_.datasetName(1);
//       std::cout << DSName << std::endl;
//       const std::vector< std::vector< std::string > > & AllDSName=hltConfig_.datasetContents(); 
//       const std::vector< std::string > & TriggNames=hltConfig_.triggerNames();


//       for(std::vector< std::string >::const_iterator iName = TriggNames.begin(); iName !=TriggNames.end(); ++iName ){

// 			std::cout << (*iName) << std::endl;
//       }
//     }
//   } else {
//     std::cout<<"TauNtuple: " << " HLT config extraction failure with process name " << processName_<<std::endl;
//   }

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
  
  Muon_numberOfChambers.clear();
  Muon_Charge.clear();
  Muon_Track_idx.clear();

  Muon_hitPattern_pixelLayerwithMeas.clear();
  Muon_numberOfMatchedStations.clear();
  Muon_normChi2.clear();
  Muon_hitPattern_numberOfValidMuonHits.clear();
  Muon_innerTrack_numberofValidHits.clear();
  Muon_numberOfMatches.clear();

  //======= PFTaus ===
  PFTau_p4.clear();
  PFTau_isTightIsolation.clear();
  PFTau_isMediumIsolation.clear();
  PFTau_isLooseIsolation.clear();
  PFTau_hpsDecayMode.clear();   
  PFTau_Charge.clear();
  PFTau_Track_idx.clear();

  //======= KinFitTaus ===
  KFTau_discriminatorByKFit.clear();
  KFTau_discriminatorByQC.clear();
  KFTau_TauVis_p4.clear();
  KFTau_TauFit_p4.clear();
  KFTau_Neutrino_p4.clear();
  KFTau_MatchedHPS_idx.clear();
  KFTau_Track_idx.clear();
  KFTau_indexOfFitInfo.clear();

  KFTau_Fit_chi2.clear();	    
  KFTau_Fit_ndf.clear();	    
  KFTau_Fit_ambiguity.clear();
  KFTau_Fit_charge.clear();   
  KFTau_Fit_csum.clear();     
  KFTau_Fit_iterations.clear();
  KFTau_Fit_IndexToPrimVertexVector.clear();
  KFTau_Fit_TauPrimVtx.clear();


  KFTau_Fit_TauEnergyFraction.clear();
  KFTau_Fit_RefitVisibleMass.clear();
  KFTau_Fit_Chi2.clear();
  KFTau_Fit_PV_PV_significance.clear();
  KFTau_Fit_SV_PV_significance.clear();
  


  //=======  Electrons ===

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
   PFJet_PFJet_neutralEmEnergyFraction.clear();


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
  MCTauandProd_pdgid.clear();
  MCTauandProd_midx.clear();
  MCTauandProd_charge.clear();
  MCTau_JAK.clear();   
  MCTau_DecayBitMask.clear(); 
  }
}





