// Package:    EventCounter
// Class:      EventCounter
//
// Original Author:  Ian Nugent
//         Created:  Thu Dec  3 11:38:49 CET 2011
// $Id: EventCounter.cc,v 1.1 2011/12/03 15:27:51 inugent Exp $
//
#include "TauDataFormat/TauNtuple/interface/EventCounter.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include <DataFormats/Candidate/interface/Candidate.h>
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>
#include "TauDataFormat/TauNtuple/interface/TauDecay_CMSSW.h"
#include "TauDataFormat/TauNtuple/interface/PdtPdgMini.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
EventCounter::EventCounter(const edm::ParameterSet& iConfig):
  CounterType_(iConfig.getUntrackedParameter<std::string>("CounterType","")),
  gensrc_(iConfig.getParameter<edm::InputTag>( "gensrc" )),
  GenEventInfo_(iConfig.getParameter<edm::InputTag>("GenEventInfo")),
  DataMC_Type_(iConfig.getUntrackedParameter<std::string>("DataMCType",""))
{

  DataMCType DMT;
  DataMC_Type_idx=DMT.GetType(DataMC_Type_);

  //now do what ever initialization is needed
  nevents_weighted.clear();
  nevents.clear();
  DataMCMap.clear();
 }


EventCounter::~EventCounter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
EventCounter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   //////////////////////////
   //
   // Determine the mc/data type 
   //
   float w=1;
   int type=-999;
   if(iEvent.isRealData()){
     type=DataMCType::Data;
   }
   else{
     type=DataMC_Type_idx;
     edm::Handle<GenEventInfoProduct> GenEventInfoProduct;
     iEvent.getByLabel("generator",GenEventInfoProduct);
     w*=GenEventInfoProduct->weight();
     ///////////////////////////////////
     //
     // Check for Signal MC if there are signal taus create a new proccess it based on the 
     //
     edm::Handle<reco::GenParticleCollection> genParticles;
     iEvent.getByLabel(gensrc_, genParticles);
     for(reco::GenParticleCollection::const_iterator itr = genParticles->begin(); itr!= genParticles->end(); ++itr){
       bool found=false;
       unsigned int pdgid=abs(itr->pdgId());
       if(pdgid==PdtPdgMini::Z0 || pdgid==PdtPdgMini::W_plus || pdgid==PdtPdgMini::Higgs0 || pdgid==PdtPdgMini::Higgs_plus){
	 // flag to only select particles that has a daughter tau                                                                                                                                                                           
	 for(unsigned int i = 0; i <itr->numberOfDaughters(); i++){
	   const reco::Candidate *dau=itr->daughter(i);
	   if(abs(dau->pdgId())==PdtPdgMini::tau_minus){
	     unsigned int JAK_ID,TauBitMask;
	     TauDecay_CMSSW TauDecay;
	     TauDecay.AnalyzeTau(static_cast<const reco::GenParticle*>(dau),JAK_ID,TauBitMask);
	     unsigned int mask=TauBitMask%127;
	     if(!found){
	       found=true;
	       type*=100000000;
	       type+=JAK_ID*1000000;
	       type+=mask*10000;
	     }
	     else{
	       type+=JAK_ID*100;
	       type+=mask;
	     }
	   }
	 }
       }
       if(found)break;
     }
   }
   ///////////////////////////////////////
   //
   // Store information
   //
   int idx=-1;
   if(DataMCMap.count(type)>0){
     idx=DataMCMap[type];
   }
   else{
     idx=nevents.size();
     DataMCMap.insert(std::pair<unsigned int,unsigned int>(type,(unsigned int)idx));
     nevents.push_back(0);
     nevents_weighted.push_back(0);
   }
   if(idx<(int)nevents_weighted.size() && idx<(int)nevents.size() && idx>=0){
     nevents_weighted.at(idx)+=w;
     nevents.at(idx)+=1;
   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
EventCounter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EventCounter::endJob() 
{
  for( std::map<unsigned int,unsigned int>::iterator it=DataMCMap.begin() ; it != DataMCMap.end(); it++ ){
    unsigned int idx=(*it).second;
    TString type="";type+=(*it).first;
    if(idx<nevents_weighted.size() && idx<nevents.size()){
      std::cout <<"[EventCounter-" << CounterType_ << "]: " << type <<  " no of weighted events: " <<  nevents_weighted.at(idx) << " no. of events: " <<  nevents.at(idx) << std::endl;
    }
  }
}

// ------------ method called when starting to processes a run  ------------
void 
EventCounter::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
EventCounter::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
EventCounter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
EventCounter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
EventCounter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


