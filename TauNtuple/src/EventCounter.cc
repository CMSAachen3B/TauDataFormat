// Package:    EventCounter
// Class:      EventCounter
//
// Original Author:  Ian Nugent
//         Created:  Thu Dec  3 11:38:49 CET 2011
// $Id: EventCounter.cc,v 1.4 2013/05/21 14:54:26 inugent Exp $
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
#include "Validation/EventGenerator/interface/PdtPdgMini.h"
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
  nevents.clear();
  DataMCMap.clear();
 }


EventCounter::~EventCounter()
{
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
   unsigned int type=DataMCType::unknown;
   bool isEmbedding = (DataMC_Type_idx >= 34 && DataMC_Type_idx <= 37);
   if(iEvent.isRealData() && !isEmbedding){
     type=DataMCType::Data;
   }
   else if(isEmbedding){
	 type=DataMC_Type_idx;
	 // splitting by tau decay not yet implemented for embedding
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
     TauDecay_CMSSW myTauDecay;
     myTauDecay.CheckForSignal(type,genParticles);
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
     nevents.push_back(TH1D("","",2,-0.5,1.5));
   }
   if(idx<(int)nevents.size() && idx>=0){
     nevents.at(idx).Fill(1.0,w);
     nevents.at(idx).Fill(0.0,1.0);
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
    if(idx<nevents.size()){
      std::cout <<"[EventCounter-" << CounterType_ << "]: " << type 
		<< " no of weighted events: " <<  nevents.at(idx).GetBinContent(2) << " +/- " << nevents.at(idx).GetBinError(2) 
		<< " no. of events: "         <<  nevents.at(idx).GetBinContent(1) << " +/- " << nevents.at(idx).GetBinError(1) 
		<< std::endl;
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

//define this as a plug-in
DEFINE_FWK_MODULE(EventCounter);

