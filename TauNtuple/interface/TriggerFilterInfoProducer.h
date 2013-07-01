#ifndef TriggerFilterInfoProducer_h
#define TriggerFilterInfoProducer_h

#include <iostream>

// essentials !!!
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "TauDataFormat/TauNtuple/interface/TriggerHelper.h"

class  TriggerFilterInfoProducer : public edm::EDProducer
{
  
 public:
  explicit TriggerFilterInfoProducer( const edm::ParameterSet& ) ;
  virtual ~TriggerFilterInfoProducer();
    
  virtual void produce( edm::Event&, const edm::EventSetup&);
  virtual void beginJob() ;
  virtual void endRun( const edm::Run&, const edm::EventSetup& ) ;
  virtual void endJob() ;

 private:
  TriggerHelper MyTriggerHelper;
  
}; 

DEFINE_FWK_MODULE(TriggerFilterInfoProducer);
#endif
