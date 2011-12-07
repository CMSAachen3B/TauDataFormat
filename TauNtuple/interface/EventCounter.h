// -*- C++ -*-
//
// Package:    EventCounter
// Class:      EventCounter
// 
/**\class EventCounter EventCounter.cc SkimmingTools/EventCounter/src/EventCounter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ian Nugent
//         Created:  Thu Dec  3 11:38:49 CET 2011
// $Id: EventCounter.cc,v 1.1 2011/12/03 15:27:51 inugent Exp $
//
//
#ifndef EventsCounter_h
#define EventsCounter_h


// system include files
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <iostream>
#include <vector>
#include <string>
#include <map>
//
// class declaration
//

class EventCounter : public edm::EDAnalyzer {
   public:
      explicit EventCounter(const edm::ParameterSet&);
      ~EventCounter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
  std::string CounterType_;
  edm::InputTag gensrc_;
  edm::InputTag GenEventInfo_;
  std::string DataMC_Type_;
  std::vector<double> nevents_weighted;
  std::vector<double> nevents;
  unsigned int DataMC_Type_idx;
  std::map<unsigned int,unsigned int>    DataMCMap;
};


//define this as a plug-in
DEFINE_FWK_MODULE(EventCounter);

#endif
