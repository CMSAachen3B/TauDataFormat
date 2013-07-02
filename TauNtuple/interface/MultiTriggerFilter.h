// -*- C++ -*-
//
// Package:    MultiTriggerFilter
// Class:      MultiTriggerFilter
// 
/*
 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ian M. Nugent
//

#ifndef MultiTriggerFilter_h
#define MultiTriggerFilter_h

// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include <FWCore/Common/interface/TriggerNames.h>
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "TauDataFormat/TauNtuple/interface/TauNtuple.h"

//
// class declaration
//

class MultiTriggerFilter : public edm::EDFilter {
  public:
    explicit MultiTriggerFilter(const edm::ParameterSet&);
    ~MultiTriggerFilter();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginJob() ;
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;
      
    virtual bool beginRun(edm::Run&, edm::EventSetup const&);
    virtual bool endRun(edm::Run&, edm::EventSetup const&);
    virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
    virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

    // ----------member data ---------------------------
    edm::Event * iEvent_;
    std::string processName_; // process name of (HLT) process for which to get HLT configuration
    HLTConfigProvider hltConfig_;/// The instance of the HLTConfigProvider as a data member
    edm::InputTag  TriggerResults_;
    std::vector<std::string> useTriggers_;
};

#endif
