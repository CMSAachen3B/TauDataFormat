#include "TauDataFormat/TauNtuple/interface/MultiTriggerFilter.h"
#include "TString.h"

MultiTriggerFilter::MultiTriggerFilter(const edm::ParameterSet& iConfig) :
		TriggerResults_(iConfig.getParameter<edm::InputTag>("TriggerResults")), useTriggers_(iConfig.getParameter<std::vector<std::string> >("useTriggers")) {

}

MultiTriggerFilter::~MultiTriggerFilter() {

}

//
// member functions
//

// ------------ method called on each new Event  ------------
bool MultiTriggerFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
	TauNtuple::MyTriggerInfoNames.clear();
	edm::Handle<edm::TriggerResults> triggerResults;
	iEvent.getByLabel(TriggerResults_, triggerResults);
	const edm::TriggerNames& trigNames = iEvent.triggerNames(*triggerResults);
	bool passed = false;
	if (triggerResults.isValid()) {
		int ntrigs = triggerResults->size();
		for (int itrig = 0; itrig < ntrigs; ++itrig) {
			TString trigname = trigNames.triggerName(itrig);
			for (unsigned int i = 0; i < useTriggers_.size(); i++) {
				if (trigname.Contains(useTriggers_.at(i))) {
					TauNtuple::MyTriggerInfoNames.push_back(trigname.Data());
					//std::cout<<trigname<<"  "<<triggerResults->accept(itrig)<<std::endl;
					if (triggerResults->accept(itrig))
						passed = true;
				}
			}
		}
	}
	return passed;
}

// ------------ method called once each job just before starting event loop  ------------
void MultiTriggerFilter::beginJob() {

}

// ------------ method called once each job just after ending the event loop  ------------
void MultiTriggerFilter::endJob() {

}

// ------------ method called when starting to processes a run  ------------
bool MultiTriggerFilter::beginRun(edm::Run& iRun, edm::EventSetup const& iSetup) {
	std::string processName = "HLT";
	bool changed(true);
	if (hltConfig_.init(iRun, iSetup, processName, changed)) {
		// if init returns TRUE, initialisation has succeeded!
		printf(" %d \n", changed);
		if (changed) {
			const std::string & DSName = hltConfig_.datasetName(1);
			std::cout << DSName << std::endl;
			//const std::vector< std::vector< std::string > > & AllDSName=hltConfig_.datasetContents();
			const std::vector<std::string> & TriggNames = hltConfig_.triggerNames();

			for (std::vector<std::string>::const_iterator iName = TriggNames.begin(); iName != TriggNames.end(); ++iName) {

				std::cout << (*iName) << std::endl;
			}
			// The HLT config has actually changed wrt the previous Run, hence rebook your
			// histograms or do anything else dependent on the revised HLT config
			std::cout << "HLTSelector:: tableName of process " << processName << " is " << hltConfig_.tableName() << std::endl;
//  printf("succes \n");
		}
	} else {
		// if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
		// with the file and/or code and needs to be investigated!
		// In this case, all access methods will return empty values!
	}

	return true;
}

// ------------ method called when ending the processing of a run  ------------
bool MultiTriggerFilter::endRun(edm::Run&, edm::EventSetup const&) {
	return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool MultiTriggerFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {
	return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool MultiTriggerFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {
	return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MultiTriggerFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MultiTriggerFilter);

