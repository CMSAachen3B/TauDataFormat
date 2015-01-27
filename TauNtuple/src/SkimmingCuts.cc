#include "TauDataFormat/TauNtuple/interface/SkimmingCuts.h"

SkimmingCuts::SkimmingCuts(const edm::ParameterSet& iConfig) :
		preselection_(iConfig.getUntrackedParameter<std::string>("preselection")) {
}

SkimmingCuts::~SkimmingCuts() {
}

bool SkimmingCuts::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
	cnt_++;
	bool pass = false;
	iEvent_ = &iEvent;
	bool AcceptMuon = MuonCuts(iEvent, iSetup);
	bool AcceptElectron = ElectronCuts(iEvent, iSetup);

	if (preselection_ == "Mu") {
		if (AcceptMuon) {
			pass = true;
			cntFound_++;
		}
	} else if (preselection_ == "Ele") {
		if (AcceptElectron) {
			pass = true;
			cntFound_++;
		}
	} else if (preselection_ == "MuOrElePre") {
		if (AcceptMuon || AcceptElectron) {
			pass = true;
			cntFound_++;
		}
	} else if (preselection_ == "DoubleMu") {
		if (DoubleMu(iEvent, iSetup)) {
			pass = true;
			cntFound_++;
		}
	} else if (preselection_ == "DoubleEle") {
		if (DoubleEle(iEvent, iSetup)) {
			pass = true;
			cntFound_++;
		}
	} else if (preselection_ == "EMuTvariable") {
		if ((AcceptMuon && AcceptElectron)
				|| DoubleMu(iEvent, iSetup)
				|| DoubleEle(iEvent, iSetup)
				) {
			pass = true;
			cntFound_++;
		}
	} else if (preselection_ == "MuTau") {
		if(AcceptMuon && PFTauCuts(iEvent, iSetup)) {
			pass = true;
			cntFound_++;
		}
	} else {
		std::cout << "No preselection given. Will use the following instead: one muon + one tau" << std::endl;
		if (AcceptMuon && PFTauCuts(iEvent, iSetup)) {
			pass = true;
			cntFound_++;
		}
	}

	return pass;
}

bool SkimmingCuts::MuonCuts(edm::Event& iEvent, const edm::EventSetup& iSetup) {
	edm::Handle<reco::MuonCollection> muonCollection;
	iEvent_->getByLabel(TauNtuple::muonsTag_, muonCollection);

	for (unsigned int iMuon = 0; iMuon < muonCollection->size(); iMuon++) {
		reco::MuonRef RefMuon(muonCollection, iMuon);
		if (RefMuon.isNonnull()) {
			if (RefMuon->p4().pt() > 8. && fabs(RefMuon->p4().eta()) < 2.5) {
				return true;
			}
		}
	}
	return false;
}

bool SkimmingCuts::ElectronCuts(edm::Event& iEvent, const edm::EventSetup& iSetup) {
	edm::Handle<reco::GsfElectronCollection> electronCollection;
	iEvent_->getByLabel(TauNtuple::PFElectronTag_, electronCollection);

	for (unsigned int iElectron = 0; iElectron < electronCollection->size();iElectron++) {
		reco::GsfElectronRef RefElectron(electronCollection, iElectron);
		if (RefElectron.isNonnull()) {
			reco::SuperClusterRef refSuperCluster = RefElectron->superCluster();
			if (RefElectron->p4().Et() > 8.
					&& fabs(refSuperCluster->eta()) < 2.5
					&& RefElectron->gsfTrack()->trackerExpectedHitsInner().numberOfHits()<= 1
					) {
				return true;
			}
		}
	}
	return false;
}

bool SkimmingCuts::PFTauCuts(edm::Event& iEvent, const edm::EventSetup& iSetup) {
	edm::Handle<std::vector<reco::PFTau> > PFTaus;
	iEvent.getByLabel(TauNtuple::hpsTauProducer_, PFTaus);

	edm::Handle<reco::PFTauDiscriminator> HPSByDecayModeFinding;
	iEvent.getByLabel(edm::InputTag("hpsPFTauDiscriminationByDecayModeFinding", "","TauNtupleProcess"), HPSByDecayModeFinding);

	for (unsigned int iPFTau = 0; iPFTau < PFTaus->size(); iPFTau++) {
		reco::PFTauRef PFTauCand(PFTaus, iPFTau);
		if (PFTauCand.isNonnull()) {
			if (PFTauCand->p4().pt() > 18.
					&& fabs(PFTauCand->p4().eta()) < 2.4
					&& (*HPSByDecayModeFinding)[PFTauCand]
				) {
				return true;
			}
		}
	}
	return false;
}

bool SkimmingCuts::DoubleMu(edm::Event& iEvent, const edm::EventSetup& iSetup) {
	unsigned int mus(0);
	edm::Handle<reco::MuonCollection> muonCollection;
	iEvent_->getByLabel(TauNtuple::muonsTag_, muonCollection);

	for (unsigned int iMuon = 0; iMuon < muonCollection->size(); iMuon++) {
		reco::MuonRef RefMuon(muonCollection, iMuon);
		if (RefMuon.isNonnull()) {
			if (RefMuon->p4().pt() > 8. && fabs(RefMuon->p4().eta()) < 2.5) {
				mus++;
			}
		}
		if (mus > 1)
			return true;
	}
	return false;
}

bool SkimmingCuts::DoubleEle(edm::Event& iEvent, const edm::EventSetup& iSetup) {
	unsigned int es(0);
	edm::Handle<reco::GsfElectronCollection> electronCollection;
	iEvent_->getByLabel(TauNtuple::PFElectronTag_, electronCollection);

	for (unsigned int iElectron = 0; iElectron < electronCollection->size();iElectron++) {
		reco::GsfElectronRef RefElectron(electronCollection, iElectron);
		if (RefElectron.isNonnull()) {
			reco::SuperClusterRef refSuperCluster = RefElectron->superCluster();
			if (RefElectron->p4().pt() > 8.
					&& fabs(refSuperCluster->eta()) < 2.5
					&& RefElectron->gsfTrack()->trackerExpectedHitsInner().numberOfHits()<= 1
				) {
				es++;
			}
		}
		if (es > 1)
			return true;
	}
	return false;
}

void SkimmingCuts::beginJob() {
	cnt_ = 0;
	cntFound_ = 0;
}

void SkimmingCuts::endJob() {
	float ratio = 0.0;
	if (cnt_ != 0)
		ratio = (float) cntFound_ / cnt_;
	std::cout << "[SkimmingCuts]-->  " << "NEvents: " << cnt_
			<< "   NEventsPass: " << cntFound_ << "   Efficiency: "
			<< ratio * 100.0 << "%" << std::endl;
}

void SkimmingCuts::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(SkimmingCuts);
