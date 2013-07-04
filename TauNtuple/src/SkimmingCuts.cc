#include "TauDataFormat/TauNtuple/interface/SkimmingCuts.h"

SkimmingCuts::SkimmingCuts(const edm::ParameterSet& iConfig):
  doMuonOnly_( iConfig.getParameter<bool>("doMuonOnly") )
{
}


SkimmingCuts::~SkimmingCuts(){}

bool SkimmingCuts::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){
  cnt_++;
  bool pass = false;
  iEvent_=&iEvent;
  bool AcceptMuon = MuonsCuts(iEvent, iSetup);
  if(doMuonOnly_){
    if(AcceptMuon){
      pass = true;
      cntFound_++;
    }
    return pass;
  }
  bool AcceptPFTau = PFTausCuts(iEvent, iSetup);
  bool AcceptElectron = ElectronCuts(iEvent, iSetup);

    //----------- This Blos is for private production to analyse muon tau decay with KF
  if(AcceptMuon && (AcceptElectron || AcceptPFTau)){
    pass = true;
    cntFound_++;
  }
  return pass;
}

bool SkimmingCuts::MuonsCuts(edm::Event& iEvent, const edm::EventSetup& iSetup){
  edm::Handle< reco::MuonCollection > muonCollection;
  iEvent_->getByLabel(TauNtuple::muonsTag_,  muonCollection);
  for(unsigned int iMuon = 0; iMuon< muonCollection->size(); iMuon++){
    reco::MuonRef RefMuon(muonCollection, iMuon);
    if(TauNtuple::isGoodMuon(RefMuon))return true;
  }
  return false;
}

bool SkimmingCuts::ElectronCuts(edm::Event& iEvent, const edm::EventSetup& iSetup){
  edm::Handle< reco::GsfElectronCollection > electronCollection;
  iEvent_->getByLabel(TauNtuple::PFElectronTag_, electronCollection);
  
  for(unsigned int iElectron=0; iElectron<electronCollection->size(); iElectron++){
    reco::GsfElectronRef RefElectron(electronCollection, iElectron);
    if(TauNtuple::isGoodElectron(RefElectron))return true;
  }
  return false;
}

bool SkimmingCuts::PFTausCuts(edm::Event& iEvent, const edm::EventSetup& iSetup){
  edm::Handle<std::vector<reco::PFTau> > PFTaus;
  iEvent.getByLabel(TauNtuple::hpsTauProducer_, PFTaus);

  edm::Handle<reco::PFTauDiscriminator> HPSByDecayModeFinding;
  iEvent.getByLabel("hpsPFTauDiscriminationByDecayModeFinding", HPSByDecayModeFinding);

  edm::Handle<reco::PFTauDiscriminator> HPSPFTauDiscriminationByMediumIsolationMVA;
  iEvent.getByLabel("hpsPFTauDiscriminationByMediumIsolationMVA", HPSPFTauDiscriminationByMediumIsolationMVA);

  for ( unsigned int iPFTau = 0; iPFTau < PFTaus->size(); iPFTau++ ) {
    reco::PFTauRef PFTauCand(PFTaus, iPFTau);
    if(TauNtuple::isGoodTau(PFTauCand,HPSPFTauDiscriminationByMediumIsolationMVA,HPSByDecayModeFinding))return true;
  }
  return false;
}

void SkimmingCuts::beginJob(){
  cnt_ =0;
  cntFound_ =0;
}

void SkimmingCuts::endJob() {
  float ratio = 0.0;
  if(cnt_!=0) ratio=(float)cntFound_/cnt_;
  std::cout<<"[SkimmingCuts]-->  "<<"NEvents: " << cnt_ <<"   NEventsPass: "<< cntFound_ <<"   Efficiency: "<< ratio*100.0 <<"%"<<std::endl;
}

void SkimmingCuts::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(SkimmingCuts);
