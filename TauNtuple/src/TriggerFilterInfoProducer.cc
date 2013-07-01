#include "TauDataFormat/TauNtuple/interface/TriggerFilterInfoProducer.h"

TriggerFilterInfoProducer::TriggerFilterInfoProducer( const edm::ParameterSet& pset ) :
  MyTriggerHelper()
{
  produces<std::vector<std::string> >("TriggerFilterInfoList").setBranchAlias("TriggerFilterInfoList");
}

TriggerFilterInfoProducer::~TriggerFilterInfoProducer(){
}

void TriggerFilterInfoProducer::beginJob(){
  return;
}

void TriggerFilterInfoProducer::produce( edm::Event& e, const edm::EventSetup& iSetup){
  std::auto_ptr<std::vector<std::string> > TriggerFilterInfoList(new std::vector<std::string>());
  std::vector<std::string> TriggerList=MyTriggerHelper.GetTriggerList();
  for(unsigned int i=0; i<TriggerList.size();i++){
    TriggerFilterInfoList->push_back(TriggerList.at(i));
  }
  e.put(TriggerFilterInfoList,"TriggerFilterInfoList");  
  return ;
}  

void TriggerFilterInfoProducer::endRun( const edm::Run& r, const edm::EventSetup& iSetup){
  return;
}


void TriggerFilterInfoProducer::endJob(){
  return;
}

