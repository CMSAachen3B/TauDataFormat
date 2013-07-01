#include "TauDataFormat/TauNtuple/interface/TriggerHelper.h"

std::vector<std::string> TriggerHelper::TriggerList;

TriggerHelper::TriggerHelper(){
  Reset();
}

TriggerHelper::~TriggerHelper(){
}


void TriggerHelper::Reset(){
  TriggerList.clear();
}
void TriggerHelper::AddTrigger(std::string name){
  TriggerList.push_back(name);
}

std::vector<std::string> TriggerHelper::GetTriggerList(){
  return  TriggerList;
}

