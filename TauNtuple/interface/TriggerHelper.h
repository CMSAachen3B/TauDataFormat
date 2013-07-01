#ifndef TriggerHelper_h
#define TriggerHelper_h

#include <iostream>
#include <vector>
#include <string>

class TriggerHelper{
 public:

  TriggerHelper();
  ~TriggerHelper();

  void Reset();
  void AddTrigger(std::string name);
  std::vector<std::string> GetTriggerList();

 private:
  static std::vector<std::string> TriggerList;

};
#endif
