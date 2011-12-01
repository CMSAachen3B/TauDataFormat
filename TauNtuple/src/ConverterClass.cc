#ifndef ConverterClass_h
#define ConverterClass_h


#include "TLorentzVector.h"
#include "iostream"
#include "TObject.h"
#include "TMath.h"
#include <TROOT.h>
#include <vector>
#include "TVector3.h"
#include "Math/SMatrix.h"


class ConverterClass {

public:  
  ConverterClass();
  int MultiplyByTwo(int a);
  void ConvertMomentaToTLorentz(std::vector<TLorentzVector> Momenta, std::vector<std::vector<double> > &inputVector);
  void ConvertVertexToTVector3(std::vector<TVector3 > VtxVector, std::vector<std::vector<double> > inputVector);
  void CreateCovarianceMatrix(std::vector<ROOT::Math::SMatrix<double,3,3,ROOT::Math::MatRepSym<double,3> > > & matrix, std::vector<std::vector<float> > inputVector);
};

#endif

ConverterClass::ConverterClass()
{}

int ConverterClass::MultiplyByTwo(int a){
  return a*2;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void ConverterClass::ConvertMomentaToTLorentz(std::vector<TLorentzVector> Momenta, std::vector<std::vector<double> > &inputVector){
// converts std::vector<std::vector<double> > to std::vector<TLorentzVector> Momenta
void
ConverterClass::ConvertMomentaToTLorentz(std::vector<TLorentzVector> Momenta, std::vector<std::vector<double> > &inputVector){
  for(std::vector<std::vector<double> >::const_iterator it = inputVector.begin(); it!= inputVector.end(); ++it){
    TLorentzVector temp;
    temp.SetE((*it).at(0));
    temp.SetPx((*it).at(1));
    temp.SetPy((*it).at(2));
    temp.SetPz((*it).at(3));
    Momenta.push_back(temp);
    if((*it).size()!=3) std::cout<<"Warning Momenta vector has size not equal to 4!!! " <<std::endl;
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void ConverterClass::ConvertVertexToTVector3(std::vector<TVector3 > VtxVector, std::vector<std::vector<double> > inputVector){
// converts std::vector<std::vector<double> > to std::vector<TVector3 >
void 
ConverterClass::ConvertVertexToTVector3(std::vector<TVector3 > VtxVector, std::vector<std::vector<double> > inputVector){
  for(std::vector<std::vector<double> >::const_iterator it = inputVector.begin(); it!= inputVector.end(); ++it){
    TVector3 temp;
    temp.SetX((*it).at(0));
    temp.SetY((*it).at(1));
    temp.SetZ((*it).at(2));
    VtxVector.push_back(temp);
    if((*it).size()!=3) std::cout<<"Warning Vertex vector has size not equal to 3!!! " <<std::endl;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void ConverterClass::CreateCovarianceMatrix(std::vector<ROOT::Math::SMatrix<double,3,3,ROOT::Math::MatRepSym<double,3> > > & matrix, std::vector<std::vector<float> > inputVector)
// converts std::vector<std::vector<float> > to std::vector<ROOT::Math::SMatrix<double,3,3,ROOT::Math::MatRepSym<double,3> > >
void 
ConverterClass::CreateCovarianceMatrix(std::vector<ROOT::Math::SMatrix<double,3,3,ROOT::Math::MatRepSym<double,3> > > & matrix, std::vector<std::vector<float> > inputVector){
  ROOT::Math::SVector<double,6> svector6;
  for(std::vector<std::vector<float> >::const_iterator it = inputVector.begin(); it!= inputVector.end(); ++it){
    unsigned int index =0;
    for(std::vector<float>::const_iterator itt = (*it).begin(); itt!=(*it).end(); ++itt,index++){svector6(index) = (*itt);}
    ROOT::Math::SMatrix<double,3,3,ROOT::Math::MatRepSym<double,3> > temp(svector6);
    matrix.push_back(temp);
  }
}


