#include "TH1F.h"
#include "TFile.h"



void MakeHistosFor_Reweighting(){


  //=============== Combine Data hists

  TFile *f_178098_180252   =  TFile::Open("Cert_178098-180252_7TeV_PromptReco_Collisions11_JSON.pileupTruth_v2_finebin.root");
  TFile *f_177878_179431   =  TFile::Open("Cert_177878-179431_7TeV_PromptReco_Collisions11_JSON.pileupTruth_v2_finebin.root");
  TFile *f_177718_178078   =  TFile::Open("Cert_177718_178078_7TeV_PromptReco_Collisons11_JSON.pileupTruth_v2_finebin.root");
  TFile *f_175832_177515   =  TFile::Open("Cert_175832-177515_PromptReco_JSON.pileupTruth_v2_finebin.root");
  TFile *f_172620_173692   =  TFile::Open("Cert_172620-173692_PromptReco_JSON.pileupTruth_v2_finebin.root");
  TFile *f_170249_172619   =  TFile::Open("Cert_170249-172619_7TeV_ReReco5Aug_Collisions11_JSON_v2.pileupTruth_v2_finebin.root");
  TFile *f_165088_167913   =  TFile::Open("Cert_165088-167913_7TeV_PromptReco_JSON.pileupTruth_v2_finebin.root");
  TFile *f_160404_163869   =  TFile::Open("Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v3.pileupTruth_v2_finebin.root");

  TFile *f_DataPileUpTrue   =  TFile::Open("DataPileUp_e4c192b65fa65af19498a24cd3e2a095_true.root");






  TFile *f_out = new TFile("Lumi_160404_180252_andMC_Flat_Tail.root","RECREATE");

  TH1D *h_178098_180252 =  (TH1D *)f_178098_180252->Get("pileup");
  TH1D *h_177878_179431 =  (TH1D *)f_177878_179431->Get("pileup");
  TH1D *h_177718_178078 =  (TH1D *)f_177718_178078->Get("pileup");
  TH1D *h_175832_177515 =  (TH1D *)f_175832_177515->Get("pileup");
  TH1D *h_172620_173692 =  (TH1D *)f_172620_173692->Get("pileup");
  TH1D *h_170249_172619 =  (TH1D *)f_170249_172619->Get("pileup");
  TH1D *h_165088_167913 =  (TH1D *)f_165088_167913->Get("pileup");
  TH1D *h_160404_163869 =  (TH1D *)f_160404_163869->Get("pileup");
  TH1D *h_DataPileUpTrue = (TH1D *)f_DataPileUpTrue->Get("pileup");


  h_178098_180252->SetName("h_178098_180252");
  h_177878_179431->SetName("h_177878_179431");
  h_177718_178078->SetName("h_177718_178078");
  h_175832_177515->SetName("h_175832_177515");
  h_172620_173692->SetName("h_172620_173692");
  h_170249_172619->SetName("h_170249_172619");
  h_165088_167913->SetName("h_165088_167913");
  h_160404_163869->SetName("h_160404_163869");
  h_DataPileUpTrue->SetName("h_DataPileUpTrue");

  h_178098_180252->SetTitle("h_178098_180252");
  h_177878_179431->SetTitle("h_177878_179431");
  h_177718_178078->SetTitle("h_177718_178078");
  h_175832_177515->SetTitle("h_175832_177515");
  h_172620_173692->SetTitle("h_172620_173692");
  h_170249_172619->SetTitle("h_170249_172619");
  h_165088_167913->SetTitle("h_165088_167913");
  h_160404_163869->SetTitle("h_160404_163869");
  h_DataPileUpTrue->SetTitle("h_DataPileUpTrue");


  printf("h_178098_180252->%d \n",h_178098_180252->GetNbinsX());
  printf("h_177878_179431->%d \n",h_177878_179431->GetNbinsX());
  printf("h_177718_178078->%d \n",h_177718_178078->GetNbinsX());
  printf("h_175832_177515->%d \n",h_175832_177515->GetNbinsX());
  printf("h_172620_173692->%d \n",h_172620_173692->GetNbinsX());
  printf("h_170249_172619->%d \n",h_170249_172619->GetNbinsX());
  printf("h_165088_167913->%d \n",h_165088_167913->GetNbinsX());
  printf("h_160404_163869->%d \n",h_160404_163869->GetNbinsX());
  printf("h_DataPileUpTrue->%d \n",h_DataPileUpTrue->GetNbinsX());


  int nBins = h_178098_180252->GetNbinsX();

  TH1D* h_160404_180252_all = new TH1D("h_160404_180252_all","h_160404_180252_all",1000,0,25);
  for(int iBin = 1; iBin < nBins; iBin++){
    h_160404_180252_all->SetBinContent(iBin,h_178098_180252->GetBinContent(iBin)+
				       h_177878_179431->GetBinContent(iBin)+
				       h_177718_178078->GetBinContent(iBin)+
				       h_175832_177515->GetBinContent(iBin)+
				       h_172620_173692->GetBinContent(iBin)+
				       h_170249_172619->GetBinContent(iBin)+
				       h_165088_167913->GetBinContent(iBin)+
				       h_160404_163869->GetBinContent(iBin) );

  }




  //==================== MC hsit 

  TH1D* MC_FLAT_PLUS_TAIL_PU = new TH1D("MC_FLAT_PLUS_TAIL_PU","MC_FLAT_PLUS_TAIL_PU",25,0,25);
  TH1D* MC_Fall11_PU = new TH1D("MC_Fall11_PU","MC_Fall11_PU",50,0,50);
  Double_t probdistFlat10[25] = {
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0630151648,
    0.0526654164,
    0.0402754482,
    0.0292988928,
    0.0194384503,
    0.0122016783,
    0.007207042,
    0.004003637,
    0.0020278322,
    0.0010739954,
    0.0004595759,
    0.0002229748,
    0.0001028162,
    4.58337152809607E-05
  };
  Double_t probdistFall11[50] = {
    0.003388501,
    0.010357558,
    0.024724258,
    0.042348605,
    0.058279812,
    0.068851751,
    0.072914824,
    0.071579609,
    0.066811668,
    0.060672356,
    0.054528356,
    0.04919354,
    0.044886042,
    0.041341896,
    0.0384679,
    0.035871463,
    0.03341952,
    0.030915649,
    0.028395374,
    0.025798107,
    0.023237445,
    0.020602754,
    0.0180688,
    0.015559693,
    0.013211063,
    0.010964293,
    0.008920993,
    0.007080504,
    0.005499239,
    0.004187022,
    0.003096474,
    0.002237361,
    0.001566428,
    0.001074149,
    0.000721755,
    0.000470838,
    0.00030268,
    0.000184665,
    0.000112883,
    6.74043E-05,
    3.82178E-05,
    2.22847E-05,
    1.20933E-05,
    6.96173E-06,
    3.4689E-06,
    1.96172E-06,
    8.49283E-07,
    5.02393E-07,
    2.15311E-07,
    9.56938E-08

  }




  int i=0;
 for (i=1;i<26;i++) {
    MC_FLAT_PLUS_TAIL_PU->SetBinContent(i,probdistFlat10[i-1]);
  }

 for (i=1;i<50;i++) {
    MC_Fall11_PU->SetBinContent(i,probdistFall11[i-1]);
  }


   MC_FLAT_PLUS_TAIL_PU->Write();
   MC_Fall11_PU->Write();


   h_178098_180252 ->Write();
   h_177878_179431 ->Write();
   h_177718_178078 ->Write();
   h_175832_177515 ->Write();
   h_172620_173692 ->Write();
   h_170249_172619 ->Write();
   h_165088_167913 ->Write();
   h_160404_163869 ->Write();
   h_160404_180252_all ->Write();
   h_DataPileUpTrue ->Write();

   f_out->Write();
   f_out->Close();
  

}


