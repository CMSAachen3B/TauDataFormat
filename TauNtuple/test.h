//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Nov 18 10:09:54 2011 by ROOT version 5.27/06b
// from TTree t/t
// found on file: output.root
//////////////////////////////////////////////////////////

#ifndef test_h
#define test_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
   const Int_t kMaxMuon_p4 = 14;
   const Int_t kMaxPFTau_p4 = 95;
   const Int_t kMaxKFTau_TauVis_p4 = 95;
   const Int_t kMaxKFTau_TauFit_p4 = 95;
   const Int_t kMaxKFTau_Neutrino_p4 = 95;
   const Int_t kMaxPFJet_p4 = 95;
   const Int_t kMaxTrack_p4 = 989;

class test {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   vector<float>   *Vtx_chi2;
   vector<float>   *Vtx_nTrk;
   vector<float>   *Vtx_ndof;
   vector<float>   *Vtx_chi2;
   vector<float>   *Vtx_x;
   vector<float>   *Vtx_y;
   vector<float>   *Vtx_z;
   vector<vector<float> > *Vtx_Cov;
   vector<vector<int> > *Vtx_Track_idx;
   Int_t           Muon_p4_;
   Double_t        Muon_p4_fCoordinates_fX[kMaxMuon_p4];   //[Muon_p4_]
   Double_t        Muon_p4_fCoordinates_fY[kMaxMuon_p4];   //[Muon_p4_]
   Double_t        Muon_p4_fCoordinates_fZ[kMaxMuon_p4];   //[Muon_p4_]
   Double_t        Muon_p4_fCoordinates_fT[kMaxMuon_p4];   //[Muon_p4_]
   vector<vector<double> > *Muon_Poca;
   vector<bool>    *Muon_isGlobalMuon;
   vector<bool>    *Muon_isStandAloneMuon;
   vector<bool>    *Muon_isTrackerMuon;
   vector<bool>    *Muon_isCaloMuon;
   vector<bool>    *Muon_isIsolationValid;
   vector<bool>    *Muon_isQualityValid;
   vector<bool>    *Muon_isTimeValid;
   vector<float>   *Muon_emEt03;
   vector<float>   *Muon_emVetoEt03;
   vector<float>   *Muon_hadEt03;
   vector<float>   *Muon_hadVetoEt03;
   vector<float>   *Muon_nJets03;
   vector<float>   *Muon_nTracks03;
   vector<float>   *Muon_sumPt03;
   vector<float>   *Muon_trackerVetoPt03;
   vector<float>   *Muon_emEt05;
   vector<float>   *Muon_emVetoEt05;
   vector<float>   *Muon_hadEt05;
   vector<float>   *Muon_hadVetoEt05;
   vector<float>   *Muon_nJets05;
   vector<float>   *Muon_nTracks05;
   vector<float>   *Muon_sumPt05;
   vector<float>   *Muon_trackerVetoPt05;
   vector<bool>    *Muon_isIsolationValid;
   vector<unsigned int> *Muon_Track_idx;
   Int_t           PFTau_p4_;
   Double_t        PFTau_p4_fCoordinates_fX[kMaxPFTau_p4];   //[PFTau_p4_]
   Double_t        PFTau_p4_fCoordinates_fY[kMaxPFTau_p4];   //[PFTau_p4_]
   Double_t        PFTau_p4_fCoordinates_fZ[kMaxPFTau_p4];   //[PFTau_p4_]
   Double_t        PFTau_p4_fCoordinates_fT[kMaxPFTau_p4];   //[PFTau_p4_]
   vector<bool>    *PFTau_isTightIsolation;
   vector<bool>    *PFTau_isMediumIsolation;
   vector<bool>    *PFTau_isLooseIsolation;
   vector<int>     *PFTau_hpsDecayMode;
   vector<int>     *PFTau_Charge;
   vector<vector<int> > *PFTau_Track_idx;
   vector<bool>    *KFTau_discriminatorByKFit;
   vector<bool>    *KFTau_discriminatorByQC;
   vector<int>     *KFTau_nKinTaus;
   Int_t           KFTau_TauVis_p4_;
   Double_t        KFTau_TauVis_p4_fCoordinates_fX[kMaxKFTau_TauVis_p4];   //[KFTau_TauVis_p4_]
   Double_t        KFTau_TauVis_p4_fCoordinates_fY[kMaxKFTau_TauVis_p4];   //[KFTau_TauVis_p4_]
   Double_t        KFTau_TauVis_p4_fCoordinates_fZ[kMaxKFTau_TauVis_p4];   //[KFTau_TauVis_p4_]
   Double_t        KFTau_TauVis_p4_fCoordinates_fT[kMaxKFTau_TauVis_p4];   //[KFTau_TauVis_p4_]
   Int_t           KFTau_TauFit_p4_;
   Double_t        KFTau_TauFit_p4_fCoordinates_fX[kMaxKFTau_TauFit_p4];   //[KFTau_TauFit_p4_]
   Double_t        KFTau_TauFit_p4_fCoordinates_fY[kMaxKFTau_TauFit_p4];   //[KFTau_TauFit_p4_]
   Double_t        KFTau_TauFit_p4_fCoordinates_fZ[kMaxKFTau_TauFit_p4];   //[KFTau_TauFit_p4_]
   Double_t        KFTau_TauFit_p4_fCoordinates_fT[kMaxKFTau_TauFit_p4];   //[KFTau_TauFit_p4_]
   Int_t           KFTau_Neutrino_p4_;
   Double_t        KFTau_Neutrino_p4_fCoordinates_fX[kMaxKFTau_Neutrino_p4];   //[KFTau_Neutrino_p4_]
   Double_t        KFTau_Neutrino_p4_fCoordinates_fY[kMaxKFTau_Neutrino_p4];   //[KFTau_Neutrino_p4_]
   Double_t        KFTau_Neutrino_p4_fCoordinates_fZ[kMaxKFTau_Neutrino_p4];   //[KFTau_Neutrino_p4_]
   Double_t        KFTau_Neutrino_p4_fCoordinates_fT[kMaxKFTau_Neutrino_p4];   //[KFTau_Neutrino_p4_]
   vector<unsigned int> *KFTau_MatchedHPS_idx;
   vector<vector<int> > *KFTau_Track_idx;
   Int_t           PFJet_p4_;
   Double_t        PFJet_p4_fCoordinates_fX[kMaxPFJet_p4];   //[PFJet_p4_]
   Double_t        PFJet_p4_fCoordinates_fY[kMaxPFJet_p4];   //[PFJet_p4_]
   Double_t        PFJet_p4_fCoordinates_fZ[kMaxPFJet_p4];   //[PFJet_p4_]
   Double_t        PFJet_p4_fCoordinates_fT[kMaxPFJet_p4];   //[PFJet_p4_]
   vector<float>   *PFJet_chargedEmEnergy;
   vector<float>   *PFJet_chargedHadronEnergy;
   vector<float>   *PFJet_chargedHadronMultiplicity;
   vector<float>   *PFJet_chargedMuEnergy;
   vector<float>   *PFJet_chargedMultiplicity;
   vector<float>   *PFJet_electronEnergy;
   vector<float>   *PFJet_electronMultiplicity;
   vector<float>   *PFJet_HFEMEnergy;
   vector<float>   *PFJet_HFEMMultiplicity;
   vector<float>   *PFJet_HFHadronEnergy;
   vector<float>   *PFJet_HFHadronMultiplicity;
   vector<float>   *PFJet_muonEnergy;
   vector<float>   *PFJet_muonMultiplicity;
   vector<float>   *PFJet_neutralEmEnergy;
   vector<float>   *PFJet_neutralHadronEnergy;
   vector<float>   *PFJet_neutralHadronMultiplicity;
   vector<float>   *PFJet_photonEnergy;
   vector<float>   *PFJet_photonMultiplicity;
   vector<float>   *PFJet_jetArea;
   vector<float>   *PFJet_maxDistance;
   vector<int>     *PFJet_nConstituents;
   vector<float>   *PFJet_pileup;
   vector<float>   *PFJet_etaetaMoment;
   vector<float>   *PFJet_etaphiMoment;
   vector<vector<int> > *PFJet_Track_idx;
   vector<unsigned int> *PFJet_MatchedHPS_idx;
   Double_t        MET_et;
   Double_t        MET_phi;
   Double_t        MET_sumET;
   UInt_t          Event_EventNumber;
   UInt_t          Event_RunNumber;
   Int_t           Event_bunchCrossing;
   Int_t           Event_orbitNumber;
   UInt_t          Event_luminosityBlock;
   Bool_t          Event_isRealData;
   Int_t           Track_p4_;
   Double_t        Track_p4_fCoordinates_fX[kMaxTrack_p4];   //[Track_p4_]
   Double_t        Track_p4_fCoordinates_fY[kMaxTrack_p4];   //[Track_p4_]
   Double_t        Track_p4_fCoordinates_fZ[kMaxTrack_p4];   //[Track_p4_]
   Double_t        Track_p4_fCoordinates_fT[kMaxTrack_p4];   //[Track_p4_]
   vector<vector<double> > *Track_Poca;
   vector<int>     *Track_charge;
   vector<float>   *Track_chi2;
   vector<float>   *Track_ndof;
   vector<unsigned short> *Track_numberOfLostHits;
   vector<unsigned short> *Track_numberOfValidHits;
   vector<unsigned int> *Track_qualityMask;

   // List of branches
   TBranch        *b_Vtx_chi2;   //!
   TBranch        *b_Vtx_nTrk;   //!
   TBranch        *b_Vtx_ndof;   //!
   TBranch        *b_Vtx_chi2;   //!
   TBranch        *b_Vtx_x;   //!
   TBranch        *b_Vtx_y;   //!
   TBranch        *b_Vtx_z;   //!
   TBranch        *b_Vtx_Cov;   //!
   TBranch        *b_Vtx_Track_idx;   //!
   TBranch        *b_Muon_p4_;   //!
   TBranch        *b_Muon_p4_fCoordinates_fX;   //!
   TBranch        *b_Muon_p4_fCoordinates_fY;   //!
   TBranch        *b_Muon_p4_fCoordinates_fZ;   //!
   TBranch        *b_Muon_p4_fCoordinates_fT;   //!
   TBranch        *b_Muon_Poca;   //!
   TBranch        *b_Muon_isGlobalMuon;   //!
   TBranch        *b_Muon_isStandAloneMuon;   //!
   TBranch        *b_Muon_isTrackerMuon;   //!
   TBranch        *b_Muon_isCaloMuon;   //!
   TBranch        *b_Muon_isIsolationValid;   //!
   TBranch        *b_Muon_isQualityValid;   //!
   TBranch        *b_Muon_isTimeValid;   //!
   TBranch        *b_Muon_emEt03;   //!
   TBranch        *b_Muon_emVetoEt03;   //!
   TBranch        *b_Muon_hadEt03;   //!
   TBranch        *b_Muon_hadVetoEt03;   //!
   TBranch        *b_Muon_nJets03;   //!
   TBranch        *b_Muon_nTracks03;   //!
   TBranch        *b_Muon_sumPt03;   //!
   TBranch        *b_Muon_trackerVetoPt03;   //!
   TBranch        *b_Muon_emEt05;   //!
   TBranch        *b_Muon_emVetoEt05;   //!
   TBranch        *b_Muon_hadEt05;   //!
   TBranch        *b_Muon_hadVetoEt05;   //!
   TBranch        *b_Muon_nJets05;   //!
   TBranch        *b_Muon_nTracks05;   //!
   TBranch        *b_Muon_sumPt05;   //!
   TBranch        *b_Muon_trackerVetoPt05;   //!
   TBranch        *b_Muon_isIsolationValid;   //!
   TBranch        *b_Muon_Track_idx;   //!
   TBranch        *b_PFTau_p4_;   //!
   TBranch        *b_PFTau_p4_fCoordinates_fX;   //!
   TBranch        *b_PFTau_p4_fCoordinates_fY;   //!
   TBranch        *b_PFTau_p4_fCoordinates_fZ;   //!
   TBranch        *b_PFTau_p4_fCoordinates_fT;   //!
   TBranch        *b_PFTau_isTightIsolation;   //!
   TBranch        *b_PFTau_isMediumIsolation;   //!
   TBranch        *b_PFTau_isLooseIsolation;   //!
   TBranch        *b_PFTau_hpsDecayMode;   //!
   TBranch        *b_PFTau_Charge;   //!
   TBranch        *b_PFTau_Track_idx;   //!
   TBranch        *b_KFTau_discriminatorByKFit;   //!
   TBranch        *b_KFTau_discriminatorByQC;   //!
   TBranch        *b_KFTau_nKinTaus;   //!
   TBranch        *b_KFTau_TauVis_p4_;   //!
   TBranch        *b_KFTau_TauVis_p4_fCoordinates_fX;   //!
   TBranch        *b_KFTau_TauVis_p4_fCoordinates_fY;   //!
   TBranch        *b_KFTau_TauVis_p4_fCoordinates_fZ;   //!
   TBranch        *b_KFTau_TauVis_p4_fCoordinates_fT;   //!
   TBranch        *b_KFTau_TauFit_p4_;   //!
   TBranch        *b_KFTau_TauFit_p4_fCoordinates_fX;   //!
   TBranch        *b_KFTau_TauFit_p4_fCoordinates_fY;   //!
   TBranch        *b_KFTau_TauFit_p4_fCoordinates_fZ;   //!
   TBranch        *b_KFTau_TauFit_p4_fCoordinates_fT;   //!
   TBranch        *b_KFTau_Neutrino_p4_;   //!
   TBranch        *b_KFTau_Neutrino_p4_fCoordinates_fX;   //!
   TBranch        *b_KFTau_Neutrino_p4_fCoordinates_fY;   //!
   TBranch        *b_KFTau_Neutrino_p4_fCoordinates_fZ;   //!
   TBranch        *b_KFTau_Neutrino_p4_fCoordinates_fT;   //!
   TBranch        *b_KFTau_MatchedHPS_idx;   //!
   TBranch        *b_KFTau_Track_idx;   //!
   TBranch        *b_PFJet_p4_;   //!
   TBranch        *b_PFJet_p4_fCoordinates_fX;   //!
   TBranch        *b_PFJet_p4_fCoordinates_fY;   //!
   TBranch        *b_PFJet_p4_fCoordinates_fZ;   //!
   TBranch        *b_PFJet_p4_fCoordinates_fT;   //!
   TBranch        *b_PFJet_chargedEmEnergy;   //!
   TBranch        *b_PFJet_chargedHadronEnergy;   //!
   TBranch        *b_PFJet_chargedHadronMultiplicity;   //!
   TBranch        *b_PFJet_chargedMuEnergy;   //!
   TBranch        *b_PFJet_chargedMultiplicity;   //!
   TBranch        *b_PFJet_electronEnergy;   //!
   TBranch        *b_PFJet_electronMultiplicity;   //!
   TBranch        *b_PFJet_HFEMEnergy;   //!
   TBranch        *b_PFJet_HFEMMultiplicity;   //!
   TBranch        *b_PFJet_HFHadronEnergy;   //!
   TBranch        *b_PFJet_HFHadronMultiplicity;   //!
   TBranch        *b_PFJet_muonEnergy;   //!
   TBranch        *b_PFJet_muonMultiplicity;   //!
   TBranch        *b_PFJet_neutralEmEnergy;   //!
   TBranch        *b_PFJet_neutralHadronEnergy;   //!
   TBranch        *b_PFJet_neutralHadronMultiplicity;   //!
   TBranch        *b_PFJet_photonEnergy;   //!
   TBranch        *b_PFJet_photonMultiplicity;   //!
   TBranch        *b_PFJet_jetArea;   //!
   TBranch        *b_PFJet_maxDistance;   //!
   TBranch        *b_PFJet_nConstituents;   //!
   TBranch        *b_PFJet_pileup;   //!
   TBranch        *b_PFJet_etaetaMoment;   //!
   TBranch        *b_PFJet_etaphiMoment;   //!
   TBranch        *b_PFJet_Track_idx;   //!
   TBranch        *b_PFJet_MatchedHPS_idx;   //!
   TBranch        *b_MET_et;   //!
   TBranch        *b_MET_phi;   //!
   TBranch        *b_MET_sumET;   //!
   TBranch        *b_Event_EventNumber;   //!
   TBranch        *b_Event_RunNumber;   //!
   TBranch        *b_Event_bunchCrossing;   //!
   TBranch        *b_Event_orbitNumber;   //!
   TBranch        *b_Event_luminosityBlock;   //!
   TBranch        *b_Event_isRealData;   //!
   TBranch        *b_Track_p4_;   //!
   TBranch        *b_Track_p4_fCoordinates_fX;   //!
   TBranch        *b_Track_p4_fCoordinates_fY;   //!
   TBranch        *b_Track_p4_fCoordinates_fZ;   //!
   TBranch        *b_Track_p4_fCoordinates_fT;   //!
   TBranch        *b_Track_Poca;   //!
   TBranch        *b_Track_charge;   //!
   TBranch        *b_Track_chi2;   //!
   TBranch        *b_Track_ndof;   //!
   TBranch        *b_Track_numberOfLostHits;   //!
   TBranch        *b_Track_numberOfValidHits;   //!
   TBranch        *b_Track_qualityMask;   //!

   test(TTree *tree=0);
   virtual ~test();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef test_cxx
test::test(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("output.root");
      if (!f) {
         f = new TFile("output.root");
      }
      tree = (TTree*)gDirectory->Get("t");

   }
   Init(tree);
}

test::~test()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t test::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t test::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void test::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Vtx_chi2 = 0;
   Vtx_nTrk = 0;
   Vtx_ndof = 0;
   Vtx_chi2 = 0;
   Vtx_x = 0;
   Vtx_y = 0;
   Vtx_z = 0;
   Vtx_Cov = 0;
   Vtx_Track_idx = 0;
   Muon_Poca = 0;
   Muon_isGlobalMuon = 0;
   Muon_isStandAloneMuon = 0;
   Muon_isTrackerMuon = 0;
   Muon_isCaloMuon = 0;
   Muon_isIsolationValid = 0;
   Muon_isQualityValid = 0;
   Muon_isTimeValid = 0;
   Muon_emEt03 = 0;
   Muon_emVetoEt03 = 0;
   Muon_hadEt03 = 0;
   Muon_hadVetoEt03 = 0;
   Muon_nJets03 = 0;
   Muon_nTracks03 = 0;
   Muon_sumPt03 = 0;
   Muon_trackerVetoPt03 = 0;
   Muon_emEt05 = 0;
   Muon_emVetoEt05 = 0;
   Muon_hadEt05 = 0;
   Muon_hadVetoEt05 = 0;
   Muon_nJets05 = 0;
   Muon_nTracks05 = 0;
   Muon_sumPt05 = 0;
   Muon_trackerVetoPt05 = 0;
   Muon_isIsolationValid = 0;
   Muon_Track_idx = 0;
   PFTau_isTightIsolation = 0;
   PFTau_isMediumIsolation = 0;
   PFTau_isLooseIsolation = 0;
   PFTau_hpsDecayMode = 0;
   PFTau_Charge = 0;
   PFTau_Track_idx = 0;
   KFTau_discriminatorByKFit = 0;
   KFTau_discriminatorByQC = 0;
   KFTau_nKinTaus = 0;
   KFTau_MatchedHPS_idx = 0;
   KFTau_Track_idx = 0;
   PFJet_chargedEmEnergy = 0;
   PFJet_chargedHadronEnergy = 0;
   PFJet_chargedHadronMultiplicity = 0;
   PFJet_chargedMuEnergy = 0;
   PFJet_chargedMultiplicity = 0;
   PFJet_electronEnergy = 0;
   PFJet_electronMultiplicity = 0;
   PFJet_HFEMEnergy = 0;
   PFJet_HFEMMultiplicity = 0;
   PFJet_HFHadronEnergy = 0;
   PFJet_HFHadronMultiplicity = 0;
   PFJet_muonEnergy = 0;
   PFJet_muonMultiplicity = 0;
   PFJet_neutralEmEnergy = 0;
   PFJet_neutralHadronEnergy = 0;
   PFJet_neutralHadronMultiplicity = 0;
   PFJet_photonEnergy = 0;
   PFJet_photonMultiplicity = 0;
   PFJet_jetArea = 0;
   PFJet_maxDistance = 0;
   PFJet_nConstituents = 0;
   PFJet_pileup = 0;
   PFJet_etaetaMoment = 0;
   PFJet_etaphiMoment = 0;
   PFJet_Track_idx = 0;
   PFJet_MatchedHPS_idx = 0;
   Track_Poca = 0;
   Track_charge = 0;
   Track_chi2 = 0;
   Track_ndof = 0;
   Track_numberOfLostHits = 0;
   Track_numberOfValidHits = 0;
   Track_qualityMask = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Vtx_chi2", &Vtx_chi2, &b_Vtx_chi2);
   fChain->SetBranchAddress("Vtx_nTrk", &Vtx_nTrk, &b_Vtx_nTrk);
   fChain->SetBranchAddress("Vtx_ndof", &Vtx_ndof, &b_Vtx_ndof);
   fChain->SetBranchAddress("Vtx_chi2", &Vtx_chi2, &b_Vtx_chi2);
   fChain->SetBranchAddress("Vtx_x", &Vtx_x, &b_Vtx_x);
   fChain->SetBranchAddress("Vtx_y", &Vtx_y, &b_Vtx_y);
   fChain->SetBranchAddress("Vtx_z", &Vtx_z, &b_Vtx_z);
   fChain->SetBranchAddress("Vtx_Cov", &Vtx_Cov, &b_Vtx_Cov);
   fChain->SetBranchAddress("Vtx_Track_idx", &Vtx_Track_idx, &b_Vtx_Track_idx);
   fChain->SetBranchAddress("Muon_p4", &Muon_p4_, &b_Muon_p4_);
   fChain->SetBranchAddress("Muon_p4.fCoordinates.fX", Muon_p4_fCoordinates_fX, &b_Muon_p4_fCoordinates_fX);
   fChain->SetBranchAddress("Muon_p4.fCoordinates.fY", Muon_p4_fCoordinates_fY, &b_Muon_p4_fCoordinates_fY);
   fChain->SetBranchAddress("Muon_p4.fCoordinates.fZ", Muon_p4_fCoordinates_fZ, &b_Muon_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("Muon_p4.fCoordinates.fT", Muon_p4_fCoordinates_fT, &b_Muon_p4_fCoordinates_fT);
   fChain->SetBranchAddress("Muon_Poca", &Muon_Poca, &b_Muon_Poca);
   fChain->SetBranchAddress("Muon_isGlobalMuon", &Muon_isGlobalMuon, &b_Muon_isGlobalMuon);
   fChain->SetBranchAddress("Muon_isStandAloneMuon", &Muon_isStandAloneMuon, &b_Muon_isStandAloneMuon);
   fChain->SetBranchAddress("Muon_isTrackerMuon", &Muon_isTrackerMuon, &b_Muon_isTrackerMuon);
   fChain->SetBranchAddress("Muon_isCaloMuon", &Muon_isCaloMuon, &b_Muon_isCaloMuon);
   fChain->SetBranchAddress("Muon_isIsolationValid", &Muon_isIsolationValid, &b_Muon_isIsolationValid);
   fChain->SetBranchAddress("Muon_isQualityValid", &Muon_isQualityValid, &b_Muon_isQualityValid);
   fChain->SetBranchAddress("Muon_isTimeValid", &Muon_isTimeValid, &b_Muon_isTimeValid);
   fChain->SetBranchAddress("Muon_emEt03", &Muon_emEt03, &b_Muon_emEt03);
   fChain->SetBranchAddress("Muon_emVetoEt03", &Muon_emVetoEt03, &b_Muon_emVetoEt03);
   fChain->SetBranchAddress("Muon_hadEt03", &Muon_hadEt03, &b_Muon_hadEt03);
   fChain->SetBranchAddress("Muon_hadVetoEt03", &Muon_hadVetoEt03, &b_Muon_hadVetoEt03);
   fChain->SetBranchAddress("Muon_nJets03", &Muon_nJets03, &b_Muon_nJets03);
   fChain->SetBranchAddress("Muon_nTracks03", &Muon_nTracks03, &b_Muon_nTracks03);
   fChain->SetBranchAddress("Muon_sumPt03", &Muon_sumPt03, &b_Muon_sumPt03);
   fChain->SetBranchAddress("Muon_trackerVetoPt03", &Muon_trackerVetoPt03, &b_Muon_trackerVetoPt03);
   fChain->SetBranchAddress("Muon_emEt05", &Muon_emEt05, &b_Muon_emEt05);
   fChain->SetBranchAddress("Muon_emVetoEt05", &Muon_emVetoEt05, &b_Muon_emVetoEt05);
   fChain->SetBranchAddress("Muon_hadEt05", &Muon_hadEt05, &b_Muon_hadEt05);
   fChain->SetBranchAddress("Muon_hadVetoEt05", &Muon_hadVetoEt05, &b_Muon_hadVetoEt05);
   fChain->SetBranchAddress("Muon_nJets05", &Muon_nJets05, &b_Muon_nJets05);
   fChain->SetBranchAddress("Muon_nTracks05", &Muon_nTracks05, &b_Muon_nTracks05);
   fChain->SetBranchAddress("Muon_sumPt05", &Muon_sumPt05, &b_Muon_sumPt05);
   fChain->SetBranchAddress("Muon_trackerVetoPt05", &Muon_trackerVetoPt05, &b_Muon_trackerVetoPt05);
   fChain->SetBranchAddress("Muon_isIsolationValid", &Muon_isIsolationValid, &b_Muon_isIsolationValid);
   fChain->SetBranchAddress("Muon_Track_idx", &Muon_Track_idx, &b_Muon_Track_idx);
   fChain->SetBranchAddress("PFTau_p4", &PFTau_p4_, &b_PFTau_p4_);
   fChain->SetBranchAddress("PFTau_p4.fCoordinates.fX", PFTau_p4_fCoordinates_fX, &b_PFTau_p4_fCoordinates_fX);
   fChain->SetBranchAddress("PFTau_p4.fCoordinates.fY", PFTau_p4_fCoordinates_fY, &b_PFTau_p4_fCoordinates_fY);
   fChain->SetBranchAddress("PFTau_p4.fCoordinates.fZ", PFTau_p4_fCoordinates_fZ, &b_PFTau_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("PFTau_p4.fCoordinates.fT", PFTau_p4_fCoordinates_fT, &b_PFTau_p4_fCoordinates_fT);
   fChain->SetBranchAddress("PFTau_isTightIsolation", &PFTau_isTightIsolation, &b_PFTau_isTightIsolation);
   fChain->SetBranchAddress("PFTau_isMediumIsolation", &PFTau_isMediumIsolation, &b_PFTau_isMediumIsolation);
   fChain->SetBranchAddress("PFTau_isLooseIsolation", &PFTau_isLooseIsolation, &b_PFTau_isLooseIsolation);
   fChain->SetBranchAddress("PFTau_hpsDecayMode", &PFTau_hpsDecayMode, &b_PFTau_hpsDecayMode);
   fChain->SetBranchAddress("PFTau_Charge", &PFTau_Charge, &b_PFTau_Charge);
   fChain->SetBranchAddress("PFTau_Track_idx", &PFTau_Track_idx, &b_PFTau_Track_idx);
   fChain->SetBranchAddress("KFTau_discriminatorByKFit", &KFTau_discriminatorByKFit, &b_KFTau_discriminatorByKFit);
   fChain->SetBranchAddress("KFTau_discriminatorByQC", &KFTau_discriminatorByQC, &b_KFTau_discriminatorByQC);
   fChain->SetBranchAddress("KFTau_nKinTaus", &KFTau_nKinTaus, &b_KFTau_nKinTaus);
   fChain->SetBranchAddress("KFTau_TauVis_p4", &KFTau_TauVis_p4_, &b_KFTau_TauVis_p4_);
   fChain->SetBranchAddress("KFTau_TauVis_p4.fCoordinates.fX", KFTau_TauVis_p4_fCoordinates_fX, &b_KFTau_TauVis_p4_fCoordinates_fX);
   fChain->SetBranchAddress("KFTau_TauVis_p4.fCoordinates.fY", KFTau_TauVis_p4_fCoordinates_fY, &b_KFTau_TauVis_p4_fCoordinates_fY);
   fChain->SetBranchAddress("KFTau_TauVis_p4.fCoordinates.fZ", KFTau_TauVis_p4_fCoordinates_fZ, &b_KFTau_TauVis_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("KFTau_TauVis_p4.fCoordinates.fT", KFTau_TauVis_p4_fCoordinates_fT, &b_KFTau_TauVis_p4_fCoordinates_fT);
   fChain->SetBranchAddress("KFTau_TauFit_p4", &KFTau_TauFit_p4_, &b_KFTau_TauFit_p4_);
   fChain->SetBranchAddress("KFTau_TauFit_p4.fCoordinates.fX", KFTau_TauFit_p4_fCoordinates_fX, &b_KFTau_TauFit_p4_fCoordinates_fX);
   fChain->SetBranchAddress("KFTau_TauFit_p4.fCoordinates.fY", KFTau_TauFit_p4_fCoordinates_fY, &b_KFTau_TauFit_p4_fCoordinates_fY);
   fChain->SetBranchAddress("KFTau_TauFit_p4.fCoordinates.fZ", KFTau_TauFit_p4_fCoordinates_fZ, &b_KFTau_TauFit_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("KFTau_TauFit_p4.fCoordinates.fT", KFTau_TauFit_p4_fCoordinates_fT, &b_KFTau_TauFit_p4_fCoordinates_fT);
   fChain->SetBranchAddress("KFTau_Neutrino_p4", &KFTau_Neutrino_p4_, &b_KFTau_Neutrino_p4_);
   fChain->SetBranchAddress("KFTau_Neutrino_p4.fCoordinates.fX", KFTau_Neutrino_p4_fCoordinates_fX, &b_KFTau_Neutrino_p4_fCoordinates_fX);
   fChain->SetBranchAddress("KFTau_Neutrino_p4.fCoordinates.fY", KFTau_Neutrino_p4_fCoordinates_fY, &b_KFTau_Neutrino_p4_fCoordinates_fY);
   fChain->SetBranchAddress("KFTau_Neutrino_p4.fCoordinates.fZ", KFTau_Neutrino_p4_fCoordinates_fZ, &b_KFTau_Neutrino_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("KFTau_Neutrino_p4.fCoordinates.fT", KFTau_Neutrino_p4_fCoordinates_fT, &b_KFTau_Neutrino_p4_fCoordinates_fT);
   fChain->SetBranchAddress("KFTau_MatchedHPS_idx", &KFTau_MatchedHPS_idx, &b_KFTau_MatchedHPS_idx);
   fChain->SetBranchAddress("KFTau_Track_idx", &KFTau_Track_idx, &b_KFTau_Track_idx);
   fChain->SetBranchAddress("PFJet_p4", &PFJet_p4_, &b_PFJet_p4_);
   fChain->SetBranchAddress("PFJet_p4.fCoordinates.fX", PFJet_p4_fCoordinates_fX, &b_PFJet_p4_fCoordinates_fX);
   fChain->SetBranchAddress("PFJet_p4.fCoordinates.fY", PFJet_p4_fCoordinates_fY, &b_PFJet_p4_fCoordinates_fY);
   fChain->SetBranchAddress("PFJet_p4.fCoordinates.fZ", PFJet_p4_fCoordinates_fZ, &b_PFJet_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("PFJet_p4.fCoordinates.fT", PFJet_p4_fCoordinates_fT, &b_PFJet_p4_fCoordinates_fT);
   fChain->SetBranchAddress("PFJet_chargedEmEnergy", &PFJet_chargedEmEnergy, &b_PFJet_chargedEmEnergy);
   fChain->SetBranchAddress("PFJet_chargedHadronEnergy", &PFJet_chargedHadronEnergy, &b_PFJet_chargedHadronEnergy);
   fChain->SetBranchAddress("PFJet_chargedHadronMultiplicity", &PFJet_chargedHadronMultiplicity, &b_PFJet_chargedHadronMultiplicity);
   fChain->SetBranchAddress("PFJet_chargedMuEnergy", &PFJet_chargedMuEnergy, &b_PFJet_chargedMuEnergy);
   fChain->SetBranchAddress("PFJet_chargedMultiplicity", &PFJet_chargedMultiplicity, &b_PFJet_chargedMultiplicity);
   fChain->SetBranchAddress("PFJet_electronEnergy", &PFJet_electronEnergy, &b_PFJet_electronEnergy);
   fChain->SetBranchAddress("PFJet_electronMultiplicity", &PFJet_electronMultiplicity, &b_PFJet_electronMultiplicity);
   fChain->SetBranchAddress("PFJet_HFEMEnergy", &PFJet_HFEMEnergy, &b_PFJet_HFEMEnergy);
   fChain->SetBranchAddress("PFJet_HFEMMultiplicity", &PFJet_HFEMMultiplicity, &b_PFJet_HFEMMultiplicity);
   fChain->SetBranchAddress("PFJet_HFHadronEnergy", &PFJet_HFHadronEnergy, &b_PFJet_HFHadronEnergy);
   fChain->SetBranchAddress("PFJet_HFHadronMultiplicity", &PFJet_HFHadronMultiplicity, &b_PFJet_HFHadronMultiplicity);
   fChain->SetBranchAddress("PFJet_muonEnergy", &PFJet_muonEnergy, &b_PFJet_muonEnergy);
   fChain->SetBranchAddress("PFJet_muonMultiplicity", &PFJet_muonMultiplicity, &b_PFJet_muonMultiplicity);
   fChain->SetBranchAddress("PFJet_neutralEmEnergy", &PFJet_neutralEmEnergy, &b_PFJet_neutralEmEnergy);
   fChain->SetBranchAddress("PFJet_neutralHadronEnergy", &PFJet_neutralHadronEnergy, &b_PFJet_neutralHadronEnergy);
   fChain->SetBranchAddress("PFJet_neutralHadronMultiplicity", &PFJet_neutralHadronMultiplicity, &b_PFJet_neutralHadronMultiplicity);
   fChain->SetBranchAddress("PFJet_photonEnergy", &PFJet_photonEnergy, &b_PFJet_photonEnergy);
   fChain->SetBranchAddress("PFJet_photonMultiplicity", &PFJet_photonMultiplicity, &b_PFJet_photonMultiplicity);
   fChain->SetBranchAddress("PFJet_jetArea", &PFJet_jetArea, &b_PFJet_jetArea);
   fChain->SetBranchAddress("PFJet_maxDistance", &PFJet_maxDistance, &b_PFJet_maxDistance);
   fChain->SetBranchAddress("PFJet_nConstituents", &PFJet_nConstituents, &b_PFJet_nConstituents);
   fChain->SetBranchAddress("PFJet_pileup", &PFJet_pileup, &b_PFJet_pileup);
   fChain->SetBranchAddress("PFJet_etaetaMoment", &PFJet_etaetaMoment, &b_PFJet_etaetaMoment);
   fChain->SetBranchAddress("PFJet_etaphiMoment", &PFJet_etaphiMoment, &b_PFJet_etaphiMoment);
   fChain->SetBranchAddress("PFJet_Track_idx", &PFJet_Track_idx, &b_PFJet_Track_idx);
   fChain->SetBranchAddress("PFJet_MatchedHPS_idx", &PFJet_MatchedHPS_idx, &b_PFJet_MatchedHPS_idx);
   fChain->SetBranchAddress("MET_et", &MET_et, &b_MET_et);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   fChain->SetBranchAddress("MET_sumET", &MET_sumET, &b_MET_sumET);
   fChain->SetBranchAddress("Event_EventNumber", &Event_EventNumber, &b_Event_EventNumber);
   fChain->SetBranchAddress("Event_RunNumber", &Event_RunNumber, &b_Event_RunNumber);
   fChain->SetBranchAddress("Event_bunchCrossing", &Event_bunchCrossing, &b_Event_bunchCrossing);
   fChain->SetBranchAddress("Event_orbitNumber", &Event_orbitNumber, &b_Event_orbitNumber);
   fChain->SetBranchAddress("Event_luminosityBlock", &Event_luminosityBlock, &b_Event_luminosityBlock);
   fChain->SetBranchAddress("Event_isRealData", &Event_isRealData, &b_Event_isRealData);
   fChain->SetBranchAddress("Track_p4", &Track_p4_, &b_Track_p4_);
   fChain->SetBranchAddress("Track_p4.fCoordinates.fX", Track_p4_fCoordinates_fX, &b_Track_p4_fCoordinates_fX);
   fChain->SetBranchAddress("Track_p4.fCoordinates.fY", Track_p4_fCoordinates_fY, &b_Track_p4_fCoordinates_fY);
   fChain->SetBranchAddress("Track_p4.fCoordinates.fZ", Track_p4_fCoordinates_fZ, &b_Track_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("Track_p4.fCoordinates.fT", Track_p4_fCoordinates_fT, &b_Track_p4_fCoordinates_fT);
   fChain->SetBranchAddress("Track_Poca", &Track_Poca, &b_Track_Poca);
   fChain->SetBranchAddress("Track_charge", &Track_charge, &b_Track_charge);
   fChain->SetBranchAddress("Track_chi2", &Track_chi2, &b_Track_chi2);
   fChain->SetBranchAddress("Track_ndof", &Track_ndof, &b_Track_ndof);
   fChain->SetBranchAddress("Track_numberOfLostHits", &Track_numberOfLostHits, &b_Track_numberOfLostHits);
   fChain->SetBranchAddress("Track_numberOfValidHits", &Track_numberOfValidHits, &b_Track_numberOfValidHits);
   fChain->SetBranchAddress("Track_qualityMask", &Track_qualityMask, &b_Track_qualityMask);
   Notify();
}

Bool_t test::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void test::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t test::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef test_cxx
