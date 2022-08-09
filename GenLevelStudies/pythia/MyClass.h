//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr 20 13:41:06 2022 by ROOT version 6.24/06
// from TTree AnalysisTree/AnalysisTree
// found on file: ../pythia/WW_DoubleBosonDecay.root
//////////////////////////////////////////////////////////

#ifndef MyClass_h
#define MyClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class MyClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        WeakBosonPlus_pT;
   Double_t        WeakBosonPlus_p;
   Double_t        WeakBosonPlus_eta;
   Double_t        WeakBosonPlus_energy;
   Double_t        WeakBosonPlus_phi;
   Double_t        WeakBosonPlus_m0;
   Int_t           WeakBosonPlus_id;
   Int_t           WeakBosonPlus_charge;
   Int_t           WeakBosonPlus_status;
   Double_t        LeptonPlus_pT;
   Double_t        LeptonPlus_p;
   Double_t        LeptonPlus_eta;
   Double_t        LeptonPlus_energy;
   Double_t        LeptonPlus_phi;
   Double_t        LeptonPlus_m0;
   Int_t           LeptonPlus_id;
   Int_t           LeptonPlus_charge;
   Int_t           LeptonPlus_status;
   Double_t        NeutrinoPlus_pT;
   Double_t        NeutrinoPlus_p;
   Double_t        NeutrinoPlus_eta;
   Double_t        NeutrinoPlus_energy;
   Double_t        NeutrinoPlus_phi;
   Double_t        NeutrinoPlus_m0;
   Int_t           NeutrinoPlus_id;
   Int_t           NeutrinoPlus_charge;
   Int_t           NeutrinoPlus_status;
   Double_t        WeakBosonMinus_pT;
   Double_t        WeakBosonMinus_p;
   Double_t        WeakBosonMinus_eta;
   Double_t        WeakBosonMinus_energy;
   Double_t        WeakBosonMinus_phi;
   Double_t        WeakBosonMinus_m0;
   Int_t           WeakBosonMinus_id;
   Int_t           WeakBosonMinus_charge;
   Int_t           WeakBosonMinus_status;
   Double_t        LeptonMinus_pT;
   Double_t        LeptonMinus_p;
   Double_t        LeptonMinus_eta;
   Double_t        LeptonMinus_energy;
   Double_t        LeptonMinus_phi;
   Double_t        LeptonMinus_m0;
   Int_t           LeptonMinus_id;
   Int_t           LeptonMinus_charge;
   Int_t           LeptonMinus_status;
   Double_t        NeutrinoMinus_pT;
   Double_t        NeutrinoMinus_p;
   Double_t        NeutrinoMinus_eta;
   Double_t        NeutrinoMinus_energy;
   Double_t        NeutrinoMinus_phi;
   Double_t        NeutrinoMinus_m0;
   Int_t           NeutrinoMinus_id;
   Int_t           NeutrinoMinus_charge;
   Int_t           NeutrinoMinus_status;

   // List of branches
   TBranch        *b_WeakBosonPlus;   //!
   TBranch        *b_LeptonPlus;   //!
   TBranch        *b_NeutrinoPlus;   //!
   TBranch        *b_WeakBosonMinus;   //!
   TBranch        *b_LeptonMinus;   //!
   TBranch        *b_NeutrinoMinus;   //!

   MyClass(TTree *tree=0);
   virtual ~MyClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MyClass_cxx
MyClass::MyClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../pythia/WW_DoubleBosonDecay.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../pythia/WW_DoubleBosonDecay.root");
      }
      f->GetObject("AnalysisTree",tree);

   }
   Init(tree);
}

MyClass::~MyClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MyClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MyClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MyClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("WeakBosonPlus", &WeakBosonPlus_pT, &b_WeakBosonPlus);
   fChain->SetBranchAddress("LeptonPlus", &LeptonPlus_pT, &b_LeptonPlus);
   fChain->SetBranchAddress("NeutrinoPlus", &NeutrinoPlus_pT, &b_NeutrinoPlus);
   fChain->SetBranchAddress("WeakBosonMinus", &WeakBosonMinus_pT, &b_WeakBosonMinus);
   fChain->SetBranchAddress("LeptonMinus", &LeptonMinus_pT, &b_LeptonMinus);
   fChain->SetBranchAddress("NeutrinoMinus", &NeutrinoMinus_pT, &b_NeutrinoMinus);
   Notify();
}

Bool_t MyClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MyClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MyClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MyClass_cxx
