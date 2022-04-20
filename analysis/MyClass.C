#define MyClass_cxx
#include "MyClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void MyClass::Loop()
{
//   In a ROOT session, you can do:
//      root> .L MyClass.C
//      root> MyClass t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   // Define an output file.
   TFile out_file("WW_DoubleBosonDecayOutput.root","RECREATE");

   // Define some histograms
   TH1F* plepton_eta_hist = new TH1F("Positive Lepton Eta", "Positive Lepton Eta",100, 2, 5);
   TH1F* plepton_pt_hist = new TH1F("Positive Lepton pT", "Positive Lepton pT",100, 0, 200);
   TH1F* mlepton_eta_hist = new TH1F("Minus Lepton Eta", "Minus Lepton Eta",100, 2, 5);
   TH1F* mlepton_pt_hist = new TH1F("Minus Lepton pT", "Minus Lepton pT",100, 0, 200);


   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      // Gabe's Code
      if (2<LeptonPlus_eta && LeptonPlus_eta<5 && 2<LeptonMinus_eta && LeptonMinus_eta<5){
         plepton_eta_hist->Fill(LeptonPlus_eta);
         plepton_pt_hist->Fill(LeptonPlus_pT);
         mlepton_eta_hist->Fill(LeptonMinus_eta);
         mlepton_pt_hist->Fill(LeptonMinus_pT);
      }

   } // Entry Loop Ends

   // Found estimated number of events in detector to be:
   // Original number of events in simulation: 10000
   // n = L * sigma = 1.665E15*7.528E-13 = 1,253
   // Rescaling the histograms
   double scaling_factor = 1253 / 10000;
   plepton_eta_hist->Scale(scaling_factor);
   plepton_pt_hist->Scale(scaling_factor);
   mlepton_eta_hist->Scale(scaling_factor);
   mlepton_pt_hist->Scale(scaling_factor);

   // Writing to the file
   plepton_eta_hist->Write();
   plepton_pt_hist->Write();
   mlepton_eta_hist->Write();
   mlepton_pt_hist->Write();
   out_file.Close();
}
