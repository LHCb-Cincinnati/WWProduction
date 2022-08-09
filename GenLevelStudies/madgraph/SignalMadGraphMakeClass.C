#define SignalMadGraphMakeClass_cxx
#include "SignalMadGraphMakeClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

// More General Imports
#include <math.h>
#include <algorithm>

void SignalMadGraphMakeClass::Loop()
{
//   In a ROOT session, you can do:
//      root> .L SignalMadGraphMakeClass.C
//      root> SignalMadGraphMakeClass t
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

   // Useful Quantities
   // Output File
   TFile ofile("test.root","RECREATE");
   
   // Set Default Indices.
   int lepton_index = 11;
   int antilepton_index = 11;
   int neutrino_index = 11;
   int antineutrino_index = 11;

   // Scaling Histograms
   double luminosity = 1.0; // Unit luminosity measured in fb^-1
   double xsection = 5045.0; // Measured in fb
   double n_gen = 100000.0; // Number of events generated
   double scale_factor = (luminosity * xsection)/n_gen;

   // Quantities to calculate
   ROOT::Math::PxPyPzEVector missing_vector;
   ROOT::Math::PxPyPzEVector lepton_pair_vector;
   double dilepton_invariant_mass;
   double leading_lepton_pt;
   double trailing_lepton_pt;
   double delta_phi;
   double pT_miss_proj;

   // Defining Vectors
   ROOT::Math::PxPyPzEVector antilepton_vector;
   ROOT::Math::PxPyPzEVector lepton_vector;
   ROOT::Math::PxPyPzEVector neutrino_vector;
   ROOT::Math::PxPyPzEVector antineutrino_vector;

   // Tree/Histogram/File Definitions
   TH1F* lepton_pair_invariant_mass_hist = new TH1F("Lepton Pair Invariant Mass", "Lepton Pair Invariant Mass", 100, 0, 200);
   TH1F* lepton_pair_energy_hist = new TH1F("Lepton Pair Transverse Energy", "Lepton Pair Transverse Energy", 100, 0, 150);
   TH1F* lepton_pair_pT_hist = new TH1F("Lepton Pair pT", "Lepton Pair pT", 100, 0, 150);
   TH1F* leading_lepton_pT_hist = new TH1F("Leading Lepton pT", "Leading Lepton pT", 100, 0, 150);
   TH1F* trailing_lepton_pT_hist = new TH1F("Trailing Lepton pT", "Trailing Lepton pT", 100, 0, 150);
   TH1F* missing_energy_hist = new TH1F("Missing Transverse Energy", "Missing Transverse Energy", 100, 0, 150);
   TH1F* missing_pT_hist = new TH1F("Missing pT", "Missing pT", 100, 0, 150);
   TH1F* delta_phi_hist = new TH1F("Delta Phi", "Delta Phi", 100, 0, 3.14);
   TH1F* missing_pT_proj_hist = new TH1F("Missing pT proj", "Missing pT proj", 100, 0, 150);
   TH1F* test_lepton_pair_invariant_mass_hist = new TH1F("Test Lepton Pair Invariant Mass", "Test Lepton Pair Invariant Mass", 100, 0, 200);


   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      // Finding appropriate particle indices.
      for(int i=0; i < *Event_Nparticles; i++){
         int PID = Particle_PID[i];
         if (PID == 13 or PID == 11)  lepton_index = i;
         else if (PID == -13 or PID == -11)  antilepton_index = i;
         else if (PID == 12 or PID == 14)  neutrino_index = i;
         else if (PID == -12 or PID == -14)  antineutrino_index = i;
      }

      // Eta Cut
      if ((Particle_Eta[lepton_index] > 2) &&  (Particle_Eta[antilepton_index] > 2) && (Particle_Eta[lepton_index] < 5) &&  (Particle_Eta[antilepton_index] < 5)){

         // Filling Vectors
         antilepton_vector = ROOT::Math::PxPyPzEVector(Particle_Px[antilepton_index], Particle_Py[antilepton_index], Particle_Pz[antilepton_index], Particle_E[antilepton_index]);
         lepton_vector = ROOT::Math::PxPyPzEVector(Particle_Px[lepton_index], Particle_Py[lepton_index], Particle_Pz[lepton_index], Particle_E[lepton_index]);
         neutrino_vector = ROOT::Math::PxPyPzEVector(Particle_Px[neutrino_index], Particle_Py[neutrino_index], Particle_Pz[neutrino_index], Particle_E[neutrino_index]);
         antineutrino_vector = ROOT::Math::PxPyPzEVector(Particle_Px[antineutrino_index], Particle_Py[antineutrino_index], Particle_Pz[antineutrino_index], Particle_E[antineutrino_index]);

         // Calculating Quantities
         // Missing Vector
         missing_vector = neutrino_vector + antineutrino_vector;

         // Lepton Pair Vector
         lepton_pair_vector = lepton_vector + antilepton_vector;

         // DiLepton Invariant Mass
         dilepton_invariant_mass = lepton_pair_vector.M();

         // Leading and Trailing Leptons
         if (lepton_vector.Pt() > antilepton_vector.Pt()){
            leading_lepton_pt = lepton_vector.Pt();
            trailing_lepton_pt = antilepton_vector.Pt();
         }
         else{
            leading_lepton_pt = antilepton_vector.Pt();
            trailing_lepton_pt = lepton_vector.Pt();
         }

         // Delta Phi
         if (abs(ROOT::Math::VectorUtil::DeltaPhi(lepton_vector, missing_vector)) < abs(ROOT::Math::VectorUtil::DeltaPhi(antilepton_vector, missing_vector))){
            delta_phi = abs(ROOT::Math::VectorUtil::DeltaPhi(lepton_vector, missing_vector));
         }
         else{
            delta_phi = abs(ROOT::Math::VectorUtil::DeltaPhi(antilepton_vector, missing_vector));
         }

         // pT Miss proj
         pT_miss_proj = missing_vector.Pt()*sin(delta_phi);

         if (*Event_Nparticles != 4) test_lepton_pair_invariant_mass_hist->Fill(dilepton_invariant_mass);

         // Filling Histograms
         lepton_pair_invariant_mass_hist->Fill(dilepton_invariant_mass);
         lepton_pair_energy_hist->Fill(lepton_pair_vector.Et());
         lepton_pair_pT_hist->Fill(lepton_pair_vector.Pt());
         leading_lepton_pT_hist->Fill(leading_lepton_pt);
         trailing_lepton_pT_hist->Fill(trailing_lepton_pt);
         missing_energy_hist->Fill(missing_vector.Et());
         missing_pT_hist->Fill(missing_vector.Pt());
         delta_phi_hist->Fill(delta_phi);
         missing_pT_proj_hist->Fill(pT_miss_proj);
      } // End of eta if loop


   } // End of Event Loop
   
   lepton_pair_invariant_mass_hist->Scale(scale_factor);
   lepton_pair_energy_hist->Scale(scale_factor);
   lepton_pair_pT_hist->Scale(scale_factor);
   leading_lepton_pT_hist->Scale(scale_factor);
   trailing_lepton_pT_hist->Scale(scale_factor);
   missing_energy_hist->Scale(scale_factor);
   missing_pT_hist->Scale(scale_factor);
   delta_phi_hist->Scale(scale_factor);
   missing_pT_proj_hist->Scale(scale_factor);

   // Writing to File
   lepton_pair_invariant_mass_hist->Write();
   lepton_pair_energy_hist->Write();
   lepton_pair_pT_hist->Write();
   leading_lepton_pT_hist->Write();
   trailing_lepton_pT_hist->Write();
   missing_energy_hist->Write();
   missing_pT_hist->Write();
   delta_phi_hist->Write();
   missing_pT_proj_hist->Write();
   test_lepton_pair_invariant_mass_hist->Write();

   // Closing File
   ofile.Close();
}
