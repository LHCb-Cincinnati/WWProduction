// Program to simulate W Boson Decay at the LHC

// Stdlib header file for input and output.
#include <iostream>
#include <math.h>

// Header file to access Pythia 8 program elements.
#include "/Applications/pythia8307/include/Pythia8/Pythia.h"

// ROOT stuff
#include "TH1.h"
#include "TTree.h"
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TFile.h"

using namespace Pythia8;

int main(int argc, char* argv[]) {

  // Create the ROOT application environment.
  TApplication theApp("hist", &argc, argv);

  // Start Pythia
  Pythia pythia;

  // Read in commands from external file.
  if (argc != 2)
    pythia.readFile("WeakSingleBosonDecay.cmnd");
  else {
    // Check that the provided input name corresponds to an existing file.
    ifstream is(argv[1]);
    if (!is) {
      cerr << " Command-line file " << argv[1] << " was not found. \n"
           << " Program stopped! " << endl;
      return 1;
    }
    pythia.readFile(argv[1]);
  }
  pythia.init(); 


  // Extract information from file
  int nEvent = pythia.mode("Main:numberOfEvents");

  // Create file on which histogram(s) can be saved.
  TFile* outFile = new TFile("WeakBosonSingleDecay.root", "RECREATE");

  // Particle Structure
  struct ParticleStruct{
    double pT;
    double p;
    double eta;
    double energy;
    double phi;
    int id;
    int charge;
    int status;
  };
  ParticleStruct w_struct;
  ParticleStruct muon_struct;
  ParticleStruct neutrino_struct;

  // ROOT objects
  TTree *Tree = new TTree("Tree","Tree");
  Tree->Branch("WBoson",&w_struct,"pT/D:p/D:eta/D:energy/D:phi/D:id/I:charge/I:status/I");
  Tree->Branch("Muon",&muon_struct,"pT/D:p/D:eta/D:energy/D:phi/D:id/I:charge/I:status/I");
  Tree->Branch("Neutrino",&neutrino_struct,"pT/D:p/D:eta/D:energy/D:phi/D:id/I:charge/I:status/I");

  // Begin event loop. Generate event; skip if generation aborted.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;
    int iW = 0;
    // pythia.event.list(); // Used to print out each event.

    // Find last W boson
    for (int i = 0; i < pythia.event.size(); ++i){
      if (pythia.event[i].id() == 24 || pythia.event[i].id() == -24)
        iW = i;
    }

    // Get indices of daughter particles from W decay.
    int imuon = pythia.event[iW].daughter1();
    int ineutrino = pythia.event[iW].daughter2();
    // Loop through daughter particles.
    // vector<int> daughter_list = pythia.event[iW].daughterList();
    // for (int i = 0; i < daughter_list.size(); i++){
    //   cout << i+1 << ": " << pythia.event[daughter_list[i]].id() << endl;
    // }

    // Fill out the structs with event data.
    double px = pythia.event[iW].px();
    double py = pythia.event[iW].py();
    double pz = pythia.event[iW].pz();
    w_struct.p = sqrt(px*px + py*py + pz*pz);
    w_struct.id = pythia.event[iW].id();
    w_struct.status = pythia.event[iW].status();
    w_struct.charge = pythia.event[iW].charge();
    w_struct.pT = pythia.event[iW].pT();
    w_struct.eta = pythia.event[iW].eta();
    w_struct.energy = pythia.event[iW].e();
    w_struct.phi = pythia.event[iW].phi();

    // Fill the tree with the new data
    Tree->Fill();
  }

  // Statistics on event generation.
  pythia.stat();

  // Save tree on file and close file.
  Tree->Write();
  delete outFile;

  // Done.
  return 0;
}
