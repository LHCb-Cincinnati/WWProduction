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

// Particle Structure
// A structure to contain the data for a particle.
struct ParticleStruct{
  double pT;
  double p;
  double eta;
  double energy;
  double phi;
  double m0;
  int id;
  int charge;
  int status;
};

// A function to fill the ParticleStruct structure given an event and index.
int struct_fill(ParticleStruct* particle_structref, Pythia8::Event event, int index){
  double px = event[index].px();
  double py = event[index].py();
  double pz = event[index].pz();
  particle_structref->p = sqrt(px*px + py*py + pz*pz);
  particle_structref->id = event[index].id();
  particle_structref->status = event[index].status();
  particle_structref->charge = event[index].charge();
  particle_structref->pT = event[index].pT();
  particle_structref->eta = event[index].eta();
  particle_structref->energy = event[index].e();
  particle_structref->phi = event[index].phi();
  particle_structref->m0 = event[index].m0();
  return(0);
}

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
  TFile* outFile = new TFile("WeakSingleBosonDecay.root", "RECREATE");

  ParticleStruct w_struct;
  ParticleStruct muon_struct;
  ParticleStruct neutrino_struct;

  // ROOT objects
  TTree *Tree = new TTree("Tree","Tree");
  Tree->Branch("WBoson",&w_struct,"pT/D:p/D:eta/D:energy/D:phi/D:m0/D:id/I:charge/I:status/I");
  Tree->Branch("Muon",&muon_struct,"pT/D:p/D:eta/D:energy/D:phi/D:m0/D:id/I:charge/I:status/I");
  Tree->Branch("Neutrino",&neutrino_struct,"pT/D:p/D:eta/D:energy/D:phi/D:m0/D:id/I:charge/I:status/I");

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

    // Fill out the structs with event data.
    struct_fill(&w_struct, pythia.event, iW);
    struct_fill(&muon_struct, pythia.event, imuon);
    struct_fill(&neutrino_struct, pythia.event, ineutrino);

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
