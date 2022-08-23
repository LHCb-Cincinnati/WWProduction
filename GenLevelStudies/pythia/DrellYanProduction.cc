// Program to simulate W Boson Decay at the LHC

// Stdlib header file for input and output.
#include <iostream>
#include <math.h>
#include <fstream>

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
  double px;
  double py;
  double pz;
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
  particle_structref->px = event[index].px();
  particle_structref->py = event[index].py();
  particle_structref->pz = event[index].pz();
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

// A function to compute the invariant mass given two four vectors.
double compute_invariant_mass(Pythia8::Vec4 v1, Pythia8::Vec4 v2){
  Pythia8::Vec4 v_new = v1+v2;
  return(v_new.mCalc());
}

int main(int argc, char* argv[]) {

  // Redirect std out to a log file.
  string fname = "DrellYanProduction.log";
  std::ofstream outstream(fname);
  std::streambuf* filebuf = outstream.rdbuf();
  std::streambuf* coutbuf = cout.rdbuf();
  cout.rdbuf(filebuf);

  // Create the ROOT application environment.
  TApplication theApp("hist", &argc, argv);

  // Start Pythia
  Pythia pythia;

  // Read in commands from external file.
  if (argc != 2)
    pythia.readFile("DrellYanProduction.cmnd");
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

  // Initialize Pythia
  pythia.init(); 

  // Extract information from file
  int nEvent = pythia.mode("Main:numberOfEvents");

  // Create file on which histogram(s) can be saved.
  TFile* outFile = new TFile("DrellYanProduction.root", "RECREATE");

  // Some Definitions
  int ilminus, ilplus;
  double dilepton_invariant_mass;
  ParticleStruct lepton_plus_struct;
  ParticleStruct lepton_minus_struct;
  TTree *Tree = new TTree("AnalysisTree","AnalysisTree");
  Tree->Branch("LeptonPlus",&lepton_plus_struct,"px/D:py/D:pz/D:pT/D:p/D:eta/D:energy/D:phi/D:m0/D:id/I:charge/I:status/I");
  Tree->Branch("LeptonMinus",&lepton_minus_struct,"px/D:py/D:pz/D:pT/D:p/D:eta/D:energy/D:phi/D:m0/D:id/I:charge/I:status/I");
  Tree->Branch("DiLeptonInvariantMass",&dilepton_invariant_mass);

  // Begin event loop. Generate event; skip if generation aborted.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;
    // Debugging Help
    // pythia.event.list(); // Used to print out each event.
    // cout << "Event: " << iEvent << endl;

    // Find the muons in the event
    for (int i = 0; i < pythia.event.size(); ++i){
      if (pythia.event[i].isFinal()){
        if (pythia.event[i].id() == -13) ilplus = i;
        else if (pythia.event[i].id() == 13) ilminus = i;
      }
    } 

    // Fill out the structs with event data.
    struct_fill(&lepton_plus_struct, pythia.event, ilplus);
    struct_fill(&lepton_minus_struct, pythia.event, ilminus);
    dilepton_invariant_mass = compute_invariant_mass(pythia.event[ilminus].p(), pythia.event[ilplus].p());

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
