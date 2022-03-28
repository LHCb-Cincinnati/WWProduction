// Program to simulate W Boson Decay at the LHC

// Stdlib header file for input and output.
#include <iostream>

// Header file to access Pythia 8 program elements.
#include "/Applications/pythia8307/include/Pythia8/Pythia.h"

// ROOT, for histogramming.
#include "TH1.h"
#include "TTree.h"

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT, for saving file.
#include "TFile.h"

using namespace Pythia8;

int main(int argc, char* argv[]) {

  // Create the ROOT application environment.
  TApplication theApp("hist", &argc, argv);

  // Start Pythia
  Pythia pythia;

  // Read in commands from external file.
  if (argc != 2)
    pythia.readFile("WeakBosonDecay.cmnd");
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
  Hist pTZ("dN/dpTZ", 100, 0., 100.);


  // Extract information from file
  int nEvent = pythia.mode("Main:numberOfEvents");

  // Create file on which histogram(s) can be saved.
  TFile* outFile = new TFile("WeakBosonDecay.root", "RECREATE");

  // Structure
  struct Particle{
    int id;
    int charge;
    int status;
    double pT;
    double p;
    double eta;
    double energy;
    double phi;
  };
  Particle test;

  // ROOT objects
  TH1F *WpT_hist = new TH1F("WpT","WpT", 100, -100, 100);
  TTree *Tree = new TTree("Tree","Tree");
  Tree->Branch("Test",&test,"id/I:charge/I:status/I:pT/D:p/D:eta/D:energy/D:phi/D");

  // Begin event loop. Generate event; skip if generation aborted.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;
    int iW = 0;

    // Find number of all final charged particles.
    for (int i = 0; i < pythia.event.size(); ++i){
      if (pythia.event[i].id() == 24 || pythia.event[i].id() == -24)
        iW = i;
    }
    // Fill charged multiplicity in histogram. End event loop.
    // WpT_hist->Fill(pythia.event[iW].pT());
    double px = pythia.event[iW].px();
    double py = pythia.event[iW].py();
    double pz = pythia.event[iW].pz();
    test.p = px*px + py*py + pz*pz;
    test.id = pythia.event[iW].id();
    test.status = pythia.event[iW].status();
    test.charge = pythia.event[iW].charge();
    test.pT = pythia.event[iW].pT();
    test.eta = pythia.event[iW].eta();
    test.energy = pythia.event[iW].e();
    test.phi = pythia.event[iW].phi();

    Tree->Fill();
    pTZ.fill(pythia.event[iW].pT());
  }

  // Statistics on event generation.
  pythia.stat();
  cout << pTZ;

  // Show histogram. Possibility to close it.
  // WpT_hist->Draw();
  // std::cout << "\nDouble click on the histogram window to quit.\n";
  // gPad->WaitPrimitive();

  // Save histogram on file and close file.
  Tree->Write();
  delete outFile;

  // Done.
  return 0;
}