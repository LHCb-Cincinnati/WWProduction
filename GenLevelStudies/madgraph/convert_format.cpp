#include <array>
#include <vector>

int convert_format() {
    // Open File
    auto ifile = TFile::Open("/data/home/ganowak/MG5_aMC_v2_9_3/ZZ_LO/Events/run_02/unweighted_events.root","READ");

    // Get Input Tree
    TTree *itree = (TTree*)ifile->Get("tree");

    // Define Variables
    const int est_pid_length = 12;
    Int_t mother_array[est_pid_length], pid_array[est_pid_length];
    Double_t px_array[est_pid_length], py_array[est_pid_length], pz_array[est_pid_length], e_array[est_pid_length], pt_array[est_pid_length];
    int ilplus, ilminus;
    std::array<Double_t, 5> lplus_array;
    std::array<Double_t, 5> lminus_array;
    double lminus_pT, lplus_pT;

    // Define Output File
    auto ofile = TFile::Open("madgraph/ZZ_MG5_LO_CMTS09_mu10.root","RECREATE");
    TTree *otree = new TTree("Tree","Tree");
    otree->Branch("TargetLepton", &lminus_array,"px/D:py/D:pz/D:e/D:pid/D");
    otree->Branch("TargetAntiLepton", &lplus_array,"px/D:py/D:pz/D:e/D:pid/D");

    // Set Branch Addresses
    itree->SetBranchAddress("Particle.PID", &pid_array);
    itree->SetBranchAddress("Particle.Mother1", &mother_array);
    itree->SetBranchAddress("Particle.Px", &px_array);
    itree->SetBranchAddress("Particle.Py", &py_array);
    itree->SetBranchAddress("Particle.Pz", &pz_array);
    itree->SetBranchAddress("Particle.E", &e_array);
    itree->SetBranchAddress("Particle.PT", &pt_array);

    // Event Loop
    for (int i = 0; i<itree->GetEntries(); i++){
    // for (int i = 0; i<10; i++){
        itree->GetEntry(i);
        lminus_pT = 0;
        lplus_pT= 0;
        for (int j = 0; j<est_pid_length; j++){
            if (pid_array[j] == 11 || pid_array[j] == 13) {
                if (pt_array[j] > lminus_pT){
                    lminus_pT = pt_array[j];
                    lminus_array = {px_array[j], py_array[j], pz_array[j], e_array[j], (Double_t) pid_array[j]};
                }
            }
            else if (pid_array[j] == -11 || pid_array[j] == -13) {
                if (pt_array[j] > lplus_pT){
                    lplus_pT = pt_array[j];
                    lplus_array = {px_array[j], py_array[j], pz_array[j], e_array[j], (Double_t) pid_array[j]};
                }
            }
        }
        otree->Fill();
    }
    ofile->Write();
    ifile->Close();
    ofile->Close();
    return(0);
//  && pid_array[mother_array[j]] == 24
}
