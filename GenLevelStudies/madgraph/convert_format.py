import  numpy as np
import pdb
import sys
import ROOT
# Personal Packages 
sys.path.append(".") # Not great form.
import AnalysisTools as at

# Parse Inputs
parser = at.Parser(sys.argv[1:])
args = parser.args
file_name =  args.input_files[0].name
ofile_name = args.output

# Open Input File
ifile = ROOT.TFile.Open(file_name, "READ")
# itree = ifile.Get("tree")
itree= ifile.Get("tree")

# Create output arrays
lminus_array = np.array([0]*5, dtype=np.float64)
lplus_array = np.array([0]*5, dtype=np.float64)

# Create Output File
ofile = ROOT.TFile.Open(ofile_name, "RECREATE")
otree = ROOT.TTree("Tree", "Tree")
otree.Branch("TargetLepton", lminus_array, "px/D:py/D:pz/D:e/D:pid/D")
otree.Branch("TargetAntiLepton", lplus_array, "px/D:py/D:pz/D:e/D:pid/D")

# ROOT Config
PxPyPzEVector = ROOT.Math.LorentzVector(ROOT.Math.PxPyPzE4D('double'))

for entry in itree:
    lminus_pT = 0
    lplus_pT = 0
    pid_array = getattr(itree, "Particle.PID")
    for index, pid in enumerate(pid_array):
        particle_vec = PxPyPzEVector(
            getattr(itree, "Particle.Px")[index],
            getattr(itree, "Particle.Py")[index],
            getattr(itree, "Particle.Pz")[index],
            getattr(itree, "Particle.E")[index]
        )
        # pdb.set_trace()
        # print(pid, particle_vec.Eta())
        if (((pid == 13) | (pid == 11)) & (particle_vec.Eta()<5) & (particle_vec.Eta()>2)):
            if (particle_vec.Pt() > lminus_pT):
                lminus_pT = particle_vec.Pt()
                lminus_array[:] = ([
                    np.double(getattr(itree, "Particle.Px")[index]),
                    np.double(getattr(itree, "Particle.Py")[index]),
                    np.double(getattr(itree, "Particle.Pz")[index]),
                    np.double(getattr(itree, "Particle.E")[index]),
                    np.double(pid)
                ])
        elif (((pid == -13) | (pid == -11)) & (particle_vec.Eta()<5) & (particle_vec.Eta()>2)):
            if (particle_vec.Pt() > lplus_pT):
                lplus_pT = particle_vec.Pt()
                lplus_array[:] = ([
                    np.double(getattr(itree, "Particle.Px")[index]),
                    np.double(getattr(itree, "Particle.Py")[index]),
                    np.double(getattr(itree, "Particle.Pz")[index]),
                    np.double(getattr(itree, "Particle.E")[index]),
                    np.double(pid)
                ])
    if ((lminus_pT == 0) | (lplus_pT == 0)):
        continue
    else: 
        otree.Fill()

otree.Print()
ifile.Close()
ofile.Write()
ofile.Close()
