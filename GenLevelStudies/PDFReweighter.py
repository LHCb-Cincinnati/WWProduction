import ROOT
import numpy as np
import os
import lhapdf

# Set lhapdf path
os.environ["LHAPDF_DATA_PATH"] = "/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current/"

class PDFReweighter:
    def __init__(self, pdf_set_from="NNPDF30_lo_as_0130", pdf_set_to="CT18", pdf_set_to_num=0):
        """Initialize the PDF sets for reweighting."""
        self.pdf_from = lhapdf.getPDFSet(pdf_set_from)
        self.pdf_from_member = self.pdf_from.mkPDF(0)  # Central member of the source PDF set
        self.pdf_set_to_num = pdf_set_to_num 
        self.pdf_to = lhapdf.getPDFSet(pdf_set_to)
        self.pdf_to_member = self.pdf_to.mkPDF(self.pdf_set_to_num)  # Central member of the target PDF set

    def calculate_reweighting_factor(self, x1, x2, id1, id2, sqrt_s=13e3):
        """
        Calculate the reweighting factor between two PDF sets.

        Args:
            x1: Momentum fraction of parton 1.
            x2: Momentum fraction of parton 2.
            id1: PDG ID of parton 1.
            id2: PDG ID of parton 2.
            sqrt_s: Center-of-mass energy in GeV (default: 13 TeV).

        Returns:
            Reweighting factor.
        """
        q2 = x1 * x2 * (sqrt_s ** 2)  # Compute q^2

        pdf_from_1 = self.pdf_from_member.xfxQ2(id1, x1, q2)
        pdf_from_2 = self.pdf_from_member.xfxQ2(id2, x2, q2)

        pdf_to_1 = self.pdf_to_member.xfxQ2(id1, x1, q2)
        pdf_to_2 = self.pdf_to_member.xfxQ2(id2, x2, q2)

        # Ensure PDFs are not zero to avoid division errors
        if pdf_from_1 * pdf_from_2 == 0:
            return 0

        # Reweighting factor
        reweighting_factor = (pdf_to_1 * pdf_to_2) / (pdf_from_1 * pdf_from_2)
        return reweighting_factor

    def add_reweight_to_root_file(self, input_file, output_file):
        """
        Adds reweighting factors as a branch in a ROOT file.

        Args:
            input_file: Path to the input ROOT file.
            output_file: Path to save the updated ROOT file.
        """
        # Open ROOT file
        infile = ROOT.TFile.Open(input_file, "READ")
        tree = infile.Get("Events")

        # Create a new ROOT file to store updated events
        outfile = ROOT.TFile(output_file, "RECREATE")
        new_tree = tree.CloneTree(0)

        # Add a branch for the reweighting factor
        reweight_factor = np.zeros(1, dtype=float)
        reweight_branch = new_tree.Branch("reweight_factor", reweight_factor, "reweight_factor/D")

        # Loop over events and calculate reweighting factors
        for event in tree:
            reweight_factor[0] = self.calculate_reweighting_factor(
                x1=event.x1, x2=event.x2, id1=event.id1, id2=event.id2, sqrt_s=13e3
            )
            new_tree.Fill()

        # Write the new tree to the output file
        outfile.Write()
        outfile.Close()
        infile.Close()

# Example usage:
# Initialize the reweighter with specific PDF sets
# reweighter = PDFReweighter(pdf_set_from="NNPDF30_lo_as_0130", pdf_set_to="CT18")

# Example: Calculate reweighting factor for specific kinematics
# x1, x2 = 0.1, 0.2
# id1, id2 = 1, -1  # PDG IDs for quark and antiquark
# sqrt_s = 13e3  # Center-of-mass energy in GeV
# reweight_factor = reweighter.calculate_reweighting_factor(x1, x2, id1, id2, sqrt_s)
# print(f"Reweighting factor: {reweight_factor}")

# Example: Add reweighting factors to a ROOT file
# reweighter.add_reweight_to_root_file("input.root", "output_with_weights.root")
