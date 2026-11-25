! @file MG5LOShower.cmd
! @brief Shower a single unmatched MG5 LO LHE sample with Pythia 8.
!
! Compile the linked executables from the repository root:
!   conda activate hep
!   make -C GenLevelStudies/pythia -f Makefile.pythia_gen_lhe all
! Run, for example:
!   GenLevelStudies/pythia/pythia_gen_lhe_wz_lo input.lhe output.root --seed 1

! General settings.
Main:numberOfEvents          = 1000000
Main:timesAllowErrors        = 10
Init:showChangedSettings     = on
Init:showChangedParticleData = on
Next:numberShowInfo          = 1
Next:numberShowProcess       = 1
Next:numberShowEvent         = 1

! LHE input; the executable replaces the placeholder path.
Beams:frameType              = 4
Beams:LHEF                   = dummy.lhe

! These samples contain one unmatched LO multiplicity. Keep the Pythia
! JetMatching integration installed, but disable its acceptance veto.
JetMatching:merge            = off
JetMatching:setMad           = off
Check:epTolErr               = 1e-2
