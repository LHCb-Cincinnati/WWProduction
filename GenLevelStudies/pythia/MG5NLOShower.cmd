! MG5NLOShower.cmd
! Shower MG5_aMC@NLO LHE events with Pythia8.
! The Python driver overwrites Beams:LHEF and may overwrite
! Main:numberOfEvents from command-line arguments.

! General settings.
Main:numberOfEvents          = 10000
Main:timesAllowErrors        = 10
Init:showChangedSettings     = on
Init:showChangedParticleData = on
Next:numberShowInfo          = 1
Next:numberShowProcess       = 1
Next:numberShowEvent         = 1

! LHE input.
Beams:frameType              = 4
Beams:LHEF                   = dummy.lhe

! FxFx defaults for MG5_aMC@NLO samples with ickkw = 3.
! Check qCut, qCutME, and nJetMax against the MG5 run card before production.
JetMatching:merge            = on
JetMatching:scheme           = 1
JetMatching:setMad           = off
JetMatching:qCut             = 20.0
JetMatching:coneRadius       = 1.0
JetMatching:etaJetMax        = 10.0
SpaceShower:MEcorrections    = off
SpaceShower:pTmaxMatch       = 1
SpaceShower:pTmaxFudge       = 1
TimeShower:pTmaxMatch        = 1
TimeShower:pTmaxFudge        = 1
TimeShower:MEcorrections     = off
TimeShower:globalRecoil      = on
TimeShower:weightGluonToQuark = 1
TimeShower:limitPTmaxGlobal  = on
TimeShower:nMaxGlobalRecoil  = 1
TimeShower:globalRecoilMode  = on
TimeShower:nMaxGlobalBranch  = 1

! Be more forgiving with finite-precision LHE momentum mismatches.
Check:epTolErr               = 1e-2

! Single-subrun FxFx setup for an LHE file containing the merged sample.
LHEFInputs:nSubruns          = 1
Main:subrun                  = 0
JetMatching:doFxFx           = on
JetMatching:qCutME           = 10.0
JetMatching:nJetMax          = 1
