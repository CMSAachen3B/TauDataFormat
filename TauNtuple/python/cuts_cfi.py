import FWCore.ParameterSet.Config as cms

#PreselectionCuts = cms.EDFilter('SkimmingCuts',doMuonOnly=cms.bool(False))
#MuonPreselectionCuts = cms.EDFilter('SkimmingCuts',doMuonOnly=cms.bool(True))
PreselectionCuts = cms.EDFilter('SkimmingCuts',preselection=cms.untracked.string(""))
MuonPreselectionCuts = cms.EDFilter('SkimmingCuts',preselection=cms.untracked.string("Mu"))
ElePreselectionCuts = cms.EDFilter('SkimmingCuts',preselection=cms.untracked.string("Ele"))
MuJetPreselectionCuts = cms.EDFilter('SkimmingCuts',preselection=cms.untracked.string("MuJet"))
DoubleMuPreselectionCuts = cms.EDFilter('SkimmingCuts',preselection=cms.untracked.string("DoubleMu"))
DoubleElePreselectionCuts = cms.EDFilter('SkimmingCuts',preselection=cms.untracked.string("DoubleEle"))
