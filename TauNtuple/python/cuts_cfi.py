import FWCore.ParameterSet.Config as cms

#PreselectionCuts = cms.EDFilter('SkimmingCuts',doMuonOnly=cms.bool(False))
#MuonPreselectionCuts = cms.EDFilter('SkimmingCuts',doMuonOnly=cms.bool(True))
NoPreselectionCuts = cms.EDFilter('SkimmingCuts',preselection=cms.untracked.string(""))
MuonPreselectionCuts = cms.EDFilter('SkimmingCuts',preselection=cms.untracked.string("Mu"))
ElePreselectionCuts = cms.EDFilter('SkimmingCuts',preselection=cms.untracked.string("Ele"))
DoubleMuPreselectionCuts = cms.EDFilter('SkimmingCuts',preselection=cms.untracked.string("DoubleMu"))
DoubleElePreselectionCuts = cms.EDFilter('SkimmingCuts',preselection=cms.untracked.string("DoubleEle"))
MuOrElePreselectionCuts = cms.EDFilter('SkimmingCuts',preselection=cms.untracked.string("MuOrElePre"))
EMuTvariablePreselectionCuts = cms.EDFilter('SkimmingCuts',preselection=cms.untracked.string("EMuTvariable"))
MuTauPreselectionCuts = cms.EDFilter('SkimmingCuts',preselection=cms.untracked.string("MuTau"))
