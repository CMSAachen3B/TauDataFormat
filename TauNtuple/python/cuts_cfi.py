import FWCore.ParameterSet.Config as cms

PreselectionCuts = cms.EDFilter('SkimmingCuts',doMuonOnly=cms.bool(False))
MuonPreselectionCuts = cms.EDFilter('SkimmingCuts',doMuonOnly=cms.bool(True))
