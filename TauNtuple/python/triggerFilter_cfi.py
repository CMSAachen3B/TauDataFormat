import FWCore.ParameterSet.Config as cms

MultiTrigFilter = cms.EDFilter("MultiTriggerFilter",
                               TriggerResults = cms.InputTag("TriggerResults","","HLT"),
#                               useTriggers = cms.vstring("HLT_QuadJet40_IsoPFTau40","HLT_QuadJet40_IsoPFTau40","HLT_Mu20","HLT_Mu24","HLT_Mu30","HLT_IsoMu12_LooseIsoPFTau10","HLT_Mu15_LooseIsoPFTau20")
                               useTriggers = cms.vstring("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v","HLT_IsoMu18_eta2p1_LooseIsoPFTau20_v","HLT_IsoMu24_eta2p1_v")
                               )
