import FWCore.ParameterSet.Config as cms

EvntCounterA = cms.EDAnalyzer('EventCounter',
                              CounterType    = cms.untracked.string("AllEvents"),
                              gensrc         = cms.InputTag('genParticles'),
                              GenEventInfo   = cms.InputTag('generator'),
                              DataMCType     = cms.untracked.string("dy_tautau")
                              )


EvntCounterB = cms.EDAnalyzer('EventCounter',
                              CounterType    = cms.untracked.string("BeforeTauNtuple"),
                              gensrc         = cms.InputTag('genParticles'),
                              GenEventInfo   = cms.InputTag('generator'),
                              DataMCType     = cms.untracked.string("dy_tautau")
                              )
