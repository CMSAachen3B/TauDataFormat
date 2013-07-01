import FWCore.ParameterSet.Config as cms

PreselectionCuts = cms.EDFilter('SkimmingCuts',
                                hpsTauProducer=cms.InputTag("hpsPFTauProducer"),
                                kinematicTausAdvanced = cms.InputTag("KinematicTauProducer","KinematicFitTau"),                                
                                MuonPtCut = cms.double(18.0),
                                MuonIsGlobal = cms.bool(True),
                                ElectronPtCut = cms.double(20.0),
                                ElectronEtaCut = cms.double(2.1),
                                NMuons      = cms.double(1.0),
                                MuonEtaCut = cms.double(2.1),
                                PFTauPtCut =cms.double(18.0),
                                PFTauEtaCut=cms.double(2.1),
                                doMuonOnly=cms.bool(False)
)

MuonPreselectionCuts = cms.EDFilter('SkimmingCuts',
                                    hpsTauProducer=cms.InputTag("hpsPFTauProducer"),
                                    kinematicTausAdvanced = cms.InputTag("KinematicTauProducer","KinematicFitTau"),
                                    MuonPtCut = cms.double(18.0),
                                    MuonIsGlobal = cms.bool(True),
                                    ElectronPtCut = cms.double(20.0),
                                    ElectronEtaCut = cms.double(2.1),
                                    NMuons      = cms.double(1.0),
                                    MuonEtaCut = cms.double(2.1),
                                    PFTauPtCut =cms.double(18.0),
                                    PFTauEtaCut=cms.double(2.1),
                                    doMuonOnly=cms.bool(True)
                                    )



ControlSample_Skim = cms.EDFilter('ControlSample_SkimCuts',
                                  muonsTag   = cms.InputTag("muons"),
                                  pfjetsTag  = cms.InputTag("ak5PFJets"),
                                  MuonPtCut  = cms.double(24.0),
                                  MuonEtaCut = cms.double(2.1),
                                  JetEtaCut  = cms.double(2.4),
                                  JetPtCut  = cms.double(10),
                                  dphiMuJet  = cms.double(2.5),
                                  )

    
