import FWCore.ParameterSet.Config as cms

Scale = '1.0'       # The "Scale" argument should be set to 1.0 if you want to reweight to exactly the 
                   #input data distribution. If you are evaluating systematic errors, you can use this factor
                   #to scale the distribution before the weights are calculated (i.e., correctly) by putting in a non-unity argument. 


NtupleMaker = cms.EDProducer('TauNtuple',
                             primVtx = cms.InputTag("offlinePrimaryVertices"),
                             hpsTauProducer=cms.InputTag("hpsPFTauProducer"),
                             beamSpot = cms.InputTag("offlineBeamSpot"),
                             hpsPFTauDiscriminationByTightIsolation =cms.InputTag("hpsPFTauDiscriminationByTightIsolation"), 
                             hpsPFTauDiscriminationByMediumIsolation=cms.InputTag("hpsPFTauDiscriminationByMediumIsolation"),
                             hpsPFTauDiscriminationByLooseIsolation =cms.InputTag("hpsPFTauDiscriminationByLooseIsolation"),
                             hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr=cms.InputTag("hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr"), 
                             hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr =cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr"), 
                             hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr=cms.InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr"),
                             hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr =cms.InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr"),
                             hpsPFTauDiscriminationAgainstElectronLoose=cms.InputTag("hpsPFTauDiscriminationByLooseElectronRejection"),
                             hpsPFTauDiscriminationAgainstElectronMedium=cms.InputTag("hpsPFTauDiscriminationByMediumElectronRejection"),
                             hpsPFTauDiscriminationAgainstElectronTight=cms.InputTag("hpsPFTauDiscriminationByTightElectronRejection"),
                             hpsPFTauDiscriminationAgainstMuonLoose=cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection"),
                             hpsPFTauDiscriminationAgainstMuonMedium=cms.InputTag("hpsPFTauDiscriminationByMediumMuonRejection"),
                             hpsPFTauDiscriminationAgainstMuonTight=cms.InputTag("hpsPFTauDiscriminationByTightMuonRejection"),
                             hpsPFTauDiscriminationByDecayModeFinding =cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),
                             pfMet = cms.InputTag('pfMet'),
                             muons   = cms.InputTag('muons'),
#                             kinematicTaus  = cms.InputTag("KinematicTauBasicProducer"),
                             kinematicTausAdvanced = cms.InputTag("KinematicTauProducer","KinematicFitTau"),
#                             kinematicTausAdvanced = cms.InputTag("KinematicTauProducer","SelectedKinematicDecays"),
#                             tauPrimaryVtx  = cms.InputTag("ThreeProngInputSelectorStep2","primVtx"),
                             tauPrimaryVtx  = cms.InputTag("offlinePrimaryVertices"),
                             pfjets         = cms.InputTag("ak5PFJets"),
                             pfelectrons =  cms.InputTag("gsfElectrons"),#"pfElectrons"),#,"electrons","RECO"),
                             generalTracks  = cms.InputTag("generalTracks"),
                             gensrc         = cms.InputTag('genParticles'),
                             GenEventInfo   = cms.InputTag('generator'),
                             jetFlavour     = cms.InputTag('patJetFlavourAssociation'),
                             discriminators = cms.vstring("PFRecoTauDiscriminationByKinematicFit","PFRecoTauDiscriminationByKinematicFitQuality"),
                             do_MCComplete  = cms.untracked.bool(True),
                             do_MCSummary   = cms.untracked.bool(True),
                             ScaleFactor    = cms.untracked.string(Scale),
                             PUInputHistoMC    = cms.untracked.string("MC_Summer12"),
                             PUInputHistoData  = cms.untracked.string("h_190456_20868"),
                             PUOutputFile  = cms.untracked.string("Weight3D.root"),
                             PUInputFile = cms.untracked.string("$CMSSW_BASE/src/data/Lumi_190456_208686MC_PU_S10_andData.root"),
                             TriggerProcessName = cms.untracked.string("HLT"),
                             TriggerInfoName = cms.InputTag("TrigFilterInfo","TriggerFilterInfoList"),
                             TriggerEvent = cms.InputTag("hltTriggerSummaryAOD","","HLT"), 
                             TriggerResults = cms.InputTag("TriggerResults","","HLT"),
                             L1GtTriggerMenuLite = cms.InputTag("l1GtTriggerMenuLite"),
                             l1TriggerNames = cms.vstring("L1_DoubleTau44er","L1_DoubleJetC64"),
                             RhoIsolAllInputTag   =  cms.InputTag('kt6PFJets', 'rho', 'RECO'),
                             TriggerJetMatchingdr = cms.untracked.double(0.3),
                             TriggerMuonMatchingdr = cms.untracked.double(0.3),
                             TriggerElectronMatchingdr = cms.untracked.double(0.3),
                             TriggerTauMatchingdr = cms.untracked.double(0.3),
                             BTagAlgorithm = cms.untracked.string("combinedSecondaryVertexBJetTags")

                             )                                                                   
 

