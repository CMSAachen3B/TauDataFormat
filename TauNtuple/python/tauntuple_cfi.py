import FWCore.ParameterSet.Config as cms

Scale = '1.0'       # The "Scale" argument should be set to 1.0 if you want to reweight to exactly the 
                   #input data distribution. If you are evaluating systematic errors, you can use this factor
                   #to scale the distribution before the weights are calculated (i.e., correctly) by putting in a non-unity argument. 


NtupleMaker = cms.EDProducer('TauNtuple',
                             primVtx = cms.InputTag("offlinePrimaryVertices"),
                             hpsTauProducer=cms.InputTag("hpsPFTauProducer"),
                             hpsPFTauDiscriminationByTightIsolation =cms.InputTag("hpsPFTauDiscriminationByTightIsolation"), 
                             hpsPFTauDiscriminationByMediumIsolation=cms.InputTag("hpsPFTauDiscriminationByMediumIsolation"),
                             hpsPFTauDiscriminationByLooseIsolation =cms.InputTag("hpsPFTauDiscriminationByLooseIsolation"),
                             pfMet = cms.InputTag('pfMet'),
                             muons   = cms.InputTag('muons'),
                             kinematicTaus  = cms.InputTag("KinematicTauBasicProducer"),
                             kinematicTausAdvanced = cms.InputTag("KinematicTauProducer","SelectedKinematicDecays"),
                             tauPrimaryVtx  =cms.InputTag("ThreeProngInputSelector","primVtx"),
                             pfjets         = cms.InputTag("ak5PFJets"),
                             generalTracks  = cms.InputTag("generalTracks"),
                             gensrc         = cms.InputTag('genParticles'),
                             GenEventInfo   = cms.InputTag('generator'),
                             discriminators = cms.vstring("PFRecoTauDiscriminationByKinematicFit","PFRecoTauDiscriminationByKinematicFitQuality"),
                             do_MCComplete  = cms.untracked.bool(False),
                             do_MCSummary   = cms.untracked.bool(True),
                             ScaleFactor    = cms.untracked.string(Scale),
                             PUInputHistoMC    = cms.untracked.string("MC_FLAT_PLUS_TAIL_PU"),               
                             PUInputHistoData  = cms.untracked.string("h_160404_180252_all"),
                             #PUInputFile = cms.untracked.string("src/data/Lumi_160404_180252_andMC_Flat_Tail.root")  # if run on the GRID
                             PUInputFile = cms.untracked.string("Lumi_160404_180252_andMC_Flat_Tail.root")            # if run on the local PC
                             
)                                                                   
 
