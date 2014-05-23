import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")#https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions#Global_Tags_for_Monte_Carlo_Prod
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr = cms.untracked.PSet(
    threshold = cms.untracked.string('INFO'),
    FwkReport = cms.untracked.PSet(limit = cms.untracked.int32(0)),
    DEBUG = cms.untracked.PSet(limit = cms.untracked.int32(1))
    )


process.load('Configuration/StandardSequences/Services_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
#process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

# global tag
#process.GlobalTag.globaltag = 'START44_V9B::All'
#process.GlobalTag.globaltag = 'START52_V9::All'
#process.GlobalTag.globaltag = 'START53_V7A::All'
#process.GlobalTag.globaltag = 'FT_R_53_V6C::All'


process.GlobalTag.globaltag = 'FT_R_53_V18::All'


process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load('HLTrigger.Configuration.HLT_GRun_cff')

#process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")


process.load("Configuration.StandardSequences.EndOfProcess_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")

######################################################

process.MessageLogger.categories.append('KinematicTau')
process.MessageLogger.categories.append('ThreeProngInputSelector_Step2')
process.MessageLogger.debugModules = cms.untracked.vstring('KinematicTau','ThreeProngInputSelector_Step2')
process.MessageLogger.cerr = cms.untracked.PSet(
     threshold = cms.untracked.string('INFO'),
     FwkReport = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
     DEBUG = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
     KinematicTau = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
     ThreeProngInputSelector_Step2 = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
)



numberOfEvents = 5000



#inPath1 = '/user/cherepanov/data_CMSSW_4_2_0/DYtoTauTau/TauSkim/dea7f64cc9262673a0c76075e9ca9700/'  #Ztau tau      number oif files = 454  Version 3 DY
#f1=open('/user/cherepanov/data_CMSSW_4_2_0/DYtoTauTau/TauSkim/flis_list')



#fdata=open('/.automount/home/home__home2/institut_3b/cherepanov/work/CMSSW_4_2_0/src/Ztautau/Muon3prong/DataFileList.dat')  # nFiles_data = 138

listOfFiles=[]
listOfFiles_data=[]


nFiles1 = 454  
nFiles_data=138

############# take set of unmerged files ###########
#for nf in range(1,nFiles1):
#    string=f1.readline()
#    listOfFiles.append('file://'+inPath1+string[:-1] )
############# take set of unmerged files ###########
############ take set of unmerged Data files ###########
#for nf in range(1,nFiles_data):
#    string_data=fdata.readline()
#    listOfFiles_data.append('file://'+string_data[:-1] )
############ take set of unmerged Data files ###########


process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/data/Run2012A/TauPlusX/AOD/13Jul2012-v1/0000/187AF3FD-BED9-E111-A974-00259073E3AC.root')
#                                'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/data/Run2012A/TauPlusX/AOD/13Jul2012-v1/0000/76C17B7B-98D7-E111-B500-20CF3027A639.root')
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/C09EAADE-3D97-E111-B32E-001E67396DF1.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/D2E4D132-3D97-E111-88B2-003048673FC0.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/32696749-3E97-E111-BB99-001A647894A0.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/02F49D49-3E97-E111-9AF6-003048670A0C.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/521DD948-9F96-E111-B544-003048D45FD8.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/EA57A1DE-CD96-E111-8D2F-0025B3E063A8.root',
    'file://dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0000/760F45FA-3D97-E111-A86A-002590200974.root')

                            
                            )



process.myProducerLabel = cms.EDProducer('TauNtuple')

#process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string('myOutputFile.root')
#)
process.load("TauDataFormat.TauNtuple.tauntuple_cfi")
#process.load("TauDataFormat.TauNtuple.eventCounter_cfi")

############################ KinematicFit 3pi decay mode; 4pi come soon  ###############################
#process.load("CommonTools.PrimVtxSelector.PrimVtxSelector_cfi")
#process.load("RecoTauTag.KinematicTau.InputTrackSelector_cfi")
#process.load("RecoTauTag.KinematicTau.ThreeProngInputSelector_cff")
#process.load("RecoTauTag.KinematicTau.kinematictau_cfi")
#process.load("RecoTauTag.KinematicTau.kinematictauAdvanced_cfi")
############################ KinematicFit 3pi decay mode; 4pi come soon  ###############################
process.load("TriggerFilter.Filter.triggerFilter_cfi")



process.EvntCounterA = cms.EDAnalyzer('EventCounter',
                                      CounterType    = cms.untracked.string("AllEvents"),
                                      gensrc         = cms.InputTag('genParticles'),
                                      GenEventInfo   = cms.InputTag('generator'),
                                      DataMCType     = cms.untracked.string("data")
                                      )


process.EvntCounterB = cms.EDAnalyzer('EventCounter',
                                      CounterType    = cms.untracked.string("BeforeTauNtuple"),
                                      gensrc         = cms.InputTag('genParticles'),
                                      GenEventInfo   = cms.InputTag('generator'),
                                      DataMCType     = cms.untracked.string("data")
                                      )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(numberOfEvents) )


process.load("RecoTauTag.KinematicTau.KinematicFitSequences_cff")
process.load("RecoTauTag.KinematicTau.KinematicTauPostProcessing_cfi")


process.load('Configuration.StandardSequences.EDMtoMEAtJobEnd_cff')
process.load("Validation.Configuration.postValidation_cff")

process.schedule = cms.Schedule()

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.EDMtoMEAtJobEnd_cff')
process.load('RecoTauTag.KinematicTau.KinematicTauPostProcessing_cfi')
process.load('RecoTauTag.KinematicTau.Tau_JAKID_Filter_cfi')
process.load('DQMServices.Components.MEtoEDMConverter_cfi')

process.load("TriggerFilter.Filter.triggerFilter_cfi")
process.load("HLTrigger.HLTfilters.triggerResultsFilter_cfi")
process.load("SkimmingTools.EventCounter.countInput_cfi")
process.load("SkimmingTools.EventCounter.countTriggerPassed_cfi")
process.load("SkimmingTools.EventCounter.countKinFitPassed_cfi")
process.load("TauDataFormat.TauNtuple.eventCounter_cfi")

process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
process.load("SkimmingTools.SkimmingCuts.cuts_cfi")
process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")
process.pfType1CorrectedMet.applyType0Corrections = cms.bool(False)
process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag(
    cms.InputTag('pfMETcorrType0'),
    cms.InputTag('pfJetMETcorr', 'type1')        
)

process.dqmSaver.convention = 'Offline'
process.dqmSaver.saveByRun = cms.untracked.int32(-1)
process.dqmSaver.saveAtJobEnd = cms.untracked.bool(True)
process.dqmSaver.forceRunNumber = cms.untracked.int32(1)
process.dqmSaver.workflow = "/KinematicFitSequencewithDQM/VAL/RECO"
#process.DQMStore.verbose=1

#process.endjob_step = cms.Path(process.endOfProcess)
#process.KinFitSkim  = cms.Path(process.PFTau*process.KinematicFitSequencewithDQM*process.MEtoEDMConverter*process.EDMtoME*process.dqmSaver)
#process.p  = cms.Path(process.TauJAKIDFilter*process.PFTau*process.EvntCounterA*process.TrigFilter*process.TrigFilterInfo*process.KinematicFitSequence*process.MEtoEDMConverter*process.EDMtoME*process.dqmSaver*process.EvntCounterB*process.NtupleMaker)



#process.p  = cms.Path(process.EvntCounterA*process.TrigFilter*process.TrigFilterInfo*process.KinematicFitSequence*process.EvntCounterB*process.recoTauClassicHPSSequence*process.NtupleMaker)

#process.p  = cms.Path(process.EvntCounterA*process.TrigFilter*process.TrigFilterInfo*process.KinematicFitSequencewithSkim*process.CountKinFitPassedEvents*process.recoTauClassicHPSSequence*process.NtupleMaker)
#process.p  = cms.Path(process.EvntCounterA*process.TrigFilter*process.TrigFilterInfo*process.KinematicFitSequencewithSkim*process.CountKinFitPassedEvents*process.EvntCounterB*process.recoTauClassicHPSSequence*process.NtupleMaker)

#process.p  = cms.Path(process.EvntCounterA*process.TrigFilter*process.TrigFilterInfo*process.CountTriggerPassedEvents*process.KinematicFitSequencewithSkim*process.CountKinFitPassedEvents*process.EvntCounterB*process.recoTauClassicHPSSequence*process.NtupleMaker)


#process.p  = cms.Path(process.TauJAKIDFilter*process.EvntCounterA*process.CountInputEvents*process.MultiTrigFilter*process.TrigFilterInfo*process.CountTriggerPassedEvents*process.KinematicFitSequence*process.CountKinFitPassedEvents*process.PreselectionCuts*process.EvntCounterB*process.recoTauClassicHPSSequence*process.type0PFMEtCorrection*process.producePFMETCorrections*process.NtupleMaker)

process.p  = cms.Path(process.EvntCounterA*process.CountInputEvents*process.MultiTrigFilter*process.TrigFilterInfo*process.CountTriggerPassedEvents*process.KinematicFitSequence*process.CountKinFitPassedEvents*process.EvntCounterB*process.recoTauClassicHPSSequence*process.PreselectionCuts*process.type0PFMEtCorrection*process.producePFMETCorrections*process.NtupleMaker)





#process.p  = cms.Path(process.EvntCounterA*process.CountInputEvents*process.TrigFilter*process.TrigFilterInfo*process.CountTriggerPassedEvents*process.EvntCounterB*process.recoTauClassicHPSSequence*process.type0PFMEtCorrection*process.producePFMETCorrections*process.NtupleMaker) # without fit

#process.p  = cms.Path(process.EvntCounterA*process.CountInputEvents*process.TrigFilter*process.TrigFilterInfo*process.KinematicFitSequencewithSkim*process.MEtoEDMConverter*process.EDMtoME*process.EvntCounterB*process.recoTauClassicHPSSequence*process.NtupleMaker)


#process.p  = cms.Path(process.EvntCounterA*process.TrigFilter*process.TrigFilterInfo*process.KinematicFitSequencewithSkim*process.EvntCounterB*process.recoTauClassicHPSSequence*process.NtupleMaker)

#process.p  = cms.Path(process.EvntCounterA*process.TrigFilter)#*process.TrigFilterInfo*process.KinematicFitSequencewithSkim*process.EvntCounterB*process.recoTauClassicHPSSequence*process.NtupleMaker)


#process.p  = cms.Path(process.EvntCounterA*process.CountInputEvents*process.TrigFilter*process.TrigFilterInfo*process.CountTriggerPassedEvents*process.KinematicFitSequencewithSkim*process.CountKinFitPassedEvents*process.MEtoEDMConverter*process.EDMtoME*process.EvntCounterB*process.recoTauClassicHPSSequence*process.NtupleMaker)


##process.KinFitSkim  = cms.Path(process.TauJAKIDFilter*process.PFTau*process.KinematicFitSequence*process.MEtoEDMConverter*process.EDMtoME*process.dqmSaver)
#process.KinFitSkim  = cms.Path(process.KinematicFitSequence*process.MEtoEDMConverter*process.EDMtoME*process.dqmSaver)
#process.KinFitSkim  = cms.Path(process.KinematicFitSequence*process.MEtoEDMConverter*process.EDMtoME*process.dqmSaver)
process.schedule = cms.Schedule(process.p)#,process.endjob_step)


#process.p = cms.Path(process.EvntCounterA*process.TrigFilter*process.TrigFilterInfo*process.PrimVtxSelector*process.InputTrackSelector*process.ThreeProngInputSelector*process.KinematicTauBasicProducer*process.KinematicTauProducer*process.EvntCounterB*process.NtupleMaker)

#process.p = cms.Path(process.NtupleMaker)
