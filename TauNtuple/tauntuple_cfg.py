import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")#https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions#Global_Tags_for_Monte_Carlo_Prod
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr = cms.untracked.PSet(
    threshold = cms.untracked.string('INFO'),
	FwkReport = cms.untracked.PSet(limit = cms.untracked.int32(0)),
	DEBUG = cms.untracked.PSet(limit = cms.untracked.int32(-1))
    )

process.GlobalTag.globaltag = 'START42_V11::All'
#process.GlobalTag.globaltag = 'GR_P_V22::All'

numberOfEvents = 100



inPath1 = '/user/cherepanov/data_CMSSW_4_2_0/DYtoTauTau/TauSkim/dea7f64cc9262673a0c76075e9ca9700/'  #Ztau tau      number oif files = 454  Version 3 DY
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
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#    listOfFiles_data)                 
#    listOfFiles)
    'file:///user/inugent/DYtoTauTau_TauQualitySkim.root')
#    'file:///user/inugent/TauXPromptReco172_798.root')
#    'file:/afs/cern.ch/cms/Tutorials/TWIKI_DATA/RelValTTbar_RECO_424.root')
    
)

process.myProducerLabel = cms.EDProducer('TauNtuple')

#process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string('myOutputFile.root')
#)

############################ KinematicFit 3pi decay mode; 4pi come soon  ###############################
process.load("CommonTools.PrimVtxSelector.PrimVtxSelector_cfi")
process.load("RecoTauTag.KinematicTau.InputTrackSelector_cfi")
process.load("RecoTauTag.KinematicTau.ThreeProngInputSelector_cfi")
process.load("RecoTauTag.KinematicTau.kinematictau_cfi")
process.load("RecoTauTag.KinematicTau.kinematictauAdvanced_cfi")
############################ KinematicFit 3pi decay mode; 4pi come soon  ###############################
process.load("TauDataFormat.TauNtuple.tauntuple_cfi")

#EvntCounter = cms.EDAnalyser('EventCounter',
#                             CounterType    = cms.untracked.string("WW"),
#                             gensrc         = cms.InputTag('genParticles'),
#                             GenEventInfo   = cms.InputTag('generator'),
#                             DataMCType     = cms.untracked.string("DataMCType")
#                             )



process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(numberOfEvents) )

process.p = cms.Path(process.PrimVtxSelector*process.InputTrackSelector*process.ThreeProngInputSelector*process.KinematicTauBasicProducer*process.DetailedProducer*process.NtupleMaker)

#process.p = cms.Path(process.NtupleMaker)

