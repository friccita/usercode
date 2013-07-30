import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'file:/data1/friccita/NMSSM_ZH_skim_PFTausNoMu.root'
    #'file:/data1/yohay/NMSSMHiggs_ZH_skim_PFTausNoMu.root'
    ),
##                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
    
)

#process.load("BoostedTauAnalysis.TauAnalyzer.HPSPFTauProducer_cfi")
process.load("BoostedTauAnalysis.TauAnalyzer.tauanalyzer_cfi")
process.TauAnalyzer.tauTag = cms.InputTag("hpsPFTauProducer", "", "OWNPARTICLES")
process.TauAnalyzer.outFileName = cms.string('/afs/cern.ch/user/f/friccita/CMSSW_5_2_5/src/BoostedTauAnalysis/TauAnalyzer/test/EffPlots_NMSSMHiggs_ZHaa_LooseIsoCombined_cleanPFJets.root')

process.TauAnalyzer.momPDGID = cms.int32(23)
process.TauAnalyzer.genMuPTMin = cms.double(20.0)
process.TauAnalyzer.effVsEtaPTMin = cms.double(15.0)

process.p = cms.Path(process.TauAnalyzer)
#process.p = cms.Path(process.HPSPFTauProducer*process.TauAnalyzer)
