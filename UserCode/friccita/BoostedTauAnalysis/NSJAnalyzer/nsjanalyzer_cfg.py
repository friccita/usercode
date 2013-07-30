import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("Nsubjettiness")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("BoostedTauAnalysis.NSJAnalyzer.nsjanalyzer_cfi")
process.load("BoostedTauAnalysis.WmunuFilter.WmunuFilter_cff")
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.GeometryExtended_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.GlobalTag.globaltag = "START52_V9::All"

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'file:/data1/yohay/NMSSMHiggs_gg_skim_1000Files.root'
#    'file:/data1/friccita/NMSSMHiggs_WH_CHS.root'
    'file:/data1/yohay/NMSSMHiggs_WH_files1-250_24Sep12.root',
    'file:/data1/yohay/NMSSMHiggs_WH_files251-500_24Sep12.root',
    'file:/data1/yohay/NMSSMHiggs_WH_files501-750_24Sep12.root',
    'file:/data1/yohay/NMSSMHiggs_WH_files751-1000_24Sep12.root'
#         'root://eoscms//eos/cms/store/user/friccita/WToMuNu_skim/Summer12_WJetsToLNu_WMuNuSkim_1000_3_TDu.root'
 )
)

process.FastJetAnalyzer_alltaus.Sample = cms.string("Wh1")
#process.FastJetAnalyzer_alltaus.jetSrc = cms.InputTag("pfJets", "", "PF2PAT")
process.FastJetAnalyzer_alltaus.outFileName = cms.string('/data1/friccita/NSJdatasets_10242012/ak5_groomed/FastJetAnalysis_signalWH_alltaus_pileup_pruned_z0p1_D0p5_latest.root')
process.FastJetAnalyzer_alltaus.outputTextFile = cms.string('JetConstituentInfo_sig.txt')

#mylist = FileUtils.loadListFromFile ('wjets.txt')
#readFiles = cms.untracked.vstring( *mylist)
#process.source = cms.Source("PoolSource",
#                            fileNames = readFiles
##                           skipEvents = cms.untracked.uint32(6763)
#                           )

process.p = cms.Path(process.WmunuFilter*process.FastJetAnalyzer_alltaus)
