import FWCore.ParameterSet.Config as cms

NSJAnalyzer = cms.EDAnalyzer('NSJAnalyzer'
                             )

#FastJetAnalyzer = cms.EDAnalyzer('FastJetAnalyzer',
#                                 momPDGID = cms.uint32(36),
#                                 genMuTauPTMin = cms.double(0.0), #GeV
#                                 genMuPTMin = cms.double(20.0) #GeV
#                                 )


#FastJetAnalyzerBkg = cms.EDAnalyzer('FastJetAnalyzerBkg',
#                                    outFileName = cms.string('/data1/friccita/NSJdatasets_10242012/ak5_groomed/FastJetAnalysis_WJets_bkg_pileup_pruned.root'),
#                                    nFilt = cms.int32(3), #orig 3
#                                    rFilt = cms.double(0.3), #orig 0.3
#                                    trimPtFracMin = cms.double(0.03), #orig 0.03
#                                    zCut = cms.double(0.1), #orig 0.1
#                                    RcutFactor = cms.double(0.5) #orig 0.5
#                                 )


FastJetAnalyzer_alltaus = cms.EDAnalyzer('FastJetAnalyzer_alltaus',
                                         Sample = cms.string('WJets'),
                                         jetSrc = cms.InputTag("ak5PFJets"),
                                         outFileName = cms.string('/data1/friccita/NSJdatasets_10242012/ak5_groomed/FastJetAnalysis_signalWH_alltaus_pileup_pruned.root'),
                                         outputTextFile = cms.string('JetConstituentInfo_signal.txt'),
                                         nFilt = cms.int32(3), #orig 3
                                         rFilt = cms.double(0.3), #orig 0.3
                                         trimPtFracMin = cms.double(0.03), #orig 0.03
                                         zCut = cms.double(0.1), #orig 0.1
                                         RcutFactor = cms.double(0.5) #orig 0.5
	        )
