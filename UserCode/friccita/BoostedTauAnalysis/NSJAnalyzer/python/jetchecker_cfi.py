import FWCore.ParameterSet.Config as cms

JetChecker = cms.EDAnalyzer('JetChecker',
                            outFileName = cms.string('/data1/friccita/JetChecker.root'),
                            JetCollectionBefore = cms.InputTag("ak5PFJets", "", "HLT"),
                            JetCollectionAfter = cms.InputTag("pfJets", "", "PF2PAT")
                            )
