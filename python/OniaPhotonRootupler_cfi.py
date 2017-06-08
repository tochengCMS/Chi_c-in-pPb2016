import FWCore.ParameterSet.Config as cms

rootuple = cms.EDAnalyzer('OniaPhotonRootupler',
                          centralitySrc = cms.InputTag("pACentrality"),
                          centralityBinSrc = cms.InputTag("centralityBin","HFtowersPlusTrunc"),
                          chi_cand = cms.InputTag("ChiProducer","ChiCandidates"),
                          ups_cand = cms.InputTag("onia2MuMuPAT"),
                          refit1S  = cms.InputTag("ChiKinFitter","ChiCandidatesFit"),
                          refit2S  = cms.InputTag("ChiKinFitter2S","ChiCandidatesFit2S"),
                          primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                          TriggerResults  = cms.InputTag("TriggerResults", "", "HLT"),
                          isMC = cms.bool(True)
                          )
