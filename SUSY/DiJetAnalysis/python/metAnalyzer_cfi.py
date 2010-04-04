import FWCore.ParameterSet.Config as cms


metAnalyzer = cms.untracked.PSet(
    doMCMET   = cms.untracked.bool(False),
    doCaloMET = cms.untracked.bool(True),
    doPfMET   = cms.untracked.bool(True),
    doTcMET   = cms.untracked.bool(True),
    
    genMETTag = cms.untracked.InputTag("genMetCalo"),
    metTag    = cms.untracked.InputTag("patMETsAK5"),
    pfmetTag  = cms.untracked.InputTag("patMETsPF"),
    tcmetTag  = cms.untracked.InputTag("patMETsTC"),

    debugMET = cms.untracked.int32(0)
    )
