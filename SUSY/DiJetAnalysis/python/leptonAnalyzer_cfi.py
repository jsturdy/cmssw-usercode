import FWCore.ParameterSet.Config as cms


leptonAnalyzer = cms.untracked.PSet(
    elecMaxEta = cms.untracked.double(2.4),
    elecMaxEt  = cms.untracked.double(15.),
    elecMinEt  = cms.untracked.double(5.),
    elecRelIso = cms.untracked.double(0.5),

    muonMaxEta = cms.untracked.double(2.),
    muonMaxEt  = cms.untracked.double(10.),
    muonMinEt  = cms.untracked.double(5.),
    muonRelIso = cms.untracked.double(0.1),
    
    #tauMaxEta = cms.untracked.double(2.),
    #tauMaxEt  = cms.untracked.double(10.),
    #tauMinEt  = cms.untracked.double(5.),
    #tauRelIso = cms.untracked.double(0.1),
    
    doMCLeps  = cms.untracked.bool(True),
    debugLeps = cms.untracked.int32(0),
    genLepTag = cms.untracked.InputTag("genParticles"),
    
    elecTag   = cms.untracked.InputTag("cleanPatElectrons"),
    pfelecTag = cms.untracked.InputTag("particleFlow_Electrons"),
    muonTag   = cms.untracked.InputTag("cleanPatMuons"),
    #do not work!
    pfmuonTag = cms.untracked.InputTag("particleFlow_Muons")
    #tauTag   = cms.untracked.InputTag("cleanPatTaus"),
    #do not work!
    #pftauTag = cms.untracked.InputTag("pfPatTaus")
)
                 