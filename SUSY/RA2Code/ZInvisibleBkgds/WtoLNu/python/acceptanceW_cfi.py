import FWCore.ParameterSet.Config as cms

acceptanceW = cms.EDProducer("AcceptanceW",
    debug    = cms.bool(True),
    debugString    = cms.string("acceptanceW"),
    genLabel = cms.InputTag("genParticles"),
    genMuonLabel     = cms.InputTag("genParticles"),
    genElectronLabel = cms.InputTag("genParticles"),
    genStatus = cms.int32(3),

    jetLabel = cms.InputTag("patJetsAK5PFPt50Eta25"),

    muonLabel  = cms.InputTag("patMuonsPF"),
    muonMinPt  = cms.double(20.),
    muonMaxEta = cms.double(2.1),

    electronLabel  = cms.InputTag("patElectronsPF"),
    electronMinPt  = cms.double(20.),
    electronMaxEta = cms.double(2.4)
)
