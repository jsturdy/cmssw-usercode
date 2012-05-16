import FWCore.ParameterSet.Config as cms

wtolnu_hists = cms.EDProducer('WtoLNu',
    debug    = cms.bool(False),
    #puWeightLabel = cms.InputTag("patJetsAK5PFPt50Eta25"),
    vertexLabel = cms.InputTag("goodVertices"),
    jetLabel = cms.InputTag("patJetsAK5PFPt30"),
    htLabel = cms.InputTag("htPFNoElecIDIso"),
    mhtLabel = cms.InputTag("mhtPFNoElecIDIso"),
    metLabel = cms.InputTag("patMETsPF"),
    #metLabel = cms.InputTag("pfMet"),

    muonLabel  = cms.InputTag("patMuonsPFIDIso"),
    electronLabel  = cms.InputTag("patElectronsPFIDIso"),
)
