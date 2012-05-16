
from PhysicsTools.PatAlgos.cleaningLayer1.jetCleaner_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.jetCountFilter_cfi import *

patJetsAK5PFPt30NoPFMuon = cleanPatJets.clone()
patJetsAK5PFPt30NoPFMuon.src = cms.InputTag('patJetsAK5PFPt30')
patJetsAK5PFPt30NoPFMuon.checkOverlaps.muons.src               = cms.InputTag('patMuonsPF')
patJetsAK5PFPt30NoPFMuon.checkOverlaps.muons.algorithm         = cms.string('byDeltaR')
patJetsAK5PFPt30NoPFMuon.checkOverlaps.muons.preselection      = cms.string('')
patJetsAK5PFPt30NoPFMuon.checkOverlaps.muons.deltaR            = cms.double(0.4)
patJetsAK5PFPt30NoPFMuon.checkOverlaps.muons.pairCut           = cms.string('')
patJetsAK5PFPt30NoPFMuon.checkOverlaps.muons.requireNoOverlaps = cms.bool(True)

patJetsAK5PFPt30NoPFMuonID                         = patJetsAK5PFPt30NoPFMuon.clone()
patJetsAK5PFPt30NoPFMuonID.checkOverlaps.muons.src = cms.InputTag('patMuonsPFID')

patJetsAK5PFPt30NoPFMuonIDIso                         = patJetsAK5PFPt30NoPFMuon.clone()
patJetsAK5PFPt30NoPFMuonIDIso.checkOverlaps.muons.src = cms.InputTag('patMuonsPFIDIso')

#####
patJetsAK5PFPt50Eta25NoPFMuon     = patJetsAK5PFPt30NoPFMuon.clone()
patJetsAK5PFPt50Eta25NoPFMuon.src = cms.InputTag('patJetsAK5PFPt50Eta25')

patJetsAK5PFPt50Eta25NoPFMuonID                         = patJetsAK5PFPt50Eta25NoPFMuon.clone()
patJetsAK5PFPt50Eta25NoPFMuonID.checkOverlaps.muons.src = cms.InputTag('patMuonsPFID')

patJetsAK5PFPt50Eta25NoPFMuonIDIso                         = patJetsAK5PFPt50Eta25NoPFMuon.clone()
patJetsAK5PFPt50Eta25NoPFMuonIDIso.checkOverlaps.muons.src = cms.InputTag('patMuonsPFIDIso')


muonCleanedPFJets = cms.Sequence(
    patJetsAK5PFPt30NoPFMuon           *
    patJetsAK5PFPt30NoPFMuonID         *
    patJetsAK5PFPt30NoPFMuonIDIso      *
    patJetsAK5PFPt50Eta25NoPFMuon      *
    patJetsAK5PFPt50Eta25NoPFMuonID    *
    patJetsAK5PFPt50Eta25NoPFMuonIDIso
)
# PFJets - filters

countJetsAK5PFPt50Eta25DiJets           = countPatJets.clone()
countJetsAK5PFPt50Eta25DiJets.src       = cms.InputTag('patJetsAK5PFPt50Eta25')
countJetsAK5PFPt50Eta25DiJets.minNumber = cms.uint32(2)

countJetsAK5PFPt50Eta25NoPFMuon           = countPatJets.clone()
countJetsAK5PFPt50Eta25NoPFMuon.src       = cms.InputTag('patJetsAK5PFPt50Eta25NoPFMuon')
countJetsAK5PFPt50Eta25NoPFMuon.minNumber = cms.uint32(3)

countJetsAK5PFPt50Eta25NoPFMuonID     = countJetsAK5PFPt50Eta25NoPFMuon.clone()
countJetsAK5PFPt50Eta25NoPFMuonID.src = cms.InputTag('patJetsAK5PFPt50Eta25NoPFMuonID')

countJetsAK5PFPt50Eta25NoPFMuonIDIso     = countJetsAK5PFPt50Eta25NoPFMuon.clone()
countJetsAK5PFPt50Eta25NoPFMuonIDIso.src = cms.InputTag('patJetsAK5PFPt50Eta25NoPFMuonIDIso')



### electron cleaned jets (same basic cleaning as against muons)
patJetsAK5PFPt30NoPFElectron = patJetsAK5PFPt30NoPFMuon.clone()
patJetsAK5PFPt30NoPFElectron.checkOverlaps.muons.src                 = cms.InputTag('patElectronsPF')

patJetsAK5PFPt30NoPFElectronID = patJetsAK5PFPt30NoPFElectron.clone()
patJetsAK5PFPt30NoPFElectronID.checkOverlaps.muons.src               = cms.InputTag('patElectronsPFID')

patJetsAK5PFPt30NoPFElectronIDIso = patJetsAK5PFPt30NoPFElectron.clone()
patJetsAK5PFPt30NoPFElectronIDIso.checkOverlaps.muons.src            = cms.InputTag('patElectronsPFIDIso')

#####
patJetsAK5PFPt50Eta25NoPFElectron = patJetsAK5PFPt30NoPFElectron.clone()
patJetsAK5PFPt50Eta25NoPFElectron.src = cms.InputTag('patJetsAK5PFPt50Eta25')

patJetsAK5PFPt50Eta25NoPFElectronID = patJetsAK5PFPt50Eta25NoPFElectron.clone()
patJetsAK5PFPt50Eta25NoPFElectronID.checkOverlaps.muons.src               = cms.InputTag('patElectronsPFID')

patJetsAK5PFPt50Eta25NoPFElectronIDIso = patJetsAK5PFPt50Eta25NoPFElectron.clone()
patJetsAK5PFPt50Eta25NoPFElectronIDIso.checkOverlaps.muons.src            = cms.InputTag('patElectronsPFIDIso')


# PFJets - filters

countJetsAK5PFPt50Eta25NoPFElectron           = countPatJets.clone()
countJetsAK5PFPt50Eta25NoPFElectron.src       = cms.InputTag('patJetsAK5PFPt50Eta25NoPFElectron')
countJetsAK5PFPt50Eta25NoPFElectron.minNumber = cms.uint32(3)

countJetsAK5PFPt50Eta25NoPFElectronID     = countJetsAK5PFPt50Eta25NoPFElectron.clone()
countJetsAK5PFPt50Eta25NoPFElectronID.src = cms.InputTag('patJetsAK5PFPt50Eta25NoPFElectronID')

countJetsAK5PFPt50Eta25NoPFElectronIDIso     = countJetsAK5PFPt50Eta25NoPFElectron.clone()
countJetsAK5PFPt50Eta25NoPFElectronIDIso.src = cms.InputTag('patJetsAK5PFPt50Eta25NoPFElectronIDIso')

electronCleanedPFJets = cms.Sequence(
    patJetsAK5PFPt30NoPFElectron           *
    patJetsAK5PFPt30NoPFElectronID         *
    patJetsAK5PFPt30NoPFElectronIDIso      *
    patJetsAK5PFPt50Eta25NoPFElectron      *
    patJetsAK5PFPt50Eta25NoPFElectronID    *
    patJetsAK5PFPt50Eta25NoPFElectronIDIso
)

######
# b-tagged jets
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
patSSVHEMBJetsAK5PF = selectedPatJets.clone()
patSSVHEMBJetsAK5PF.src = cms.InputTag('patJetsAK5PF')
patSSVHEMBJetsAK5PF.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsAK5PF = selectedPatJets.clone()
patSSVHPTBJetsAK5PF.src = cms.InputTag('patJetsAK5PF')
patSSVHPTBJetsAK5PF.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

patSSVHEMBJetsAK5PFPt30 = selectedPatJets.clone()
patSSVHEMBJetsAK5PFPt30.src = cms.InputTag('patJetsAK5PFPt30')
patSSVHEMBJetsAK5PFPt30.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsAK5PFPt30 = selectedPatJets.clone()
patSSVHPTBJetsAK5PFPt30.src = cms.InputTag('patJetsAK5PFPt30')
patSSVHPTBJetsAK5PFPt30.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

patSSVHEMBJetsAK5PFPt50Eta25 = selectedPatJets.clone()
patSSVHEMBJetsAK5PFPt50Eta25.src = cms.InputTag('patJetsAK5PFPt50Eta25')
patSSVHEMBJetsAK5PFPt50Eta25.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsAK5PFPt50Eta25 = selectedPatJets.clone()
patSSVHPTBJetsAK5PFPt50Eta25.src = cms.InputTag('patJetsAK5PFPt50Eta25')
patSSVHPTBJetsAK5PFPt50Eta25.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

### count the b-jets
countSSVHEMBJetsAK5PF = countPatJets.clone()
countSSVHEMBJetsAK5PF.src = cms.InputTag('patSSVHEMBJetsAK5PF')
countSSVHEMBJetsAK5PF.minNumber = cms.uint32(1)

countSSVHPTBJetsAK5PF = countPatJets.clone()
countSSVHPTBJetsAK5PF.src = cms.InputTag('patSSVHPTBJetsAK5PF')
countSSVHPTBJetsAK5PF.minNumber = cms.uint32(1)

countSSVHEMBJetsAK5PFPt30 = countPatJets.clone()
countSSVHEMBJetsAK5PFPt30.src = cms.InputTag('patSSVHEMBJetsAK5PFPt30')
countSSVHEMBJetsAK5PFPt30.minNumber = cms.uint32(1)

countSSVHPTBJetsAK5PFPt30 = countPatJets.clone()
countSSVHPTBJetsAK5PFPt30.src = cms.InputTag('patSSVHPTBJetsAK5PFPt30')
countSSVHPTBJetsAK5PFPt30.minNumber = cms.uint32(1)

countSSVHEMBJetsAK5PFPt50Eta25 = countPatJets.clone()
countSSVHEMBJetsAK5PFPt50Eta25.src = cms.InputTag('patSSVHEMBJetsAK5PFPt50Eta25')
countSSVHEMBJetsAK5PFPt50Eta25.minNumber = cms.uint32(1)

countSSVHPTBJetsAK5PFPt50Eta25 = countPatJets.clone()
countSSVHPTBJetsAK5PFPt50Eta25.src = cms.InputTag('patSSVHPTBJetsAK5PFPt50Eta25')
countSSVHPTBJetsAK5PFPt50Eta25.minNumber = cms.uint32(1)

zinvBVeto = cms.Sequence(
    ~countSSVHEMBJetsAK5PF
)
zinvBVetoPt30 = cms.Sequence(
    ~countSSVHEMBJetsAK5PFPt30
)
zinvBVetoPt50Eta25 = cms.Sequence(
    ~countSSVHEMBJetsAK5PFPt50Eta25
)
### create the jet collections
zinvBJets = cms.Sequence(
    patSSVHEMBJetsAK5PF*
    patSSVHPTBJetsAK5PF*
    patSSVHEMBJetsAK5PFPt30*
    patSSVHPTBJetsAK5PFPt30*
    patSSVHEMBJetsAK5PFPt50Eta25*
    patSSVHPTBJetsAK5PFPt50Eta25
)
