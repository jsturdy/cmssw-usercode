# b-tagged jets
#from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
from ZInvisibleBkgds.Photons.specialJetSelector_cff import selectedRA2PatJets

patSSVHEMBJetsAK5PF = selectedRA2PatJets.clone()
patSSVHEMBJetsAK5PF.src = cms.InputTag('patJetsAK5PF')
patSSVHEMBJetsAK5PF.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsAK5PF = selectedRA2PatJets.clone()
patSSVHPTBJetsAK5PF.src = cms.InputTag('patJetsAK5PF')
patSSVHPTBJetsAK5PF.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

patSSVHEMBJetsAK5PFPt30 = selectedRA2PatJets.clone()
patSSVHEMBJetsAK5PFPt30.src = cms.InputTag('patJetsAK5PFPt30')
patSSVHEMBJetsAK5PFPt30.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsAK5PFPt30 = selectedRA2PatJets.clone()
patSSVHPTBJetsAK5PFPt30.src = cms.InputTag('patJetsAK5PFPt30')
patSSVHPTBJetsAK5PFPt30.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

patSSVHEMBJetsAK5PFPt50Eta25 = selectedRA2PatJets.clone()
patSSVHEMBJetsAK5PFPt50Eta25.src = cms.InputTag('patJetsAK5PFPt50Eta25')
patSSVHEMBJetsAK5PFPt50Eta25.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsAK5PFPt50Eta25 = selectedRA2PatJets.clone()
patSSVHPTBJetsAK5PFPt50Eta25.src = cms.InputTag('patJetsAK5PFPt50Eta25')
patSSVHPTBJetsAK5PFPt50Eta25.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

patCSVJetsAK5PF = selectedRA2PatJets.clone()
patCSVJetsAK5PF.src = cms.InputTag('patJetsAK5PF')
patCSVJetsAK5PF.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > 0.898')

patCSVJetsAK5PFPt30Eta24 = selectedRA2PatJets.clone()
patCSVJetsAK5PFPt30Eta24.src = cms.InputTag('patCSVJetsAK5PF')
patCSVJetsAK5PFPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVJetsAK5PFPt30Eta50 = selectedRA2PatJets.clone()
patCSVJetsAK5PFPt30Eta50.src = cms.InputTag('patCSVJetsAK5PF')
patCSVJetsAK5PFPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVJetsAK5PFPt50Eta25 = selectedRA2PatJets.clone()
patCSVJetsAK5PFPt50Eta25.src = cms.InputTag('patCSVJetsAK5PF')
patCSVJetsAK5PFPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

patCSVMVAJetsAK5PF = selectedRA2PatJets.clone()
patCSVMVAJetsAK5PF.src = cms.InputTag('patJetsAK5PF')
patCSVMVAJetsAK5PF.cut = cms.string('bDiscriminator("combinedSecondaryVertexMVABJetTags") > 0.898')

patCSVMVAJetsAK5PFPt30Eta24 = selectedRA2PatJets.clone()
patCSVMVAJetsAK5PFPt30Eta24.src = cms.InputTag('patCSVMVAJetsAK5PF')
patCSVMVAJetsAK5PFPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMVAJetsAK5PFPt30Eta50 = selectedRA2PatJets.clone()
patCSVMVAJetsAK5PFPt30Eta50.src = cms.InputTag('patCSVMVAJetsAK5PF')
patCSVMVAJetsAK5PFPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMVAJetsAK5PFPt50Eta25 = selectedRA2PatJets.clone()
patCSVMVAJetsAK5PFPt50Eta25.src = cms.InputTag('patCSVMVAJetsAK5PF')
patCSVMVAJetsAK5PFPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

###patJetsPF
patSSVHEMBJetsPF = selectedRA2PatJets.clone()
patSSVHEMBJetsPF.src = cms.InputTag('patJetsPF')
patSSVHEMBJetsPF.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsPF = selectedRA2PatJets.clone()
patSSVHPTBJetsPF.src = cms.InputTag('patJetsPF')
patSSVHPTBJetsPF.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

patSSVHEMBJetsPFPt30 = selectedRA2PatJets.clone()
patSSVHEMBJetsPFPt30.src = cms.InputTag('patJetsPFPt30')
patSSVHEMBJetsPFPt30.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsPFPt30 = selectedRA2PatJets.clone()
patSSVHPTBJetsPFPt30.src = cms.InputTag('patJetsPFPt30')
patSSVHPTBJetsPFPt30.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

patSSVHEMBJetsPFPt50Eta25 = selectedRA2PatJets.clone()
patSSVHEMBJetsPFPt50Eta25.src = cms.InputTag('patJetsPFPt50Eta25')
patSSVHEMBJetsPFPt50Eta25.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsPFPt50Eta25 = selectedRA2PatJets.clone()
patSSVHPTBJetsPFPt50Eta25.src = cms.InputTag('patJetsPFPt50Eta25')
patSSVHPTBJetsPFPt50Eta25.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

patCSVJetsPF = selectedRA2PatJets.clone()
patCSVJetsPF.src = cms.InputTag('patJetsPF')
patCSVJetsPF.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > 0.898')

patCSVJetsPFPt30Eta24 = selectedRA2PatJets.clone()
patCSVJetsPFPt30Eta24.src = cms.InputTag('patCSVJetsPF')
patCSVJetsPFPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVJetsPFPt30Eta50 = selectedRA2PatJets.clone()
patCSVJetsPFPt30Eta50.src = cms.InputTag('patCSVJetsPF')
patCSVJetsPFPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVJetsPFPt50Eta25 = selectedRA2PatJets.clone()
patCSVJetsPFPt50Eta25.src = cms.InputTag('patCSVJetsPF')
patCSVJetsPFPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

patCSVMVAJetsPF = selectedRA2PatJets.clone()
patCSVMVAJetsPF.src = cms.InputTag('patJetsPF')
patCSVMVAJetsPF.cut = cms.string('bDiscriminator("combinedSecondaryVertexMVABJetTags") > 0.898')

patCSVMVAJetsPFPt30Eta24 = selectedRA2PatJets.clone()
patCSVMVAJetsPFPt30Eta24.src = cms.InputTag('patCSVMVAJetsPF')
patCSVMVAJetsPFPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMVAJetsPFPt30Eta50 = selectedRA2PatJets.clone()
patCSVMVAJetsPFPt30Eta50.src = cms.InputTag('patCSVMVAJetsPF')
patCSVMVAJetsPFPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMVAJetsPFPt50Eta25 = selectedRA2PatJets.clone()
patCSVMVAJetsPFPt50Eta25.src = cms.InputTag('patCSVMVAJetsPF')
patCSVMVAJetsPFPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

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

###patJetsPF
countSSVHEMBJetsPF = countPatJets.clone()
countSSVHEMBJetsPF.src = cms.InputTag('patSSVHEMBJetsPF')
countSSVHEMBJetsPF.minNumber = cms.uint32(1)

countSSVHPTBJetsPF = countPatJets.clone()
countSSVHPTBJetsPF.src = cms.InputTag('patSSVHPTBJetsPF')
countSSVHPTBJetsPF.minNumber = cms.uint32(1)

countSSVHEMBJetsPFPt30 = countPatJets.clone()
countSSVHEMBJetsPFPt30.src = cms.InputTag('patSSVHEMBJetsPFPt30')
countSSVHEMBJetsPFPt30.minNumber = cms.uint32(1)

countSSVHPTBJetsPFPt30 = countPatJets.clone()
countSSVHPTBJetsPFPt30.src = cms.InputTag('patSSVHPTBJetsPFPt30')
countSSVHPTBJetsPFPt30.minNumber = cms.uint32(1)

countSSVHEMBJetsPFPt50Eta25 = countPatJets.clone()
countSSVHEMBJetsPFPt50Eta25.src = cms.InputTag('patSSVHEMBJetsPFPt50Eta25')
countSSVHEMBJetsPFPt50Eta25.minNumber = cms.uint32(1)

countSSVHPTBJetsPFPt50Eta25 = countPatJets.clone()
countSSVHPTBJetsPFPt50Eta25.src = cms.InputTag('patSSVHPTBJetsPFPt50Eta25')
countSSVHPTBJetsPFPt50Eta25.minNumber = cms.uint32(1)


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
zinvBJetsAK5PF = cms.Sequence(
    #  patSSVHEMBJetsAK5PF
    #* patSSVHPTBJetsAK5PF
    #* patSSVHEMBJetsAK5PFPt30
    #* patSSVHPTBJetsAK5PFPt30
    #* patSSVHEMBJetsAK5PFPt50Eta25
    #* patSSVHPTBJetsAK5PFPt50Eta25
      patCSVJetsAK5PF
    * patCSVJetsAK5PFPt30Eta24
    #* patCSVJetsAK5PFPt30Eta50
    #* patCSVJetsAK5PFPt50Eta25
    #* patCSVMVAJetsAK5PF
    #* patCSVMVAJetsAK5PFPt30Eta24
    #* patCSVMVAJetsAK5PFPt30Eta50
    #* patCSVMVAJetsAK5PFPt50Eta25
)

zinvBJetsPF = cms.Sequence(
    #  patSSVHEMBJetsPF
    #* patSSVHPTBJetsPF
    #* patSSVHEMBJetsPFPt30
    #* patSSVHPTBJetsPFPt30
    #* patSSVHEMBJetsPFPt50Eta25
    #* patSSVHPTBJetsPFPt50Eta25
      patCSVJetsPF
    * patCSVJetsPFPt30Eta24
    #* patCSVJetsPFPt30Eta50
    #* patCSVJetsPFPt50Eta25
    #* patCSVMVAJetsPF
    #* patCSVMVAJetsPFPt30Eta24
    #* patCSVMVAJetsPFPt30Eta50
    #* patCSVMVAJetsPFPt50Eta25
)
