
# muon filters
from SandBox.Skims.RA2Photons_cff import *
from PhysicsTools.PatAlgos.selectionLayer1.photonCountFilter_cfi import *

photonidcut    = cms.string('pt > 70.0 && abs(eta) < 2.5 && hadronicOverEm < 0.05  && hasPixelSeed == 0 && ((isEB && sigmaIetaIeta < 0.01) || (isEE && sigmaIetaIeta < 0.03))')
photonisocut   = cms.string('')
photonpfisocut = cms.string('')

patPhotonsPFID = patPhotonsID.clone()

patPhotonsPFIDIso = patPhotonsID.clone()
patPhotonsPFIDIso.src = cms.InputTag('patPhotonsPFID')
patPhotonsPFIDIso.cut = photonpfisocut

patPhotonsIDIso = patPhotonsID.clone()
patPhotonsIDIso.src = cms.InputTag('patPhotonsID')
patPhotonsIDIso.cut = photonisocut

photonsPF = cms.Sequence(
    patPhotonsPFID+
    patPhotonsPFIDIso
)

photons = cms.Sequence(
    patPhotonsID+
    patPhotonsIDIso
)

countPhotonsPFIDGamma = countPatPhotons.clone()
countPhotonsPFIDGamma.src = cms.InputTag('patPhotonsPFID')
countPhotonsPFIDGamma.minNumber = cms.uint32(1)
countPhotonsPFIDGamma.maxNumber = cms.uint32(1)

countPhotonsPFIDIsoGamma = countPatPhotons.clone()
countPhotonsPFIDIsoGamma.src = cms.InputTag('patPhotonsPFIDIso')
countPhotonsPFIDIsoGamma.minNumber = cms.uint32(1)
countPhotonsPFIDIsoGamma.maxNumber = cms.uint32(1)

countPhotonsIDGamma = countPatPhotons.clone()
countPhotonsIDGamma.src = cms.InputTag('patPhotonsID')
countPhotonsIDGamma.minNumber = cms.uint32(1)
countPhotonsIDGamma.maxNumber = cms.uint32(1)

countPhotonsIDIsoGamma = countPatPhotons.clone()
countPhotonsIDIsoGamma.src = cms.InputTag('patPhotonsIDIso')
countPhotonsIDIsoGamma.minNumber = cms.uint32(1)
countPhotonsIDIsoGamma.maxNumber = cms.uint32(1)

countPFGamma = cms.Sequence(
    countPhotonsPFIDIsoGamma
)

countGamma = cms.Sequence(
    countPhotonsIDIsoGamma
)
