
from PhysicsTools.PatAlgos.selectionLayer1.photonSelector_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.photonCountFilter_cfi import *

photonIDCutTight  = cms.string('et > 50.0 && (abs(eta) < 1.4442 || (abs(eta) > 1.566 && abs(eta) < 2.5)) && '
                               'hadronicOverEm < 0.5  && hasPixelSeed == 0 && '
                               'hadTowOverEm < userFloat("hadTowOverEmTightCut") && '
                               'sigmaIetaIeta < userFloat("showerShapeTightCut")'
                              )

photonIDCutMedium = cms.string('et > 50.0 && (abs(eta) < 1.4442 || (abs(eta) > 1.566 && abs(eta) < 2.5)) && '
                               'hadronicOverEm < 0.5  && hasPixelSeed == 0 && '
                               'hadTowOverEm < userFloat("hadTowOverEmMediumCut") && '
                               'sigmaIetaIeta < userFloat("showerShapeMediumCut")'
                              )

photonIDCutLoose  = cms.string('et > 50.0 && (abs(eta) < 1.4442 || (abs(eta) > 1.566 && abs(eta) < 2.5)) && '
                               'hadronicOverEm < 0.5  && hasPixelSeed == 0 && '
                               'hadTowOverEm < userFloat("hadTowOverEmLooseCut") && '
                               'sigmaIetaIeta < userFloat("showerShapeLooseCut")'
                              )

photoncombiso03cut = cms.string('trkSumPtSolidConeDR03 + ecalRecHitSumEtConeDR03 + '
                                'userFloat("hcalIsoConeDR03") - 3.141593*0.3*0.3*userFloat("rho25")< 5.'
                               )
photoncombiso04cut = cms.string('trkSumPtSolidConeDR04 + ecalRecHitSumEtConeDR04 + '
                                'userFloat("hcalIsoConeDR04") - 3.141593*0.4*0.4*userFloat("rho25")< 5.'
                               )

photonISOCutTight  = cms.string('userFloat("pfChargedPURel") < userFloat("pfChargedTightCut") && '
                                'userFloat("pfNeutralPURel") < userFloat("pfNeutralTightCut") && '
                                'userFloat("pfGammaPURel")   < userFloat("pfGammaTightCut")'
                               )

photonISOCutMedium = cms.string('userFloat("pfChargedPURel") < userFloat("pfChargedMediumCut") && '
                                'userFloat("pfNeutralPURel") < userFloat("pfNeutralMediumCut") && '
                                'userFloat("pfGammaPURel")   < userFloat("pfGammaMediumCut")'
                               )

photonISOCutLoose  = cms.string('userFloat("pfChargedPURel") < userFloat("pfChargedLooseCut") && '
                                'userFloat("pfNeutralPURel") < userFloat("pfNeutralLooseCut") && '
                                'userFloat("pfGammaPURel")   < userFloat("pfGammaLooseCut")'
                               )


patPhotonsID = cms.EDFilter(
  "PATPhotonSelector",
   src = cms.InputTag('patPhotonsUserData'),
   cut = photonIDCutTight,
   filter = cms.bool(False),
)
patPhotonsIDPFIso = patPhotonsID.clone(src = cms.InputTag('patPhotonsID'))
patPhotonsIDPFIso.cut = photonISOCutTight

patPhotonsIDIso = patPhotonsID.clone(
   src = cms.InputTag('patPhotonsID'),
   cut = photoncombiso03cut,
)

patPhotonRefsID = cms.EDFilter(
  "PATPhotonRefSelector",
   src = cms.InputTag('patPhotonsUserData'),
   cut = photonIDCutTight,
   filter = cms.bool(False)
)

zinvPhotons = cms.Sequence(
    patPhotonsID
  * patPhotonsIDPFIso
)

countPhotonsID = countPatPhotons.clone()
countPhotonsID.src = cms.InputTag("patPhotonsID")
countPhotonsID.minNumber = cms.uint32(1)

countPhotonsIDPFIso = countPatPhotons.clone()
countPhotonsIDPFIso.src = cms.InputTag("patPhotonsIDPFIso")
countPhotonsIDPFIso.minNumber = cms.uint32(1)

countPhotonRefsID = countPatPhotons.clone()
countPhotonRefsID.src = cms.InputTag("patPhotonsRefsID")
countPhotonRefsID.minNumber = cms.uint32(1)