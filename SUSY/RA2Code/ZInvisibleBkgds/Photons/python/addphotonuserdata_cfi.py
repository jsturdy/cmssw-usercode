import FWCore.ParameterSet.Config as cms

addphotonuserbasic = cms.EDProducer("AddPhotonUserData",
    debug          = cms.bool(False),
    debugString    = cms.string("addphotonuserdatabasic"),
    photonLabel    = cms.InputTag("patPhotons"),
    floatLabels    = cms.VInputTag(),
    floatNames     = cms.vstring(),
    embedConversionInfo = cms.bool(True),
    gsfElectronLabel = cms.InputTag("gsfElectrons"),
    conversionsLabel = cms.InputTag("conversions"),
    beamspotLabel    = cms.InputTag("offlineBeamSpot"),
    useAlternateIsolations = cms.bool(True),
    vetoConeSize     = cms.double(0.3),
    candidateLabel   = cms.InputTag("particleFlow"),
    vertexLabel      = cms.InputTag("offlinePrimaryVertices"),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(
            'hcalIsoConeDR03_2012','hcalIsoConeDR04_2012',
            'pfChargedEA','pfNeutralEA','pfGammaEA',
            'hadTowOverEmTightCut' ,'showerShapeTightCut' ,'pfChargedTightCut' ,'pfNeutralTightCut' ,'pfGammaTightCut',
            'hadTowOverEmMediumCut','showerShapeMediumCut','pfChargedMediumCut','pfNeutralMediumCut','pfGammaMediumCut',
            'hadTowOverEmLooseCut' ,'showerShapeLooseCut' ,'pfChargedLooseCut' ,'pfNeutralLooseCut' ,'pfGammaLooseCut',
            'hadTowOverEmVeryLooseCut' ,'showerShapeVeryLooseCut' ,'pfChargedVeryLooseCut' ,'pfNeutralVeryLooseCut' ,'pfGammaVeryLooseCut'
        ),
        userFunctions = cms.vstring(
            """hcalTowerSumEtConeDR03 +
               (hadronicOverEm - hadTowOverEm)*superCluster.energy/cosh(superCluster.eta)""", #hcalIsoConeDR03_2012
            """hcalTowerSumEtConeDR04 +
               (hadronicOverEm - hadTowOverEm)*superCluster.energy/cosh(superCluster.eta)""", #hcalIsoConeDR04_2012
            """?0.0   <= abs(superCluster.eta) < 1.0   ? 0.012 :
               ?1.0   <= abs(superCluster.eta) < 1.479 ? 0.010 :
               ?1.479 <= abs(superCluster.eta) < 2.0   ? 0.014 :
               ?2.0   <= abs(superCluster.eta) < 2.2   ? 0.012 :
               ?2.2   <= abs(superCluster.eta) < 2.3   ? 0.016 :
               ?2.3   <= abs(superCluster.eta) < 2.4   ? 0.020 :
                                                         0.012""", #chargedEA
            """?0.0   <= abs(superCluster.eta) < 1.0   ? 0.030 :
               ?1.0   <= abs(superCluster.eta) < 1.479 ? 0.057 :
               ?1.479 <= abs(superCluster.eta) < 2.0   ? 0.039 :
               ?2.0   <= abs(superCluster.eta) < 2.2   ? 0.015 :
               ?2.2   <= abs(superCluster.eta) < 2.3   ? 0.024 :
               ?2.3   <= abs(superCluster.eta) < 2.4   ? 0.039 :
                                                         0.072""", #neutralEA
            """?0.0   <= abs(superCluster.eta) < 1.0   ? 0.148 :
               ?1.0   <= abs(superCluster.eta) < 1.479 ? 0.130 :
               ?1.479 <= abs(superCluster.eta) < 2.0   ? 0.112 :
               ?2.0   <= abs(superCluster.eta) < 2.2   ? 0.216 :
               ?2.2   <= abs(superCluster.eta) < 2.3   ? 0.262 :
               ?2.3   <= abs(superCluster.eta) < 2.4   ? 0.260 :
                                                         0.266""", #gammaEA

            ####Cut values:singleTowerHOverE,showerShape,pfCharged,pfNeutral,pfGamma
            ###Tight
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.05 :0.05 """,#singleTowerHOverE
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.011:0.031""",#showerShape
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.7 :0.5 """,  #pfCharged
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.4 + 0.04*pt :1.5 + 0.04*pt """,  #pfNeutral
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.5 + 0.005*pt :1.0 + 0.005*pt """,#pfGamma
            ###Medium
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.05 :0.05 """,#singleTowerHOverE
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.011:0.033""",#showerShape      
            """?0.0   <= abs(superCluster.eta) < 1.4442? 1.5 :1.2 """,  #pfCharged
            """?0.0   <= abs(superCluster.eta) < 1.4442? 1.0 + 0.04*pt :1.5 + 0.04*pt """,  #pfNeutral
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.7 + 0.005*pt :1.0 + 0.005*pt """,#pfGamma
            ###Loose
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.05 :0.05 """,#singleTowerHOverE
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.012:0.034""",#showerShape      
            """?0.0   <= abs(superCluster.eta) < 1.4442? 2.6 :2.3 """,  #pfCharged
            """?0.0   <= abs(superCluster.eta) < 1.4442? 3.5 + 0.04*pt :2.9 + 0.04*pt """,#pfNeutral
            """?0.0   <= abs(superCluster.eta) < 1.4442? 1.3 + 0.005*pt : 1e10 """,   #pfGamma
            ###VeryLoose
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.05 :0.05 """,#singleTowerHOverE
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.015:0.037""",#showerShape
            """?0.0   <= abs(superCluster.eta) < 1.4442? 5*(0.7) :5*(0.5) """,  #pfCharged
            """?0.0   <= abs(superCluster.eta) < 1.4442? 5*(0.4 + 0.04*pt) :5*(1.5 + 0.04*pt) """,  #pfNeutral
            """?0.0   <= abs(superCluster.eta) < 1.4442? 5*(0.5 + 0.005*pt):5*(1.0 + 0.005*pt) """,#pfGamma
        )
    )
)

addphotonuserdata1 = cms.EDProducer("AddPhotonUserData",
    debug          = cms.bool(False),
    debugString    = cms.string("addphotonuserdata1"),
    photonLabel    = cms.InputTag("patPhotonsAlt"),
    floatLabels    = cms.VInputTag(cms.InputTag("kt6PFJets","rho")),
    floatNames     = cms.vstring("rho25"),
    embedConversionInfo = cms.bool(True),
    gsfElectronLabel = cms.InputTag("gsfElectrons"),
    conversionsLabel = cms.InputTag("conversions"),
    beamspotLabel    = cms.InputTag("offlineBeamSpot"),
    useAlternateIsolations = cms.bool(True),
    vetoConeSize     = cms.double(0.3),
    candidateLabel   = cms.InputTag("particleFlow"),
    vertexLabel      = cms.InputTag("offlinePrimaryVertices"),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(
            'hcalIsoConeDR03_2012','hcalIsoConeDR04_2012',
            'pfChargedEA','pfNeutralEA','pfGammaEA',
            'hadTowOverEmTightCut' ,'showerShapeTightCut' ,'pfChargedTightCut' ,'pfNeutralTightCut' ,'pfGammaTightCut',
            'hadTowOverEmMediumCut','showerShapeMediumCut','pfChargedMediumCut','pfNeutralMediumCut','pfGammaMediumCut',
            'hadTowOverEmLooseCut' ,'showerShapeLooseCut' ,'pfChargedLooseCut' ,'pfNeutralLooseCut' ,'pfGammaLooseCut',
            'hadTowOverEmVeryLooseCut' ,'showerShapeVeryLooseCut' ,'pfChargedVeryLooseCut' ,'pfNeutralVeryLooseCut' ,'pfGammaVeryLooseCut',
            'combIsoR03','combIsoR04',
        ),
        userFunctions = cms.vstring(
            """hcalTowerSumEtConeDR03 +
               (hadronicOverEm - hadTowOverEm)*superCluster.energy/cosh(superCluster.eta)""", #hcalIsoConeDR03_2012
            """hcalTowerSumEtConeDR04 +
               (hadronicOverEm - hadTowOverEm)*superCluster.energy/cosh(superCluster.eta)""", #hcalIsoConeDR04_2012
            """?0.0   <= abs(superCluster.eta) < 1.0   ? 0.012 :
               ?1.0   <= abs(superCluster.eta) < 1.479 ? 0.010 :
               ?1.479 <= abs(superCluster.eta) < 2.0   ? 0.014 :
               ?2.0   <= abs(superCluster.eta) < 2.2   ? 0.012 :
               ?2.2   <= abs(superCluster.eta) < 2.3   ? 0.016 :
               ?2.3   <= abs(superCluster.eta) < 2.4   ? 0.020 :
                                                         0.012""", #chargedEA
            """?0.0   <= abs(superCluster.eta) < 1.0   ? 0.030 :
               ?1.0   <= abs(superCluster.eta) < 1.479 ? 0.057 :
               ?1.479 <= abs(superCluster.eta) < 2.0   ? 0.039 :
               ?2.0   <= abs(superCluster.eta) < 2.2   ? 0.015 :
               ?2.2   <= abs(superCluster.eta) < 2.3   ? 0.024 :
               ?2.3   <= abs(superCluster.eta) < 2.4   ? 0.039 :
                                                         0.072""", #neutralEA
            """?0.0   <= abs(superCluster.eta) < 1.0   ? 0.148 :
               ?1.0   <= abs(superCluster.eta) < 1.479 ? 0.130 :
               ?1.479 <= abs(superCluster.eta) < 2.0   ? 0.112 :
               ?2.0   <= abs(superCluster.eta) < 2.2   ? 0.216 :
               ?2.2   <= abs(superCluster.eta) < 2.3   ? 0.262 :
               ?2.3   <= abs(superCluster.eta) < 2.4   ? 0.260 :
                                                         0.266""", #gammaEA

            ####Cut values:singleTowerHOverE,showerShape,pfCharged,pfNeutral,pfGamma
            ###Tight
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.05 :0.05 """,#singleTowerHOverE
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.011:0.031""",#showerShape
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.7 :0.5 """,  #pfCharged
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.4 + 0.04*pt :1.5 + 0.04*pt """,  #pfNeutral
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.5 + 0.005*pt :1.0 + 0.005*pt """,#pfGamma
            ###Medium
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.05 :0.05 """,#singleTowerHOverE
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.011:0.033""",#showerShape      
            """?0.0   <= abs(superCluster.eta) < 1.4442? 1.5 :1.2 """,  #pfCharged
            """?0.0   <= abs(superCluster.eta) < 1.4442? 1.0 + 0.04*pt :1.5 + 0.04*pt """,  #pfNeutral
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.7 + 0.005*pt :1.0 + 0.005*pt """,#pfGamma
            ###Loose
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.05 :0.05 """,#singleTowerHOverE
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.012:0.034""",#showerShape      
            """?0.0   <= abs(superCluster.eta) < 1.4442? 2.6 :2.3 """,  #pfCharged
            """?0.0   <= abs(superCluster.eta) < 1.4442? 3.5 + 0.04*pt :2.9 + 0.04*pt """,#pfNeutral
            """?0.0   <= abs(superCluster.eta) < 1.4442? 1.3 + 0.005*pt : 1e10 """,   #pfGamma
            ###VeryLoose
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.05 :0.05 """,#singleTowerHOverE
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.015:0.037""",#showerShape
            """?0.0   <= abs(superCluster.eta) < 1.4442? 5*(0.7) :5*(0.5) """,  #pfCharged
            """?0.0   <= abs(superCluster.eta) < 1.4442? 5*(0.4 + 0.04*pt) :5*(1.5 + 0.04*pt) """,  #pfNeutral
            """?0.0   <= abs(superCluster.eta) < 1.4442? 5*(0.5 + 0.005*pt):5*(1.0 + 0.005*pt) """,#pfGamma
            'trkSumPtSolidConeDR03 + ecalRecHitSumEtConeDR03 + userFloat("hcalIsoConeDR03_2012")',
            'trkSumPtSolidConeDR04 + ecalRecHitSumEtConeDR04 + userFloat("hcalIsoConeDR04_2012")'
        )
    )
)
addphotonuserdata2 = cms.EDProducer("AddPhotonUserData",
    debug          = cms.bool(False),
    debugString    = cms.string("addphotonuserdata2"),
    photonLabel    = cms.InputTag("patPhotonsUser1"),
    floatLabels    = cms.VInputTag(),
    floatNames     = cms.vstring(),
    embedConversionInfo = cms.bool(False),
    gsfElectronLabel = cms.InputTag("gsfElectrons"),
    conversionsLabel = cms.InputTag("conversions"),
    beamspotLabel    = cms.InputTag("offlineBeamSpot"),
    useAlternateIsolations = cms.bool(False),
    vetoConeSize     = cms.double(0.3),
    candidateLabel   = cms.InputTag("particleFlow"),
    vertexLabel      = cms.InputTag("offlinePrimaryVertices"),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring('pfChargedPUSub','pfNeutralPUSub','pfGammaPUSub',
                                         'pfChargedPU','pfNeutralPU','pfGammaPU',
                                         'combIsoR03PU','combIsoR04PU'),
        userFunctions = cms.vstring(
            'userFloat("pfChargedEA")*userFloat("rho25")', #chargedSub
            'userFloat("pfNeutralEA")*userFloat("rho25")', #neutralSub
            'userFloat("pfGammaEA")  *userFloat("rho25")', #gammaSub
            'max((userFloat("pfChargedIsoAlt") - userFloat("pfChargedEA")*userFloat("rho25")),0.)',
            'max((userFloat("pfNeutralIsoAlt") - userFloat("pfNeutralEA")*userFloat("rho25")),0.)',
            'max((userFloat("pfGammaIsoAlt")   - userFloat("pfGammaEA")  *userFloat("rho25")),0.)',
            'trkSumPtSolidConeDR03 + ecalRecHitSumEtConeDR03 + userFloat("hcalIsoConeDR03_2012") - 3.141593*0.3*0.3*userFloat("rho25")',
            'trkSumPtSolidConeDR04 + ecalRecHitSumEtConeDR04 + userFloat("hcalIsoConeDR04_2012") - 3.141593*0.4*0.4*userFloat("rho25")'
        )
    )
)
