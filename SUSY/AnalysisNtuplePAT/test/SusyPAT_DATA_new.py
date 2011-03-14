import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")
#process = cms.Process("NTUPLE")

#-- VarParsing ----------------------------------------------------------------
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('standard')
# setup any defaults you want
options.register ('patout',
    'PATtuple_V9_DATA.root', # default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.string,          # string, int, or float
    "Name of PATtuple file"
)
options.register ('outputfile',
    'SUSYPAT_V9_DATA.root', # default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.string,          # string, int, or float
    "Name of PATtuple file"
)
options.register ('globtag',
    'GR_R_39X_V6', # default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.string,          # string, int, or float
    "Name of GlobalTag"
)
options.register ('hltname',
    'HLT', # default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.string,          # string, int, or float
    "Name of TriggerResults"
)

options.parseArguments()
options._tagOrder =[]

print("outputfile:"+options.outputfile)
print("patout:"+options.patout)
print("globtag:"+options.globtag)
print("hltname:"+options.hltname)


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/0217C878-BB0D-E011-B7D2-003048678B14.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/02FE72AD-DB0D-E011-814A-0018F3D096A4.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/0C274620-C70D-E011-A6C7-001A92971BD6.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/160C265E-B50D-E011-B18D-0026189437E8.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/1855D63E-ED0D-E011-B5E0-001A928116EE.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/1C3A372C-CF0D-E011-A0D6-003048D15E2C.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/1C854310-DA0D-E011-A8BF-002618FDA210.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/1E180735-E40D-E011-AB2B-0026189438DE.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/24D5AA5A-D00D-E011-B833-001A92971B96.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/28576202-E60D-E011-BD95-002618943870.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/2C2733AA-E90D-E011-90DF-002618943870.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/2C4F4445-CF0D-E011-A473-0018F3D096BA.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/30FC8258-D10D-E011-BC95-002618943950.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/32B08C1E-DD0D-E011-8063-00304867904E.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/44C67412-C50D-E011-8EC5-003048678D52.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/4644E549-010E-E011-8214-00261894387E.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/4A05893F-CB0D-E011-9C34-003048678F92.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/4A204987-F80D-E011-9DAB-00304867C0C4.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/4C13F7A9-BF0D-E011-A34C-003048678DA2.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/4C37A366-D30D-E011-9FDA-0026189438B0.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/5A6A6C22-C20D-E011-9258-00304867BED8.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/5CBF9D4E-D10D-E011-8D7F-001A9281170C.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/62625828-D00D-E011-9C5A-002618943950.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/64A78527-F30D-E011-A7E3-0018F3D096F6.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/68F3C4F3-D80D-E011-8891-00304867C1BA.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/7247BFA2-EF0D-E011-9237-003048678F02.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/743647B1-DB0D-E011-A8BC-001A92971B96.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/822B117A-D20D-E011-962E-002618943885.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/863D7935-DA0D-E011-89B0-003048679214.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/8CCEAF45-010E-E011-B15C-002354EF3BE3.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/8CF873F0-AE0D-E011-A380-0026189438B0.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/92C336BE-D90D-E011-A381-00304867C1BA.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/94F8DFB5-C00D-E011-BC08-002618943864.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/9AD0BD9C-BE0D-E011-BF31-0018F3D0970A.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/A467A738-CF0D-E011-82F6-003048678FF4.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/A615C775-B30D-E011-9AC6-0030486792B6.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/B27FCB43-E30D-E011-909F-003048678FDE.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/BE26312B-C50D-E011-B796-0030486790B0.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/C2685434-DA0D-E011-B6E5-0018F3D09692.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/C6DFC096-FE0D-E011-8E83-001A92971B7E.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/C83B269D-020E-E011-A42F-003048D15DDA.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/CAA64DC3-AC0D-E011-AAF3-003048D15DB6.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/DC696035-CB0D-E011-950B-003048678D52.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/E6EAEE81-E20D-E011-B91F-001A92811702.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/F64BF2B1-E80D-E011-A7DF-002354EF3BE0.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/F6888514-F70D-E011-9488-001A92971B48.root',
    '/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0000/F6B2C77C-B40D-E011-B2A7-0018F3D0962C.root'
    #'file:/tmp/sturdy/SUSYPAT_V9_DATA.root'
    ),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)
process.JPTAntiKt5JetExtender = cms.EDProducer("JetExtender",
    jets = cms.InputTag("ak5CaloJets"),
    jet2TracksAtCALO = cms.InputTag("JPTAntiKt5JetTracksAssociatorAtCaloFace"),
    jet2TracksAtVX = cms.InputTag("JPTAntiKt5JetTracksAssociatorAtVertex"),
    coneSize = cms.double(0.5)
)


process.JPTAntiKt5JetTracksAssociatorAtCaloFace = cms.EDProducer("JetTracksAssociatorAtCaloFace",
    trackQuality = cms.string('goodIterative'),
    jets = cms.InputTag("TCTauJetPlusTrackZSPCorJetAntiKt5"),
    tracks = cms.InputTag("generalTracks"),
    coneSize = cms.double(0.5),
    extrapolations = cms.InputTag("trackExtrapolator")
)


process.JPTAntiKt5JetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
    jets = cms.InputTag("TCTauJetPlusTrackZSPCorJetAntiKt5"),
    tracks = cms.InputTag("generalTracks"),
    coneSize = cms.double(0.5)
)


process.JPTCaloRecoTauProducer = cms.EDProducer("CaloRecoTauProducer",
    LeadTrack_minPt = cms.double(1.0),
    MatchingConeSize_min = cms.double(0.0),
    ECALSignalConeSizeFormula = cms.string('0.15'),
    TrackerIsolConeMetric = cms.string('DR'),
    TrackerSignalConeMetric = cms.string('DR'),
    EBRecHitsSource = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    IsolationTrack_minPt = cms.double(1.0),
    ECALSignalConeSize_min = cms.double(0.0),
    ECALRecHit_minEt = cms.double(0.5),
    MatchingConeMetric = cms.string('DR'),
    TrackerSignalConeSizeFormula = cms.string('0.07'),
    MatchingConeSizeFormula = cms.string('0.10'),
    TrackerIsolConeSize_min = cms.double(0.0),
    TrackerIsolConeSize_max = cms.double(0.6),
    TrackerSignalConeSize_max = cms.double(0.6),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    ESRecHitsSource = cms.InputTag("ecalPreshowerRecHit","EcalRecHitsES"),
    TrackerSignalConeSize_min = cms.double(0.0),
    JetPtMin = cms.double(0.0),
    AreaMetric_recoElements_maxabsEta = cms.double(2.5),
    ECALIsolConeMetric = cms.string('DR'),
    ECALIsolConeSizeFormula = cms.string('0.50'),
    ECALIsolConeSize_max = cms.double(0.6),
    EERecHitsSource = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    IsolationTrack_minHits = cms.uint32(0),
    ECALSignalConeMetric = cms.string('DR'),
    TrackLeadTrack_maxDZ = cms.double(0.2),
    Track_minPt = cms.double(0.5),
    TrackerIsolConeSizeFormula = cms.string('0.50'),
    ECALSignalConeSize_max = cms.double(0.6),
    ECALIsolConeSize_min = cms.double(0.0),
    UseTrackLeadTrackDZconstraint = cms.bool(True),
    smearedPVsigmaY = cms.double(0.0015),
    smearedPVsigmaX = cms.double(0.0015),
    smearedPVsigmaZ = cms.double(0.005),
    CaloRecoTauTagInfoProducer = cms.InputTag("caloRecoTauTagInfoProducer"),
    MatchingConeSize_max = cms.double(0.6)
)


process.JPTeidTight = cms.EDProducer("EleIdCutBasedExtProducer",
    electronQuality = cms.string('tight'),
    classbasedtightEleIDCutsV02 = cms.PSet(
        cutisohcal = cms.vdouble(10.9, 7.01, 8.75, 3.51, 7.75, 
            1.62, 11.6, 9.9, 4.97, 5.33, 
            3.18, 2.32, 0.164, 5.46, 12.0, 
            0.00604, 4.1, 0.000628),
        cutmishits = cms.vdouble(5.5, 1.5, 0.5, 1.5, 2.5, 
            0.5, 3.5, 5.5, 0.5, 0.5, 
            0.5, 0.5, 0.5, 1.5, 0.5, 
            0.5, 0.5, 0.5),
        cuthoe = cms.vdouble(0.0871, 0.0289, 0.0783, 0.0946, 0.0245, 
            0.0363, 0.0671, 0.048, 0.0614, 0.0924, 
            0.0158, 0.049, 0.0382, 0.0915, 0.0451, 
            0.0452, 0.00196, 0.0043),
        cutdeta = cms.vdouble(0.00915, 0.00302, 0.0061, 0.0135, 0.00565, 
            0.00793, 0.0102, 0.00266, 0.0106, 0.00903, 
            0.00766, 0.00723, 0.0116, 0.00203, 0.00659, 
            0.0148, 0.00555, 0.0128),
        cuteopin = cms.vdouble(0.878, 0.859, 0.874, 0.944, 0.737, 
            0.773, 0.86, 0.967, 0.917, 0.812, 
            0.915, 1.01, 0.847, 0.953, 0.979, 
            0.841, 0.771, 1.09),
        cutip = cms.vdouble(0.0239, 0.027, 0.0768, 0.0231, 0.178, 
            0.0957, 0.0102, 0.0168, 0.043, 0.0166, 
            0.0594, 0.0308, 2.1, 0.00527, 3.17, 
            4.91, 0.769, 5.9),
        cutisotk = cms.vdouble(6.53, 4.6, 6.0, 8.63, 3.11, 
            7.77, 5.42, 4.81, 4.06, 6.47, 
            2.8, 3.45, 5.29, 5.18, 15.4, 
            5.38, 4.47, 0.0347),
        cutsee = cms.vdouble(0.0131, 0.0106, 0.0115, 0.0306, 0.028, 
            0.0293, 0.0131, 0.0106, 0.0115, 0.0317, 
            0.029, 0.0289, 0.0142, 0.0106, 0.0103, 
            0.035, 0.0296, 0.0333),
        cutdphi = cms.vdouble(0.0369, 0.0307, 0.117, 0.0475, 0.0216, 
            0.117, 0.0372, 0.0246, 0.0426, 0.0612, 
            0.0142, 0.039, 0.0737, 0.0566, 0.0359, 
            0.0187, 0.012, 0.0358),
        cutisoecal = cms.vdouble(20.0, 27.2, 4.48, 13.5, 4.56, 
            3.19, 12.2, 13.1, 7.42, 7.67, 
            4.12, 4.85, 10.1, 12.4, 11.1, 
            11.0, 10.6, 13.4)
    ),
    classbasedtightEleIDCutsV03 = cms.PSet(
        cutdetain = cms.vdouble(0.00811, 0.00341, 0.00633, 0.0103, 0.00667, 
            0.01, 0.0106, 0.0145, 0.0163, 0.0076, 
            0.00259, 0.00511, 0.00941, 0.0043, 0.00857, 
            0.012, 0.0169, 0.00172, 0.00861, 0.00362, 
            0.00601, 0.00925, 0.00489, 0.00832, 0.0119, 
            0.0169, 0.000996),
        cutiso_sum = cms.vdouble(11.8, 8.31, 6.26, 6.18, 3.28, 
            4.38, 4.17, 5.4, 1.57, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0213, 0.0422, 0.0632, 0.0361, 0.073, 
            0.126, 0.171, 0.119, 0.0372, 0.0131, 
            0.0146, 0.0564, 0.0152, 0.0222, 0.0268, 
            0.0314, 0.0884, 0.00374, 0.00852, 0.00761, 
            0.0143, 0.0106, 0.0127, 0.0119, 0.0123, 
            0.0235, 0.00363),
        cuthoe = cms.vdouble(0.0783, 0.0387, 0.105, 0.118, 0.0227, 
            0.062, 0.13, 2.47, 0.38, 0.0888, 
            0.0503, 0.0955, 0.0741, 0.015, 0.03, 
            0.589, 1.13, 0.612, 0.0494, 0.0461, 
            0.0292, 0.0369, 0.0113, 0.0145, 0.124, 
            2.05, 0.61),
        cutfmishits = cms.vdouble(2.5, 1.5, 1.5, 1.5, 1.5, 
            0.5, 2.5, 0.5, 0.5, 2.5, 
            1.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5, -0.5, 2.5, 1.5, 
            0.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(13.7, 11.6, 7.14, 9.98, 3.52, 
            4.87, 6.24, 7.96, 2.53, 11.2, 
            11.9, 7.88, 8.16, 5.58, 5.03, 
            11.4, 8.15, 5.79, 10.4, 11.1, 
            10.4, 7.47, 5.08, 5.9, 11.8, 
            14.1, 11.7),
        cutdcotdist = cms.vdouble(0.0393, 0.0256, 0.00691, 0.0394, 0.0386, 
            0.039, 0.0325, 0.0384, 0.0382, 0.0245, 
            0.000281, 5.46e-05, 0.0342, 0.0232, 0.00107, 
            0.0178, 0.0193, 0.000758, 0.000108, 0.0248, 
            0.000458, 0.0129, 0.00119, 0.0182, 4.53e-05, 
            0.0189, 0.000928),
        cutsee = cms.vdouble(0.0143, 0.0105, 0.0123, 0.0324, 0.0307, 
            0.0301, 0.0109, 0.027, 0.0292, 0.0133, 
            0.0104, 0.0116, 0.0332, 0.0296, 0.031, 
            0.00981, 0.0307, 0.072, 0.0149, 0.0105, 
            0.011, 0.0342, 0.0307, 0.0303, 0.00954, 
            0.0265, 0.0101),
        cuteseedopcor = cms.vdouble(0.784, 0.366, 0.57, 0.911, 0.298, 
            0.645, 0.51, 0.497, 0.932, 0.835, 
            0.968, 0.969, 0.923, 0.898, 0.98, 
            0.63, 0.971, 1.0, 0.515, 0.963, 
            0.986, 0.823, 0.879, 1.01, 0.931, 
            0.937, 1.05),
        cutdphiin = cms.vdouble(0.0404, 0.0499, 0.263, 0.042, 0.0484, 
            0.241, 0.242, 0.231, 0.286, 0.0552, 
            0.0338, 0.154, 0.0623, 0.0183, 0.0392, 
            0.0547, 0.0588, 0.00654, 0.042, 0.0217, 
            0.0885, 0.0445, 0.0141, 0.0234, 0.065, 
            0.0258, 0.0346),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 13.7, 13.2, 
            13.6, 14.2, 14.1, 13.9, 12.9, 
            14.9, 17.7)
    ),
    classbasedtightEleIDCutsV00 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.032, 0.016, 0.0525, 0.09, 0.025, 
            0.035, 0.065, 0.092),
        hOverE = cms.vdouble(0.05, 0.042, 0.045, 0.0, 0.055, 
            0.037, 0.05, 0.0),
        sigmaEtaEta = cms.vdouble(0.0125, 0.011, 0.01, 0.0, 0.0265, 
            0.0252, 0.026, 0.0),
        deltaEtaIn = cms.vdouble(0.0055, 0.003, 0.0065, 0.0, 0.006, 
            0.0055, 0.0075, 0.0),
        eSeedOverPin = cms.vdouble(0.24, 0.94, 0.11, 0.0, 0.32, 
            0.83, 0.0, 0.0)
    ),
    classbasedtightEleIDCutsV01 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.0225, 0.0114, 0.0234, 0.039, 0.0215, 
            0.0095, 0.0148, 0.0167),
        hOverE = cms.vdouble(0.056, 0.0221, 0.037, 0.0, 0.0268, 
            0.0102, 0.0104, 0.0),
        sigmaEtaEta = cms.vdouble(0.0095, 0.0094, 0.0094, 0.0, 0.026, 
            0.0257, 0.0246, 0.0),
        deltaEtaIn = cms.vdouble(0.0043, 0.00282, 0.0036, 0.0, 0.0066, 
            0.0049, 0.0041, 0.0),
        eSeedOverPin = cms.vdouble(0.32, 0.94, 0.221, 0.0, 0.74, 
            0.89, 0.66, 0.0)
    ),
    classbasedtightEleIDCutsV06 = cms.PSet(
        cutdetain = cms.vdouble(0.0116, 0.00449, 0.00938, 0.0184, 0.00678, 
            0.0109, 0.0252, 0.0268, 0.0139),
        cutiso_sum = cms.vdouble(15.5, 12.2, 12.2, 11.7, 7.16, 
            9.71, 8.66, 11.9, 2.98),
        cutip_gsf = cms.vdouble(0.0131, 0.0586, 0.0839, 0.0366, 0.452, 
            0.204, 0.0913, 0.0802, 0.0731),
        cutip_gsfl = cms.vdouble(0.0119, 0.0527, 0.0471, 0.0212, 0.233, 
            0.267, 0.109, 0.122, 0.0479),
        cuthoe = cms.vdouble(0.215, 0.0608, 0.147, 0.369, 0.0349, 
            0.102, 0.52, 0.422, 0.404),
        cutiso_sumoetl = cms.vdouble(6.21, 6.81, 5.3, 5.39, 2.73, 
            4.73, 4.84, 3.46, 3.73),
        cutfmishits = cms.vdouble(1.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 1.5, 2.5, 0.5),
        cuthoel = cms.vdouble(0.228, 0.0836, 0.143, 0.37, 0.0392, 
            0.0979, 0.3, 0.381, 0.339),
        cutdphiin = cms.vdouble(0.0897, 0.0993, 0.295, 0.0979, 0.151, 
            0.252, 0.341, 0.308, 0.328),
        cutseel = cms.vdouble(0.0132, 0.0117, 0.0112, 0.0387, 0.0281, 
            0.0287, 0.00987, 0.0296, 0.0544),
        cutiso_sumoet = cms.vdouble(11.9, 7.81, 6.28, 8.92, 4.65, 
            5.49, 9.36, 8.84, 5.94),
        cutsee = cms.vdouble(0.0145, 0.0116, 0.012, 0.039, 0.0297, 
            0.0311, 0.00987, 0.0347, 0.0917),
        cuteseedopcor = cms.vdouble(0.637, 0.943, 0.742, 0.748, 0.763, 
            0.631, 0.214, 0.873, 0.473),
        cutdphiinl = cms.vdouble(0.061, 0.14, 0.286, 0.0921, 0.197, 
            0.24, 0.333, 0.303, 0.258),
        cutdetainl = cms.vdouble(0.00816, 0.00401, 0.0081, 0.019, 0.00588, 
            0.00893, 0.0171, 0.0434, 0.0143)
    ),
    classbasedtightEleIDCutsV04 = cms.PSet(
        cutdetain = cms.vdouble(0.00811, 0.00341, 0.00633, 0.0103, 0.00667, 
            0.01, 0.0106, 0.0145, 0.0163, 0.0076, 
            0.00259, 0.00511, 0.00941, 0.0043, 0.00857, 
            0.012, 0.0169, 0.00172, 0.00861, 0.00362, 
            0.00601, 0.00925, 0.00489, 0.00832, 0.0119, 
            0.0169, 0.000996),
        cutiso_sum = cms.vdouble(11.8, 8.31, 6.26, 6.18, 3.28, 
            4.38, 4.17, 5.4, 1.57, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0213, 0.0422, 0.0632, 0.0361, 0.073, 
            0.126, 0.171, 0.119, 0.0372, 0.0131, 
            0.0146, 0.0564, 0.0152, 0.0222, 0.0268, 
            0.0314, 0.0884, 0.00374, 0.00852, 0.00761, 
            0.0143, 0.0106, 0.0127, 0.0119, 0.0123, 
            0.0235, 0.00363),
        cuthoe = cms.vdouble(0.0783, 0.0387, 0.105, 0.118, 0.0227, 
            0.062, 0.13, 2.47, 0.38, 0.0888, 
            0.0503, 0.0955, 0.0741, 0.015, 0.03, 
            0.589, 1.13, 0.612, 0.0494, 0.0461, 
            0.0292, 0.0369, 0.0113, 0.0145, 0.124, 
            2.05, 0.61),
        cutfmishits = cms.vdouble(2.5, 1.5, 1.5, 1.5, 1.5, 
            0.5, 2.5, 0.5, 0.5, 2.5, 
            1.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5, -0.5, 2.5, 1.5, 
            0.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(13.7, 11.6, 7.14, 9.98, 3.52, 
            4.87, 6.24, 7.96, 2.53, 11.2, 
            11.9, 7.88, 8.16, 5.58, 5.03, 
            11.4, 8.15, 5.79, 10.4, 11.1, 
            10.4, 7.47, 5.08, 5.9, 11.8, 
            14.1, 11.7),
        cutdcotdist = cms.vdouble(0.0393, 0.0256, 0.00691, 0.0394, 0.0386, 
            0.039, 0.0325, 0.0384, 0.0382, 0.0245, 
            0.000281, 5.46e-05, 0.0342, 0.0232, 0.00107, 
            0.0178, 0.0193, 0.000758, 0.000108, 0.0248, 
            0.000458, 0.0129, 0.00119, 0.0182, 4.53e-05, 
            0.0189, 0.000928),
        cutsee = cms.vdouble(0.0143, 0.0105, 0.0123, 0.0324, 0.0307, 
            0.0301, 0.0109, 0.027, 0.0292, 0.0133, 
            0.0104, 0.0116, 0.0332, 0.0296, 0.031, 
            0.00981, 0.0307, 0.072, 0.0149, 0.0105, 
            0.011, 0.0342, 0.0307, 0.0303, 0.00954, 
            0.0265, 0.0101),
        cuteseedopcor = cms.vdouble(0.784, 0.366, 0.57, 0.911, 0.298, 
            0.645, 0.51, 0.497, 0.932, 0.835, 
            0.968, 0.969, 0.923, 0.898, 0.98, 
            0.63, 0.971, 1.0, 0.515, 0.963, 
            0.986, 0.823, 0.879, 1.01, 0.931, 
            0.937, 1.05),
        cutdphiin = cms.vdouble(0.0404, 0.0499, 0.263, 0.042, 0.0484, 
            0.241, 0.242, 0.231, 0.286, 0.0552, 
            0.0338, 0.154, 0.0623, 0.0183, 0.0392, 
            0.0547, 0.0588, 0.00654, 0.042, 0.0217, 
            0.0885, 0.0445, 0.0141, 0.0234, 0.065, 
            0.0258, 0.0346),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 13.7, 13.2, 
            13.6, 14.2, 14.1, 13.9, 12.9, 
            14.9, 17.7)
    ),
    electronIDType = cms.string('classbased'),
    robusttightEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    electronVersion = cms.string(''),
    robusttightEleIDCutsV00 = cms.PSet(
        barrel = cms.vdouble(0.015, 0.0092, 0.02, 0.0025, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.018, 0.025, 0.02, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusttightEleIDCutsV01 = cms.PSet(
        barrel = cms.vdouble(0.01, 0.0099, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.01, 0.028, 0.02, 0.0066, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusttightEleIDCutsV02 = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusttightEleIDCutsV03 = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    reducedEndcapRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    verticesCollection = cms.InputTag("offlinePrimaryVerticesWithBS"),
    classbasedlooseEleIDCuts = cms.PSet(
        cutdetain = cms.vdouble(0.0137, 0.00678, 0.0241, 0.0187, 0.0161, 
            0.0224, 0.0252, 0.0308, 0.0273),
        cutiso_sum = cms.vdouble(33.0, 17.0, 17.9, 18.8, 8.55, 
            12.5, 17.6, 18.5, 2.98),
        cutip_gsf = cms.vdouble(0.0551, 0.0765, 0.143, 0.0874, 0.594, 
            0.37, 0.0913, 1.15, 0.231),
        cutip_gsfl = cms.vdouble(0.0186, 0.0759, 0.138, 0.0473, 0.62, 
            0.304, 0.109, 0.775, 0.0479),
        cuthoe = cms.vdouble(0.247, 0.137, 0.147, 0.371, 0.0588, 
            0.147, 0.52, 0.452, 0.404),
        cutiso_sumoetl = cms.vdouble(11.3, 9.05, 9.07, 9.94, 5.25, 
            6.15, 10.7, 10.8, 4.4),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 4.5, 3.5, 3.5),
        cuthoel = cms.vdouble(0.236, 0.126, 0.147, 0.375, 0.0392, 
            0.145, 0.365, 0.383, 0.384),
        cutdphiin = cms.vdouble(0.0897, 0.262, 0.353, 0.116, 0.357, 
            0.319, 0.342, 0.404, 0.336),
        cutseel = cms.vdouble(0.0164, 0.0118, 0.015, 0.0523, 0.0326, 
            0.0456, 0.0185, 0.0589, 0.0544),
        cutiso_sumoet = cms.vdouble(34.5, 12.7, 12.1, 19.9, 6.35, 
            8.85, 14.0, 10.5, 9.74),
        cutsee = cms.vdouble(0.0176, 0.0125, 0.0181, 0.0415, 0.0364, 
            0.0418, 0.0146, 0.0678, 0.133),
        cuteseedopcor = cms.vdouble(0.63, 0.82, 0.401, 0.718, 0.4, 
            0.458, 0.15, 0.664, 0.373),
        cutdphiinl = cms.vdouble(0.0747, 0.25, 0.356, 0.0956, 0.347, 
            0.326, 0.333, 0.647, 0.289),
        cutdetainl = cms.vdouble(0.0124, 0.00503, 0.0257, 0.0228, 0.0118, 
            0.0178, 0.0188, 0.14, 0.024)
    ),
    reducedBarrelRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    robusthighenergyEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV01 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.0275, 0.09, 0.007, -1, 
            -1, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV00 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.011, 0.09, 0.005, -1, 
            -1, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.1, 0.0275, 0.09, 0.007, -1, 
            -1, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV03 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV02 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    classbasedlooseEleIDCutsV00 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.05, 0.025, 0.053, 0.09, 0.07, 
            0.03, 0.092, 0.092),
        hOverE = cms.vdouble(0.115, 0.1, 0.055, 0.0, 0.145, 
            0.12, 0.15, 0.0),
        sigmaEtaEta = cms.vdouble(0.014, 0.012, 0.0115, 0.0, 0.0275, 
            0.0265, 0.0265, 0.0),
        deltaEtaIn = cms.vdouble(0.009, 0.0045, 0.0085, 0.0, 0.0105, 
            0.0068, 0.01, 0.0),
        eSeedOverPin = cms.vdouble(0.11, 0.91, 0.11, 0.0, 0.0, 
            0.85, 0.0, 0.0)
    ),
    classbasedlooseEleIDCutsV01 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.053, 0.0189, 0.059, 0.099, 0.0278, 
            0.0157, 0.042, 0.08),
        hOverE = cms.vdouble(0.076, 0.033, 0.07, 0.0, 0.083, 
            0.0148, 0.033, 0.0),
        sigmaEtaEta = cms.vdouble(0.0101, 0.0095, 0.0097, 0.0, 0.0271, 
            0.0267, 0.0259, 0.0),
        deltaEtaIn = cms.vdouble(0.0078, 0.00259, 0.0062, 0.0, 0.0078, 
            0.0061, 0.0061, 0.0),
        eSeedOverPin = cms.vdouble(0.3, 0.92, 0.211, 0.0, 0.42, 
            0.88, 0.68, 0.0)
    ),
    classbasedlooseEleIDCutsV02 = cms.PSet(
        cutisohcal = cms.vdouble(13.5, 9.93, 7.56, 14.8, 8.1, 
            10.8, 42.7, 20.1, 9.11, 10.4, 
            6.89, 5.59, 8.53, 9.59, 24.2, 
            2.78, 8.67, 0.288),
        cutmishits = cms.vdouble(5.5, 1.5, 5.5, 2.5, 2.5, 
            2.5, 3.5, 5.5, 0.5, 1.5, 
            2.5, 0.5, 1.5, 1.5, 0.5, 
            0.5, 0.5, 0.5),
        cuthoe = cms.vdouble(0.0887, 0.0934, 0.0949, 0.0986, 0.0431, 
            0.0878, 0.097, 0.0509, 0.098, 0.0991, 
            0.0321, 0.0928, 0.0663, 0.0717, 0.0966, 
            0.0758, 0.0149, 0.0131),
        cutdeta = cms.vdouble(0.00958, 0.00406, 0.0122, 0.0137, 0.00837, 
            0.0127, 0.011, 0.00336, 0.00977, 0.015, 
            0.00675, 0.0109, 0.014, 0.00508, 0.0109, 
            0.0146, 0.00506, 0.0127),
        cuteopin = cms.vdouble(0.878, 0.802, 0.814, 0.942, 0.735, 
            0.774, 0.829, 0.909, 0.829, 0.813, 
            0.86, 0.897, 0.817, 0.831, 0.818, 
            0.861, 0.787, 0.789),
        cutip = cms.vdouble(0.0246, 0.076, 0.0966, 0.0885, 0.441, 
            0.205, 0.0292, 0.0293, 0.0619, 0.0251, 
            0.159, 0.0815, 7.29, 0.0106, 5.76, 
            6.89, 1.27, 5.89),
        cutisotk = cms.vdouble(24.3, 8.45, 14.4, 27.8, 6.02, 
            10.5, 14.1, 10.2, 14.5, 19.1, 
            6.1, 14.1, 8.59, 8.33, 8.3, 
            8.93, 8.6, 16.0),
        cutsee = cms.vdouble(0.0172, 0.0115, 0.0143, 0.0344, 0.0295, 
            0.0304, 0.0145, 0.0108, 0.0128, 0.0347, 
            0.0307, 0.0316, 0.018, 0.011, 0.0132, 
            0.0349, 0.031, 0.0327),
        cutdphi = cms.vdouble(0.0372, 0.114, 0.118, 0.0488, 0.117, 
            0.119, 0.0606, 0.0548, 0.117, 0.07, 
            0.0355, 0.117, 0.088, 0.045, 0.118, 
            0.0919, 0.0236, 0.0515),
        cutisoecal = cms.vdouble(33.4, 28.1, 7.32, 27.4, 7.33, 
            21.7, 93.8, 102.0, 12.1, 26.0, 
            8.91, 10.0, 16.1, 31.3, 16.9, 
            15.4, 13.3, 37.7)
    ),
    classbasedlooseEleIDCutsV03 = cms.PSet(
        cutdetain = cms.vdouble(0.00989, 0.00484, 0.0146, 0.0146, 0.00902, 
            0.0172, 0.0137, 0.0477, 0.0275, 0.00967, 
            0.00377, 0.00924, 0.013, 0.00666, 0.0123, 
            0.0125, 0.0228, 0.0112, 0.0106, 0.0038, 
            0.00897, 0.0139, 0.00667, 0.0122, 0.0122, 
            0.0193, 0.00239),
        cutiso_sum = cms.vdouble(31.5, 10.3, 8.8, 11.0, 6.13, 
            6.94, 7.52, 9.0, 3.5, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0431, 0.0767, 0.139, 0.101, 0.149, 
            0.154, 0.932, 0.15, 0.124, 0.0238, 
            0.0467, 0.0759, 0.0369, 0.147, 0.0986, 
            0.0626, 0.195, 0.116, 0.0122, 0.0125, 
            0.0693, 0.0162, 0.089, 0.0673, 0.0467, 
            0.0651, 0.0221),
        cuthoe = cms.vdouble(0.166, 0.0771, 0.144, 0.37, 0.0497, 
            0.139, 0.401, 2.68, 0.516, 0.234, 
            0.0556, 0.144, 0.368, 0.031, 0.12, 
            0.602, 2.01, 1.05, 0.104, 0.063, 
            0.0565, 0.38, 0.0192, 0.0294, 0.537, 
            4.65, 1.87),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 2.5, 2.5, 1.5, 2.5, 
            1.5, 1.5, 1.5, 1.5, 0.5, 
            2.5, 2.5, 0.5, 2.5, 1.5, 
            0.5, 1.5, 1.5, 0.5, 2.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(28.9, 15.3, 12.0, 18.3, 7.17, 
            9.42, 11.0, 9.81, 3.94, 22.7, 
            15.9, 12.3, 17.0, 7.58, 8.89, 
            15.2, 12.7, 6.17, 20.8, 21.2, 
            17.2, 15.5, 9.37, 10.6, 19.8, 
            22.1, 15.6),
        cutdcotdist = cms.vdouble(0.0393, 0.0392, 0.0397, 0.0394, 0.0393, 
            0.039, 0.0378, 0.0388, 0.0382, 0.0385, 
            0.0167, 0.00325, 0.0394, 0.0387, 0.0388, 
            0.0227, 0.0258, 0.0127, 0.0298, 0.03, 
            0.00946, 0.039, 0.0231, 0.0278, 0.00162, 
            0.0367, 0.0199),
        cutsee = cms.vdouble(0.0175, 0.0127, 0.0177, 0.0373, 0.0314, 
            0.0329, 0.0157, 0.0409, 0.14, 0.0169, 
            0.0106, 0.0142, 0.0363, 0.0322, 0.0354, 
            0.0117, 0.0372, 28.2, 0.0171, 0.0113, 
            0.014, 0.0403, 0.0323, 0.0411, 0.0104, 
            0.0436, 0.0114),
        cuteseedopcor = cms.vdouble(0.78, 0.302, 0.483, 0.904, 0.168, 
            0.645, 0.108, 0.284, 0.324, 0.591, 
            0.286, 0.488, 0.813, 0.791, 0.672, 
            0.398, 0.834, 0.878, 0.515, 0.937, 
            0.806, 0.816, 0.85, 0.507, 0.367, 
            0.83, 0.648),
        cutdphiin = cms.vdouble(0.041, 0.275, 0.365, 0.047, 0.273, 
            0.296, 0.329, 0.465, 0.627, 0.0581, 
            0.0954, 0.327, 0.0702, 0.0582, 0.279, 
            0.117, 0.318, 0.246, 0.0821, 0.052, 
            0.292, 0.116, 0.0435, 0.312, 0.118, 
            0.296, 0.0459),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 12.0, 12.0, 
            12.0, 12.0, 12.0, 12.0, 12.0, 
            12.0, 12.5)
    ),
    classbasedlooseEleIDCutsV04 = cms.PSet(
        cutdetain = cms.vdouble(0.00989, 0.00484, 0.0146, 0.0146, 0.00902, 
            0.0172, 0.0137, 0.0477, 0.0275, 0.00967, 
            0.00377, 0.00924, 0.013, 0.00666, 0.0123, 
            0.0125, 0.0228, 0.0112, 0.0106, 0.0038, 
            0.00897, 0.0139, 0.00667, 0.0122, 0.0122, 
            0.0193, 0.00239),
        cutiso_sum = cms.vdouble(31.5, 10.3, 8.8, 11.0, 6.13, 
            6.94, 7.52, 9.0, 3.5, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0431, 0.0767, 0.139, 0.101, 0.149, 
            0.154, 0.932, 0.15, 0.124, 0.0238, 
            0.0467, 0.0759, 0.0369, 0.147, 0.0986, 
            0.0626, 0.195, 0.116, 0.0122, 0.0125, 
            0.0693, 0.0162, 0.089, 0.0673, 0.0467, 
            0.0651, 0.0221),
        cuthoe = cms.vdouble(0.166, 0.0771, 0.144, 0.37, 0.0497, 
            0.139, 0.401, 2.68, 0.516, 0.234, 
            0.0556, 0.144, 0.368, 0.031, 0.12, 
            0.602, 2.01, 1.05, 0.104, 0.063, 
            0.0565, 0.38, 0.0192, 0.0294, 0.537, 
            4.65, 1.87),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 2.5, 2.5, 1.5, 2.5, 
            1.5, 1.5, 1.5, 1.5, 0.5, 
            2.5, 2.5, 0.5, 2.5, 1.5, 
            0.5, 1.5, 1.5, 0.5, 2.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(28.9, 15.3, 12.0, 18.3, 7.17, 
            9.42, 11.0, 9.81, 3.94, 22.7, 
            15.9, 12.3, 17.0, 7.58, 8.89, 
            15.2, 12.7, 6.17, 20.8, 21.2, 
            17.2, 15.5, 9.37, 10.6, 19.8, 
            22.1, 15.6),
        cutdcotdist = cms.vdouble(0.0393, 0.0392, 0.0397, 0.0394, 0.0393, 
            0.039, 0.0378, 0.0388, 0.0382, 0.0385, 
            0.0167, 0.00325, 0.0394, 0.0387, 0.0388, 
            0.0227, 0.0258, 0.0127, 0.0298, 0.03, 
            0.00946, 0.039, 0.0231, 0.0278, 0.00162, 
            0.0367, 0.0199),
        cutsee = cms.vdouble(0.0175, 0.0127, 0.0177, 0.0373, 0.0314, 
            0.0329, 0.0157, 0.0409, 0.14, 0.0169, 
            0.0106, 0.0142, 0.0363, 0.0322, 0.0354, 
            0.0117, 0.0372, 28.2, 0.0171, 0.0113, 
            0.014, 0.0403, 0.0323, 0.0411, 0.0104, 
            0.0436, 0.0114),
        cuteseedopcor = cms.vdouble(0.78, 0.302, 0.483, 0.904, 0.168, 
            0.645, 0.108, 0.284, 0.324, 0.591, 
            0.286, 0.488, 0.813, 0.791, 0.672, 
            0.398, 0.834, 0.878, 0.515, 0.937, 
            0.806, 0.816, 0.85, 0.507, 0.367, 
            0.83, 0.648),
        cutdphiin = cms.vdouble(0.041, 0.275, 0.365, 0.047, 0.273, 
            0.296, 0.329, 0.465, 0.627, 0.0581, 
            0.0954, 0.327, 0.0702, 0.0582, 0.279, 
            0.117, 0.318, 0.246, 0.0821, 0.052, 
            0.292, 0.116, 0.0435, 0.312, 0.118, 
            0.296, 0.0459),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 12.0, 12.0, 
            12.0, 12.0, 12.0, 12.0, 12.0, 
            12.0, 12.5)
    ),
    classbasedlooseEleIDCutsV06 = cms.PSet(
        cutdetain = cms.vdouble(0.0137, 0.00678, 0.0241, 0.0187, 0.0161, 
            0.0224, 0.0252, 0.0308, 0.0273),
        cutiso_sum = cms.vdouble(33.0, 17.0, 17.9, 18.8, 8.55, 
            12.5, 17.6, 18.5, 2.98),
        cutip_gsf = cms.vdouble(0.0551, 0.0765, 0.143, 0.0874, 0.594, 
            0.37, 0.0913, 1.15, 0.231),
        cutip_gsfl = cms.vdouble(0.0186, 0.0759, 0.138, 0.0473, 0.62, 
            0.304, 0.109, 0.775, 0.0479),
        cuthoe = cms.vdouble(0.247, 0.137, 0.147, 0.371, 0.0588, 
            0.147, 0.52, 0.452, 0.404),
        cutiso_sumoetl = cms.vdouble(11.3, 9.05, 9.07, 9.94, 5.25, 
            6.15, 10.7, 10.8, 4.4),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 4.5, 3.5, 3.5),
        cuthoel = cms.vdouble(0.236, 0.126, 0.147, 0.375, 0.0392, 
            0.145, 0.365, 0.383, 0.384),
        cutdphiin = cms.vdouble(0.0897, 0.262, 0.353, 0.116, 0.357, 
            0.319, 0.342, 0.404, 0.336),
        cutseel = cms.vdouble(0.0164, 0.0118, 0.015, 0.0523, 0.0326, 
            0.0456, 0.0185, 0.0589, 0.0544),
        cutiso_sumoet = cms.vdouble(34.5, 12.7, 12.1, 19.9, 6.35, 
            8.85, 14.0, 10.5, 9.74),
        cutsee = cms.vdouble(0.0176, 0.0125, 0.0181, 0.0415, 0.0364, 
            0.0418, 0.0146, 0.0678, 0.133),
        cuteseedopcor = cms.vdouble(0.63, 0.82, 0.401, 0.718, 0.4, 
            0.458, 0.15, 0.664, 0.373),
        cutdphiinl = cms.vdouble(0.0747, 0.25, 0.356, 0.0956, 0.347, 
            0.326, 0.333, 0.647, 0.289),
        cutdetainl = cms.vdouble(0.0124, 0.00503, 0.0257, 0.0228, 0.0118, 
            0.0178, 0.0188, 0.14, 0.024)
    ),
    src = cms.InputTag("gsfElectrons"),
    robusttightEleIDCuts = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    classbasedtightEleIDCuts = cms.PSet(
        cutdetain = cms.vdouble(0.0116, 0.00449, 0.00938, 0.0184, 0.00678, 
            0.0109, 0.0252, 0.0268, 0.0139),
        cutiso_sum = cms.vdouble(15.5, 12.2, 12.2, 11.7, 7.16, 
            9.71, 8.66, 11.9, 2.98),
        cutip_gsf = cms.vdouble(0.0131, 0.0586, 0.0839, 0.0366, 0.452, 
            0.204, 0.0913, 0.0802, 0.0731),
        cutip_gsfl = cms.vdouble(0.0119, 0.0527, 0.0471, 0.0212, 0.233, 
            0.267, 0.109, 0.122, 0.0479),
        cuthoe = cms.vdouble(0.215, 0.0608, 0.147, 0.369, 0.0349, 
            0.102, 0.52, 0.422, 0.404),
        cutiso_sumoetl = cms.vdouble(6.21, 6.81, 5.3, 5.39, 2.73, 
            4.73, 4.84, 3.46, 3.73),
        cutfmishits = cms.vdouble(1.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 1.5, 2.5, 0.5),
        cuthoel = cms.vdouble(0.228, 0.0836, 0.143, 0.37, 0.0392, 
            0.0979, 0.3, 0.381, 0.339),
        cutdphiin = cms.vdouble(0.0897, 0.0993, 0.295, 0.0979, 0.151, 
            0.252, 0.341, 0.308, 0.328),
        cutseel = cms.vdouble(0.0132, 0.0117, 0.0112, 0.0387, 0.0281, 
            0.0287, 0.00987, 0.0296, 0.0544),
        cutiso_sumoet = cms.vdouble(11.9, 7.81, 6.28, 8.92, 4.65, 
            5.49, 9.36, 8.84, 5.94),
        cutsee = cms.vdouble(0.0145, 0.0116, 0.012, 0.039, 0.0297, 
            0.0311, 0.00987, 0.0347, 0.0917),
        cuteseedopcor = cms.vdouble(0.637, 0.943, 0.742, 0.748, 0.763, 
            0.631, 0.214, 0.873, 0.473),
        cutdphiinl = cms.vdouble(0.061, 0.14, 0.286, 0.0921, 0.197, 
            0.24, 0.333, 0.303, 0.258),
        cutdetainl = cms.vdouble(0.00816, 0.00401, 0.0081, 0.019, 0.00588, 
            0.00893, 0.0171, 0.0434, 0.0143)
    ),
    algorithm = cms.string('eIDCB'),
    robusthighenergyEleIDCuts = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCuts = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV02 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV03 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV00 = cms.PSet(
        barrel = cms.vdouble(0.115, 0.014, 0.09, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.15, 0.0275, 0.092, 0.0105, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV01 = cms.PSet(
        barrel = cms.vdouble(0.075, 0.0132, 0.058, 0.0077, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.083, 0.027, 0.042, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    additionalCategories = cms.bool(True),
    etBinning = cms.bool(True)
)


process.JetPlusTrackZSPCorJetAntiKt5 = cms.EDProducer("JetPlusTrackProducer",
    VectorialCorrection = cms.bool(True),
    ElectronIds = cms.InputTag("JPTeidTight"),
    UseMuons = cms.bool(True),
    Muons = cms.InputTag("muons"),
    UseTrackQuality = cms.bool(True),
    JetTracksAssociationAtCaloFace = cms.InputTag("ak5JetTracksAssociatorAtCaloFace"),
    LeakageMap = cms.string('CondFormats/JetMETObjects/data/CMSSW_362_TrackLeakage.txt'),
    UseOutOfConeTracks = cms.bool(True),
    UseInConeTracks = cms.bool(True),
    UseOutOfVertexTracks = cms.bool(True),
    UseResponseInVecCorr = cms.bool(False),
    ResponseMap = cms.string('CondFormats/JetMETObjects/data/CMSSW_31X_resptowers.txt'),
    EfficiencyMap = cms.string('CondFormats/JetMETObjects/data/CMSSW_362_TrackNonEff.txt'),
    UsePions = cms.bool(True),
    Electrons = cms.InputTag("gsfElectrons"),
    JetSplitMerge = cms.int32(2),
    DzVertexCut = cms.double(0.2),
    MaxJetEta = cms.double(3.0),
    Verbose = cms.bool(True),
    UseElectrons = cms.bool(True),
    JetTracksAssociationAtVertex = cms.InputTag("ak5JetTracksAssociatorAtVertex"),
    PtErrorQuality = cms.double(0.05),
    TrackQuality = cms.string('highPurity'),
    UseEfficiency = cms.bool(False),
    src = cms.InputTag("ak5CaloJets"),
    PU = cms.int32(-1),
    srcPVs = cms.InputTag("offlinePrimaryVertices"),
    FixedPU = cms.int32(0),
    alias = cms.untracked.string('JetPlusTrackZSPCorJetAntiKt5'),
    tagName = cms.vstring('ZSP_CMSSW361_Akt_05_PU0'),
    tagNameOffset = cms.vstring(),
    UseZSP = cms.bool(True)
)


process.JetPlusTrackZSPCorJetIcone5 = cms.EDProducer("JetPlusTrackProducer",
    VectorialCorrection = cms.bool(True),
    ElectronIds = cms.InputTag("JPTeidTight"),
    UseMuons = cms.bool(True),
    Muons = cms.InputTag("muons"),
    UseTrackQuality = cms.bool(True),
    JetTracksAssociationAtCaloFace = cms.InputTag("iterativeCone5JetTracksAssociatorAtCaloFace"),
    LeakageMap = cms.string('CondFormats/JetMETObjects/data/CMSSW_362_TrackLeakage.txt'),
    UseOutOfConeTracks = cms.bool(True),
    UseInConeTracks = cms.bool(True),
    UseOutOfVertexTracks = cms.bool(True),
    UseResponseInVecCorr = cms.bool(False),
    ResponseMap = cms.string('CondFormats/JetMETObjects/data/CMSSW_31X_resptowers.txt'),
    EfficiencyMap = cms.string('CondFormats/JetMETObjects/data/CMSSW_362_TrackNonEff.txt'),
    UsePions = cms.bool(True),
    Electrons = cms.InputTag("gsfElectrons"),
    JetSplitMerge = cms.int32(0),
    DzVertexCut = cms.double(0.2),
    MaxJetEta = cms.double(3.0),
    Verbose = cms.bool(True),
    UseElectrons = cms.bool(True),
    JetTracksAssociationAtVertex = cms.InputTag("iterativeCone5JetTracksAssociatorAtVertex"),
    PtErrorQuality = cms.double(0.05),
    TrackQuality = cms.string('highPurity'),
    UseEfficiency = cms.bool(False),
    src = cms.InputTag("iterativeCone5CaloJets"),
    PU = cms.int32(-1),
    srcPVs = cms.InputTag("offlinePrimaryVertices"),
    FixedPU = cms.int32(0),
    alias = cms.untracked.string('JetPlusTrackZSPCorJetIcone5'),
    tagName = cms.vstring('ZSP_CMSSW361_Akt_05_PU0'),
    tagNameOffset = cms.vstring(),
    UseZSP = cms.bool(True)
)


process.JetPlusTrackZSPCorJetSiscone5 = cms.EDProducer("JetPlusTrackProducer",
    VectorialCorrection = cms.bool(True),
    ElectronIds = cms.InputTag("JPTeidTight"),
    UseMuons = cms.bool(True),
    Muons = cms.InputTag("muons"),
    UseTrackQuality = cms.bool(True),
    JetTracksAssociationAtCaloFace = cms.InputTag("sisCone5JetTracksAssociatorAtCaloFace"),
    LeakageMap = cms.string('CondFormats/JetMETObjects/data/CMSSW_362_TrackLeakage.txt'),
    UseOutOfConeTracks = cms.bool(True),
    UseInConeTracks = cms.bool(True),
    UseOutOfVertexTracks = cms.bool(True),
    UseResponseInVecCorr = cms.bool(False),
    ResponseMap = cms.string('CondFormats/JetMETObjects/data/CMSSW_31X_resptowers.txt'),
    EfficiencyMap = cms.string('CondFormats/JetMETObjects/data/CMSSW_362_TrackNonEff.txt'),
    UsePions = cms.bool(True),
    Electrons = cms.InputTag("gsfElectrons"),
    JetSplitMerge = cms.int32(1),
    DzVertexCut = cms.double(0.2),
    MaxJetEta = cms.double(3.0),
    Verbose = cms.bool(True),
    UseElectrons = cms.bool(True),
    JetTracksAssociationAtVertex = cms.InputTag("sisCone5JetTracksAssociatorAtVertex"),
    PtErrorQuality = cms.double(0.05),
    TrackQuality = cms.string('highPurity'),
    UseEfficiency = cms.bool(False),
    src = cms.InputTag("ak5CaloJets"),
    PU = cms.int32(-1),
    srcPVs = cms.InputTag("offlinePrimaryVertices"),
    FixedPU = cms.int32(0),
    alias = cms.untracked.string('JetPlusTrackZSPCorJetSiscone5'),
    tagName = cms.vstring('ZSP_CMSSW361_Akt_05_PU0'),
    tagNameOffset = cms.vstring(),
    UseZSP = cms.bool(True)
)


process.TCTauJetPlusTrackZSPCorJetAntiKt5 = cms.EDProducer("JetPlusTrackProducer",
    VectorialCorrection = cms.bool(True),
    ElectronIds = cms.InputTag("JPTeidTight"),
    UseMuons = cms.bool(True),
    Muons = cms.InputTag("muons"),
    UseTrackQuality = cms.bool(True),
    JetTracksAssociationAtCaloFace = cms.InputTag("ak5JetTracksAssociatorAtCaloFace"),
    LeakageMap = cms.string('CondFormats/JetMETObjects/data/CMSSW_362_TrackLeakage.txt'),
    UseOutOfConeTracks = cms.bool(True),
    UseInConeTracks = cms.bool(True),
    UseOutOfVertexTracks = cms.bool(True),
    UseResponseInVecCorr = cms.bool(False),
    ResponseMap = cms.string('CondFormats/JetMETObjects/data/CMSSW_31X_resptowers.txt'),
    EfficiencyMap = cms.string('CondFormats/JetMETObjects/data/CMSSW_362_TrackNonEff.txt'),
    UsePions = cms.bool(True),
    Electrons = cms.InputTag("gsfElectrons"),
    JetSplitMerge = cms.int32(2),
    DzVertexCut = cms.double(0.2),
    MaxJetEta = cms.double(3.0),
    Verbose = cms.bool(True),
    UseElectrons = cms.bool(True),
    JetTracksAssociationAtVertex = cms.InputTag("ak5JetTracksAssociatorAtVertex"),
    PtErrorQuality = cms.double(0.05),
    TrackQuality = cms.string('highPurity'),
    UseEfficiency = cms.bool(False),
    src = cms.InputTag("ak5CaloJets"),
    PU = cms.int32(-1),
    srcPVs = cms.InputTag("offlinePrimaryVertices"),
    FixedPU = cms.int32(0),
    alias = cms.untracked.string('JetPlusTrackZSPCorJetAntiKt5'),
    tagName = cms.vstring('ZSP_CMSSW361_Akt_05_PU0'),
    tagNameOffset = cms.vstring(),
    UseZSP = cms.bool(False)
)


process.TauDecayModeCutMutliplexerPrototype = cms.EDProducer("PFTauDecayModeCutMultiplexer",
    RemapOutput = cms.bool(True),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducer"),
    computers = cms.VPSet(cms.PSet(
        applyIsolation = cms.bool(False),
        cut = cms.double(-10.0),
        computerName = cms.string('OneProngNoPiZero'),
        decayModeIndices = cms.vint32(0)
    ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(-10.0),
            computerName = cms.string('OneProngOnePiZero'),
            decayModeIndices = cms.vint32(1)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(-10.0),
            computerName = cms.string('OneProngTwoPiZero'),
            decayModeIndices = cms.vint32(2)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(-10.0),
            computerName = cms.string('ThreeProngNoPiZero'),
            decayModeIndices = cms.vint32(10)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(-10.0),
            computerName = cms.string('ThreeProngOnePiZero'),
            decayModeIndices = cms.vint32(11)
        )),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFinding")
        )
    ),
    PFTauDiscriminantToMultiplex = cms.InputTag("shrinkingConePFTauDiscriminationByTaNC"),
    PFTauDecayModeSrc = cms.InputTag("shrinkingConePFTauDecayModeIndexProducer")
)


process.ak3HiGenJets = cms.EDProducer("SubEventGenJetProducer",
    Active_Area_Repeats = cms.int32(5),
    src = cms.InputTag("hiGenParticlesForJets"),
    doAreaFastjet = cms.bool(True),
    doPVCorrection = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    doRhoFastjet = cms.bool(False),
    srcPVs = cms.InputTag(""),
    inputEtMin = cms.double(0.0),
    nSigmaPU = cms.double(1.0),
    radiusPU = cms.double(0.5),
    Rho_EtaMax = cms.double(4.5),
    GhostArea = cms.double(0.01),
    jetPtMin = cms.double(3.0),
    jetType = cms.string('GenJet'),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    jetAlgorithm = cms.string('AntiKt'),
    rParam = cms.double(0.3)
)


process.ak4HiGenJets = cms.EDProducer("SubEventGenJetProducer",
    Active_Area_Repeats = cms.int32(5),
    src = cms.InputTag("hiGenParticlesForJets"),
    doAreaFastjet = cms.bool(True),
    doPVCorrection = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    doRhoFastjet = cms.bool(False),
    srcPVs = cms.InputTag(""),
    inputEtMin = cms.double(0.0),
    nSigmaPU = cms.double(1.0),
    radiusPU = cms.double(0.5),
    Rho_EtaMax = cms.double(4.5),
    GhostArea = cms.double(0.01),
    jetPtMin = cms.double(3.0),
    jetType = cms.string('GenJet'),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    jetAlgorithm = cms.string('AntiKt'),
    rParam = cms.double(0.4)
)


process.ak5GenJets = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    src = cms.InputTag("genParticlesForJets"),
    doAreaFastjet = cms.bool(False),
    doPVCorrection = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    doRhoFastjet = cms.bool(False),
    srcPVs = cms.InputTag(""),
    inputEtMin = cms.double(0.0),
    nSigmaPU = cms.double(1.0),
    radiusPU = cms.double(0.5),
    Rho_EtaMax = cms.double(4.5),
    GhostArea = cms.double(0.01),
    jetPtMin = cms.double(3.0),
    jetType = cms.string('GenJet'),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    jetAlgorithm = cms.string('AntiKt'),
    rParam = cms.double(0.5)
)


process.ak5GenJetsNoMuNoNu = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('AntiKt'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.5),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJetsNoMuNoNu"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.ak5GenJetsNoNu = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('AntiKt'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.5),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJetsNoNu"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.ak5HiGenJets = cms.EDProducer("SubEventGenJetProducer",
    Active_Area_Repeats = cms.int32(5),
    src = cms.InputTag("hiGenParticlesForJets"),
    doAreaFastjet = cms.bool(True),
    doPVCorrection = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    doRhoFastjet = cms.bool(False),
    srcPVs = cms.InputTag(""),
    inputEtMin = cms.double(0.0),
    nSigmaPU = cms.double(1.0),
    radiusPU = cms.double(0.5),
    Rho_EtaMax = cms.double(4.5),
    GhostArea = cms.double(0.01),
    jetPtMin = cms.double(3.0),
    jetType = cms.string('GenJet'),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    jetAlgorithm = cms.string('AntiKt'),
    rParam = cms.double(0.5)
)


process.ak5JetExtender = cms.EDProducer("JetExtender",
    jets = cms.InputTag("ak5CaloJets"),
    jet2TracksAtCALO = cms.InputTag("ak5JetTracksAssociatorAtCaloFace"),
    jet2TracksAtVX = cms.InputTag("ak5JetTracksAssociatorAtVertex"),
    coneSize = cms.double(0.5)
)


process.ak5JetID = cms.EDProducer("JetIDProducer",
    eeRecHitsColl = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    hbheRecHitsColl = cms.InputTag("hbhereco"),
    rpcRecHits = cms.InputTag("rpcRecHits"),
    hoRecHitsColl = cms.InputTag("horeco"),
    ebRecHitsColl = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    hfRecHitsColl = cms.InputTag("hfreco"),
    useRecHits = cms.bool(True),
    src = cms.InputTag("ak5CaloJets")
)


process.ak5JetTracksAssociatorAtCaloFace = cms.EDProducer("JetTracksAssociatorAtCaloFace",
    trackQuality = cms.string('goodIterative'),
    tracks = cms.InputTag("generalTracks"),
    coneSize = cms.double(0.5),
    extrapolations = cms.InputTag("trackExtrapolator"),
    jets = cms.InputTag("ak5CaloJets")
)


process.ak5JetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
    tracks = cms.InputTag("generalTracks"),
    coneSize = cms.double(0.5),
    jets = cms.InputTag("ak5CaloJets")
)


process.ak5PFJetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
    jets = cms.InputTag("ak5PFJets"),
    tracks = cms.InputTag("generalTracks"),
    coneSize = cms.double(0.5)
)


process.ak7GenJets = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('AntiKt'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.7),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJets"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.ak7GenJetsNoMuNoNu = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('AntiKt'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.7),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJetsNoMuNoNu"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.ak7GenJetsNoNu = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('AntiKt'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.7),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJetsNoNu"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.ak7HiGenJets = cms.EDProducer("SubEventGenJetProducer",
    Active_Area_Repeats = cms.int32(5),
    src = cms.InputTag("hiGenParticlesForJets"),
    doAreaFastjet = cms.bool(True),
    doPVCorrection = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    doRhoFastjet = cms.bool(False),
    srcPVs = cms.InputTag(""),
    inputEtMin = cms.double(0.0),
    nSigmaPU = cms.double(1.0),
    radiusPU = cms.double(0.5),
    Rho_EtaMax = cms.double(4.5),
    GhostArea = cms.double(0.01),
    jetPtMin = cms.double(3.0),
    jetType = cms.string('GenJet'),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    jetAlgorithm = cms.string('AntiKt'),
    rParam = cms.double(0.7)
)


process.ca4GenJets = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    src = cms.InputTag("genParticlesForJets"),
    doAreaFastjet = cms.bool(False),
    doPVCorrection = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    doRhoFastjet = cms.bool(False),
    srcPVs = cms.InputTag(""),
    inputEtMin = cms.double(0.0),
    nSigmaPU = cms.double(1.0),
    radiusPU = cms.double(0.5),
    Rho_EtaMax = cms.double(4.5),
    GhostArea = cms.double(0.01),
    jetPtMin = cms.double(3.0),
    jetType = cms.string('GenJet'),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    jetAlgorithm = cms.string('CambridgeAachen'),
    rParam = cms.double(0.4)
)


process.ca4GenJetsNoMuNoNu = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('CambridgeAachen'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.4),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJetsNoMuNoNu"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.ca4GenJetsNoNu = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('CambridgeAachen'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.4),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJetsNoNu"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.ca6GenJets = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('CambridgeAachen'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.6),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJets"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.ca6GenJetsNoMuNoNu = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('CambridgeAachen'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.6),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJetsNoMuNoNu"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.ca6GenJetsNoNu = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('CambridgeAachen'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.6),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJetsNoNu"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.caloRecoTauDiscriminationAgainstElectron = cms.EDProducer("CaloRecoTauDiscriminationAgainstElectron",
    CaloTauProducer = cms.InputTag("caloRecoTauProducer"),
    leadTrack_HCAL3x3hitsEtSumOverPt_minvalue = cms.double(0.1),
    maxleadTrackHCAL3x3hottesthitDEta = cms.double(0.1),
    ApplyCut_maxleadTrackHCAL3x3hottesthitDEta = cms.bool(False),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("caloRecoTauDiscriminationByLeadingTrackFinding")
        )
    ),
    ApplyCut_leadTrackavoidsECALcrack = cms.bool(True)
)


process.caloRecoTauDiscriminationAgainstMuon = cms.EDProducer("CaloRecoTauDiscriminationAgainstMuon",
    CaloTauProducer = cms.InputTag("caloRecoTauProducer"),
    muonSource = cms.InputTag("muons"),
    segmCompCoefficient = cms.double(0.5),
    muonCompCut = cms.double(0.0),
    dRmatch = cms.double(0.5),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("caloRecoTauDiscriminationByLeadingTrackFinding")
        )
    ),
    discriminatorOption = cms.string('noSegMatch'),
    caloCompCoefficient = cms.double(0.5)
)


process.caloRecoTauDiscriminationByECALIsolation = cms.EDProducer("CaloRecoTauDiscriminationByIsolation",
    ECALisolAnnulus_maximumSumEtCut = cms.double(1.5),
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    TrackerIsolAnnulus_maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("caloRecoTauDiscriminationByLeadingTrackFinding")
        )
    ),
    ApplyDiscriminationByTrackerIsolation = cms.bool(False),
    CaloTauProducer = cms.InputTag("caloRecoTauProducer")
)


process.caloRecoTauDiscriminationByIsolation = cms.EDProducer("CaloRecoTauDiscriminationByIsolation",
    ECALisolAnnulus_maximumSumEtCut = cms.double(1.5),
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    TrackerIsolAnnulus_maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("caloRecoTauDiscriminationByLeadingTrackFinding")
        )
    ),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    CaloTauProducer = cms.InputTag("caloRecoTauProducer")
)


process.caloRecoTauDiscriminationByLeadingTrackFinding = cms.EDProducer("CaloRecoTauDiscriminationByLeadingTrackPtCut",
    MinPtLeadingTrack = cms.double(0.0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    CaloTauProducer = cms.InputTag("caloRecoTauProducer")
)


process.caloRecoTauDiscriminationByLeadingTrackPtCut = cms.EDProducer("CaloRecoTauDiscriminationByLeadingTrackPtCut",
    MinPtLeadingTrack = cms.double(5.0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("caloRecoTauDiscriminationByLeadingTrackFinding")
        )
    ),
    CaloTauProducer = cms.InputTag("caloRecoTauProducer")
)


process.caloRecoTauDiscriminationByTrackIsolation = cms.EDProducer("CaloRecoTauDiscriminationByIsolation",
    ECALisolAnnulus_maximumSumEtCut = cms.double(1.5),
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    TrackerIsolAnnulus_maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("caloRecoTauDiscriminationByLeadingTrackFinding")
        )
    ),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    CaloTauProducer = cms.InputTag("caloRecoTauProducer")
)


process.caloRecoTauProducer = cms.EDProducer("TCRecoTauProducer",
    tkminTrackerHitsn = cms.int32(5),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    EtCaloOverTrackMin = cms.double(-0.9),
    EERecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    tkmaxChi2 = cms.double(100.0),
    EtHcalOverTrackMin = cms.double(-0.3),
    CaloRecoTauProducer = cms.InputTag("JPTCaloRecoTauProducer"),
    EBRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    tkminPixelHitsn = cms.int32(0),
    MatchingConeSize = cms.double(0.1),
    TrackCollection = cms.InputTag("generalTracks"),
    HBHERecHitCollection = cms.InputTag("hbhereco"),
    EtCaloOverTrackMax = cms.double(0.0),
    TrackAssociatorParameters = cms.PSet(
        muonMaxDistanceSigmaX = cms.double(0.0),
        muonMaxDistanceSigmaY = cms.double(0.0),
        CSCSegmentCollectionLabel = cms.InputTag("cscSegments"),
        dRHcal = cms.double(9999.0),
        dRPreshowerPreselection = cms.double(0.2),
        CaloTowerCollectionLabel = cms.InputTag("towerMaker"),
        useEcal = cms.bool(True),
        dREcal = cms.double(9999.0),
        dREcalPreselection = cms.double(0.05),
        HORecHitCollectionLabel = cms.InputTag("horeco"),
        dRMuon = cms.double(9999.0),
        propagateAllDirections = cms.bool(True),
        muonMaxDistanceX = cms.double(5.0),
        muonMaxDistanceY = cms.double(5.0),
        useHO = cms.bool(True),
        trajectoryUncertaintyTolerance = cms.double(-1.0),
        usePreshower = cms.bool(False),
        DTRecSegment4DCollectionLabel = cms.InputTag("dt4DSegments"),
        EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
        dRHcalPreselection = cms.double(0.2),
        useMuon = cms.bool(True),
        useCalo = cms.bool(False),
        accountForTrajectoryChangeCalo = cms.bool(False),
        EBRecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
        dRMuonPreselection = cms.double(0.2),
        truthMatch = cms.bool(False),
        HBHERecHitCollectionLabel = cms.InputTag("hbhereco"),
        useHcal = cms.bool(True)
    ),
    SignalConeSize = cms.double(0.2),
    EcalConeSize = cms.double(0.5),
    Track_minPt = cms.double(1.0),
    EtHcalOverTrackMax = cms.double(1.0),
    HFRecHitCollection = cms.InputTag("hfreco"),
    DropRejectedJets = cms.untracked.bool(False),
    HORecHitCollection = cms.InputTag("horeco"),
    DropCaloJets = cms.untracked.bool(False),
    tkmaxipt = cms.double(0.1)
)


process.caloRecoTauTagInfoProducer = cms.EDProducer("CaloRecoTauTagInfoProducer",
    tkminTrackerHitsn = cms.int32(3),
    tkminPixelHitsn = cms.int32(0),
    ECALBasicClusterpropagTrack_matchingDRConeSize = cms.double(0.015),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    tkminPt = cms.double(0.5),
    BarrelBasicClustersSource = cms.InputTag("hybridSuperClusters","hybridBarrelBasicClusters"),
    UsePVconstraint = cms.bool(True),
    tkmaxChi2 = cms.double(100.0),
    EndcapBasicClustersSource = cms.InputTag("multi5x5BasicClusters","multi5x5EndcapBasicClusters"),
    smearedPVsigmaY = cms.double(0.0015),
    smearedPVsigmaX = cms.double(0.0015),
    ECALBasicClusterminE = cms.double(1.0),
    smearedPVsigmaZ = cms.double(0.005),
    tkQuality = cms.string('highPurity'),
    tkPVmaxDZ = cms.double(0.2),
    UseTrackQuality = cms.bool(True),
    ECALBasicClustersAroundCaloJet_DRConeSize = cms.double(0.5),
    tkmaxipt = cms.double(0.03),
    CaloJetTracksAssociatorProducer = cms.InputTag("JPTAntiKt5JetTracksAssociatorAtVertex")
)


process.cleanPatElectrons = cms.EDProducer("PATElectronCleaner",
    finalCut = cms.string(''),
    src = cms.InputTag("selectedPatElectrons"),
    checkOverlaps = cms.PSet(
        muons = cms.PSet(
            src = cms.InputTag("cleanPatMuons"),
            deltaR = cms.double(0.3),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False)
        )
    ),
    preselection = cms.string('')
)


process.cleanPatElectronsPF = cms.EDProducer("PATElectronCleaner",
    finalCut = cms.string(''),
    src = cms.InputTag("selectedPatElectronsPF"),
    checkOverlaps = cms.PSet(
        muons = cms.PSet(
            src = cms.InputTag("cleanPatMuonsPF"),
            deltaR = cms.double(0.3),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False)
        )
    ),
    preselection = cms.string('')
)


process.cleanPatElectronsTriggerMatch = cms.EDProducer("PATTriggerMatchElectronEmbedder",
    src = cms.InputTag("cleanPatElectrons"),
    matches = cms.VInputTag("patElectronMatch")
)


process.cleanPatJets = cms.EDProducer("PATJetCleaner",
    finalCut = cms.string(''),
    src = cms.InputTag("selectedPatJets"),
    checkOverlaps = cms.PSet(
        taus = cms.PSet(
            src = cms.InputTag("cleanPatTaus"),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False)
        ),
        photons = cms.PSet(
            src = cms.InputTag("cleanPatPhotons"),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False)
        ),
        electrons = cms.PSet(
            src = cms.InputTag("cleanPatElectrons"),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False)
        ),
        muons = cms.PSet(
            src = cms.InputTag("cleanPatMuons"),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False)
        ),
        tkIsoElectrons = cms.PSet(
            src = cms.InputTag("cleanPatElectrons"),
            deltaR = cms.double(0.3),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string('pt > 10 && trackIso < 3'),
            requireNoOverlaps = cms.bool(False)
        )
    ),
    preselection = cms.string('')
)


process.cleanPatJetsAK5Calo = cms.EDProducer("PATJetCleaner",
    finalCut = cms.string(''),
    src = cms.InputTag("selectedPatJets"),
    checkOverlaps = cms.PSet(
        taus = cms.PSet(
            src = cms.InputTag("cleanPatTaus"),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False)
        ),
        photons = cms.PSet(
            src = cms.InputTag("cleanPatPhotons"),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False)
        ),
        electrons = cms.PSet(
            src = cms.InputTag("cleanPatElectrons"),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False)
        ),
        muons = cms.PSet(
            src = cms.InputTag("cleanPatMuons"),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False)
        ),
        tkIsoElectrons = cms.PSet(
            src = cms.InputTag("cleanPatElectrons"),
            deltaR = cms.double(0.3),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string('pt > 10 && trackIso < 3'),
            requireNoOverlaps = cms.bool(False)
        )
    ),
    preselection = cms.string('')
)


process.cleanPatJetsAK5CaloTriggerMatch = cms.EDProducer("PATTriggerMatchJetEmbedder",
    src = cms.InputTag("cleanPatJetsAK5Calo"),
    matches = cms.VInputTag("patJetMatchAK5Calo")
)


process.cleanPatJetsAK5JPT = cms.EDProducer("PATJetCleaner",
    finalCut = cms.string(''),
    src = cms.InputTag("selectedPatJetsAK5JPT"),
    checkOverlaps = cms.PSet(
        taus = cms.PSet(
            src = cms.InputTag("cleanPatTaus"),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False)
        ),
        photons = cms.PSet(
            src = cms.InputTag("cleanPatPhotons"),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False)
        ),
        electrons = cms.PSet(
            src = cms.InputTag("cleanPatElectrons"),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False)
        ),
        muons = cms.PSet(
            src = cms.InputTag("cleanPatMuons"),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False)
        ),
        tkIsoElectrons = cms.PSet(
            src = cms.InputTag("cleanPatElectrons"),
            deltaR = cms.double(0.3),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string('pt > 10 && trackIso < 3'),
            requireNoOverlaps = cms.bool(False)
        )
    ),
    preselection = cms.string('')
)


process.cleanPatJetsAK5JPTTriggerMatch = cms.EDProducer("PATTriggerMatchJetEmbedder",
    src = cms.InputTag("cleanPatJetsAK5JPT"),
    matches = cms.VInputTag("patJetMatchAK5JPT")
)


process.cleanPatJetsAK5PF = cms.EDProducer("PATJetCleaner",
    finalCut = cms.string(''),
    src = cms.InputTag("selectedPatJetsAK5PF"),
    checkOverlaps = cms.PSet(
        taus = cms.PSet(
            src = cms.InputTag("cleanPatTaus"),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False)
        ),
        photons = cms.PSet(
            src = cms.InputTag("cleanPatPhotons"),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False)
        ),
        electrons = cms.PSet(
            src = cms.InputTag("cleanPatElectrons"),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False)
        ),
        muons = cms.PSet(
            src = cms.InputTag("cleanPatMuons"),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False)
        ),
        tkIsoElectrons = cms.PSet(
            src = cms.InputTag("cleanPatElectrons"),
            deltaR = cms.double(0.3),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string('pt > 10 && trackIso < 3'),
            requireNoOverlaps = cms.bool(False)
        )
    ),
    preselection = cms.string('')
)


process.cleanPatJetsAK5PFTriggerMatch = cms.EDProducer("PATTriggerMatchJetEmbedder",
    src = cms.InputTag("cleanPatJetsAK5PF"),
    matches = cms.VInputTag("patJetMatchAK5PF")
)


process.cleanPatJetsPF = cms.EDProducer("PATJetCleaner",
    finalCut = cms.string(''),
    src = cms.InputTag("selectedPatJetsPF"),
    checkOverlaps = cms.PSet(
        taus = cms.PSet(
            src = cms.InputTag("cleanPatTausPF"),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False)
        ),
        photons = cms.PSet(
            src = cms.InputTag("cleanPatPhotonsPF"),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False)
        ),
        electrons = cms.PSet(
            src = cms.InputTag("cleanPatElectronsPF"),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False)
        ),
        muons = cms.PSet(
            src = cms.InputTag("cleanPatMuonsPF"),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False)
        ),
        tkIsoElectrons = cms.PSet(
            src = cms.InputTag("cleanPatElectronsPF"),
            deltaR = cms.double(0.3),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string('pt > 10 && trackIso < 3'),
            requireNoOverlaps = cms.bool(False)
        )
    ),
    preselection = cms.string('')
)


process.cleanPatMuons = cms.EDProducer("PATMuonCleaner",
    finalCut = cms.string(''),
    src = cms.InputTag("selectedPatMuons"),
    checkOverlaps = cms.PSet(

    ),
    preselection = cms.string('')
)


process.cleanPatMuonsPF = cms.EDProducer("PATMuonCleaner",
    finalCut = cms.string(''),
    src = cms.InputTag("selectedPatMuonsPF"),
    checkOverlaps = cms.PSet(

    ),
    preselection = cms.string('')
)


process.cleanPatMuonsTriggerMatch = cms.EDProducer("PATTriggerMatchMuonEmbedder",
    src = cms.InputTag("cleanPatMuons"),
    matches = cms.VInputTag("patMuonMatch")
)


process.cleanPatPhotons = cms.EDProducer("PATPhotonCleaner",
    finalCut = cms.string(''),
    src = cms.InputTag("selectedPatPhotons"),
    checkOverlaps = cms.PSet(
        electrons = cms.PSet(
            src = cms.InputTag("cleanPatElectrons"),
            requireNoOverlaps = cms.bool(False),
            algorithm = cms.string('bySuperClusterSeed')
        )
    ),
    preselection = cms.string('')
)


process.cleanPatPhotonsPF = cms.EDProducer("PATPhotonCleaner",
    finalCut = cms.string(''),
    src = cms.InputTag("selectedPatPhotonsPF"),
    checkOverlaps = cms.PSet(
        electrons = cms.PSet(
            src = cms.InputTag("cleanPatElectronsPF"),
            requireNoOverlaps = cms.bool(False),
            algorithm = cms.string('bySuperClusterSeed')
        )
    ),
    preselection = cms.string('')
)


process.cleanPatPhotonsTriggerMatch = cms.EDProducer("PATTriggerMatchPhotonEmbedder",
    src = cms.InputTag("cleanPatPhotons"),
    matches = cms.VInputTag("patPhotonMatch")
)


process.cleanPatTaus = cms.EDProducer("PATTauCleaner",
    finalCut = cms.string(''),
    src = cms.InputTag("selectedPatTaus"),
    checkOverlaps = cms.PSet(
        muons = cms.PSet(
            src = cms.InputTag("cleanPatMuons"),
            deltaR = cms.double(0.3),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False)
        ),
        electrons = cms.PSet(
            src = cms.InputTag("cleanPatElectrons"),
            deltaR = cms.double(0.3),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False)
        )
    ),
    preselection = cms.string('tauID("leadingTrackFinding") > 0.5 & tauID("leadingPionPtCut") > 0.5 & tauID("byIsolationUsingLeadingPion") > 0.5 & tauID("againstMuon") > 0.5 & tauID("againstElectron") > 0.5 & (signalPFChargedHadrCands.size() = 1 | signalPFChargedHadrCands.size() = 3)')
)


process.cleanPatTausPF = cms.EDProducer("PATTauCleaner",
    finalCut = cms.string(''),
    src = cms.InputTag("selectedPatTausPF"),
    checkOverlaps = cms.PSet(
        muons = cms.PSet(
            src = cms.InputTag("cleanPatMuonsPF"),
            deltaR = cms.double(0.3),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False)
        ),
        electrons = cms.PSet(
            src = cms.InputTag("cleanPatElectronsPF"),
            deltaR = cms.double(0.3),
            pairCut = cms.string(''),
            checkRecoComponents = cms.bool(False),
            algorithm = cms.string('byDeltaR'),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False)
        )
    ),
    preselection = cms.string('tauID("leadingTrackFinding") > 0.5 & tauID("leadingPionPtCut") > 0.5 & tauID("byIsolationUsingLeadingPion") > 0.5 & tauID("againstMuon") > 0.5 & tauID("againstElectron") > 0.5 & (signalPFChargedHadrCands.size() = 1 | signalPFChargedHadrCands.size() = 3)')
)


process.cleanPatTausTriggerMatch = cms.EDProducer("PATTriggerMatchTauEmbedder",
    src = cms.InputTag("cleanPatTaus"),
    matches = cms.VInputTag("patTauMatch")
)


process.combinedMVABJetTags = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('combinedMVA'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfos"), cms.InputTag("secondaryVertexTagInfos"), cms.InputTag("softMuonTagInfos"), cms.InputTag("softElectronTagInfos"))
)


process.combinedSecondaryVertexBJetTags = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('combinedSecondaryVertex'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfos"), cms.InputTag("secondaryVertexTagInfos"))
)


process.combinedSecondaryVertexBJetTagsAK5JPT = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('combinedSecondaryVertex'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAK5JPT"), cms.InputTag("secondaryVertexTagInfosAK5JPT"))
)


process.combinedSecondaryVertexBJetTagsAK5PF = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('combinedSecondaryVertex'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAK5PF"), cms.InputTag("secondaryVertexTagInfosAK5PF"))
)


process.combinedSecondaryVertexBJetTagsAOD = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('combinedSecondaryVertex'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAOD"), cms.InputTag("secondaryVertexTagInfosAOD"))
)


process.combinedSecondaryVertexBJetTagsAODPF = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('combinedSecondaryVertex'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAODPF"), cms.InputTag("secondaryVertexTagInfosAODPF"))
)


process.combinedSecondaryVertexMVABJetTags = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('combinedSecondaryVertexMVA'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfos"), cms.InputTag("secondaryVertexTagInfos"))
)


process.combinedSecondaryVertexMVABJetTagsAK5JPT = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('combinedSecondaryVertexMVA'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAK5JPT"), cms.InputTag("secondaryVertexTagInfosAK5JPT"))
)


process.combinedSecondaryVertexMVABJetTagsAK5PF = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('combinedSecondaryVertexMVA'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAK5PF"), cms.InputTag("secondaryVertexTagInfosAK5PF"))
)


process.combinedSecondaryVertexMVABJetTagsAOD = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('combinedSecondaryVertexMVA'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAOD"), cms.InputTag("secondaryVertexTagInfosAOD"))
)


process.combinedSecondaryVertexMVABJetTagsAODPF = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('combinedSecondaryVertexMVA'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAODPF"), cms.InputTag("secondaryVertexTagInfosAODPF"))
)


process.corMetGlobalMuons = cms.EDProducer("MuonMET",
    muonMETDepositValueMapInputTag = cms.InputTag("muonMETValueMapProducer","muCorrData"),
    metTypeInputTag = cms.InputTag("CaloMET"),
    muonsInputTag = cms.InputTag("muons"),
    uncorMETInputTag = cms.InputTag("met")
)


process.eidCutBasedExt = cms.EDProducer("EleIdCutBasedExtProducer",
    electronQuality = cms.string('loose'),
    classbasedtightEleIDCutsV02 = cms.PSet(
        cutisohcal = cms.vdouble(10.9, 7.01, 8.75, 3.51, 7.75, 
            1.62, 11.6, 9.9, 4.97, 5.33, 
            3.18, 2.32, 0.164, 5.46, 12.0, 
            0.00604, 4.1, 0.000628),
        cutmishits = cms.vdouble(5.5, 1.5, 0.5, 1.5, 2.5, 
            0.5, 3.5, 5.5, 0.5, 0.5, 
            0.5, 0.5, 0.5, 1.5, 0.5, 
            0.5, 0.5, 0.5),
        cuthoe = cms.vdouble(0.0871, 0.0289, 0.0783, 0.0946, 0.0245, 
            0.0363, 0.0671, 0.048, 0.0614, 0.0924, 
            0.0158, 0.049, 0.0382, 0.0915, 0.0451, 
            0.0452, 0.00196, 0.0043),
        cutdeta = cms.vdouble(0.00915, 0.00302, 0.0061, 0.0135, 0.00565, 
            0.00793, 0.0102, 0.00266, 0.0106, 0.00903, 
            0.00766, 0.00723, 0.0116, 0.00203, 0.00659, 
            0.0148, 0.00555, 0.0128),
        cuteopin = cms.vdouble(0.878, 0.859, 0.874, 0.944, 0.737, 
            0.773, 0.86, 0.967, 0.917, 0.812, 
            0.915, 1.01, 0.847, 0.953, 0.979, 
            0.841, 0.771, 1.09),
        cutip = cms.vdouble(0.0239, 0.027, 0.0768, 0.0231, 0.178, 
            0.0957, 0.0102, 0.0168, 0.043, 0.0166, 
            0.0594, 0.0308, 2.1, 0.00527, 3.17, 
            4.91, 0.769, 5.9),
        cutisotk = cms.vdouble(6.53, 4.6, 6.0, 8.63, 3.11, 
            7.77, 5.42, 4.81, 4.06, 6.47, 
            2.8, 3.45, 5.29, 5.18, 15.4, 
            5.38, 4.47, 0.0347),
        cutsee = cms.vdouble(0.0131, 0.0106, 0.0115, 0.0306, 0.028, 
            0.0293, 0.0131, 0.0106, 0.0115, 0.0317, 
            0.029, 0.0289, 0.0142, 0.0106, 0.0103, 
            0.035, 0.0296, 0.0333),
        cutdphi = cms.vdouble(0.0369, 0.0307, 0.117, 0.0475, 0.0216, 
            0.117, 0.0372, 0.0246, 0.0426, 0.0612, 
            0.0142, 0.039, 0.0737, 0.0566, 0.0359, 
            0.0187, 0.012, 0.0358),
        cutisoecal = cms.vdouble(20.0, 27.2, 4.48, 13.5, 4.56, 
            3.19, 12.2, 13.1, 7.42, 7.67, 
            4.12, 4.85, 10.1, 12.4, 11.1, 
            11.0, 10.6, 13.4)
    ),
    classbasedtightEleIDCutsV03 = cms.PSet(
        cutdetain = cms.vdouble(0.00811, 0.00341, 0.00633, 0.0103, 0.00667, 
            0.01, 0.0106, 0.0145, 0.0163, 0.0076, 
            0.00259, 0.00511, 0.00941, 0.0043, 0.00857, 
            0.012, 0.0169, 0.00172, 0.00861, 0.00362, 
            0.00601, 0.00925, 0.00489, 0.00832, 0.0119, 
            0.0169, 0.000996),
        cutiso_sum = cms.vdouble(11.8, 8.31, 6.26, 6.18, 3.28, 
            4.38, 4.17, 5.4, 1.57, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0213, 0.0422, 0.0632, 0.0361, 0.073, 
            0.126, 0.171, 0.119, 0.0372, 0.0131, 
            0.0146, 0.0564, 0.0152, 0.0222, 0.0268, 
            0.0314, 0.0884, 0.00374, 0.00852, 0.00761, 
            0.0143, 0.0106, 0.0127, 0.0119, 0.0123, 
            0.0235, 0.00363),
        cuthoe = cms.vdouble(0.0783, 0.0387, 0.105, 0.118, 0.0227, 
            0.062, 0.13, 2.47, 0.38, 0.0888, 
            0.0503, 0.0955, 0.0741, 0.015, 0.03, 
            0.589, 1.13, 0.612, 0.0494, 0.0461, 
            0.0292, 0.0369, 0.0113, 0.0145, 0.124, 
            2.05, 0.61),
        cutfmishits = cms.vdouble(2.5, 1.5, 1.5, 1.5, 1.5, 
            0.5, 2.5, 0.5, 0.5, 2.5, 
            1.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5, -0.5, 2.5, 1.5, 
            0.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(13.7, 11.6, 7.14, 9.98, 3.52, 
            4.87, 6.24, 7.96, 2.53, 11.2, 
            11.9, 7.88, 8.16, 5.58, 5.03, 
            11.4, 8.15, 5.79, 10.4, 11.1, 
            10.4, 7.47, 5.08, 5.9, 11.8, 
            14.1, 11.7),
        cutdcotdist = cms.vdouble(0.0393, 0.0256, 0.00691, 0.0394, 0.0386, 
            0.039, 0.0325, 0.0384, 0.0382, 0.0245, 
            0.000281, 5.46e-05, 0.0342, 0.0232, 0.00107, 
            0.0178, 0.0193, 0.000758, 0.000108, 0.0248, 
            0.000458, 0.0129, 0.00119, 0.0182, 4.53e-05, 
            0.0189, 0.000928),
        cutsee = cms.vdouble(0.0143, 0.0105, 0.0123, 0.0324, 0.0307, 
            0.0301, 0.0109, 0.027, 0.0292, 0.0133, 
            0.0104, 0.0116, 0.0332, 0.0296, 0.031, 
            0.00981, 0.0307, 0.072, 0.0149, 0.0105, 
            0.011, 0.0342, 0.0307, 0.0303, 0.00954, 
            0.0265, 0.0101),
        cuteseedopcor = cms.vdouble(0.784, 0.366, 0.57, 0.911, 0.298, 
            0.645, 0.51, 0.497, 0.932, 0.835, 
            0.968, 0.969, 0.923, 0.898, 0.98, 
            0.63, 0.971, 1.0, 0.515, 0.963, 
            0.986, 0.823, 0.879, 1.01, 0.931, 
            0.937, 1.05),
        cutdphiin = cms.vdouble(0.0404, 0.0499, 0.263, 0.042, 0.0484, 
            0.241, 0.242, 0.231, 0.286, 0.0552, 
            0.0338, 0.154, 0.0623, 0.0183, 0.0392, 
            0.0547, 0.0588, 0.00654, 0.042, 0.0217, 
            0.0885, 0.0445, 0.0141, 0.0234, 0.065, 
            0.0258, 0.0346),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 13.7, 13.2, 
            13.6, 14.2, 14.1, 13.9, 12.9, 
            14.9, 17.7)
    ),
    classbasedtightEleIDCutsV00 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.032, 0.016, 0.0525, 0.09, 0.025, 
            0.035, 0.065, 0.092),
        hOverE = cms.vdouble(0.05, 0.042, 0.045, 0.0, 0.055, 
            0.037, 0.05, 0.0),
        sigmaEtaEta = cms.vdouble(0.0125, 0.011, 0.01, 0.0, 0.0265, 
            0.0252, 0.026, 0.0),
        deltaEtaIn = cms.vdouble(0.0055, 0.003, 0.0065, 0.0, 0.006, 
            0.0055, 0.0075, 0.0),
        eSeedOverPin = cms.vdouble(0.24, 0.94, 0.11, 0.0, 0.32, 
            0.83, 0.0, 0.0)
    ),
    classbasedtightEleIDCutsV01 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.0225, 0.0114, 0.0234, 0.039, 0.0215, 
            0.0095, 0.0148, 0.0167),
        hOverE = cms.vdouble(0.056, 0.0221, 0.037, 0.0, 0.0268, 
            0.0102, 0.0104, 0.0),
        sigmaEtaEta = cms.vdouble(0.0095, 0.0094, 0.0094, 0.0, 0.026, 
            0.0257, 0.0246, 0.0),
        deltaEtaIn = cms.vdouble(0.0043, 0.00282, 0.0036, 0.0, 0.0066, 
            0.0049, 0.0041, 0.0),
        eSeedOverPin = cms.vdouble(0.32, 0.94, 0.221, 0.0, 0.74, 
            0.89, 0.66, 0.0)
    ),
    classbasedtightEleIDCutsV06 = cms.PSet(
        cutdetain = cms.vdouble(0.0116, 0.00449, 0.00938, 0.0184, 0.00678, 
            0.0109, 0.0252, 0.0268, 0.0139),
        cutiso_sum = cms.vdouble(15.5, 12.2, 12.2, 11.7, 7.16, 
            9.71, 8.66, 11.9, 2.98),
        cutip_gsf = cms.vdouble(0.0131, 0.0586, 0.0839, 0.0366, 0.452, 
            0.204, 0.0913, 0.0802, 0.0731),
        cutip_gsfl = cms.vdouble(0.0119, 0.0527, 0.0471, 0.0212, 0.233, 
            0.267, 0.109, 0.122, 0.0479),
        cuthoe = cms.vdouble(0.215, 0.0608, 0.147, 0.369, 0.0349, 
            0.102, 0.52, 0.422, 0.404),
        cutiso_sumoetl = cms.vdouble(6.21, 6.81, 5.3, 5.39, 2.73, 
            4.73, 4.84, 3.46, 3.73),
        cutfmishits = cms.vdouble(1.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 1.5, 2.5, 0.5),
        cuthoel = cms.vdouble(0.228, 0.0836, 0.143, 0.37, 0.0392, 
            0.0979, 0.3, 0.381, 0.339),
        cutdphiin = cms.vdouble(0.0897, 0.0993, 0.295, 0.0979, 0.151, 
            0.252, 0.341, 0.308, 0.328),
        cutseel = cms.vdouble(0.0132, 0.0117, 0.0112, 0.0387, 0.0281, 
            0.0287, 0.00987, 0.0296, 0.0544),
        cutiso_sumoet = cms.vdouble(11.9, 7.81, 6.28, 8.92, 4.65, 
            5.49, 9.36, 8.84, 5.94),
        cutsee = cms.vdouble(0.0145, 0.0116, 0.012, 0.039, 0.0297, 
            0.0311, 0.00987, 0.0347, 0.0917),
        cuteseedopcor = cms.vdouble(0.637, 0.943, 0.742, 0.748, 0.763, 
            0.631, 0.214, 0.873, 0.473),
        cutdphiinl = cms.vdouble(0.061, 0.14, 0.286, 0.0921, 0.197, 
            0.24, 0.333, 0.303, 0.258),
        cutdetainl = cms.vdouble(0.00816, 0.00401, 0.0081, 0.019, 0.00588, 
            0.00893, 0.0171, 0.0434, 0.0143)
    ),
    classbasedtightEleIDCutsV04 = cms.PSet(
        cutdetain = cms.vdouble(0.00811, 0.00341, 0.00633, 0.0103, 0.00667, 
            0.01, 0.0106, 0.0145, 0.0163, 0.0076, 
            0.00259, 0.00511, 0.00941, 0.0043, 0.00857, 
            0.012, 0.0169, 0.00172, 0.00861, 0.00362, 
            0.00601, 0.00925, 0.00489, 0.00832, 0.0119, 
            0.0169, 0.000996),
        cutiso_sum = cms.vdouble(11.8, 8.31, 6.26, 6.18, 3.28, 
            4.38, 4.17, 5.4, 1.57, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0213, 0.0422, 0.0632, 0.0361, 0.073, 
            0.126, 0.171, 0.119, 0.0372, 0.0131, 
            0.0146, 0.0564, 0.0152, 0.0222, 0.0268, 
            0.0314, 0.0884, 0.00374, 0.00852, 0.00761, 
            0.0143, 0.0106, 0.0127, 0.0119, 0.0123, 
            0.0235, 0.00363),
        cuthoe = cms.vdouble(0.0783, 0.0387, 0.105, 0.118, 0.0227, 
            0.062, 0.13, 2.47, 0.38, 0.0888, 
            0.0503, 0.0955, 0.0741, 0.015, 0.03, 
            0.589, 1.13, 0.612, 0.0494, 0.0461, 
            0.0292, 0.0369, 0.0113, 0.0145, 0.124, 
            2.05, 0.61),
        cutfmishits = cms.vdouble(2.5, 1.5, 1.5, 1.5, 1.5, 
            0.5, 2.5, 0.5, 0.5, 2.5, 
            1.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5, -0.5, 2.5, 1.5, 
            0.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(13.7, 11.6, 7.14, 9.98, 3.52, 
            4.87, 6.24, 7.96, 2.53, 11.2, 
            11.9, 7.88, 8.16, 5.58, 5.03, 
            11.4, 8.15, 5.79, 10.4, 11.1, 
            10.4, 7.47, 5.08, 5.9, 11.8, 
            14.1, 11.7),
        cutdcotdist = cms.vdouble(0.0393, 0.0256, 0.00691, 0.0394, 0.0386, 
            0.039, 0.0325, 0.0384, 0.0382, 0.0245, 
            0.000281, 5.46e-05, 0.0342, 0.0232, 0.00107, 
            0.0178, 0.0193, 0.000758, 0.000108, 0.0248, 
            0.000458, 0.0129, 0.00119, 0.0182, 4.53e-05, 
            0.0189, 0.000928),
        cutsee = cms.vdouble(0.0143, 0.0105, 0.0123, 0.0324, 0.0307, 
            0.0301, 0.0109, 0.027, 0.0292, 0.0133, 
            0.0104, 0.0116, 0.0332, 0.0296, 0.031, 
            0.00981, 0.0307, 0.072, 0.0149, 0.0105, 
            0.011, 0.0342, 0.0307, 0.0303, 0.00954, 
            0.0265, 0.0101),
        cuteseedopcor = cms.vdouble(0.784, 0.366, 0.57, 0.911, 0.298, 
            0.645, 0.51, 0.497, 0.932, 0.835, 
            0.968, 0.969, 0.923, 0.898, 0.98, 
            0.63, 0.971, 1.0, 0.515, 0.963, 
            0.986, 0.823, 0.879, 1.01, 0.931, 
            0.937, 1.05),
        cutdphiin = cms.vdouble(0.0404, 0.0499, 0.263, 0.042, 0.0484, 
            0.241, 0.242, 0.231, 0.286, 0.0552, 
            0.0338, 0.154, 0.0623, 0.0183, 0.0392, 
            0.0547, 0.0588, 0.00654, 0.042, 0.0217, 
            0.0885, 0.0445, 0.0141, 0.0234, 0.065, 
            0.0258, 0.0346),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 13.7, 13.2, 
            13.6, 14.2, 14.1, 13.9, 12.9, 
            14.9, 17.7)
    ),
    electronIDType = cms.string('robust'),
    robusttightEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    electronVersion = cms.string(''),
    robusttightEleIDCutsV00 = cms.PSet(
        barrel = cms.vdouble(0.015, 0.0092, 0.02, 0.0025, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.018, 0.025, 0.02, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusttightEleIDCutsV01 = cms.PSet(
        barrel = cms.vdouble(0.01, 0.0099, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.01, 0.028, 0.02, 0.0066, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusttightEleIDCutsV02 = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusttightEleIDCutsV03 = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    reducedEndcapRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    verticesCollection = cms.InputTag("offlinePrimaryVerticesWithBS"),
    classbasedlooseEleIDCuts = cms.PSet(
        cutdetain = cms.vdouble(0.0137, 0.00678, 0.0241, 0.0187, 0.0161, 
            0.0224, 0.0252, 0.0308, 0.0273),
        cutiso_sum = cms.vdouble(33.0, 17.0, 17.9, 18.8, 8.55, 
            12.5, 17.6, 18.5, 2.98),
        cutip_gsf = cms.vdouble(0.0551, 0.0765, 0.143, 0.0874, 0.594, 
            0.37, 0.0913, 1.15, 0.231),
        cutip_gsfl = cms.vdouble(0.0186, 0.0759, 0.138, 0.0473, 0.62, 
            0.304, 0.109, 0.775, 0.0479),
        cuthoe = cms.vdouble(0.247, 0.137, 0.147, 0.371, 0.0588, 
            0.147, 0.52, 0.452, 0.404),
        cutiso_sumoetl = cms.vdouble(11.3, 9.05, 9.07, 9.94, 5.25, 
            6.15, 10.7, 10.8, 4.4),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 4.5, 3.5, 3.5),
        cuthoel = cms.vdouble(0.236, 0.126, 0.147, 0.375, 0.0392, 
            0.145, 0.365, 0.383, 0.384),
        cutdphiin = cms.vdouble(0.0897, 0.262, 0.353, 0.116, 0.357, 
            0.319, 0.342, 0.404, 0.336),
        cutseel = cms.vdouble(0.0164, 0.0118, 0.015, 0.0523, 0.0326, 
            0.0456, 0.0185, 0.0589, 0.0544),
        cutiso_sumoet = cms.vdouble(34.5, 12.7, 12.1, 19.9, 6.35, 
            8.85, 14.0, 10.5, 9.74),
        cutsee = cms.vdouble(0.0176, 0.0125, 0.0181, 0.0415, 0.0364, 
            0.0418, 0.0146, 0.0678, 0.133),
        cuteseedopcor = cms.vdouble(0.63, 0.82, 0.401, 0.718, 0.4, 
            0.458, 0.15, 0.664, 0.373),
        cutdphiinl = cms.vdouble(0.0747, 0.25, 0.356, 0.0956, 0.347, 
            0.326, 0.333, 0.647, 0.289),
        cutdetainl = cms.vdouble(0.0124, 0.00503, 0.0257, 0.0228, 0.0118, 
            0.0178, 0.0188, 0.14, 0.024)
    ),
    reducedBarrelRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    robusthighenergyEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV01 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.0275, 0.09, 0.007, -1, 
            -1, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV00 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.011, 0.09, 0.005, -1, 
            -1, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.1, 0.0275, 0.09, 0.007, -1, 
            -1, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV03 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV02 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    classbasedlooseEleIDCutsV00 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.05, 0.025, 0.053, 0.09, 0.07, 
            0.03, 0.092, 0.092),
        hOverE = cms.vdouble(0.115, 0.1, 0.055, 0.0, 0.145, 
            0.12, 0.15, 0.0),
        sigmaEtaEta = cms.vdouble(0.014, 0.012, 0.0115, 0.0, 0.0275, 
            0.0265, 0.0265, 0.0),
        deltaEtaIn = cms.vdouble(0.009, 0.0045, 0.0085, 0.0, 0.0105, 
            0.0068, 0.01, 0.0),
        eSeedOverPin = cms.vdouble(0.11, 0.91, 0.11, 0.0, 0.0, 
            0.85, 0.0, 0.0)
    ),
    classbasedlooseEleIDCutsV01 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.053, 0.0189, 0.059, 0.099, 0.0278, 
            0.0157, 0.042, 0.08),
        hOverE = cms.vdouble(0.076, 0.033, 0.07, 0.0, 0.083, 
            0.0148, 0.033, 0.0),
        sigmaEtaEta = cms.vdouble(0.0101, 0.0095, 0.0097, 0.0, 0.0271, 
            0.0267, 0.0259, 0.0),
        deltaEtaIn = cms.vdouble(0.0078, 0.00259, 0.0062, 0.0, 0.0078, 
            0.0061, 0.0061, 0.0),
        eSeedOverPin = cms.vdouble(0.3, 0.92, 0.211, 0.0, 0.42, 
            0.88, 0.68, 0.0)
    ),
    classbasedlooseEleIDCutsV02 = cms.PSet(
        cutisohcal = cms.vdouble(13.5, 9.93, 7.56, 14.8, 8.1, 
            10.8, 42.7, 20.1, 9.11, 10.4, 
            6.89, 5.59, 8.53, 9.59, 24.2, 
            2.78, 8.67, 0.288),
        cutmishits = cms.vdouble(5.5, 1.5, 5.5, 2.5, 2.5, 
            2.5, 3.5, 5.5, 0.5, 1.5, 
            2.5, 0.5, 1.5, 1.5, 0.5, 
            0.5, 0.5, 0.5),
        cuthoe = cms.vdouble(0.0887, 0.0934, 0.0949, 0.0986, 0.0431, 
            0.0878, 0.097, 0.0509, 0.098, 0.0991, 
            0.0321, 0.0928, 0.0663, 0.0717, 0.0966, 
            0.0758, 0.0149, 0.0131),
        cutdeta = cms.vdouble(0.00958, 0.00406, 0.0122, 0.0137, 0.00837, 
            0.0127, 0.011, 0.00336, 0.00977, 0.015, 
            0.00675, 0.0109, 0.014, 0.00508, 0.0109, 
            0.0146, 0.00506, 0.0127),
        cuteopin = cms.vdouble(0.878, 0.802, 0.814, 0.942, 0.735, 
            0.774, 0.829, 0.909, 0.829, 0.813, 
            0.86, 0.897, 0.817, 0.831, 0.818, 
            0.861, 0.787, 0.789),
        cutip = cms.vdouble(0.0246, 0.076, 0.0966, 0.0885, 0.441, 
            0.205, 0.0292, 0.0293, 0.0619, 0.0251, 
            0.159, 0.0815, 7.29, 0.0106, 5.76, 
            6.89, 1.27, 5.89),
        cutisotk = cms.vdouble(24.3, 8.45, 14.4, 27.8, 6.02, 
            10.5, 14.1, 10.2, 14.5, 19.1, 
            6.1, 14.1, 8.59, 8.33, 8.3, 
            8.93, 8.6, 16.0),
        cutsee = cms.vdouble(0.0172, 0.0115, 0.0143, 0.0344, 0.0295, 
            0.0304, 0.0145, 0.0108, 0.0128, 0.0347, 
            0.0307, 0.0316, 0.018, 0.011, 0.0132, 
            0.0349, 0.031, 0.0327),
        cutdphi = cms.vdouble(0.0372, 0.114, 0.118, 0.0488, 0.117, 
            0.119, 0.0606, 0.0548, 0.117, 0.07, 
            0.0355, 0.117, 0.088, 0.045, 0.118, 
            0.0919, 0.0236, 0.0515),
        cutisoecal = cms.vdouble(33.4, 28.1, 7.32, 27.4, 7.33, 
            21.7, 93.8, 102.0, 12.1, 26.0, 
            8.91, 10.0, 16.1, 31.3, 16.9, 
            15.4, 13.3, 37.7)
    ),
    classbasedlooseEleIDCutsV03 = cms.PSet(
        cutdetain = cms.vdouble(0.00989, 0.00484, 0.0146, 0.0146, 0.00902, 
            0.0172, 0.0137, 0.0477, 0.0275, 0.00967, 
            0.00377, 0.00924, 0.013, 0.00666, 0.0123, 
            0.0125, 0.0228, 0.0112, 0.0106, 0.0038, 
            0.00897, 0.0139, 0.00667, 0.0122, 0.0122, 
            0.0193, 0.00239),
        cutiso_sum = cms.vdouble(31.5, 10.3, 8.8, 11.0, 6.13, 
            6.94, 7.52, 9.0, 3.5, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0431, 0.0767, 0.139, 0.101, 0.149, 
            0.154, 0.932, 0.15, 0.124, 0.0238, 
            0.0467, 0.0759, 0.0369, 0.147, 0.0986, 
            0.0626, 0.195, 0.116, 0.0122, 0.0125, 
            0.0693, 0.0162, 0.089, 0.0673, 0.0467, 
            0.0651, 0.0221),
        cuthoe = cms.vdouble(0.166, 0.0771, 0.144, 0.37, 0.0497, 
            0.139, 0.401, 2.68, 0.516, 0.234, 
            0.0556, 0.144, 0.368, 0.031, 0.12, 
            0.602, 2.01, 1.05, 0.104, 0.063, 
            0.0565, 0.38, 0.0192, 0.0294, 0.537, 
            4.65, 1.87),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 2.5, 2.5, 1.5, 2.5, 
            1.5, 1.5, 1.5, 1.5, 0.5, 
            2.5, 2.5, 0.5, 2.5, 1.5, 
            0.5, 1.5, 1.5, 0.5, 2.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(28.9, 15.3, 12.0, 18.3, 7.17, 
            9.42, 11.0, 9.81, 3.94, 22.7, 
            15.9, 12.3, 17.0, 7.58, 8.89, 
            15.2, 12.7, 6.17, 20.8, 21.2, 
            17.2, 15.5, 9.37, 10.6, 19.8, 
            22.1, 15.6),
        cutdcotdist = cms.vdouble(0.0393, 0.0392, 0.0397, 0.0394, 0.0393, 
            0.039, 0.0378, 0.0388, 0.0382, 0.0385, 
            0.0167, 0.00325, 0.0394, 0.0387, 0.0388, 
            0.0227, 0.0258, 0.0127, 0.0298, 0.03, 
            0.00946, 0.039, 0.0231, 0.0278, 0.00162, 
            0.0367, 0.0199),
        cutsee = cms.vdouble(0.0175, 0.0127, 0.0177, 0.0373, 0.0314, 
            0.0329, 0.0157, 0.0409, 0.14, 0.0169, 
            0.0106, 0.0142, 0.0363, 0.0322, 0.0354, 
            0.0117, 0.0372, 28.2, 0.0171, 0.0113, 
            0.014, 0.0403, 0.0323, 0.0411, 0.0104, 
            0.0436, 0.0114),
        cuteseedopcor = cms.vdouble(0.78, 0.302, 0.483, 0.904, 0.168, 
            0.645, 0.108, 0.284, 0.324, 0.591, 
            0.286, 0.488, 0.813, 0.791, 0.672, 
            0.398, 0.834, 0.878, 0.515, 0.937, 
            0.806, 0.816, 0.85, 0.507, 0.367, 
            0.83, 0.648),
        cutdphiin = cms.vdouble(0.041, 0.275, 0.365, 0.047, 0.273, 
            0.296, 0.329, 0.465, 0.627, 0.0581, 
            0.0954, 0.327, 0.0702, 0.0582, 0.279, 
            0.117, 0.318, 0.246, 0.0821, 0.052, 
            0.292, 0.116, 0.0435, 0.312, 0.118, 
            0.296, 0.0459),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 12.0, 12.0, 
            12.0, 12.0, 12.0, 12.0, 12.0, 
            12.0, 12.5)
    ),
    classbasedlooseEleIDCutsV04 = cms.PSet(
        cutdetain = cms.vdouble(0.00989, 0.00484, 0.0146, 0.0146, 0.00902, 
            0.0172, 0.0137, 0.0477, 0.0275, 0.00967, 
            0.00377, 0.00924, 0.013, 0.00666, 0.0123, 
            0.0125, 0.0228, 0.0112, 0.0106, 0.0038, 
            0.00897, 0.0139, 0.00667, 0.0122, 0.0122, 
            0.0193, 0.00239),
        cutiso_sum = cms.vdouble(31.5, 10.3, 8.8, 11.0, 6.13, 
            6.94, 7.52, 9.0, 3.5, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0431, 0.0767, 0.139, 0.101, 0.149, 
            0.154, 0.932, 0.15, 0.124, 0.0238, 
            0.0467, 0.0759, 0.0369, 0.147, 0.0986, 
            0.0626, 0.195, 0.116, 0.0122, 0.0125, 
            0.0693, 0.0162, 0.089, 0.0673, 0.0467, 
            0.0651, 0.0221),
        cuthoe = cms.vdouble(0.166, 0.0771, 0.144, 0.37, 0.0497, 
            0.139, 0.401, 2.68, 0.516, 0.234, 
            0.0556, 0.144, 0.368, 0.031, 0.12, 
            0.602, 2.01, 1.05, 0.104, 0.063, 
            0.0565, 0.38, 0.0192, 0.0294, 0.537, 
            4.65, 1.87),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 2.5, 2.5, 1.5, 2.5, 
            1.5, 1.5, 1.5, 1.5, 0.5, 
            2.5, 2.5, 0.5, 2.5, 1.5, 
            0.5, 1.5, 1.5, 0.5, 2.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(28.9, 15.3, 12.0, 18.3, 7.17, 
            9.42, 11.0, 9.81, 3.94, 22.7, 
            15.9, 12.3, 17.0, 7.58, 8.89, 
            15.2, 12.7, 6.17, 20.8, 21.2, 
            17.2, 15.5, 9.37, 10.6, 19.8, 
            22.1, 15.6),
        cutdcotdist = cms.vdouble(0.0393, 0.0392, 0.0397, 0.0394, 0.0393, 
            0.039, 0.0378, 0.0388, 0.0382, 0.0385, 
            0.0167, 0.00325, 0.0394, 0.0387, 0.0388, 
            0.0227, 0.0258, 0.0127, 0.0298, 0.03, 
            0.00946, 0.039, 0.0231, 0.0278, 0.00162, 
            0.0367, 0.0199),
        cutsee = cms.vdouble(0.0175, 0.0127, 0.0177, 0.0373, 0.0314, 
            0.0329, 0.0157, 0.0409, 0.14, 0.0169, 
            0.0106, 0.0142, 0.0363, 0.0322, 0.0354, 
            0.0117, 0.0372, 28.2, 0.0171, 0.0113, 
            0.014, 0.0403, 0.0323, 0.0411, 0.0104, 
            0.0436, 0.0114),
        cuteseedopcor = cms.vdouble(0.78, 0.302, 0.483, 0.904, 0.168, 
            0.645, 0.108, 0.284, 0.324, 0.591, 
            0.286, 0.488, 0.813, 0.791, 0.672, 
            0.398, 0.834, 0.878, 0.515, 0.937, 
            0.806, 0.816, 0.85, 0.507, 0.367, 
            0.83, 0.648),
        cutdphiin = cms.vdouble(0.041, 0.275, 0.365, 0.047, 0.273, 
            0.296, 0.329, 0.465, 0.627, 0.0581, 
            0.0954, 0.327, 0.0702, 0.0582, 0.279, 
            0.117, 0.318, 0.246, 0.0821, 0.052, 
            0.292, 0.116, 0.0435, 0.312, 0.118, 
            0.296, 0.0459),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 12.0, 12.0, 
            12.0, 12.0, 12.0, 12.0, 12.0, 
            12.0, 12.5)
    ),
    classbasedlooseEleIDCutsV06 = cms.PSet(
        cutdetain = cms.vdouble(0.0137, 0.00678, 0.0241, 0.0187, 0.0161, 
            0.0224, 0.0252, 0.0308, 0.0273),
        cutiso_sum = cms.vdouble(33.0, 17.0, 17.9, 18.8, 8.55, 
            12.5, 17.6, 18.5, 2.98),
        cutip_gsf = cms.vdouble(0.0551, 0.0765, 0.143, 0.0874, 0.594, 
            0.37, 0.0913, 1.15, 0.231),
        cutip_gsfl = cms.vdouble(0.0186, 0.0759, 0.138, 0.0473, 0.62, 
            0.304, 0.109, 0.775, 0.0479),
        cuthoe = cms.vdouble(0.247, 0.137, 0.147, 0.371, 0.0588, 
            0.147, 0.52, 0.452, 0.404),
        cutiso_sumoetl = cms.vdouble(11.3, 9.05, 9.07, 9.94, 5.25, 
            6.15, 10.7, 10.8, 4.4),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 4.5, 3.5, 3.5),
        cuthoel = cms.vdouble(0.236, 0.126, 0.147, 0.375, 0.0392, 
            0.145, 0.365, 0.383, 0.384),
        cutdphiin = cms.vdouble(0.0897, 0.262, 0.353, 0.116, 0.357, 
            0.319, 0.342, 0.404, 0.336),
        cutseel = cms.vdouble(0.0164, 0.0118, 0.015, 0.0523, 0.0326, 
            0.0456, 0.0185, 0.0589, 0.0544),
        cutiso_sumoet = cms.vdouble(34.5, 12.7, 12.1, 19.9, 6.35, 
            8.85, 14.0, 10.5, 9.74),
        cutsee = cms.vdouble(0.0176, 0.0125, 0.0181, 0.0415, 0.0364, 
            0.0418, 0.0146, 0.0678, 0.133),
        cuteseedopcor = cms.vdouble(0.63, 0.82, 0.401, 0.718, 0.4, 
            0.458, 0.15, 0.664, 0.373),
        cutdphiinl = cms.vdouble(0.0747, 0.25, 0.356, 0.0956, 0.347, 
            0.326, 0.333, 0.647, 0.289),
        cutdetainl = cms.vdouble(0.0124, 0.00503, 0.0257, 0.0228, 0.0118, 
            0.0178, 0.0188, 0.14, 0.024)
    ),
    src = cms.InputTag("gsfElectrons"),
    robusttightEleIDCuts = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    classbasedtightEleIDCuts = cms.PSet(
        cutdetain = cms.vdouble(0.0116, 0.00449, 0.00938, 0.0184, 0.00678, 
            0.0109, 0.0252, 0.0268, 0.0139),
        cutiso_sum = cms.vdouble(15.5, 12.2, 12.2, 11.7, 7.16, 
            9.71, 8.66, 11.9, 2.98),
        cutip_gsf = cms.vdouble(0.0131, 0.0586, 0.0839, 0.0366, 0.452, 
            0.204, 0.0913, 0.0802, 0.0731),
        cutip_gsfl = cms.vdouble(0.0119, 0.0527, 0.0471, 0.0212, 0.233, 
            0.267, 0.109, 0.122, 0.0479),
        cuthoe = cms.vdouble(0.215, 0.0608, 0.147, 0.369, 0.0349, 
            0.102, 0.52, 0.422, 0.404),
        cutiso_sumoetl = cms.vdouble(6.21, 6.81, 5.3, 5.39, 2.73, 
            4.73, 4.84, 3.46, 3.73),
        cutfmishits = cms.vdouble(1.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 1.5, 2.5, 0.5),
        cuthoel = cms.vdouble(0.228, 0.0836, 0.143, 0.37, 0.0392, 
            0.0979, 0.3, 0.381, 0.339),
        cutdphiin = cms.vdouble(0.0897, 0.0993, 0.295, 0.0979, 0.151, 
            0.252, 0.341, 0.308, 0.328),
        cutseel = cms.vdouble(0.0132, 0.0117, 0.0112, 0.0387, 0.0281, 
            0.0287, 0.00987, 0.0296, 0.0544),
        cutiso_sumoet = cms.vdouble(11.9, 7.81, 6.28, 8.92, 4.65, 
            5.49, 9.36, 8.84, 5.94),
        cutsee = cms.vdouble(0.0145, 0.0116, 0.012, 0.039, 0.0297, 
            0.0311, 0.00987, 0.0347, 0.0917),
        cuteseedopcor = cms.vdouble(0.637, 0.943, 0.742, 0.748, 0.763, 
            0.631, 0.214, 0.873, 0.473),
        cutdphiinl = cms.vdouble(0.061, 0.14, 0.286, 0.0921, 0.197, 
            0.24, 0.333, 0.303, 0.258),
        cutdetainl = cms.vdouble(0.00816, 0.00401, 0.0081, 0.019, 0.00588, 
            0.00893, 0.0171, 0.0434, 0.0143)
    ),
    algorithm = cms.string('eIDCB'),
    robusthighenergyEleIDCuts = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCuts = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV02 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV03 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV00 = cms.PSet(
        barrel = cms.vdouble(0.115, 0.014, 0.09, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.15, 0.0275, 0.092, 0.0105, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV01 = cms.PSet(
        barrel = cms.vdouble(0.075, 0.0132, 0.058, 0.0077, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.083, 0.027, 0.042, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    additionalCategories = cms.bool(True),
    etBinning = cms.bool(True)
)


process.eidLoose = cms.EDProducer("EleIdCutBasedExtProducer",
    electronQuality = cms.string('loose'),
    classbasedtightEleIDCutsV02 = cms.PSet(
        cutisohcal = cms.vdouble(10.9, 7.01, 8.75, 3.51, 7.75, 
            1.62, 11.6, 9.9, 4.97, 5.33, 
            3.18, 2.32, 0.164, 5.46, 12.0, 
            0.00604, 4.1, 0.000628),
        cutmishits = cms.vdouble(5.5, 1.5, 0.5, 1.5, 2.5, 
            0.5, 3.5, 5.5, 0.5, 0.5, 
            0.5, 0.5, 0.5, 1.5, 0.5, 
            0.5, 0.5, 0.5),
        cuthoe = cms.vdouble(0.0871, 0.0289, 0.0783, 0.0946, 0.0245, 
            0.0363, 0.0671, 0.048, 0.0614, 0.0924, 
            0.0158, 0.049, 0.0382, 0.0915, 0.0451, 
            0.0452, 0.00196, 0.0043),
        cutdeta = cms.vdouble(0.00915, 0.00302, 0.0061, 0.0135, 0.00565, 
            0.00793, 0.0102, 0.00266, 0.0106, 0.00903, 
            0.00766, 0.00723, 0.0116, 0.00203, 0.00659, 
            0.0148, 0.00555, 0.0128),
        cuteopin = cms.vdouble(0.878, 0.859, 0.874, 0.944, 0.737, 
            0.773, 0.86, 0.967, 0.917, 0.812, 
            0.915, 1.01, 0.847, 0.953, 0.979, 
            0.841, 0.771, 1.09),
        cutip = cms.vdouble(0.0239, 0.027, 0.0768, 0.0231, 0.178, 
            0.0957, 0.0102, 0.0168, 0.043, 0.0166, 
            0.0594, 0.0308, 2.1, 0.00527, 3.17, 
            4.91, 0.769, 5.9),
        cutisotk = cms.vdouble(6.53, 4.6, 6.0, 8.63, 3.11, 
            7.77, 5.42, 4.81, 4.06, 6.47, 
            2.8, 3.45, 5.29, 5.18, 15.4, 
            5.38, 4.47, 0.0347),
        cutsee = cms.vdouble(0.0131, 0.0106, 0.0115, 0.0306, 0.028, 
            0.0293, 0.0131, 0.0106, 0.0115, 0.0317, 
            0.029, 0.0289, 0.0142, 0.0106, 0.0103, 
            0.035, 0.0296, 0.0333),
        cutdphi = cms.vdouble(0.0369, 0.0307, 0.117, 0.0475, 0.0216, 
            0.117, 0.0372, 0.0246, 0.0426, 0.0612, 
            0.0142, 0.039, 0.0737, 0.0566, 0.0359, 
            0.0187, 0.012, 0.0358),
        cutisoecal = cms.vdouble(20.0, 27.2, 4.48, 13.5, 4.56, 
            3.19, 12.2, 13.1, 7.42, 7.67, 
            4.12, 4.85, 10.1, 12.4, 11.1, 
            11.0, 10.6, 13.4)
    ),
    classbasedtightEleIDCutsV03 = cms.PSet(
        cutdetain = cms.vdouble(0.00811, 0.00341, 0.00633, 0.0103, 0.00667, 
            0.01, 0.0106, 0.0145, 0.0163, 0.0076, 
            0.00259, 0.00511, 0.00941, 0.0043, 0.00857, 
            0.012, 0.0169, 0.00172, 0.00861, 0.00362, 
            0.00601, 0.00925, 0.00489, 0.00832, 0.0119, 
            0.0169, 0.000996),
        cutiso_sum = cms.vdouble(11.8, 8.31, 6.26, 6.18, 3.28, 
            4.38, 4.17, 5.4, 1.57, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0213, 0.0422, 0.0632, 0.0361, 0.073, 
            0.126, 0.171, 0.119, 0.0372, 0.0131, 
            0.0146, 0.0564, 0.0152, 0.0222, 0.0268, 
            0.0314, 0.0884, 0.00374, 0.00852, 0.00761, 
            0.0143, 0.0106, 0.0127, 0.0119, 0.0123, 
            0.0235, 0.00363),
        cuthoe = cms.vdouble(0.0783, 0.0387, 0.105, 0.118, 0.0227, 
            0.062, 0.13, 2.47, 0.38, 0.0888, 
            0.0503, 0.0955, 0.0741, 0.015, 0.03, 
            0.589, 1.13, 0.612, 0.0494, 0.0461, 
            0.0292, 0.0369, 0.0113, 0.0145, 0.124, 
            2.05, 0.61),
        cutfmishits = cms.vdouble(2.5, 1.5, 1.5, 1.5, 1.5, 
            0.5, 2.5, 0.5, 0.5, 2.5, 
            1.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5, -0.5, 2.5, 1.5, 
            0.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(13.7, 11.6, 7.14, 9.98, 3.52, 
            4.87, 6.24, 7.96, 2.53, 11.2, 
            11.9, 7.88, 8.16, 5.58, 5.03, 
            11.4, 8.15, 5.79, 10.4, 11.1, 
            10.4, 7.47, 5.08, 5.9, 11.8, 
            14.1, 11.7),
        cutdcotdist = cms.vdouble(0.0393, 0.0256, 0.00691, 0.0394, 0.0386, 
            0.039, 0.0325, 0.0384, 0.0382, 0.0245, 
            0.000281, 5.46e-05, 0.0342, 0.0232, 0.00107, 
            0.0178, 0.0193, 0.000758, 0.000108, 0.0248, 
            0.000458, 0.0129, 0.00119, 0.0182, 4.53e-05, 
            0.0189, 0.000928),
        cutsee = cms.vdouble(0.0143, 0.0105, 0.0123, 0.0324, 0.0307, 
            0.0301, 0.0109, 0.027, 0.0292, 0.0133, 
            0.0104, 0.0116, 0.0332, 0.0296, 0.031, 
            0.00981, 0.0307, 0.072, 0.0149, 0.0105, 
            0.011, 0.0342, 0.0307, 0.0303, 0.00954, 
            0.0265, 0.0101),
        cuteseedopcor = cms.vdouble(0.784, 0.366, 0.57, 0.911, 0.298, 
            0.645, 0.51, 0.497, 0.932, 0.835, 
            0.968, 0.969, 0.923, 0.898, 0.98, 
            0.63, 0.971, 1.0, 0.515, 0.963, 
            0.986, 0.823, 0.879, 1.01, 0.931, 
            0.937, 1.05),
        cutdphiin = cms.vdouble(0.0404, 0.0499, 0.263, 0.042, 0.0484, 
            0.241, 0.242, 0.231, 0.286, 0.0552, 
            0.0338, 0.154, 0.0623, 0.0183, 0.0392, 
            0.0547, 0.0588, 0.00654, 0.042, 0.0217, 
            0.0885, 0.0445, 0.0141, 0.0234, 0.065, 
            0.0258, 0.0346),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 13.7, 13.2, 
            13.6, 14.2, 14.1, 13.9, 12.9, 
            14.9, 17.7)
    ),
    classbasedtightEleIDCutsV00 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.032, 0.016, 0.0525, 0.09, 0.025, 
            0.035, 0.065, 0.092),
        hOverE = cms.vdouble(0.05, 0.042, 0.045, 0.0, 0.055, 
            0.037, 0.05, 0.0),
        sigmaEtaEta = cms.vdouble(0.0125, 0.011, 0.01, 0.0, 0.0265, 
            0.0252, 0.026, 0.0),
        deltaEtaIn = cms.vdouble(0.0055, 0.003, 0.0065, 0.0, 0.006, 
            0.0055, 0.0075, 0.0),
        eSeedOverPin = cms.vdouble(0.24, 0.94, 0.11, 0.0, 0.32, 
            0.83, 0.0, 0.0)
    ),
    classbasedtightEleIDCutsV01 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.0225, 0.0114, 0.0234, 0.039, 0.0215, 
            0.0095, 0.0148, 0.0167),
        hOverE = cms.vdouble(0.056, 0.0221, 0.037, 0.0, 0.0268, 
            0.0102, 0.0104, 0.0),
        sigmaEtaEta = cms.vdouble(0.0095, 0.0094, 0.0094, 0.0, 0.026, 
            0.0257, 0.0246, 0.0),
        deltaEtaIn = cms.vdouble(0.0043, 0.00282, 0.0036, 0.0, 0.0066, 
            0.0049, 0.0041, 0.0),
        eSeedOverPin = cms.vdouble(0.32, 0.94, 0.221, 0.0, 0.74, 
            0.89, 0.66, 0.0)
    ),
    classbasedtightEleIDCutsV06 = cms.PSet(
        cutdetain = cms.vdouble(0.0116, 0.00449, 0.00938, 0.0184, 0.00678, 
            0.0109, 0.0252, 0.0268, 0.0139),
        cutiso_sum = cms.vdouble(15.5, 12.2, 12.2, 11.7, 7.16, 
            9.71, 8.66, 11.9, 2.98),
        cutip_gsf = cms.vdouble(0.0131, 0.0586, 0.0839, 0.0366, 0.452, 
            0.204, 0.0913, 0.0802, 0.0731),
        cutip_gsfl = cms.vdouble(0.0119, 0.0527, 0.0471, 0.0212, 0.233, 
            0.267, 0.109, 0.122, 0.0479),
        cuthoe = cms.vdouble(0.215, 0.0608, 0.147, 0.369, 0.0349, 
            0.102, 0.52, 0.422, 0.404),
        cutiso_sumoetl = cms.vdouble(6.21, 6.81, 5.3, 5.39, 2.73, 
            4.73, 4.84, 3.46, 3.73),
        cutfmishits = cms.vdouble(1.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 1.5, 2.5, 0.5),
        cuthoel = cms.vdouble(0.228, 0.0836, 0.143, 0.37, 0.0392, 
            0.0979, 0.3, 0.381, 0.339),
        cutdphiin = cms.vdouble(0.0897, 0.0993, 0.295, 0.0979, 0.151, 
            0.252, 0.341, 0.308, 0.328),
        cutseel = cms.vdouble(0.0132, 0.0117, 0.0112, 0.0387, 0.0281, 
            0.0287, 0.00987, 0.0296, 0.0544),
        cutiso_sumoet = cms.vdouble(11.9, 7.81, 6.28, 8.92, 4.65, 
            5.49, 9.36, 8.84, 5.94),
        cutsee = cms.vdouble(0.0145, 0.0116, 0.012, 0.039, 0.0297, 
            0.0311, 0.00987, 0.0347, 0.0917),
        cuteseedopcor = cms.vdouble(0.637, 0.943, 0.742, 0.748, 0.763, 
            0.631, 0.214, 0.873, 0.473),
        cutdphiinl = cms.vdouble(0.061, 0.14, 0.286, 0.0921, 0.197, 
            0.24, 0.333, 0.303, 0.258),
        cutdetainl = cms.vdouble(0.00816, 0.00401, 0.0081, 0.019, 0.00588, 
            0.00893, 0.0171, 0.0434, 0.0143)
    ),
    classbasedtightEleIDCutsV04 = cms.PSet(
        cutdetain = cms.vdouble(0.00811, 0.00341, 0.00633, 0.0103, 0.00667, 
            0.01, 0.0106, 0.0145, 0.0163, 0.0076, 
            0.00259, 0.00511, 0.00941, 0.0043, 0.00857, 
            0.012, 0.0169, 0.00172, 0.00861, 0.00362, 
            0.00601, 0.00925, 0.00489, 0.00832, 0.0119, 
            0.0169, 0.000996),
        cutiso_sum = cms.vdouble(11.8, 8.31, 6.26, 6.18, 3.28, 
            4.38, 4.17, 5.4, 1.57, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0213, 0.0422, 0.0632, 0.0361, 0.073, 
            0.126, 0.171, 0.119, 0.0372, 0.0131, 
            0.0146, 0.0564, 0.0152, 0.0222, 0.0268, 
            0.0314, 0.0884, 0.00374, 0.00852, 0.00761, 
            0.0143, 0.0106, 0.0127, 0.0119, 0.0123, 
            0.0235, 0.00363),
        cuthoe = cms.vdouble(0.0783, 0.0387, 0.105, 0.118, 0.0227, 
            0.062, 0.13, 2.47, 0.38, 0.0888, 
            0.0503, 0.0955, 0.0741, 0.015, 0.03, 
            0.589, 1.13, 0.612, 0.0494, 0.0461, 
            0.0292, 0.0369, 0.0113, 0.0145, 0.124, 
            2.05, 0.61),
        cutfmishits = cms.vdouble(2.5, 1.5, 1.5, 1.5, 1.5, 
            0.5, 2.5, 0.5, 0.5, 2.5, 
            1.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5, -0.5, 2.5, 1.5, 
            0.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(13.7, 11.6, 7.14, 9.98, 3.52, 
            4.87, 6.24, 7.96, 2.53, 11.2, 
            11.9, 7.88, 8.16, 5.58, 5.03, 
            11.4, 8.15, 5.79, 10.4, 11.1, 
            10.4, 7.47, 5.08, 5.9, 11.8, 
            14.1, 11.7),
        cutdcotdist = cms.vdouble(0.0393, 0.0256, 0.00691, 0.0394, 0.0386, 
            0.039, 0.0325, 0.0384, 0.0382, 0.0245, 
            0.000281, 5.46e-05, 0.0342, 0.0232, 0.00107, 
            0.0178, 0.0193, 0.000758, 0.000108, 0.0248, 
            0.000458, 0.0129, 0.00119, 0.0182, 4.53e-05, 
            0.0189, 0.000928),
        cutsee = cms.vdouble(0.0143, 0.0105, 0.0123, 0.0324, 0.0307, 
            0.0301, 0.0109, 0.027, 0.0292, 0.0133, 
            0.0104, 0.0116, 0.0332, 0.0296, 0.031, 
            0.00981, 0.0307, 0.072, 0.0149, 0.0105, 
            0.011, 0.0342, 0.0307, 0.0303, 0.00954, 
            0.0265, 0.0101),
        cuteseedopcor = cms.vdouble(0.784, 0.366, 0.57, 0.911, 0.298, 
            0.645, 0.51, 0.497, 0.932, 0.835, 
            0.968, 0.969, 0.923, 0.898, 0.98, 
            0.63, 0.971, 1.0, 0.515, 0.963, 
            0.986, 0.823, 0.879, 1.01, 0.931, 
            0.937, 1.05),
        cutdphiin = cms.vdouble(0.0404, 0.0499, 0.263, 0.042, 0.0484, 
            0.241, 0.242, 0.231, 0.286, 0.0552, 
            0.0338, 0.154, 0.0623, 0.0183, 0.0392, 
            0.0547, 0.0588, 0.00654, 0.042, 0.0217, 
            0.0885, 0.0445, 0.0141, 0.0234, 0.065, 
            0.0258, 0.0346),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 13.7, 13.2, 
            13.6, 14.2, 14.1, 13.9, 12.9, 
            14.9, 17.7)
    ),
    electronIDType = cms.string('classbased'),
    robusttightEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    electronVersion = cms.string(''),
    robusttightEleIDCutsV00 = cms.PSet(
        barrel = cms.vdouble(0.015, 0.0092, 0.02, 0.0025, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.018, 0.025, 0.02, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusttightEleIDCutsV01 = cms.PSet(
        barrel = cms.vdouble(0.01, 0.0099, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.01, 0.028, 0.02, 0.0066, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    reducedBarrelRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    robusttightEleIDCutsV03 = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    reducedEndcapRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    verticesCollection = cms.InputTag("offlinePrimaryVerticesWithBS"),
    classbasedlooseEleIDCuts = cms.PSet(
        cutdetain = cms.vdouble(0.0137, 0.00678, 0.0241, 0.0187, 0.0161, 
            0.0224, 0.0252, 0.0308, 0.0273),
        cutiso_sum = cms.vdouble(33.0, 17.0, 17.9, 18.8, 8.55, 
            12.5, 17.6, 18.5, 2.98),
        cutip_gsf = cms.vdouble(0.0551, 0.0765, 0.143, 0.0874, 0.594, 
            0.37, 0.0913, 1.15, 0.231),
        cutip_gsfl = cms.vdouble(0.0186, 0.0759, 0.138, 0.0473, 0.62, 
            0.304, 0.109, 0.775, 0.0479),
        cuthoe = cms.vdouble(0.247, 0.137, 0.147, 0.371, 0.0588, 
            0.147, 0.52, 0.452, 0.404),
        cutiso_sumoetl = cms.vdouble(11.3, 9.05, 9.07, 9.94, 5.25, 
            6.15, 10.7, 10.8, 4.4),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 4.5, 3.5, 3.5),
        cuthoel = cms.vdouble(0.236, 0.126, 0.147, 0.375, 0.0392, 
            0.145, 0.365, 0.383, 0.384),
        cutdphiin = cms.vdouble(0.0897, 0.262, 0.353, 0.116, 0.357, 
            0.319, 0.342, 0.404, 0.336),
        cutseel = cms.vdouble(0.0164, 0.0118, 0.015, 0.0523, 0.0326, 
            0.0456, 0.0185, 0.0589, 0.0544),
        cutiso_sumoet = cms.vdouble(34.5, 12.7, 12.1, 19.9, 6.35, 
            8.85, 14.0, 10.5, 9.74),
        cutsee = cms.vdouble(0.0176, 0.0125, 0.0181, 0.0415, 0.0364, 
            0.0418, 0.0146, 0.0678, 0.133),
        cuteseedopcor = cms.vdouble(0.63, 0.82, 0.401, 0.718, 0.4, 
            0.458, 0.15, 0.664, 0.373),
        cutdphiinl = cms.vdouble(0.0747, 0.25, 0.356, 0.0956, 0.347, 
            0.326, 0.333, 0.647, 0.289),
        cutdetainl = cms.vdouble(0.0124, 0.00503, 0.0257, 0.0228, 0.0118, 
            0.0178, 0.0188, 0.14, 0.024)
    ),
    robusttightEleIDCutsV02 = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV01 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.0275, 0.09, 0.007, -1, 
            -1, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV00 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.011, 0.09, 0.005, -1, 
            -1, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.1, 0.0275, 0.09, 0.007, -1, 
            -1, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV03 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV02 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    classbasedlooseEleIDCutsV00 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.05, 0.025, 0.053, 0.09, 0.07, 
            0.03, 0.092, 0.092),
        hOverE = cms.vdouble(0.115, 0.1, 0.055, 0.0, 0.145, 
            0.12, 0.15, 0.0),
        sigmaEtaEta = cms.vdouble(0.014, 0.012, 0.0115, 0.0, 0.0275, 
            0.0265, 0.0265, 0.0),
        deltaEtaIn = cms.vdouble(0.009, 0.0045, 0.0085, 0.0, 0.0105, 
            0.0068, 0.01, 0.0),
        eSeedOverPin = cms.vdouble(0.11, 0.91, 0.11, 0.0, 0.0, 
            0.85, 0.0, 0.0)
    ),
    classbasedlooseEleIDCutsV01 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.053, 0.0189, 0.059, 0.099, 0.0278, 
            0.0157, 0.042, 0.08),
        hOverE = cms.vdouble(0.076, 0.033, 0.07, 0.0, 0.083, 
            0.0148, 0.033, 0.0),
        sigmaEtaEta = cms.vdouble(0.0101, 0.0095, 0.0097, 0.0, 0.0271, 
            0.0267, 0.0259, 0.0),
        deltaEtaIn = cms.vdouble(0.0078, 0.00259, 0.0062, 0.0, 0.0078, 
            0.0061, 0.0061, 0.0),
        eSeedOverPin = cms.vdouble(0.3, 0.92, 0.211, 0.0, 0.42, 
            0.88, 0.68, 0.0)
    ),
    classbasedlooseEleIDCutsV02 = cms.PSet(
        cutisohcal = cms.vdouble(13.5, 9.93, 7.56, 14.8, 8.1, 
            10.8, 42.7, 20.1, 9.11, 10.4, 
            6.89, 5.59, 8.53, 9.59, 24.2, 
            2.78, 8.67, 0.288),
        cutmishits = cms.vdouble(5.5, 1.5, 5.5, 2.5, 2.5, 
            2.5, 3.5, 5.5, 0.5, 1.5, 
            2.5, 0.5, 1.5, 1.5, 0.5, 
            0.5, 0.5, 0.5),
        cuthoe = cms.vdouble(0.0887, 0.0934, 0.0949, 0.0986, 0.0431, 
            0.0878, 0.097, 0.0509, 0.098, 0.0991, 
            0.0321, 0.0928, 0.0663, 0.0717, 0.0966, 
            0.0758, 0.0149, 0.0131),
        cutdeta = cms.vdouble(0.00958, 0.00406, 0.0122, 0.0137, 0.00837, 
            0.0127, 0.011, 0.00336, 0.00977, 0.015, 
            0.00675, 0.0109, 0.014, 0.00508, 0.0109, 
            0.0146, 0.00506, 0.0127),
        cuteopin = cms.vdouble(0.878, 0.802, 0.814, 0.942, 0.735, 
            0.774, 0.829, 0.909, 0.829, 0.813, 
            0.86, 0.897, 0.817, 0.831, 0.818, 
            0.861, 0.787, 0.789),
        cutip = cms.vdouble(0.0246, 0.076, 0.0966, 0.0885, 0.441, 
            0.205, 0.0292, 0.0293, 0.0619, 0.0251, 
            0.159, 0.0815, 7.29, 0.0106, 5.76, 
            6.89, 1.27, 5.89),
        cutisotk = cms.vdouble(24.3, 8.45, 14.4, 27.8, 6.02, 
            10.5, 14.1, 10.2, 14.5, 19.1, 
            6.1, 14.1, 8.59, 8.33, 8.3, 
            8.93, 8.6, 16.0),
        cutsee = cms.vdouble(0.0172, 0.0115, 0.0143, 0.0344, 0.0295, 
            0.0304, 0.0145, 0.0108, 0.0128, 0.0347, 
            0.0307, 0.0316, 0.018, 0.011, 0.0132, 
            0.0349, 0.031, 0.0327),
        cutdphi = cms.vdouble(0.0372, 0.114, 0.118, 0.0488, 0.117, 
            0.119, 0.0606, 0.0548, 0.117, 0.07, 
            0.0355, 0.117, 0.088, 0.045, 0.118, 
            0.0919, 0.0236, 0.0515),
        cutisoecal = cms.vdouble(33.4, 28.1, 7.32, 27.4, 7.33, 
            21.7, 93.8, 102.0, 12.1, 26.0, 
            8.91, 10.0, 16.1, 31.3, 16.9, 
            15.4, 13.3, 37.7)
    ),
    classbasedlooseEleIDCutsV03 = cms.PSet(
        cutdetain = cms.vdouble(0.00989, 0.00484, 0.0146, 0.0146, 0.00902, 
            0.0172, 0.0137, 0.0477, 0.0275, 0.00967, 
            0.00377, 0.00924, 0.013, 0.00666, 0.0123, 
            0.0125, 0.0228, 0.0112, 0.0106, 0.0038, 
            0.00897, 0.0139, 0.00667, 0.0122, 0.0122, 
            0.0193, 0.00239),
        cutiso_sum = cms.vdouble(31.5, 10.3, 8.8, 11.0, 6.13, 
            6.94, 7.52, 9.0, 3.5, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0431, 0.0767, 0.139, 0.101, 0.149, 
            0.154, 0.932, 0.15, 0.124, 0.0238, 
            0.0467, 0.0759, 0.0369, 0.147, 0.0986, 
            0.0626, 0.195, 0.116, 0.0122, 0.0125, 
            0.0693, 0.0162, 0.089, 0.0673, 0.0467, 
            0.0651, 0.0221),
        cuthoe = cms.vdouble(0.166, 0.0771, 0.144, 0.37, 0.0497, 
            0.139, 0.401, 2.68, 0.516, 0.234, 
            0.0556, 0.144, 0.368, 0.031, 0.12, 
            0.602, 2.01, 1.05, 0.104, 0.063, 
            0.0565, 0.38, 0.0192, 0.0294, 0.537, 
            4.65, 1.87),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 2.5, 2.5, 1.5, 2.5, 
            1.5, 1.5, 1.5, 1.5, 0.5, 
            2.5, 2.5, 0.5, 2.5, 1.5, 
            0.5, 1.5, 1.5, 0.5, 2.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(28.9, 15.3, 12.0, 18.3, 7.17, 
            9.42, 11.0, 9.81, 3.94, 22.7, 
            15.9, 12.3, 17.0, 7.58, 8.89, 
            15.2, 12.7, 6.17, 20.8, 21.2, 
            17.2, 15.5, 9.37, 10.6, 19.8, 
            22.1, 15.6),
        cutdcotdist = cms.vdouble(0.0393, 0.0392, 0.0397, 0.0394, 0.0393, 
            0.039, 0.0378, 0.0388, 0.0382, 0.0385, 
            0.0167, 0.00325, 0.0394, 0.0387, 0.0388, 
            0.0227, 0.0258, 0.0127, 0.0298, 0.03, 
            0.00946, 0.039, 0.0231, 0.0278, 0.00162, 
            0.0367, 0.0199),
        cutsee = cms.vdouble(0.0175, 0.0127, 0.0177, 0.0373, 0.0314, 
            0.0329, 0.0157, 0.0409, 0.14, 0.0169, 
            0.0106, 0.0142, 0.0363, 0.0322, 0.0354, 
            0.0117, 0.0372, 28.2, 0.0171, 0.0113, 
            0.014, 0.0403, 0.0323, 0.0411, 0.0104, 
            0.0436, 0.0114),
        cuteseedopcor = cms.vdouble(0.78, 0.302, 0.483, 0.904, 0.168, 
            0.645, 0.108, 0.284, 0.324, 0.591, 
            0.286, 0.488, 0.813, 0.791, 0.672, 
            0.398, 0.834, 0.878, 0.515, 0.937, 
            0.806, 0.816, 0.85, 0.507, 0.367, 
            0.83, 0.648),
        cutdphiin = cms.vdouble(0.041, 0.275, 0.365, 0.047, 0.273, 
            0.296, 0.329, 0.465, 0.627, 0.0581, 
            0.0954, 0.327, 0.0702, 0.0582, 0.279, 
            0.117, 0.318, 0.246, 0.0821, 0.052, 
            0.292, 0.116, 0.0435, 0.312, 0.118, 
            0.296, 0.0459),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 12.0, 12.0, 
            12.0, 12.0, 12.0, 12.0, 12.0, 
            12.0, 12.5)
    ),
    classbasedlooseEleIDCutsV04 = cms.PSet(
        cutdetain = cms.vdouble(0.00989, 0.00484, 0.0146, 0.0146, 0.00902, 
            0.0172, 0.0137, 0.0477, 0.0275, 0.00967, 
            0.00377, 0.00924, 0.013, 0.00666, 0.0123, 
            0.0125, 0.0228, 0.0112, 0.0106, 0.0038, 
            0.00897, 0.0139, 0.00667, 0.0122, 0.0122, 
            0.0193, 0.00239),
        cutiso_sum = cms.vdouble(31.5, 10.3, 8.8, 11.0, 6.13, 
            6.94, 7.52, 9.0, 3.5, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0431, 0.0767, 0.139, 0.101, 0.149, 
            0.154, 0.932, 0.15, 0.124, 0.0238, 
            0.0467, 0.0759, 0.0369, 0.147, 0.0986, 
            0.0626, 0.195, 0.116, 0.0122, 0.0125, 
            0.0693, 0.0162, 0.089, 0.0673, 0.0467, 
            0.0651, 0.0221),
        cuthoe = cms.vdouble(0.166, 0.0771, 0.144, 0.37, 0.0497, 
            0.139, 0.401, 2.68, 0.516, 0.234, 
            0.0556, 0.144, 0.368, 0.031, 0.12, 
            0.602, 2.01, 1.05, 0.104, 0.063, 
            0.0565, 0.38, 0.0192, 0.0294, 0.537, 
            4.65, 1.87),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 2.5, 2.5, 1.5, 2.5, 
            1.5, 1.5, 1.5, 1.5, 0.5, 
            2.5, 2.5, 0.5, 2.5, 1.5, 
            0.5, 1.5, 1.5, 0.5, 2.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(28.9, 15.3, 12.0, 18.3, 7.17, 
            9.42, 11.0, 9.81, 3.94, 22.7, 
            15.9, 12.3, 17.0, 7.58, 8.89, 
            15.2, 12.7, 6.17, 20.8, 21.2, 
            17.2, 15.5, 9.37, 10.6, 19.8, 
            22.1, 15.6),
        cutdcotdist = cms.vdouble(0.0393, 0.0392, 0.0397, 0.0394, 0.0393, 
            0.039, 0.0378, 0.0388, 0.0382, 0.0385, 
            0.0167, 0.00325, 0.0394, 0.0387, 0.0388, 
            0.0227, 0.0258, 0.0127, 0.0298, 0.03, 
            0.00946, 0.039, 0.0231, 0.0278, 0.00162, 
            0.0367, 0.0199),
        cutsee = cms.vdouble(0.0175, 0.0127, 0.0177, 0.0373, 0.0314, 
            0.0329, 0.0157, 0.0409, 0.14, 0.0169, 
            0.0106, 0.0142, 0.0363, 0.0322, 0.0354, 
            0.0117, 0.0372, 28.2, 0.0171, 0.0113, 
            0.014, 0.0403, 0.0323, 0.0411, 0.0104, 
            0.0436, 0.0114),
        cuteseedopcor = cms.vdouble(0.78, 0.302, 0.483, 0.904, 0.168, 
            0.645, 0.108, 0.284, 0.324, 0.591, 
            0.286, 0.488, 0.813, 0.791, 0.672, 
            0.398, 0.834, 0.878, 0.515, 0.937, 
            0.806, 0.816, 0.85, 0.507, 0.367, 
            0.83, 0.648),
        cutdphiin = cms.vdouble(0.041, 0.275, 0.365, 0.047, 0.273, 
            0.296, 0.329, 0.465, 0.627, 0.0581, 
            0.0954, 0.327, 0.0702, 0.0582, 0.279, 
            0.117, 0.318, 0.246, 0.0821, 0.052, 
            0.292, 0.116, 0.0435, 0.312, 0.118, 
            0.296, 0.0459),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 12.0, 12.0, 
            12.0, 12.0, 12.0, 12.0, 12.0, 
            12.0, 12.5)
    ),
    classbasedlooseEleIDCutsV06 = cms.PSet(
        cutdetain = cms.vdouble(0.0137, 0.00678, 0.0241, 0.0187, 0.0161, 
            0.0224, 0.0252, 0.0308, 0.0273),
        cutiso_sum = cms.vdouble(33.0, 17.0, 17.9, 18.8, 8.55, 
            12.5, 17.6, 18.5, 2.98),
        cutip_gsf = cms.vdouble(0.0551, 0.0765, 0.143, 0.0874, 0.594, 
            0.37, 0.0913, 1.15, 0.231),
        cutip_gsfl = cms.vdouble(0.0186, 0.0759, 0.138, 0.0473, 0.62, 
            0.304, 0.109, 0.775, 0.0479),
        cuthoe = cms.vdouble(0.247, 0.137, 0.147, 0.371, 0.0588, 
            0.147, 0.52, 0.452, 0.404),
        cutiso_sumoetl = cms.vdouble(11.3, 9.05, 9.07, 9.94, 5.25, 
            6.15, 10.7, 10.8, 4.4),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 4.5, 3.5, 3.5),
        cuthoel = cms.vdouble(0.236, 0.126, 0.147, 0.375, 0.0392, 
            0.145, 0.365, 0.383, 0.384),
        cutdphiin = cms.vdouble(0.0897, 0.262, 0.353, 0.116, 0.357, 
            0.319, 0.342, 0.404, 0.336),
        cutseel = cms.vdouble(0.0164, 0.0118, 0.015, 0.0523, 0.0326, 
            0.0456, 0.0185, 0.0589, 0.0544),
        cutiso_sumoet = cms.vdouble(34.5, 12.7, 12.1, 19.9, 6.35, 
            8.85, 14.0, 10.5, 9.74),
        cutsee = cms.vdouble(0.0176, 0.0125, 0.0181, 0.0415, 0.0364, 
            0.0418, 0.0146, 0.0678, 0.133),
        cuteseedopcor = cms.vdouble(0.63, 0.82, 0.401, 0.718, 0.4, 
            0.458, 0.15, 0.664, 0.373),
        cutdphiinl = cms.vdouble(0.0747, 0.25, 0.356, 0.0956, 0.347, 
            0.326, 0.333, 0.647, 0.289),
        cutdetainl = cms.vdouble(0.0124, 0.00503, 0.0257, 0.0228, 0.0118, 
            0.0178, 0.0188, 0.14, 0.024)
    ),
    src = cms.InputTag("gsfElectrons"),
    robusttightEleIDCuts = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    classbasedtightEleIDCuts = cms.PSet(
        cutdetain = cms.vdouble(0.0116, 0.00449, 0.00938, 0.0184, 0.00678, 
            0.0109, 0.0252, 0.0268, 0.0139),
        cutiso_sum = cms.vdouble(15.5, 12.2, 12.2, 11.7, 7.16, 
            9.71, 8.66, 11.9, 2.98),
        cutip_gsf = cms.vdouble(0.0131, 0.0586, 0.0839, 0.0366, 0.452, 
            0.204, 0.0913, 0.0802, 0.0731),
        cutip_gsfl = cms.vdouble(0.0119, 0.0527, 0.0471, 0.0212, 0.233, 
            0.267, 0.109, 0.122, 0.0479),
        cuthoe = cms.vdouble(0.215, 0.0608, 0.147, 0.369, 0.0349, 
            0.102, 0.52, 0.422, 0.404),
        cutiso_sumoetl = cms.vdouble(6.21, 6.81, 5.3, 5.39, 2.73, 
            4.73, 4.84, 3.46, 3.73),
        cutfmishits = cms.vdouble(1.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 1.5, 2.5, 0.5),
        cuthoel = cms.vdouble(0.228, 0.0836, 0.143, 0.37, 0.0392, 
            0.0979, 0.3, 0.381, 0.339),
        cutdphiin = cms.vdouble(0.0897, 0.0993, 0.295, 0.0979, 0.151, 
            0.252, 0.341, 0.308, 0.328),
        cutseel = cms.vdouble(0.0132, 0.0117, 0.0112, 0.0387, 0.0281, 
            0.0287, 0.00987, 0.0296, 0.0544),
        cutiso_sumoet = cms.vdouble(11.9, 7.81, 6.28, 8.92, 4.65, 
            5.49, 9.36, 8.84, 5.94),
        cutsee = cms.vdouble(0.0145, 0.0116, 0.012, 0.039, 0.0297, 
            0.0311, 0.00987, 0.0347, 0.0917),
        cuteseedopcor = cms.vdouble(0.637, 0.943, 0.742, 0.748, 0.763, 
            0.631, 0.214, 0.873, 0.473),
        cutdphiinl = cms.vdouble(0.061, 0.14, 0.286, 0.0921, 0.197, 
            0.24, 0.333, 0.303, 0.258),
        cutdetainl = cms.vdouble(0.00816, 0.00401, 0.0081, 0.019, 0.00588, 
            0.00893, 0.0171, 0.0434, 0.0143)
    ),
    algorithm = cms.string('eIDCB'),
    robusthighenergyEleIDCuts = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCuts = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV02 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV03 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV00 = cms.PSet(
        barrel = cms.vdouble(0.115, 0.014, 0.09, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.15, 0.0275, 0.092, 0.0105, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV01 = cms.PSet(
        barrel = cms.vdouble(0.075, 0.0132, 0.058, 0.0077, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.083, 0.027, 0.042, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    additionalCategories = cms.bool(True),
    etBinning = cms.bool(True)
)


process.eidRobustHighEnergy = cms.EDProducer("EleIdCutBasedExtProducer",
    electronQuality = cms.string('highenergy'),
    classbasedtightEleIDCutsV02 = cms.PSet(
        cutisohcal = cms.vdouble(10.9, 7.01, 8.75, 3.51, 7.75, 
            1.62, 11.6, 9.9, 4.97, 5.33, 
            3.18, 2.32, 0.164, 5.46, 12.0, 
            0.00604, 4.1, 0.000628),
        cutmishits = cms.vdouble(5.5, 1.5, 0.5, 1.5, 2.5, 
            0.5, 3.5, 5.5, 0.5, 0.5, 
            0.5, 0.5, 0.5, 1.5, 0.5, 
            0.5, 0.5, 0.5),
        cuthoe = cms.vdouble(0.0871, 0.0289, 0.0783, 0.0946, 0.0245, 
            0.0363, 0.0671, 0.048, 0.0614, 0.0924, 
            0.0158, 0.049, 0.0382, 0.0915, 0.0451, 
            0.0452, 0.00196, 0.0043),
        cutdeta = cms.vdouble(0.00915, 0.00302, 0.0061, 0.0135, 0.00565, 
            0.00793, 0.0102, 0.00266, 0.0106, 0.00903, 
            0.00766, 0.00723, 0.0116, 0.00203, 0.00659, 
            0.0148, 0.00555, 0.0128),
        cuteopin = cms.vdouble(0.878, 0.859, 0.874, 0.944, 0.737, 
            0.773, 0.86, 0.967, 0.917, 0.812, 
            0.915, 1.01, 0.847, 0.953, 0.979, 
            0.841, 0.771, 1.09),
        cutip = cms.vdouble(0.0239, 0.027, 0.0768, 0.0231, 0.178, 
            0.0957, 0.0102, 0.0168, 0.043, 0.0166, 
            0.0594, 0.0308, 2.1, 0.00527, 3.17, 
            4.91, 0.769, 5.9),
        cutisotk = cms.vdouble(6.53, 4.6, 6.0, 8.63, 3.11, 
            7.77, 5.42, 4.81, 4.06, 6.47, 
            2.8, 3.45, 5.29, 5.18, 15.4, 
            5.38, 4.47, 0.0347),
        cutsee = cms.vdouble(0.0131, 0.0106, 0.0115, 0.0306, 0.028, 
            0.0293, 0.0131, 0.0106, 0.0115, 0.0317, 
            0.029, 0.0289, 0.0142, 0.0106, 0.0103, 
            0.035, 0.0296, 0.0333),
        cutdphi = cms.vdouble(0.0369, 0.0307, 0.117, 0.0475, 0.0216, 
            0.117, 0.0372, 0.0246, 0.0426, 0.0612, 
            0.0142, 0.039, 0.0737, 0.0566, 0.0359, 
            0.0187, 0.012, 0.0358),
        cutisoecal = cms.vdouble(20.0, 27.2, 4.48, 13.5, 4.56, 
            3.19, 12.2, 13.1, 7.42, 7.67, 
            4.12, 4.85, 10.1, 12.4, 11.1, 
            11.0, 10.6, 13.4)
    ),
    classbasedtightEleIDCutsV03 = cms.PSet(
        cutdetain = cms.vdouble(0.00811, 0.00341, 0.00633, 0.0103, 0.00667, 
            0.01, 0.0106, 0.0145, 0.0163, 0.0076, 
            0.00259, 0.00511, 0.00941, 0.0043, 0.00857, 
            0.012, 0.0169, 0.00172, 0.00861, 0.00362, 
            0.00601, 0.00925, 0.00489, 0.00832, 0.0119, 
            0.0169, 0.000996),
        cutiso_sum = cms.vdouble(11.8, 8.31, 6.26, 6.18, 3.28, 
            4.38, 4.17, 5.4, 1.57, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0213, 0.0422, 0.0632, 0.0361, 0.073, 
            0.126, 0.171, 0.119, 0.0372, 0.0131, 
            0.0146, 0.0564, 0.0152, 0.0222, 0.0268, 
            0.0314, 0.0884, 0.00374, 0.00852, 0.00761, 
            0.0143, 0.0106, 0.0127, 0.0119, 0.0123, 
            0.0235, 0.00363),
        cuthoe = cms.vdouble(0.0783, 0.0387, 0.105, 0.118, 0.0227, 
            0.062, 0.13, 2.47, 0.38, 0.0888, 
            0.0503, 0.0955, 0.0741, 0.015, 0.03, 
            0.589, 1.13, 0.612, 0.0494, 0.0461, 
            0.0292, 0.0369, 0.0113, 0.0145, 0.124, 
            2.05, 0.61),
        cutfmishits = cms.vdouble(2.5, 1.5, 1.5, 1.5, 1.5, 
            0.5, 2.5, 0.5, 0.5, 2.5, 
            1.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5, -0.5, 2.5, 1.5, 
            0.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(13.7, 11.6, 7.14, 9.98, 3.52, 
            4.87, 6.24, 7.96, 2.53, 11.2, 
            11.9, 7.88, 8.16, 5.58, 5.03, 
            11.4, 8.15, 5.79, 10.4, 11.1, 
            10.4, 7.47, 5.08, 5.9, 11.8, 
            14.1, 11.7),
        cutdcotdist = cms.vdouble(0.0393, 0.0256, 0.00691, 0.0394, 0.0386, 
            0.039, 0.0325, 0.0384, 0.0382, 0.0245, 
            0.000281, 5.46e-05, 0.0342, 0.0232, 0.00107, 
            0.0178, 0.0193, 0.000758, 0.000108, 0.0248, 
            0.000458, 0.0129, 0.00119, 0.0182, 4.53e-05, 
            0.0189, 0.000928),
        cutsee = cms.vdouble(0.0143, 0.0105, 0.0123, 0.0324, 0.0307, 
            0.0301, 0.0109, 0.027, 0.0292, 0.0133, 
            0.0104, 0.0116, 0.0332, 0.0296, 0.031, 
            0.00981, 0.0307, 0.072, 0.0149, 0.0105, 
            0.011, 0.0342, 0.0307, 0.0303, 0.00954, 
            0.0265, 0.0101),
        cuteseedopcor = cms.vdouble(0.784, 0.366, 0.57, 0.911, 0.298, 
            0.645, 0.51, 0.497, 0.932, 0.835, 
            0.968, 0.969, 0.923, 0.898, 0.98, 
            0.63, 0.971, 1.0, 0.515, 0.963, 
            0.986, 0.823, 0.879, 1.01, 0.931, 
            0.937, 1.05),
        cutdphiin = cms.vdouble(0.0404, 0.0499, 0.263, 0.042, 0.0484, 
            0.241, 0.242, 0.231, 0.286, 0.0552, 
            0.0338, 0.154, 0.0623, 0.0183, 0.0392, 
            0.0547, 0.0588, 0.00654, 0.042, 0.0217, 
            0.0885, 0.0445, 0.0141, 0.0234, 0.065, 
            0.0258, 0.0346),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 13.7, 13.2, 
            13.6, 14.2, 14.1, 13.9, 12.9, 
            14.9, 17.7)
    ),
    classbasedtightEleIDCutsV00 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.032, 0.016, 0.0525, 0.09, 0.025, 
            0.035, 0.065, 0.092),
        hOverE = cms.vdouble(0.05, 0.042, 0.045, 0.0, 0.055, 
            0.037, 0.05, 0.0),
        sigmaEtaEta = cms.vdouble(0.0125, 0.011, 0.01, 0.0, 0.0265, 
            0.0252, 0.026, 0.0),
        deltaEtaIn = cms.vdouble(0.0055, 0.003, 0.0065, 0.0, 0.006, 
            0.0055, 0.0075, 0.0),
        eSeedOverPin = cms.vdouble(0.24, 0.94, 0.11, 0.0, 0.32, 
            0.83, 0.0, 0.0)
    ),
    classbasedtightEleIDCutsV01 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.0225, 0.0114, 0.0234, 0.039, 0.0215, 
            0.0095, 0.0148, 0.0167),
        hOverE = cms.vdouble(0.056, 0.0221, 0.037, 0.0, 0.0268, 
            0.0102, 0.0104, 0.0),
        sigmaEtaEta = cms.vdouble(0.0095, 0.0094, 0.0094, 0.0, 0.026, 
            0.0257, 0.0246, 0.0),
        deltaEtaIn = cms.vdouble(0.0043, 0.00282, 0.0036, 0.0, 0.0066, 
            0.0049, 0.0041, 0.0),
        eSeedOverPin = cms.vdouble(0.32, 0.94, 0.221, 0.0, 0.74, 
            0.89, 0.66, 0.0)
    ),
    classbasedtightEleIDCutsV06 = cms.PSet(
        cutdetain = cms.vdouble(0.0116, 0.00449, 0.00938, 0.0184, 0.00678, 
            0.0109, 0.0252, 0.0268, 0.0139),
        cutiso_sum = cms.vdouble(15.5, 12.2, 12.2, 11.7, 7.16, 
            9.71, 8.66, 11.9, 2.98),
        cutip_gsf = cms.vdouble(0.0131, 0.0586, 0.0839, 0.0366, 0.452, 
            0.204, 0.0913, 0.0802, 0.0731),
        cutip_gsfl = cms.vdouble(0.0119, 0.0527, 0.0471, 0.0212, 0.233, 
            0.267, 0.109, 0.122, 0.0479),
        cuthoe = cms.vdouble(0.215, 0.0608, 0.147, 0.369, 0.0349, 
            0.102, 0.52, 0.422, 0.404),
        cutiso_sumoetl = cms.vdouble(6.21, 6.81, 5.3, 5.39, 2.73, 
            4.73, 4.84, 3.46, 3.73),
        cutfmishits = cms.vdouble(1.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 1.5, 2.5, 0.5),
        cuthoel = cms.vdouble(0.228, 0.0836, 0.143, 0.37, 0.0392, 
            0.0979, 0.3, 0.381, 0.339),
        cutdphiin = cms.vdouble(0.0897, 0.0993, 0.295, 0.0979, 0.151, 
            0.252, 0.341, 0.308, 0.328),
        cutseel = cms.vdouble(0.0132, 0.0117, 0.0112, 0.0387, 0.0281, 
            0.0287, 0.00987, 0.0296, 0.0544),
        cutiso_sumoet = cms.vdouble(11.9, 7.81, 6.28, 8.92, 4.65, 
            5.49, 9.36, 8.84, 5.94),
        cutsee = cms.vdouble(0.0145, 0.0116, 0.012, 0.039, 0.0297, 
            0.0311, 0.00987, 0.0347, 0.0917),
        cuteseedopcor = cms.vdouble(0.637, 0.943, 0.742, 0.748, 0.763, 
            0.631, 0.214, 0.873, 0.473),
        cutdphiinl = cms.vdouble(0.061, 0.14, 0.286, 0.0921, 0.197, 
            0.24, 0.333, 0.303, 0.258),
        cutdetainl = cms.vdouble(0.00816, 0.00401, 0.0081, 0.019, 0.00588, 
            0.00893, 0.0171, 0.0434, 0.0143)
    ),
    classbasedtightEleIDCutsV04 = cms.PSet(
        cutdetain = cms.vdouble(0.00811, 0.00341, 0.00633, 0.0103, 0.00667, 
            0.01, 0.0106, 0.0145, 0.0163, 0.0076, 
            0.00259, 0.00511, 0.00941, 0.0043, 0.00857, 
            0.012, 0.0169, 0.00172, 0.00861, 0.00362, 
            0.00601, 0.00925, 0.00489, 0.00832, 0.0119, 
            0.0169, 0.000996),
        cutiso_sum = cms.vdouble(11.8, 8.31, 6.26, 6.18, 3.28, 
            4.38, 4.17, 5.4, 1.57, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0213, 0.0422, 0.0632, 0.0361, 0.073, 
            0.126, 0.171, 0.119, 0.0372, 0.0131, 
            0.0146, 0.0564, 0.0152, 0.0222, 0.0268, 
            0.0314, 0.0884, 0.00374, 0.00852, 0.00761, 
            0.0143, 0.0106, 0.0127, 0.0119, 0.0123, 
            0.0235, 0.00363),
        cuthoe = cms.vdouble(0.0783, 0.0387, 0.105, 0.118, 0.0227, 
            0.062, 0.13, 2.47, 0.38, 0.0888, 
            0.0503, 0.0955, 0.0741, 0.015, 0.03, 
            0.589, 1.13, 0.612, 0.0494, 0.0461, 
            0.0292, 0.0369, 0.0113, 0.0145, 0.124, 
            2.05, 0.61),
        cutfmishits = cms.vdouble(2.5, 1.5, 1.5, 1.5, 1.5, 
            0.5, 2.5, 0.5, 0.5, 2.5, 
            1.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5, -0.5, 2.5, 1.5, 
            0.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(13.7, 11.6, 7.14, 9.98, 3.52, 
            4.87, 6.24, 7.96, 2.53, 11.2, 
            11.9, 7.88, 8.16, 5.58, 5.03, 
            11.4, 8.15, 5.79, 10.4, 11.1, 
            10.4, 7.47, 5.08, 5.9, 11.8, 
            14.1, 11.7),
        cutdcotdist = cms.vdouble(0.0393, 0.0256, 0.00691, 0.0394, 0.0386, 
            0.039, 0.0325, 0.0384, 0.0382, 0.0245, 
            0.000281, 5.46e-05, 0.0342, 0.0232, 0.00107, 
            0.0178, 0.0193, 0.000758, 0.000108, 0.0248, 
            0.000458, 0.0129, 0.00119, 0.0182, 4.53e-05, 
            0.0189, 0.000928),
        cutsee = cms.vdouble(0.0143, 0.0105, 0.0123, 0.0324, 0.0307, 
            0.0301, 0.0109, 0.027, 0.0292, 0.0133, 
            0.0104, 0.0116, 0.0332, 0.0296, 0.031, 
            0.00981, 0.0307, 0.072, 0.0149, 0.0105, 
            0.011, 0.0342, 0.0307, 0.0303, 0.00954, 
            0.0265, 0.0101),
        cuteseedopcor = cms.vdouble(0.784, 0.366, 0.57, 0.911, 0.298, 
            0.645, 0.51, 0.497, 0.932, 0.835, 
            0.968, 0.969, 0.923, 0.898, 0.98, 
            0.63, 0.971, 1.0, 0.515, 0.963, 
            0.986, 0.823, 0.879, 1.01, 0.931, 
            0.937, 1.05),
        cutdphiin = cms.vdouble(0.0404, 0.0499, 0.263, 0.042, 0.0484, 
            0.241, 0.242, 0.231, 0.286, 0.0552, 
            0.0338, 0.154, 0.0623, 0.0183, 0.0392, 
            0.0547, 0.0588, 0.00654, 0.042, 0.0217, 
            0.0885, 0.0445, 0.0141, 0.0234, 0.065, 
            0.0258, 0.0346),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 13.7, 13.2, 
            13.6, 14.2, 14.1, 13.9, 12.9, 
            14.9, 17.7)
    ),
    electronIDType = cms.string('robust'),
    robusttightEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    electronVersion = cms.string(''),
    robusttightEleIDCutsV00 = cms.PSet(
        barrel = cms.vdouble(0.015, 0.0092, 0.02, 0.0025, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.018, 0.025, 0.02, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusttightEleIDCutsV01 = cms.PSet(
        barrel = cms.vdouble(0.01, 0.0099, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.01, 0.028, 0.02, 0.0066, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    reducedBarrelRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    robusttightEleIDCutsV03 = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    reducedEndcapRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    verticesCollection = cms.InputTag("offlinePrimaryVerticesWithBS"),
    classbasedlooseEleIDCuts = cms.PSet(
        cutdetain = cms.vdouble(0.0137, 0.00678, 0.0241, 0.0187, 0.0161, 
            0.0224, 0.0252, 0.0308, 0.0273),
        cutiso_sum = cms.vdouble(33.0, 17.0, 17.9, 18.8, 8.55, 
            12.5, 17.6, 18.5, 2.98),
        cutip_gsf = cms.vdouble(0.0551, 0.0765, 0.143, 0.0874, 0.594, 
            0.37, 0.0913, 1.15, 0.231),
        cutip_gsfl = cms.vdouble(0.0186, 0.0759, 0.138, 0.0473, 0.62, 
            0.304, 0.109, 0.775, 0.0479),
        cuthoe = cms.vdouble(0.247, 0.137, 0.147, 0.371, 0.0588, 
            0.147, 0.52, 0.452, 0.404),
        cutiso_sumoetl = cms.vdouble(11.3, 9.05, 9.07, 9.94, 5.25, 
            6.15, 10.7, 10.8, 4.4),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 4.5, 3.5, 3.5),
        cuthoel = cms.vdouble(0.236, 0.126, 0.147, 0.375, 0.0392, 
            0.145, 0.365, 0.383, 0.384),
        cutdphiin = cms.vdouble(0.0897, 0.262, 0.353, 0.116, 0.357, 
            0.319, 0.342, 0.404, 0.336),
        cutseel = cms.vdouble(0.0164, 0.0118, 0.015, 0.0523, 0.0326, 
            0.0456, 0.0185, 0.0589, 0.0544),
        cutiso_sumoet = cms.vdouble(34.5, 12.7, 12.1, 19.9, 6.35, 
            8.85, 14.0, 10.5, 9.74),
        cutsee = cms.vdouble(0.0176, 0.0125, 0.0181, 0.0415, 0.0364, 
            0.0418, 0.0146, 0.0678, 0.133),
        cuteseedopcor = cms.vdouble(0.63, 0.82, 0.401, 0.718, 0.4, 
            0.458, 0.15, 0.664, 0.373),
        cutdphiinl = cms.vdouble(0.0747, 0.25, 0.356, 0.0956, 0.347, 
            0.326, 0.333, 0.647, 0.289),
        cutdetainl = cms.vdouble(0.0124, 0.00503, 0.0257, 0.0228, 0.0118, 
            0.0178, 0.0188, 0.14, 0.024)
    ),
    robusttightEleIDCutsV02 = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV01 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.0275, 0.09, 0.007, -1, 
            -1, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV00 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.011, 0.09, 0.005, -1, 
            -1, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.1, 0.0275, 0.09, 0.007, -1, 
            -1, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV03 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV02 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    classbasedlooseEleIDCutsV00 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.05, 0.025, 0.053, 0.09, 0.07, 
            0.03, 0.092, 0.092),
        hOverE = cms.vdouble(0.115, 0.1, 0.055, 0.0, 0.145, 
            0.12, 0.15, 0.0),
        sigmaEtaEta = cms.vdouble(0.014, 0.012, 0.0115, 0.0, 0.0275, 
            0.0265, 0.0265, 0.0),
        deltaEtaIn = cms.vdouble(0.009, 0.0045, 0.0085, 0.0, 0.0105, 
            0.0068, 0.01, 0.0),
        eSeedOverPin = cms.vdouble(0.11, 0.91, 0.11, 0.0, 0.0, 
            0.85, 0.0, 0.0)
    ),
    classbasedlooseEleIDCutsV01 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.053, 0.0189, 0.059, 0.099, 0.0278, 
            0.0157, 0.042, 0.08),
        hOverE = cms.vdouble(0.076, 0.033, 0.07, 0.0, 0.083, 
            0.0148, 0.033, 0.0),
        sigmaEtaEta = cms.vdouble(0.0101, 0.0095, 0.0097, 0.0, 0.0271, 
            0.0267, 0.0259, 0.0),
        deltaEtaIn = cms.vdouble(0.0078, 0.00259, 0.0062, 0.0, 0.0078, 
            0.0061, 0.0061, 0.0),
        eSeedOverPin = cms.vdouble(0.3, 0.92, 0.211, 0.0, 0.42, 
            0.88, 0.68, 0.0)
    ),
    classbasedlooseEleIDCutsV02 = cms.PSet(
        cutisohcal = cms.vdouble(13.5, 9.93, 7.56, 14.8, 8.1, 
            10.8, 42.7, 20.1, 9.11, 10.4, 
            6.89, 5.59, 8.53, 9.59, 24.2, 
            2.78, 8.67, 0.288),
        cutmishits = cms.vdouble(5.5, 1.5, 5.5, 2.5, 2.5, 
            2.5, 3.5, 5.5, 0.5, 1.5, 
            2.5, 0.5, 1.5, 1.5, 0.5, 
            0.5, 0.5, 0.5),
        cuthoe = cms.vdouble(0.0887, 0.0934, 0.0949, 0.0986, 0.0431, 
            0.0878, 0.097, 0.0509, 0.098, 0.0991, 
            0.0321, 0.0928, 0.0663, 0.0717, 0.0966, 
            0.0758, 0.0149, 0.0131),
        cutdeta = cms.vdouble(0.00958, 0.00406, 0.0122, 0.0137, 0.00837, 
            0.0127, 0.011, 0.00336, 0.00977, 0.015, 
            0.00675, 0.0109, 0.014, 0.00508, 0.0109, 
            0.0146, 0.00506, 0.0127),
        cuteopin = cms.vdouble(0.878, 0.802, 0.814, 0.942, 0.735, 
            0.774, 0.829, 0.909, 0.829, 0.813, 
            0.86, 0.897, 0.817, 0.831, 0.818, 
            0.861, 0.787, 0.789),
        cutip = cms.vdouble(0.0246, 0.076, 0.0966, 0.0885, 0.441, 
            0.205, 0.0292, 0.0293, 0.0619, 0.0251, 
            0.159, 0.0815, 7.29, 0.0106, 5.76, 
            6.89, 1.27, 5.89),
        cutisotk = cms.vdouble(24.3, 8.45, 14.4, 27.8, 6.02, 
            10.5, 14.1, 10.2, 14.5, 19.1, 
            6.1, 14.1, 8.59, 8.33, 8.3, 
            8.93, 8.6, 16.0),
        cutsee = cms.vdouble(0.0172, 0.0115, 0.0143, 0.0344, 0.0295, 
            0.0304, 0.0145, 0.0108, 0.0128, 0.0347, 
            0.0307, 0.0316, 0.018, 0.011, 0.0132, 
            0.0349, 0.031, 0.0327),
        cutdphi = cms.vdouble(0.0372, 0.114, 0.118, 0.0488, 0.117, 
            0.119, 0.0606, 0.0548, 0.117, 0.07, 
            0.0355, 0.117, 0.088, 0.045, 0.118, 
            0.0919, 0.0236, 0.0515),
        cutisoecal = cms.vdouble(33.4, 28.1, 7.32, 27.4, 7.33, 
            21.7, 93.8, 102.0, 12.1, 26.0, 
            8.91, 10.0, 16.1, 31.3, 16.9, 
            15.4, 13.3, 37.7)
    ),
    classbasedlooseEleIDCutsV03 = cms.PSet(
        cutdetain = cms.vdouble(0.00989, 0.00484, 0.0146, 0.0146, 0.00902, 
            0.0172, 0.0137, 0.0477, 0.0275, 0.00967, 
            0.00377, 0.00924, 0.013, 0.00666, 0.0123, 
            0.0125, 0.0228, 0.0112, 0.0106, 0.0038, 
            0.00897, 0.0139, 0.00667, 0.0122, 0.0122, 
            0.0193, 0.00239),
        cutiso_sum = cms.vdouble(31.5, 10.3, 8.8, 11.0, 6.13, 
            6.94, 7.52, 9.0, 3.5, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0431, 0.0767, 0.139, 0.101, 0.149, 
            0.154, 0.932, 0.15, 0.124, 0.0238, 
            0.0467, 0.0759, 0.0369, 0.147, 0.0986, 
            0.0626, 0.195, 0.116, 0.0122, 0.0125, 
            0.0693, 0.0162, 0.089, 0.0673, 0.0467, 
            0.0651, 0.0221),
        cuthoe = cms.vdouble(0.166, 0.0771, 0.144, 0.37, 0.0497, 
            0.139, 0.401, 2.68, 0.516, 0.234, 
            0.0556, 0.144, 0.368, 0.031, 0.12, 
            0.602, 2.01, 1.05, 0.104, 0.063, 
            0.0565, 0.38, 0.0192, 0.0294, 0.537, 
            4.65, 1.87),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 2.5, 2.5, 1.5, 2.5, 
            1.5, 1.5, 1.5, 1.5, 0.5, 
            2.5, 2.5, 0.5, 2.5, 1.5, 
            0.5, 1.5, 1.5, 0.5, 2.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(28.9, 15.3, 12.0, 18.3, 7.17, 
            9.42, 11.0, 9.81, 3.94, 22.7, 
            15.9, 12.3, 17.0, 7.58, 8.89, 
            15.2, 12.7, 6.17, 20.8, 21.2, 
            17.2, 15.5, 9.37, 10.6, 19.8, 
            22.1, 15.6),
        cutdcotdist = cms.vdouble(0.0393, 0.0392, 0.0397, 0.0394, 0.0393, 
            0.039, 0.0378, 0.0388, 0.0382, 0.0385, 
            0.0167, 0.00325, 0.0394, 0.0387, 0.0388, 
            0.0227, 0.0258, 0.0127, 0.0298, 0.03, 
            0.00946, 0.039, 0.0231, 0.0278, 0.00162, 
            0.0367, 0.0199),
        cutsee = cms.vdouble(0.0175, 0.0127, 0.0177, 0.0373, 0.0314, 
            0.0329, 0.0157, 0.0409, 0.14, 0.0169, 
            0.0106, 0.0142, 0.0363, 0.0322, 0.0354, 
            0.0117, 0.0372, 28.2, 0.0171, 0.0113, 
            0.014, 0.0403, 0.0323, 0.0411, 0.0104, 
            0.0436, 0.0114),
        cuteseedopcor = cms.vdouble(0.78, 0.302, 0.483, 0.904, 0.168, 
            0.645, 0.108, 0.284, 0.324, 0.591, 
            0.286, 0.488, 0.813, 0.791, 0.672, 
            0.398, 0.834, 0.878, 0.515, 0.937, 
            0.806, 0.816, 0.85, 0.507, 0.367, 
            0.83, 0.648),
        cutdphiin = cms.vdouble(0.041, 0.275, 0.365, 0.047, 0.273, 
            0.296, 0.329, 0.465, 0.627, 0.0581, 
            0.0954, 0.327, 0.0702, 0.0582, 0.279, 
            0.117, 0.318, 0.246, 0.0821, 0.052, 
            0.292, 0.116, 0.0435, 0.312, 0.118, 
            0.296, 0.0459),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 12.0, 12.0, 
            12.0, 12.0, 12.0, 12.0, 12.0, 
            12.0, 12.5)
    ),
    classbasedlooseEleIDCutsV04 = cms.PSet(
        cutdetain = cms.vdouble(0.00989, 0.00484, 0.0146, 0.0146, 0.00902, 
            0.0172, 0.0137, 0.0477, 0.0275, 0.00967, 
            0.00377, 0.00924, 0.013, 0.00666, 0.0123, 
            0.0125, 0.0228, 0.0112, 0.0106, 0.0038, 
            0.00897, 0.0139, 0.00667, 0.0122, 0.0122, 
            0.0193, 0.00239),
        cutiso_sum = cms.vdouble(31.5, 10.3, 8.8, 11.0, 6.13, 
            6.94, 7.52, 9.0, 3.5, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0431, 0.0767, 0.139, 0.101, 0.149, 
            0.154, 0.932, 0.15, 0.124, 0.0238, 
            0.0467, 0.0759, 0.0369, 0.147, 0.0986, 
            0.0626, 0.195, 0.116, 0.0122, 0.0125, 
            0.0693, 0.0162, 0.089, 0.0673, 0.0467, 
            0.0651, 0.0221),
        cuthoe = cms.vdouble(0.166, 0.0771, 0.144, 0.37, 0.0497, 
            0.139, 0.401, 2.68, 0.516, 0.234, 
            0.0556, 0.144, 0.368, 0.031, 0.12, 
            0.602, 2.01, 1.05, 0.104, 0.063, 
            0.0565, 0.38, 0.0192, 0.0294, 0.537, 
            4.65, 1.87),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 2.5, 2.5, 1.5, 2.5, 
            1.5, 1.5, 1.5, 1.5, 0.5, 
            2.5, 2.5, 0.5, 2.5, 1.5, 
            0.5, 1.5, 1.5, 0.5, 2.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(28.9, 15.3, 12.0, 18.3, 7.17, 
            9.42, 11.0, 9.81, 3.94, 22.7, 
            15.9, 12.3, 17.0, 7.58, 8.89, 
            15.2, 12.7, 6.17, 20.8, 21.2, 
            17.2, 15.5, 9.37, 10.6, 19.8, 
            22.1, 15.6),
        cutdcotdist = cms.vdouble(0.0393, 0.0392, 0.0397, 0.0394, 0.0393, 
            0.039, 0.0378, 0.0388, 0.0382, 0.0385, 
            0.0167, 0.00325, 0.0394, 0.0387, 0.0388, 
            0.0227, 0.0258, 0.0127, 0.0298, 0.03, 
            0.00946, 0.039, 0.0231, 0.0278, 0.00162, 
            0.0367, 0.0199),
        cutsee = cms.vdouble(0.0175, 0.0127, 0.0177, 0.0373, 0.0314, 
            0.0329, 0.0157, 0.0409, 0.14, 0.0169, 
            0.0106, 0.0142, 0.0363, 0.0322, 0.0354, 
            0.0117, 0.0372, 28.2, 0.0171, 0.0113, 
            0.014, 0.0403, 0.0323, 0.0411, 0.0104, 
            0.0436, 0.0114),
        cuteseedopcor = cms.vdouble(0.78, 0.302, 0.483, 0.904, 0.168, 
            0.645, 0.108, 0.284, 0.324, 0.591, 
            0.286, 0.488, 0.813, 0.791, 0.672, 
            0.398, 0.834, 0.878, 0.515, 0.937, 
            0.806, 0.816, 0.85, 0.507, 0.367, 
            0.83, 0.648),
        cutdphiin = cms.vdouble(0.041, 0.275, 0.365, 0.047, 0.273, 
            0.296, 0.329, 0.465, 0.627, 0.0581, 
            0.0954, 0.327, 0.0702, 0.0582, 0.279, 
            0.117, 0.318, 0.246, 0.0821, 0.052, 
            0.292, 0.116, 0.0435, 0.312, 0.118, 
            0.296, 0.0459),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 12.0, 12.0, 
            12.0, 12.0, 12.0, 12.0, 12.0, 
            12.0, 12.5)
    ),
    classbasedlooseEleIDCutsV06 = cms.PSet(
        cutdetain = cms.vdouble(0.0137, 0.00678, 0.0241, 0.0187, 0.0161, 
            0.0224, 0.0252, 0.0308, 0.0273),
        cutiso_sum = cms.vdouble(33.0, 17.0, 17.9, 18.8, 8.55, 
            12.5, 17.6, 18.5, 2.98),
        cutip_gsf = cms.vdouble(0.0551, 0.0765, 0.143, 0.0874, 0.594, 
            0.37, 0.0913, 1.15, 0.231),
        cutip_gsfl = cms.vdouble(0.0186, 0.0759, 0.138, 0.0473, 0.62, 
            0.304, 0.109, 0.775, 0.0479),
        cuthoe = cms.vdouble(0.247, 0.137, 0.147, 0.371, 0.0588, 
            0.147, 0.52, 0.452, 0.404),
        cutiso_sumoetl = cms.vdouble(11.3, 9.05, 9.07, 9.94, 5.25, 
            6.15, 10.7, 10.8, 4.4),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 4.5, 3.5, 3.5),
        cuthoel = cms.vdouble(0.236, 0.126, 0.147, 0.375, 0.0392, 
            0.145, 0.365, 0.383, 0.384),
        cutdphiin = cms.vdouble(0.0897, 0.262, 0.353, 0.116, 0.357, 
            0.319, 0.342, 0.404, 0.336),
        cutseel = cms.vdouble(0.0164, 0.0118, 0.015, 0.0523, 0.0326, 
            0.0456, 0.0185, 0.0589, 0.0544),
        cutiso_sumoet = cms.vdouble(34.5, 12.7, 12.1, 19.9, 6.35, 
            8.85, 14.0, 10.5, 9.74),
        cutsee = cms.vdouble(0.0176, 0.0125, 0.0181, 0.0415, 0.0364, 
            0.0418, 0.0146, 0.0678, 0.133),
        cuteseedopcor = cms.vdouble(0.63, 0.82, 0.401, 0.718, 0.4, 
            0.458, 0.15, 0.664, 0.373),
        cutdphiinl = cms.vdouble(0.0747, 0.25, 0.356, 0.0956, 0.347, 
            0.326, 0.333, 0.647, 0.289),
        cutdetainl = cms.vdouble(0.0124, 0.00503, 0.0257, 0.0228, 0.0118, 
            0.0178, 0.0188, 0.14, 0.024)
    ),
    src = cms.InputTag("gsfElectrons"),
    robusttightEleIDCuts = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    classbasedtightEleIDCuts = cms.PSet(
        cutdetain = cms.vdouble(0.0116, 0.00449, 0.00938, 0.0184, 0.00678, 
            0.0109, 0.0252, 0.0268, 0.0139),
        cutiso_sum = cms.vdouble(15.5, 12.2, 12.2, 11.7, 7.16, 
            9.71, 8.66, 11.9, 2.98),
        cutip_gsf = cms.vdouble(0.0131, 0.0586, 0.0839, 0.0366, 0.452, 
            0.204, 0.0913, 0.0802, 0.0731),
        cutip_gsfl = cms.vdouble(0.0119, 0.0527, 0.0471, 0.0212, 0.233, 
            0.267, 0.109, 0.122, 0.0479),
        cuthoe = cms.vdouble(0.215, 0.0608, 0.147, 0.369, 0.0349, 
            0.102, 0.52, 0.422, 0.404),
        cutiso_sumoetl = cms.vdouble(6.21, 6.81, 5.3, 5.39, 2.73, 
            4.73, 4.84, 3.46, 3.73),
        cutfmishits = cms.vdouble(1.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 1.5, 2.5, 0.5),
        cuthoel = cms.vdouble(0.228, 0.0836, 0.143, 0.37, 0.0392, 
            0.0979, 0.3, 0.381, 0.339),
        cutdphiin = cms.vdouble(0.0897, 0.0993, 0.295, 0.0979, 0.151, 
            0.252, 0.341, 0.308, 0.328),
        cutseel = cms.vdouble(0.0132, 0.0117, 0.0112, 0.0387, 0.0281, 
            0.0287, 0.00987, 0.0296, 0.0544),
        cutiso_sumoet = cms.vdouble(11.9, 7.81, 6.28, 8.92, 4.65, 
            5.49, 9.36, 8.84, 5.94),
        cutsee = cms.vdouble(0.0145, 0.0116, 0.012, 0.039, 0.0297, 
            0.0311, 0.00987, 0.0347, 0.0917),
        cuteseedopcor = cms.vdouble(0.637, 0.943, 0.742, 0.748, 0.763, 
            0.631, 0.214, 0.873, 0.473),
        cutdphiinl = cms.vdouble(0.061, 0.14, 0.286, 0.0921, 0.197, 
            0.24, 0.333, 0.303, 0.258),
        cutdetainl = cms.vdouble(0.00816, 0.00401, 0.0081, 0.019, 0.00588, 
            0.00893, 0.0171, 0.0434, 0.0143)
    ),
    algorithm = cms.string('eIDCB'),
    robusthighenergyEleIDCuts = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCuts = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV02 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV03 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV00 = cms.PSet(
        barrel = cms.vdouble(0.115, 0.014, 0.09, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.15, 0.0275, 0.092, 0.0105, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV01 = cms.PSet(
        barrel = cms.vdouble(0.075, 0.0132, 0.058, 0.0077, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.083, 0.027, 0.042, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    additionalCategories = cms.bool(True),
    etBinning = cms.bool(True)
)


process.eidRobustLoose = cms.EDProducer("EleIdCutBasedExtProducer",
    electronQuality = cms.string('loose'),
    classbasedtightEleIDCutsV02 = cms.PSet(
        cutisohcal = cms.vdouble(10.9, 7.01, 8.75, 3.51, 7.75, 
            1.62, 11.6, 9.9, 4.97, 5.33, 
            3.18, 2.32, 0.164, 5.46, 12.0, 
            0.00604, 4.1, 0.000628),
        cutmishits = cms.vdouble(5.5, 1.5, 0.5, 1.5, 2.5, 
            0.5, 3.5, 5.5, 0.5, 0.5, 
            0.5, 0.5, 0.5, 1.5, 0.5, 
            0.5, 0.5, 0.5),
        cuthoe = cms.vdouble(0.0871, 0.0289, 0.0783, 0.0946, 0.0245, 
            0.0363, 0.0671, 0.048, 0.0614, 0.0924, 
            0.0158, 0.049, 0.0382, 0.0915, 0.0451, 
            0.0452, 0.00196, 0.0043),
        cutdeta = cms.vdouble(0.00915, 0.00302, 0.0061, 0.0135, 0.00565, 
            0.00793, 0.0102, 0.00266, 0.0106, 0.00903, 
            0.00766, 0.00723, 0.0116, 0.00203, 0.00659, 
            0.0148, 0.00555, 0.0128),
        cuteopin = cms.vdouble(0.878, 0.859, 0.874, 0.944, 0.737, 
            0.773, 0.86, 0.967, 0.917, 0.812, 
            0.915, 1.01, 0.847, 0.953, 0.979, 
            0.841, 0.771, 1.09),
        cutip = cms.vdouble(0.0239, 0.027, 0.0768, 0.0231, 0.178, 
            0.0957, 0.0102, 0.0168, 0.043, 0.0166, 
            0.0594, 0.0308, 2.1, 0.00527, 3.17, 
            4.91, 0.769, 5.9),
        cutisotk = cms.vdouble(6.53, 4.6, 6.0, 8.63, 3.11, 
            7.77, 5.42, 4.81, 4.06, 6.47, 
            2.8, 3.45, 5.29, 5.18, 15.4, 
            5.38, 4.47, 0.0347),
        cutsee = cms.vdouble(0.0131, 0.0106, 0.0115, 0.0306, 0.028, 
            0.0293, 0.0131, 0.0106, 0.0115, 0.0317, 
            0.029, 0.0289, 0.0142, 0.0106, 0.0103, 
            0.035, 0.0296, 0.0333),
        cutdphi = cms.vdouble(0.0369, 0.0307, 0.117, 0.0475, 0.0216, 
            0.117, 0.0372, 0.0246, 0.0426, 0.0612, 
            0.0142, 0.039, 0.0737, 0.0566, 0.0359, 
            0.0187, 0.012, 0.0358),
        cutisoecal = cms.vdouble(20.0, 27.2, 4.48, 13.5, 4.56, 
            3.19, 12.2, 13.1, 7.42, 7.67, 
            4.12, 4.85, 10.1, 12.4, 11.1, 
            11.0, 10.6, 13.4)
    ),
    classbasedtightEleIDCutsV03 = cms.PSet(
        cutdetain = cms.vdouble(0.00811, 0.00341, 0.00633, 0.0103, 0.00667, 
            0.01, 0.0106, 0.0145, 0.0163, 0.0076, 
            0.00259, 0.00511, 0.00941, 0.0043, 0.00857, 
            0.012, 0.0169, 0.00172, 0.00861, 0.00362, 
            0.00601, 0.00925, 0.00489, 0.00832, 0.0119, 
            0.0169, 0.000996),
        cutiso_sum = cms.vdouble(11.8, 8.31, 6.26, 6.18, 3.28, 
            4.38, 4.17, 5.4, 1.57, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0213, 0.0422, 0.0632, 0.0361, 0.073, 
            0.126, 0.171, 0.119, 0.0372, 0.0131, 
            0.0146, 0.0564, 0.0152, 0.0222, 0.0268, 
            0.0314, 0.0884, 0.00374, 0.00852, 0.00761, 
            0.0143, 0.0106, 0.0127, 0.0119, 0.0123, 
            0.0235, 0.00363),
        cuthoe = cms.vdouble(0.0783, 0.0387, 0.105, 0.118, 0.0227, 
            0.062, 0.13, 2.47, 0.38, 0.0888, 
            0.0503, 0.0955, 0.0741, 0.015, 0.03, 
            0.589, 1.13, 0.612, 0.0494, 0.0461, 
            0.0292, 0.0369, 0.0113, 0.0145, 0.124, 
            2.05, 0.61),
        cutfmishits = cms.vdouble(2.5, 1.5, 1.5, 1.5, 1.5, 
            0.5, 2.5, 0.5, 0.5, 2.5, 
            1.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5, -0.5, 2.5, 1.5, 
            0.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(13.7, 11.6, 7.14, 9.98, 3.52, 
            4.87, 6.24, 7.96, 2.53, 11.2, 
            11.9, 7.88, 8.16, 5.58, 5.03, 
            11.4, 8.15, 5.79, 10.4, 11.1, 
            10.4, 7.47, 5.08, 5.9, 11.8, 
            14.1, 11.7),
        cutdcotdist = cms.vdouble(0.0393, 0.0256, 0.00691, 0.0394, 0.0386, 
            0.039, 0.0325, 0.0384, 0.0382, 0.0245, 
            0.000281, 5.46e-05, 0.0342, 0.0232, 0.00107, 
            0.0178, 0.0193, 0.000758, 0.000108, 0.0248, 
            0.000458, 0.0129, 0.00119, 0.0182, 4.53e-05, 
            0.0189, 0.000928),
        cutsee = cms.vdouble(0.0143, 0.0105, 0.0123, 0.0324, 0.0307, 
            0.0301, 0.0109, 0.027, 0.0292, 0.0133, 
            0.0104, 0.0116, 0.0332, 0.0296, 0.031, 
            0.00981, 0.0307, 0.072, 0.0149, 0.0105, 
            0.011, 0.0342, 0.0307, 0.0303, 0.00954, 
            0.0265, 0.0101),
        cuteseedopcor = cms.vdouble(0.784, 0.366, 0.57, 0.911, 0.298, 
            0.645, 0.51, 0.497, 0.932, 0.835, 
            0.968, 0.969, 0.923, 0.898, 0.98, 
            0.63, 0.971, 1.0, 0.515, 0.963, 
            0.986, 0.823, 0.879, 1.01, 0.931, 
            0.937, 1.05),
        cutdphiin = cms.vdouble(0.0404, 0.0499, 0.263, 0.042, 0.0484, 
            0.241, 0.242, 0.231, 0.286, 0.0552, 
            0.0338, 0.154, 0.0623, 0.0183, 0.0392, 
            0.0547, 0.0588, 0.00654, 0.042, 0.0217, 
            0.0885, 0.0445, 0.0141, 0.0234, 0.065, 
            0.0258, 0.0346),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 13.7, 13.2, 
            13.6, 14.2, 14.1, 13.9, 12.9, 
            14.9, 17.7)
    ),
    classbasedtightEleIDCutsV00 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.032, 0.016, 0.0525, 0.09, 0.025, 
            0.035, 0.065, 0.092),
        hOverE = cms.vdouble(0.05, 0.042, 0.045, 0.0, 0.055, 
            0.037, 0.05, 0.0),
        sigmaEtaEta = cms.vdouble(0.0125, 0.011, 0.01, 0.0, 0.0265, 
            0.0252, 0.026, 0.0),
        deltaEtaIn = cms.vdouble(0.0055, 0.003, 0.0065, 0.0, 0.006, 
            0.0055, 0.0075, 0.0),
        eSeedOverPin = cms.vdouble(0.24, 0.94, 0.11, 0.0, 0.32, 
            0.83, 0.0, 0.0)
    ),
    classbasedtightEleIDCutsV01 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.0225, 0.0114, 0.0234, 0.039, 0.0215, 
            0.0095, 0.0148, 0.0167),
        hOverE = cms.vdouble(0.056, 0.0221, 0.037, 0.0, 0.0268, 
            0.0102, 0.0104, 0.0),
        sigmaEtaEta = cms.vdouble(0.0095, 0.0094, 0.0094, 0.0, 0.026, 
            0.0257, 0.0246, 0.0),
        deltaEtaIn = cms.vdouble(0.0043, 0.00282, 0.0036, 0.0, 0.0066, 
            0.0049, 0.0041, 0.0),
        eSeedOverPin = cms.vdouble(0.32, 0.94, 0.221, 0.0, 0.74, 
            0.89, 0.66, 0.0)
    ),
    classbasedtightEleIDCutsV06 = cms.PSet(
        cutdetain = cms.vdouble(0.0116, 0.00449, 0.00938, 0.0184, 0.00678, 
            0.0109, 0.0252, 0.0268, 0.0139),
        cutiso_sum = cms.vdouble(15.5, 12.2, 12.2, 11.7, 7.16, 
            9.71, 8.66, 11.9, 2.98),
        cutip_gsf = cms.vdouble(0.0131, 0.0586, 0.0839, 0.0366, 0.452, 
            0.204, 0.0913, 0.0802, 0.0731),
        cutip_gsfl = cms.vdouble(0.0119, 0.0527, 0.0471, 0.0212, 0.233, 
            0.267, 0.109, 0.122, 0.0479),
        cuthoe = cms.vdouble(0.215, 0.0608, 0.147, 0.369, 0.0349, 
            0.102, 0.52, 0.422, 0.404),
        cutiso_sumoetl = cms.vdouble(6.21, 6.81, 5.3, 5.39, 2.73, 
            4.73, 4.84, 3.46, 3.73),
        cutfmishits = cms.vdouble(1.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 1.5, 2.5, 0.5),
        cuthoel = cms.vdouble(0.228, 0.0836, 0.143, 0.37, 0.0392, 
            0.0979, 0.3, 0.381, 0.339),
        cutdphiin = cms.vdouble(0.0897, 0.0993, 0.295, 0.0979, 0.151, 
            0.252, 0.341, 0.308, 0.328),
        cutseel = cms.vdouble(0.0132, 0.0117, 0.0112, 0.0387, 0.0281, 
            0.0287, 0.00987, 0.0296, 0.0544),
        cutiso_sumoet = cms.vdouble(11.9, 7.81, 6.28, 8.92, 4.65, 
            5.49, 9.36, 8.84, 5.94),
        cutsee = cms.vdouble(0.0145, 0.0116, 0.012, 0.039, 0.0297, 
            0.0311, 0.00987, 0.0347, 0.0917),
        cuteseedopcor = cms.vdouble(0.637, 0.943, 0.742, 0.748, 0.763, 
            0.631, 0.214, 0.873, 0.473),
        cutdphiinl = cms.vdouble(0.061, 0.14, 0.286, 0.0921, 0.197, 
            0.24, 0.333, 0.303, 0.258),
        cutdetainl = cms.vdouble(0.00816, 0.00401, 0.0081, 0.019, 0.00588, 
            0.00893, 0.0171, 0.0434, 0.0143)
    ),
    classbasedtightEleIDCutsV04 = cms.PSet(
        cutdetain = cms.vdouble(0.00811, 0.00341, 0.00633, 0.0103, 0.00667, 
            0.01, 0.0106, 0.0145, 0.0163, 0.0076, 
            0.00259, 0.00511, 0.00941, 0.0043, 0.00857, 
            0.012, 0.0169, 0.00172, 0.00861, 0.00362, 
            0.00601, 0.00925, 0.00489, 0.00832, 0.0119, 
            0.0169, 0.000996),
        cutiso_sum = cms.vdouble(11.8, 8.31, 6.26, 6.18, 3.28, 
            4.38, 4.17, 5.4, 1.57, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0213, 0.0422, 0.0632, 0.0361, 0.073, 
            0.126, 0.171, 0.119, 0.0372, 0.0131, 
            0.0146, 0.0564, 0.0152, 0.0222, 0.0268, 
            0.0314, 0.0884, 0.00374, 0.00852, 0.00761, 
            0.0143, 0.0106, 0.0127, 0.0119, 0.0123, 
            0.0235, 0.00363),
        cuthoe = cms.vdouble(0.0783, 0.0387, 0.105, 0.118, 0.0227, 
            0.062, 0.13, 2.47, 0.38, 0.0888, 
            0.0503, 0.0955, 0.0741, 0.015, 0.03, 
            0.589, 1.13, 0.612, 0.0494, 0.0461, 
            0.0292, 0.0369, 0.0113, 0.0145, 0.124, 
            2.05, 0.61),
        cutfmishits = cms.vdouble(2.5, 1.5, 1.5, 1.5, 1.5, 
            0.5, 2.5, 0.5, 0.5, 2.5, 
            1.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5, -0.5, 2.5, 1.5, 
            0.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(13.7, 11.6, 7.14, 9.98, 3.52, 
            4.87, 6.24, 7.96, 2.53, 11.2, 
            11.9, 7.88, 8.16, 5.58, 5.03, 
            11.4, 8.15, 5.79, 10.4, 11.1, 
            10.4, 7.47, 5.08, 5.9, 11.8, 
            14.1, 11.7),
        cutdcotdist = cms.vdouble(0.0393, 0.0256, 0.00691, 0.0394, 0.0386, 
            0.039, 0.0325, 0.0384, 0.0382, 0.0245, 
            0.000281, 5.46e-05, 0.0342, 0.0232, 0.00107, 
            0.0178, 0.0193, 0.000758, 0.000108, 0.0248, 
            0.000458, 0.0129, 0.00119, 0.0182, 4.53e-05, 
            0.0189, 0.000928),
        cutsee = cms.vdouble(0.0143, 0.0105, 0.0123, 0.0324, 0.0307, 
            0.0301, 0.0109, 0.027, 0.0292, 0.0133, 
            0.0104, 0.0116, 0.0332, 0.0296, 0.031, 
            0.00981, 0.0307, 0.072, 0.0149, 0.0105, 
            0.011, 0.0342, 0.0307, 0.0303, 0.00954, 
            0.0265, 0.0101),
        cuteseedopcor = cms.vdouble(0.784, 0.366, 0.57, 0.911, 0.298, 
            0.645, 0.51, 0.497, 0.932, 0.835, 
            0.968, 0.969, 0.923, 0.898, 0.98, 
            0.63, 0.971, 1.0, 0.515, 0.963, 
            0.986, 0.823, 0.879, 1.01, 0.931, 
            0.937, 1.05),
        cutdphiin = cms.vdouble(0.0404, 0.0499, 0.263, 0.042, 0.0484, 
            0.241, 0.242, 0.231, 0.286, 0.0552, 
            0.0338, 0.154, 0.0623, 0.0183, 0.0392, 
            0.0547, 0.0588, 0.00654, 0.042, 0.0217, 
            0.0885, 0.0445, 0.0141, 0.0234, 0.065, 
            0.0258, 0.0346),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 13.7, 13.2, 
            13.6, 14.2, 14.1, 13.9, 12.9, 
            14.9, 17.7)
    ),
    electronIDType = cms.string('robust'),
    robusttightEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    electronVersion = cms.string(''),
    robusttightEleIDCutsV00 = cms.PSet(
        barrel = cms.vdouble(0.015, 0.0092, 0.02, 0.0025, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.018, 0.025, 0.02, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusttightEleIDCutsV01 = cms.PSet(
        barrel = cms.vdouble(0.01, 0.0099, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.01, 0.028, 0.02, 0.0066, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    reducedBarrelRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    robusttightEleIDCutsV03 = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    reducedEndcapRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    verticesCollection = cms.InputTag("offlinePrimaryVerticesWithBS"),
    classbasedlooseEleIDCuts = cms.PSet(
        cutdetain = cms.vdouble(0.0137, 0.00678, 0.0241, 0.0187, 0.0161, 
            0.0224, 0.0252, 0.0308, 0.0273),
        cutiso_sum = cms.vdouble(33.0, 17.0, 17.9, 18.8, 8.55, 
            12.5, 17.6, 18.5, 2.98),
        cutip_gsf = cms.vdouble(0.0551, 0.0765, 0.143, 0.0874, 0.594, 
            0.37, 0.0913, 1.15, 0.231),
        cutip_gsfl = cms.vdouble(0.0186, 0.0759, 0.138, 0.0473, 0.62, 
            0.304, 0.109, 0.775, 0.0479),
        cuthoe = cms.vdouble(0.247, 0.137, 0.147, 0.371, 0.0588, 
            0.147, 0.52, 0.452, 0.404),
        cutiso_sumoetl = cms.vdouble(11.3, 9.05, 9.07, 9.94, 5.25, 
            6.15, 10.7, 10.8, 4.4),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 4.5, 3.5, 3.5),
        cuthoel = cms.vdouble(0.236, 0.126, 0.147, 0.375, 0.0392, 
            0.145, 0.365, 0.383, 0.384),
        cutdphiin = cms.vdouble(0.0897, 0.262, 0.353, 0.116, 0.357, 
            0.319, 0.342, 0.404, 0.336),
        cutseel = cms.vdouble(0.0164, 0.0118, 0.015, 0.0523, 0.0326, 
            0.0456, 0.0185, 0.0589, 0.0544),
        cutiso_sumoet = cms.vdouble(34.5, 12.7, 12.1, 19.9, 6.35, 
            8.85, 14.0, 10.5, 9.74),
        cutsee = cms.vdouble(0.0176, 0.0125, 0.0181, 0.0415, 0.0364, 
            0.0418, 0.0146, 0.0678, 0.133),
        cuteseedopcor = cms.vdouble(0.63, 0.82, 0.401, 0.718, 0.4, 
            0.458, 0.15, 0.664, 0.373),
        cutdphiinl = cms.vdouble(0.0747, 0.25, 0.356, 0.0956, 0.347, 
            0.326, 0.333, 0.647, 0.289),
        cutdetainl = cms.vdouble(0.0124, 0.00503, 0.0257, 0.0228, 0.0118, 
            0.0178, 0.0188, 0.14, 0.024)
    ),
    robusttightEleIDCutsV02 = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV01 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.0275, 0.09, 0.007, -1, 
            -1, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV00 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.011, 0.09, 0.005, -1, 
            -1, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.1, 0.0275, 0.09, 0.007, -1, 
            -1, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV03 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV02 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    classbasedlooseEleIDCutsV00 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.05, 0.025, 0.053, 0.09, 0.07, 
            0.03, 0.092, 0.092),
        hOverE = cms.vdouble(0.115, 0.1, 0.055, 0.0, 0.145, 
            0.12, 0.15, 0.0),
        sigmaEtaEta = cms.vdouble(0.014, 0.012, 0.0115, 0.0, 0.0275, 
            0.0265, 0.0265, 0.0),
        deltaEtaIn = cms.vdouble(0.009, 0.0045, 0.0085, 0.0, 0.0105, 
            0.0068, 0.01, 0.0),
        eSeedOverPin = cms.vdouble(0.11, 0.91, 0.11, 0.0, 0.0, 
            0.85, 0.0, 0.0)
    ),
    classbasedlooseEleIDCutsV01 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.053, 0.0189, 0.059, 0.099, 0.0278, 
            0.0157, 0.042, 0.08),
        hOverE = cms.vdouble(0.076, 0.033, 0.07, 0.0, 0.083, 
            0.0148, 0.033, 0.0),
        sigmaEtaEta = cms.vdouble(0.0101, 0.0095, 0.0097, 0.0, 0.0271, 
            0.0267, 0.0259, 0.0),
        deltaEtaIn = cms.vdouble(0.0078, 0.00259, 0.0062, 0.0, 0.0078, 
            0.0061, 0.0061, 0.0),
        eSeedOverPin = cms.vdouble(0.3, 0.92, 0.211, 0.0, 0.42, 
            0.88, 0.68, 0.0)
    ),
    classbasedlooseEleIDCutsV02 = cms.PSet(
        cutisohcal = cms.vdouble(13.5, 9.93, 7.56, 14.8, 8.1, 
            10.8, 42.7, 20.1, 9.11, 10.4, 
            6.89, 5.59, 8.53, 9.59, 24.2, 
            2.78, 8.67, 0.288),
        cutmishits = cms.vdouble(5.5, 1.5, 5.5, 2.5, 2.5, 
            2.5, 3.5, 5.5, 0.5, 1.5, 
            2.5, 0.5, 1.5, 1.5, 0.5, 
            0.5, 0.5, 0.5),
        cuthoe = cms.vdouble(0.0887, 0.0934, 0.0949, 0.0986, 0.0431, 
            0.0878, 0.097, 0.0509, 0.098, 0.0991, 
            0.0321, 0.0928, 0.0663, 0.0717, 0.0966, 
            0.0758, 0.0149, 0.0131),
        cutdeta = cms.vdouble(0.00958, 0.00406, 0.0122, 0.0137, 0.00837, 
            0.0127, 0.011, 0.00336, 0.00977, 0.015, 
            0.00675, 0.0109, 0.014, 0.00508, 0.0109, 
            0.0146, 0.00506, 0.0127),
        cuteopin = cms.vdouble(0.878, 0.802, 0.814, 0.942, 0.735, 
            0.774, 0.829, 0.909, 0.829, 0.813, 
            0.86, 0.897, 0.817, 0.831, 0.818, 
            0.861, 0.787, 0.789),
        cutip = cms.vdouble(0.0246, 0.076, 0.0966, 0.0885, 0.441, 
            0.205, 0.0292, 0.0293, 0.0619, 0.0251, 
            0.159, 0.0815, 7.29, 0.0106, 5.76, 
            6.89, 1.27, 5.89),
        cutisotk = cms.vdouble(24.3, 8.45, 14.4, 27.8, 6.02, 
            10.5, 14.1, 10.2, 14.5, 19.1, 
            6.1, 14.1, 8.59, 8.33, 8.3, 
            8.93, 8.6, 16.0),
        cutsee = cms.vdouble(0.0172, 0.0115, 0.0143, 0.0344, 0.0295, 
            0.0304, 0.0145, 0.0108, 0.0128, 0.0347, 
            0.0307, 0.0316, 0.018, 0.011, 0.0132, 
            0.0349, 0.031, 0.0327),
        cutdphi = cms.vdouble(0.0372, 0.114, 0.118, 0.0488, 0.117, 
            0.119, 0.0606, 0.0548, 0.117, 0.07, 
            0.0355, 0.117, 0.088, 0.045, 0.118, 
            0.0919, 0.0236, 0.0515),
        cutisoecal = cms.vdouble(33.4, 28.1, 7.32, 27.4, 7.33, 
            21.7, 93.8, 102.0, 12.1, 26.0, 
            8.91, 10.0, 16.1, 31.3, 16.9, 
            15.4, 13.3, 37.7)
    ),
    classbasedlooseEleIDCutsV03 = cms.PSet(
        cutdetain = cms.vdouble(0.00989, 0.00484, 0.0146, 0.0146, 0.00902, 
            0.0172, 0.0137, 0.0477, 0.0275, 0.00967, 
            0.00377, 0.00924, 0.013, 0.00666, 0.0123, 
            0.0125, 0.0228, 0.0112, 0.0106, 0.0038, 
            0.00897, 0.0139, 0.00667, 0.0122, 0.0122, 
            0.0193, 0.00239),
        cutiso_sum = cms.vdouble(31.5, 10.3, 8.8, 11.0, 6.13, 
            6.94, 7.52, 9.0, 3.5, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0431, 0.0767, 0.139, 0.101, 0.149, 
            0.154, 0.932, 0.15, 0.124, 0.0238, 
            0.0467, 0.0759, 0.0369, 0.147, 0.0986, 
            0.0626, 0.195, 0.116, 0.0122, 0.0125, 
            0.0693, 0.0162, 0.089, 0.0673, 0.0467, 
            0.0651, 0.0221),
        cuthoe = cms.vdouble(0.166, 0.0771, 0.144, 0.37, 0.0497, 
            0.139, 0.401, 2.68, 0.516, 0.234, 
            0.0556, 0.144, 0.368, 0.031, 0.12, 
            0.602, 2.01, 1.05, 0.104, 0.063, 
            0.0565, 0.38, 0.0192, 0.0294, 0.537, 
            4.65, 1.87),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 2.5, 2.5, 1.5, 2.5, 
            1.5, 1.5, 1.5, 1.5, 0.5, 
            2.5, 2.5, 0.5, 2.5, 1.5, 
            0.5, 1.5, 1.5, 0.5, 2.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(28.9, 15.3, 12.0, 18.3, 7.17, 
            9.42, 11.0, 9.81, 3.94, 22.7, 
            15.9, 12.3, 17.0, 7.58, 8.89, 
            15.2, 12.7, 6.17, 20.8, 21.2, 
            17.2, 15.5, 9.37, 10.6, 19.8, 
            22.1, 15.6),
        cutdcotdist = cms.vdouble(0.0393, 0.0392, 0.0397, 0.0394, 0.0393, 
            0.039, 0.0378, 0.0388, 0.0382, 0.0385, 
            0.0167, 0.00325, 0.0394, 0.0387, 0.0388, 
            0.0227, 0.0258, 0.0127, 0.0298, 0.03, 
            0.00946, 0.039, 0.0231, 0.0278, 0.00162, 
            0.0367, 0.0199),
        cutsee = cms.vdouble(0.0175, 0.0127, 0.0177, 0.0373, 0.0314, 
            0.0329, 0.0157, 0.0409, 0.14, 0.0169, 
            0.0106, 0.0142, 0.0363, 0.0322, 0.0354, 
            0.0117, 0.0372, 28.2, 0.0171, 0.0113, 
            0.014, 0.0403, 0.0323, 0.0411, 0.0104, 
            0.0436, 0.0114),
        cuteseedopcor = cms.vdouble(0.78, 0.302, 0.483, 0.904, 0.168, 
            0.645, 0.108, 0.284, 0.324, 0.591, 
            0.286, 0.488, 0.813, 0.791, 0.672, 
            0.398, 0.834, 0.878, 0.515, 0.937, 
            0.806, 0.816, 0.85, 0.507, 0.367, 
            0.83, 0.648),
        cutdphiin = cms.vdouble(0.041, 0.275, 0.365, 0.047, 0.273, 
            0.296, 0.329, 0.465, 0.627, 0.0581, 
            0.0954, 0.327, 0.0702, 0.0582, 0.279, 
            0.117, 0.318, 0.246, 0.0821, 0.052, 
            0.292, 0.116, 0.0435, 0.312, 0.118, 
            0.296, 0.0459),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 12.0, 12.0, 
            12.0, 12.0, 12.0, 12.0, 12.0, 
            12.0, 12.5)
    ),
    classbasedlooseEleIDCutsV04 = cms.PSet(
        cutdetain = cms.vdouble(0.00989, 0.00484, 0.0146, 0.0146, 0.00902, 
            0.0172, 0.0137, 0.0477, 0.0275, 0.00967, 
            0.00377, 0.00924, 0.013, 0.00666, 0.0123, 
            0.0125, 0.0228, 0.0112, 0.0106, 0.0038, 
            0.00897, 0.0139, 0.00667, 0.0122, 0.0122, 
            0.0193, 0.00239),
        cutiso_sum = cms.vdouble(31.5, 10.3, 8.8, 11.0, 6.13, 
            6.94, 7.52, 9.0, 3.5, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0431, 0.0767, 0.139, 0.101, 0.149, 
            0.154, 0.932, 0.15, 0.124, 0.0238, 
            0.0467, 0.0759, 0.0369, 0.147, 0.0986, 
            0.0626, 0.195, 0.116, 0.0122, 0.0125, 
            0.0693, 0.0162, 0.089, 0.0673, 0.0467, 
            0.0651, 0.0221),
        cuthoe = cms.vdouble(0.166, 0.0771, 0.144, 0.37, 0.0497, 
            0.139, 0.401, 2.68, 0.516, 0.234, 
            0.0556, 0.144, 0.368, 0.031, 0.12, 
            0.602, 2.01, 1.05, 0.104, 0.063, 
            0.0565, 0.38, 0.0192, 0.0294, 0.537, 
            4.65, 1.87),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 2.5, 2.5, 1.5, 2.5, 
            1.5, 1.5, 1.5, 1.5, 0.5, 
            2.5, 2.5, 0.5, 2.5, 1.5, 
            0.5, 1.5, 1.5, 0.5, 2.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(28.9, 15.3, 12.0, 18.3, 7.17, 
            9.42, 11.0, 9.81, 3.94, 22.7, 
            15.9, 12.3, 17.0, 7.58, 8.89, 
            15.2, 12.7, 6.17, 20.8, 21.2, 
            17.2, 15.5, 9.37, 10.6, 19.8, 
            22.1, 15.6),
        cutdcotdist = cms.vdouble(0.0393, 0.0392, 0.0397, 0.0394, 0.0393, 
            0.039, 0.0378, 0.0388, 0.0382, 0.0385, 
            0.0167, 0.00325, 0.0394, 0.0387, 0.0388, 
            0.0227, 0.0258, 0.0127, 0.0298, 0.03, 
            0.00946, 0.039, 0.0231, 0.0278, 0.00162, 
            0.0367, 0.0199),
        cutsee = cms.vdouble(0.0175, 0.0127, 0.0177, 0.0373, 0.0314, 
            0.0329, 0.0157, 0.0409, 0.14, 0.0169, 
            0.0106, 0.0142, 0.0363, 0.0322, 0.0354, 
            0.0117, 0.0372, 28.2, 0.0171, 0.0113, 
            0.014, 0.0403, 0.0323, 0.0411, 0.0104, 
            0.0436, 0.0114),
        cuteseedopcor = cms.vdouble(0.78, 0.302, 0.483, 0.904, 0.168, 
            0.645, 0.108, 0.284, 0.324, 0.591, 
            0.286, 0.488, 0.813, 0.791, 0.672, 
            0.398, 0.834, 0.878, 0.515, 0.937, 
            0.806, 0.816, 0.85, 0.507, 0.367, 
            0.83, 0.648),
        cutdphiin = cms.vdouble(0.041, 0.275, 0.365, 0.047, 0.273, 
            0.296, 0.329, 0.465, 0.627, 0.0581, 
            0.0954, 0.327, 0.0702, 0.0582, 0.279, 
            0.117, 0.318, 0.246, 0.0821, 0.052, 
            0.292, 0.116, 0.0435, 0.312, 0.118, 
            0.296, 0.0459),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 12.0, 12.0, 
            12.0, 12.0, 12.0, 12.0, 12.0, 
            12.0, 12.5)
    ),
    classbasedlooseEleIDCutsV06 = cms.PSet(
        cutdetain = cms.vdouble(0.0137, 0.00678, 0.0241, 0.0187, 0.0161, 
            0.0224, 0.0252, 0.0308, 0.0273),
        cutiso_sum = cms.vdouble(33.0, 17.0, 17.9, 18.8, 8.55, 
            12.5, 17.6, 18.5, 2.98),
        cutip_gsf = cms.vdouble(0.0551, 0.0765, 0.143, 0.0874, 0.594, 
            0.37, 0.0913, 1.15, 0.231),
        cutip_gsfl = cms.vdouble(0.0186, 0.0759, 0.138, 0.0473, 0.62, 
            0.304, 0.109, 0.775, 0.0479),
        cuthoe = cms.vdouble(0.247, 0.137, 0.147, 0.371, 0.0588, 
            0.147, 0.52, 0.452, 0.404),
        cutiso_sumoetl = cms.vdouble(11.3, 9.05, 9.07, 9.94, 5.25, 
            6.15, 10.7, 10.8, 4.4),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 4.5, 3.5, 3.5),
        cuthoel = cms.vdouble(0.236, 0.126, 0.147, 0.375, 0.0392, 
            0.145, 0.365, 0.383, 0.384),
        cutdphiin = cms.vdouble(0.0897, 0.262, 0.353, 0.116, 0.357, 
            0.319, 0.342, 0.404, 0.336),
        cutseel = cms.vdouble(0.0164, 0.0118, 0.015, 0.0523, 0.0326, 
            0.0456, 0.0185, 0.0589, 0.0544),
        cutiso_sumoet = cms.vdouble(34.5, 12.7, 12.1, 19.9, 6.35, 
            8.85, 14.0, 10.5, 9.74),
        cutsee = cms.vdouble(0.0176, 0.0125, 0.0181, 0.0415, 0.0364, 
            0.0418, 0.0146, 0.0678, 0.133),
        cuteseedopcor = cms.vdouble(0.63, 0.82, 0.401, 0.718, 0.4, 
            0.458, 0.15, 0.664, 0.373),
        cutdphiinl = cms.vdouble(0.0747, 0.25, 0.356, 0.0956, 0.347, 
            0.326, 0.333, 0.647, 0.289),
        cutdetainl = cms.vdouble(0.0124, 0.00503, 0.0257, 0.0228, 0.0118, 
            0.0178, 0.0188, 0.14, 0.024)
    ),
    src = cms.InputTag("gsfElectrons"),
    robusttightEleIDCuts = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    classbasedtightEleIDCuts = cms.PSet(
        cutdetain = cms.vdouble(0.0116, 0.00449, 0.00938, 0.0184, 0.00678, 
            0.0109, 0.0252, 0.0268, 0.0139),
        cutiso_sum = cms.vdouble(15.5, 12.2, 12.2, 11.7, 7.16, 
            9.71, 8.66, 11.9, 2.98),
        cutip_gsf = cms.vdouble(0.0131, 0.0586, 0.0839, 0.0366, 0.452, 
            0.204, 0.0913, 0.0802, 0.0731),
        cutip_gsfl = cms.vdouble(0.0119, 0.0527, 0.0471, 0.0212, 0.233, 
            0.267, 0.109, 0.122, 0.0479),
        cuthoe = cms.vdouble(0.215, 0.0608, 0.147, 0.369, 0.0349, 
            0.102, 0.52, 0.422, 0.404),
        cutiso_sumoetl = cms.vdouble(6.21, 6.81, 5.3, 5.39, 2.73, 
            4.73, 4.84, 3.46, 3.73),
        cutfmishits = cms.vdouble(1.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 1.5, 2.5, 0.5),
        cuthoel = cms.vdouble(0.228, 0.0836, 0.143, 0.37, 0.0392, 
            0.0979, 0.3, 0.381, 0.339),
        cutdphiin = cms.vdouble(0.0897, 0.0993, 0.295, 0.0979, 0.151, 
            0.252, 0.341, 0.308, 0.328),
        cutseel = cms.vdouble(0.0132, 0.0117, 0.0112, 0.0387, 0.0281, 
            0.0287, 0.00987, 0.0296, 0.0544),
        cutiso_sumoet = cms.vdouble(11.9, 7.81, 6.28, 8.92, 4.65, 
            5.49, 9.36, 8.84, 5.94),
        cutsee = cms.vdouble(0.0145, 0.0116, 0.012, 0.039, 0.0297, 
            0.0311, 0.00987, 0.0347, 0.0917),
        cuteseedopcor = cms.vdouble(0.637, 0.943, 0.742, 0.748, 0.763, 
            0.631, 0.214, 0.873, 0.473),
        cutdphiinl = cms.vdouble(0.061, 0.14, 0.286, 0.0921, 0.197, 
            0.24, 0.333, 0.303, 0.258),
        cutdetainl = cms.vdouble(0.00816, 0.00401, 0.0081, 0.019, 0.00588, 
            0.00893, 0.0171, 0.0434, 0.0143)
    ),
    algorithm = cms.string('eIDCB'),
    robusthighenergyEleIDCuts = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCuts = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV02 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV03 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV00 = cms.PSet(
        barrel = cms.vdouble(0.115, 0.014, 0.09, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.15, 0.0275, 0.092, 0.0105, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV01 = cms.PSet(
        barrel = cms.vdouble(0.075, 0.0132, 0.058, 0.0077, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.083, 0.027, 0.042, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    additionalCategories = cms.bool(True),
    etBinning = cms.bool(True)
)


process.eidRobustTight = cms.EDProducer("EleIdCutBasedExtProducer",
    electronQuality = cms.string('tight'),
    classbasedtightEleIDCutsV02 = cms.PSet(
        cutisohcal = cms.vdouble(10.9, 7.01, 8.75, 3.51, 7.75, 
            1.62, 11.6, 9.9, 4.97, 5.33, 
            3.18, 2.32, 0.164, 5.46, 12.0, 
            0.00604, 4.1, 0.000628),
        cutmishits = cms.vdouble(5.5, 1.5, 0.5, 1.5, 2.5, 
            0.5, 3.5, 5.5, 0.5, 0.5, 
            0.5, 0.5, 0.5, 1.5, 0.5, 
            0.5, 0.5, 0.5),
        cuthoe = cms.vdouble(0.0871, 0.0289, 0.0783, 0.0946, 0.0245, 
            0.0363, 0.0671, 0.048, 0.0614, 0.0924, 
            0.0158, 0.049, 0.0382, 0.0915, 0.0451, 
            0.0452, 0.00196, 0.0043),
        cutdeta = cms.vdouble(0.00915, 0.00302, 0.0061, 0.0135, 0.00565, 
            0.00793, 0.0102, 0.00266, 0.0106, 0.00903, 
            0.00766, 0.00723, 0.0116, 0.00203, 0.00659, 
            0.0148, 0.00555, 0.0128),
        cuteopin = cms.vdouble(0.878, 0.859, 0.874, 0.944, 0.737, 
            0.773, 0.86, 0.967, 0.917, 0.812, 
            0.915, 1.01, 0.847, 0.953, 0.979, 
            0.841, 0.771, 1.09),
        cutip = cms.vdouble(0.0239, 0.027, 0.0768, 0.0231, 0.178, 
            0.0957, 0.0102, 0.0168, 0.043, 0.0166, 
            0.0594, 0.0308, 2.1, 0.00527, 3.17, 
            4.91, 0.769, 5.9),
        cutisotk = cms.vdouble(6.53, 4.6, 6.0, 8.63, 3.11, 
            7.77, 5.42, 4.81, 4.06, 6.47, 
            2.8, 3.45, 5.29, 5.18, 15.4, 
            5.38, 4.47, 0.0347),
        cutsee = cms.vdouble(0.0131, 0.0106, 0.0115, 0.0306, 0.028, 
            0.0293, 0.0131, 0.0106, 0.0115, 0.0317, 
            0.029, 0.0289, 0.0142, 0.0106, 0.0103, 
            0.035, 0.0296, 0.0333),
        cutdphi = cms.vdouble(0.0369, 0.0307, 0.117, 0.0475, 0.0216, 
            0.117, 0.0372, 0.0246, 0.0426, 0.0612, 
            0.0142, 0.039, 0.0737, 0.0566, 0.0359, 
            0.0187, 0.012, 0.0358),
        cutisoecal = cms.vdouble(20.0, 27.2, 4.48, 13.5, 4.56, 
            3.19, 12.2, 13.1, 7.42, 7.67, 
            4.12, 4.85, 10.1, 12.4, 11.1, 
            11.0, 10.6, 13.4)
    ),
    classbasedtightEleIDCutsV03 = cms.PSet(
        cutdetain = cms.vdouble(0.00811, 0.00341, 0.00633, 0.0103, 0.00667, 
            0.01, 0.0106, 0.0145, 0.0163, 0.0076, 
            0.00259, 0.00511, 0.00941, 0.0043, 0.00857, 
            0.012, 0.0169, 0.00172, 0.00861, 0.00362, 
            0.00601, 0.00925, 0.00489, 0.00832, 0.0119, 
            0.0169, 0.000996),
        cutiso_sum = cms.vdouble(11.8, 8.31, 6.26, 6.18, 3.28, 
            4.38, 4.17, 5.4, 1.57, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0213, 0.0422, 0.0632, 0.0361, 0.073, 
            0.126, 0.171, 0.119, 0.0372, 0.0131, 
            0.0146, 0.0564, 0.0152, 0.0222, 0.0268, 
            0.0314, 0.0884, 0.00374, 0.00852, 0.00761, 
            0.0143, 0.0106, 0.0127, 0.0119, 0.0123, 
            0.0235, 0.00363),
        cuthoe = cms.vdouble(0.0783, 0.0387, 0.105, 0.118, 0.0227, 
            0.062, 0.13, 2.47, 0.38, 0.0888, 
            0.0503, 0.0955, 0.0741, 0.015, 0.03, 
            0.589, 1.13, 0.612, 0.0494, 0.0461, 
            0.0292, 0.0369, 0.0113, 0.0145, 0.124, 
            2.05, 0.61),
        cutfmishits = cms.vdouble(2.5, 1.5, 1.5, 1.5, 1.5, 
            0.5, 2.5, 0.5, 0.5, 2.5, 
            1.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5, -0.5, 2.5, 1.5, 
            0.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(13.7, 11.6, 7.14, 9.98, 3.52, 
            4.87, 6.24, 7.96, 2.53, 11.2, 
            11.9, 7.88, 8.16, 5.58, 5.03, 
            11.4, 8.15, 5.79, 10.4, 11.1, 
            10.4, 7.47, 5.08, 5.9, 11.8, 
            14.1, 11.7),
        cutdcotdist = cms.vdouble(0.0393, 0.0256, 0.00691, 0.0394, 0.0386, 
            0.039, 0.0325, 0.0384, 0.0382, 0.0245, 
            0.000281, 5.46e-05, 0.0342, 0.0232, 0.00107, 
            0.0178, 0.0193, 0.000758, 0.000108, 0.0248, 
            0.000458, 0.0129, 0.00119, 0.0182, 4.53e-05, 
            0.0189, 0.000928),
        cutsee = cms.vdouble(0.0143, 0.0105, 0.0123, 0.0324, 0.0307, 
            0.0301, 0.0109, 0.027, 0.0292, 0.0133, 
            0.0104, 0.0116, 0.0332, 0.0296, 0.031, 
            0.00981, 0.0307, 0.072, 0.0149, 0.0105, 
            0.011, 0.0342, 0.0307, 0.0303, 0.00954, 
            0.0265, 0.0101),
        cuteseedopcor = cms.vdouble(0.784, 0.366, 0.57, 0.911, 0.298, 
            0.645, 0.51, 0.497, 0.932, 0.835, 
            0.968, 0.969, 0.923, 0.898, 0.98, 
            0.63, 0.971, 1.0, 0.515, 0.963, 
            0.986, 0.823, 0.879, 1.01, 0.931, 
            0.937, 1.05),
        cutdphiin = cms.vdouble(0.0404, 0.0499, 0.263, 0.042, 0.0484, 
            0.241, 0.242, 0.231, 0.286, 0.0552, 
            0.0338, 0.154, 0.0623, 0.0183, 0.0392, 
            0.0547, 0.0588, 0.00654, 0.042, 0.0217, 
            0.0885, 0.0445, 0.0141, 0.0234, 0.065, 
            0.0258, 0.0346),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 13.7, 13.2, 
            13.6, 14.2, 14.1, 13.9, 12.9, 
            14.9, 17.7)
    ),
    classbasedtightEleIDCutsV00 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.032, 0.016, 0.0525, 0.09, 0.025, 
            0.035, 0.065, 0.092),
        hOverE = cms.vdouble(0.05, 0.042, 0.045, 0.0, 0.055, 
            0.037, 0.05, 0.0),
        sigmaEtaEta = cms.vdouble(0.0125, 0.011, 0.01, 0.0, 0.0265, 
            0.0252, 0.026, 0.0),
        deltaEtaIn = cms.vdouble(0.0055, 0.003, 0.0065, 0.0, 0.006, 
            0.0055, 0.0075, 0.0),
        eSeedOverPin = cms.vdouble(0.24, 0.94, 0.11, 0.0, 0.32, 
            0.83, 0.0, 0.0)
    ),
    classbasedtightEleIDCutsV01 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.0225, 0.0114, 0.0234, 0.039, 0.0215, 
            0.0095, 0.0148, 0.0167),
        hOverE = cms.vdouble(0.056, 0.0221, 0.037, 0.0, 0.0268, 
            0.0102, 0.0104, 0.0),
        sigmaEtaEta = cms.vdouble(0.0095, 0.0094, 0.0094, 0.0, 0.026, 
            0.0257, 0.0246, 0.0),
        deltaEtaIn = cms.vdouble(0.0043, 0.00282, 0.0036, 0.0, 0.0066, 
            0.0049, 0.0041, 0.0),
        eSeedOverPin = cms.vdouble(0.32, 0.94, 0.221, 0.0, 0.74, 
            0.89, 0.66, 0.0)
    ),
    classbasedtightEleIDCutsV06 = cms.PSet(
        cutdetain = cms.vdouble(0.0116, 0.00449, 0.00938, 0.0184, 0.00678, 
            0.0109, 0.0252, 0.0268, 0.0139),
        cutiso_sum = cms.vdouble(15.5, 12.2, 12.2, 11.7, 7.16, 
            9.71, 8.66, 11.9, 2.98),
        cutip_gsf = cms.vdouble(0.0131, 0.0586, 0.0839, 0.0366, 0.452, 
            0.204, 0.0913, 0.0802, 0.0731),
        cutip_gsfl = cms.vdouble(0.0119, 0.0527, 0.0471, 0.0212, 0.233, 
            0.267, 0.109, 0.122, 0.0479),
        cuthoe = cms.vdouble(0.215, 0.0608, 0.147, 0.369, 0.0349, 
            0.102, 0.52, 0.422, 0.404),
        cutiso_sumoetl = cms.vdouble(6.21, 6.81, 5.3, 5.39, 2.73, 
            4.73, 4.84, 3.46, 3.73),
        cutfmishits = cms.vdouble(1.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 1.5, 2.5, 0.5),
        cuthoel = cms.vdouble(0.228, 0.0836, 0.143, 0.37, 0.0392, 
            0.0979, 0.3, 0.381, 0.339),
        cutdphiin = cms.vdouble(0.0897, 0.0993, 0.295, 0.0979, 0.151, 
            0.252, 0.341, 0.308, 0.328),
        cutseel = cms.vdouble(0.0132, 0.0117, 0.0112, 0.0387, 0.0281, 
            0.0287, 0.00987, 0.0296, 0.0544),
        cutiso_sumoet = cms.vdouble(11.9, 7.81, 6.28, 8.92, 4.65, 
            5.49, 9.36, 8.84, 5.94),
        cutsee = cms.vdouble(0.0145, 0.0116, 0.012, 0.039, 0.0297, 
            0.0311, 0.00987, 0.0347, 0.0917),
        cuteseedopcor = cms.vdouble(0.637, 0.943, 0.742, 0.748, 0.763, 
            0.631, 0.214, 0.873, 0.473),
        cutdphiinl = cms.vdouble(0.061, 0.14, 0.286, 0.0921, 0.197, 
            0.24, 0.333, 0.303, 0.258),
        cutdetainl = cms.vdouble(0.00816, 0.00401, 0.0081, 0.019, 0.00588, 
            0.00893, 0.0171, 0.0434, 0.0143)
    ),
    classbasedtightEleIDCutsV04 = cms.PSet(
        cutdetain = cms.vdouble(0.00811, 0.00341, 0.00633, 0.0103, 0.00667, 
            0.01, 0.0106, 0.0145, 0.0163, 0.0076, 
            0.00259, 0.00511, 0.00941, 0.0043, 0.00857, 
            0.012, 0.0169, 0.00172, 0.00861, 0.00362, 
            0.00601, 0.00925, 0.00489, 0.00832, 0.0119, 
            0.0169, 0.000996),
        cutiso_sum = cms.vdouble(11.8, 8.31, 6.26, 6.18, 3.28, 
            4.38, 4.17, 5.4, 1.57, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0213, 0.0422, 0.0632, 0.0361, 0.073, 
            0.126, 0.171, 0.119, 0.0372, 0.0131, 
            0.0146, 0.0564, 0.0152, 0.0222, 0.0268, 
            0.0314, 0.0884, 0.00374, 0.00852, 0.00761, 
            0.0143, 0.0106, 0.0127, 0.0119, 0.0123, 
            0.0235, 0.00363),
        cuthoe = cms.vdouble(0.0783, 0.0387, 0.105, 0.118, 0.0227, 
            0.062, 0.13, 2.47, 0.38, 0.0888, 
            0.0503, 0.0955, 0.0741, 0.015, 0.03, 
            0.589, 1.13, 0.612, 0.0494, 0.0461, 
            0.0292, 0.0369, 0.0113, 0.0145, 0.124, 
            2.05, 0.61),
        cutfmishits = cms.vdouble(2.5, 1.5, 1.5, 1.5, 1.5, 
            0.5, 2.5, 0.5, 0.5, 2.5, 
            1.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5, -0.5, 2.5, 1.5, 
            0.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(13.7, 11.6, 7.14, 9.98, 3.52, 
            4.87, 6.24, 7.96, 2.53, 11.2, 
            11.9, 7.88, 8.16, 5.58, 5.03, 
            11.4, 8.15, 5.79, 10.4, 11.1, 
            10.4, 7.47, 5.08, 5.9, 11.8, 
            14.1, 11.7),
        cutdcotdist = cms.vdouble(0.0393, 0.0256, 0.00691, 0.0394, 0.0386, 
            0.039, 0.0325, 0.0384, 0.0382, 0.0245, 
            0.000281, 5.46e-05, 0.0342, 0.0232, 0.00107, 
            0.0178, 0.0193, 0.000758, 0.000108, 0.0248, 
            0.000458, 0.0129, 0.00119, 0.0182, 4.53e-05, 
            0.0189, 0.000928),
        cutsee = cms.vdouble(0.0143, 0.0105, 0.0123, 0.0324, 0.0307, 
            0.0301, 0.0109, 0.027, 0.0292, 0.0133, 
            0.0104, 0.0116, 0.0332, 0.0296, 0.031, 
            0.00981, 0.0307, 0.072, 0.0149, 0.0105, 
            0.011, 0.0342, 0.0307, 0.0303, 0.00954, 
            0.0265, 0.0101),
        cuteseedopcor = cms.vdouble(0.784, 0.366, 0.57, 0.911, 0.298, 
            0.645, 0.51, 0.497, 0.932, 0.835, 
            0.968, 0.969, 0.923, 0.898, 0.98, 
            0.63, 0.971, 1.0, 0.515, 0.963, 
            0.986, 0.823, 0.879, 1.01, 0.931, 
            0.937, 1.05),
        cutdphiin = cms.vdouble(0.0404, 0.0499, 0.263, 0.042, 0.0484, 
            0.241, 0.242, 0.231, 0.286, 0.0552, 
            0.0338, 0.154, 0.0623, 0.0183, 0.0392, 
            0.0547, 0.0588, 0.00654, 0.042, 0.0217, 
            0.0885, 0.0445, 0.0141, 0.0234, 0.065, 
            0.0258, 0.0346),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 13.7, 13.2, 
            13.6, 14.2, 14.1, 13.9, 12.9, 
            14.9, 17.7)
    ),
    electronIDType = cms.string('robust'),
    robusttightEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    electronVersion = cms.string(''),
    robusttightEleIDCutsV00 = cms.PSet(
        barrel = cms.vdouble(0.015, 0.0092, 0.02, 0.0025, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.018, 0.025, 0.02, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusttightEleIDCutsV01 = cms.PSet(
        barrel = cms.vdouble(0.01, 0.0099, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.01, 0.028, 0.02, 0.0066, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    reducedBarrelRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    robusttightEleIDCutsV03 = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    reducedEndcapRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    verticesCollection = cms.InputTag("offlinePrimaryVerticesWithBS"),
    classbasedlooseEleIDCuts = cms.PSet(
        cutdetain = cms.vdouble(0.0137, 0.00678, 0.0241, 0.0187, 0.0161, 
            0.0224, 0.0252, 0.0308, 0.0273),
        cutiso_sum = cms.vdouble(33.0, 17.0, 17.9, 18.8, 8.55, 
            12.5, 17.6, 18.5, 2.98),
        cutip_gsf = cms.vdouble(0.0551, 0.0765, 0.143, 0.0874, 0.594, 
            0.37, 0.0913, 1.15, 0.231),
        cutip_gsfl = cms.vdouble(0.0186, 0.0759, 0.138, 0.0473, 0.62, 
            0.304, 0.109, 0.775, 0.0479),
        cuthoe = cms.vdouble(0.247, 0.137, 0.147, 0.371, 0.0588, 
            0.147, 0.52, 0.452, 0.404),
        cutiso_sumoetl = cms.vdouble(11.3, 9.05, 9.07, 9.94, 5.25, 
            6.15, 10.7, 10.8, 4.4),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 4.5, 3.5, 3.5),
        cuthoel = cms.vdouble(0.236, 0.126, 0.147, 0.375, 0.0392, 
            0.145, 0.365, 0.383, 0.384),
        cutdphiin = cms.vdouble(0.0897, 0.262, 0.353, 0.116, 0.357, 
            0.319, 0.342, 0.404, 0.336),
        cutseel = cms.vdouble(0.0164, 0.0118, 0.015, 0.0523, 0.0326, 
            0.0456, 0.0185, 0.0589, 0.0544),
        cutiso_sumoet = cms.vdouble(34.5, 12.7, 12.1, 19.9, 6.35, 
            8.85, 14.0, 10.5, 9.74),
        cutsee = cms.vdouble(0.0176, 0.0125, 0.0181, 0.0415, 0.0364, 
            0.0418, 0.0146, 0.0678, 0.133),
        cuteseedopcor = cms.vdouble(0.63, 0.82, 0.401, 0.718, 0.4, 
            0.458, 0.15, 0.664, 0.373),
        cutdphiinl = cms.vdouble(0.0747, 0.25, 0.356, 0.0956, 0.347, 
            0.326, 0.333, 0.647, 0.289),
        cutdetainl = cms.vdouble(0.0124, 0.00503, 0.0257, 0.0228, 0.0118, 
            0.0178, 0.0188, 0.14, 0.024)
    ),
    robusttightEleIDCutsV02 = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV01 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.0275, 0.09, 0.007, -1, 
            -1, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV00 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.011, 0.09, 0.005, -1, 
            -1, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.1, 0.0275, 0.09, 0.007, -1, 
            -1, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV03 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV02 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    classbasedlooseEleIDCutsV00 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.05, 0.025, 0.053, 0.09, 0.07, 
            0.03, 0.092, 0.092),
        hOverE = cms.vdouble(0.115, 0.1, 0.055, 0.0, 0.145, 
            0.12, 0.15, 0.0),
        sigmaEtaEta = cms.vdouble(0.014, 0.012, 0.0115, 0.0, 0.0275, 
            0.0265, 0.0265, 0.0),
        deltaEtaIn = cms.vdouble(0.009, 0.0045, 0.0085, 0.0, 0.0105, 
            0.0068, 0.01, 0.0),
        eSeedOverPin = cms.vdouble(0.11, 0.91, 0.11, 0.0, 0.0, 
            0.85, 0.0, 0.0)
    ),
    classbasedlooseEleIDCutsV01 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.053, 0.0189, 0.059, 0.099, 0.0278, 
            0.0157, 0.042, 0.08),
        hOverE = cms.vdouble(0.076, 0.033, 0.07, 0.0, 0.083, 
            0.0148, 0.033, 0.0),
        sigmaEtaEta = cms.vdouble(0.0101, 0.0095, 0.0097, 0.0, 0.0271, 
            0.0267, 0.0259, 0.0),
        deltaEtaIn = cms.vdouble(0.0078, 0.00259, 0.0062, 0.0, 0.0078, 
            0.0061, 0.0061, 0.0),
        eSeedOverPin = cms.vdouble(0.3, 0.92, 0.211, 0.0, 0.42, 
            0.88, 0.68, 0.0)
    ),
    classbasedlooseEleIDCutsV02 = cms.PSet(
        cutisohcal = cms.vdouble(13.5, 9.93, 7.56, 14.8, 8.1, 
            10.8, 42.7, 20.1, 9.11, 10.4, 
            6.89, 5.59, 8.53, 9.59, 24.2, 
            2.78, 8.67, 0.288),
        cutmishits = cms.vdouble(5.5, 1.5, 5.5, 2.5, 2.5, 
            2.5, 3.5, 5.5, 0.5, 1.5, 
            2.5, 0.5, 1.5, 1.5, 0.5, 
            0.5, 0.5, 0.5),
        cuthoe = cms.vdouble(0.0887, 0.0934, 0.0949, 0.0986, 0.0431, 
            0.0878, 0.097, 0.0509, 0.098, 0.0991, 
            0.0321, 0.0928, 0.0663, 0.0717, 0.0966, 
            0.0758, 0.0149, 0.0131),
        cutdeta = cms.vdouble(0.00958, 0.00406, 0.0122, 0.0137, 0.00837, 
            0.0127, 0.011, 0.00336, 0.00977, 0.015, 
            0.00675, 0.0109, 0.014, 0.00508, 0.0109, 
            0.0146, 0.00506, 0.0127),
        cuteopin = cms.vdouble(0.878, 0.802, 0.814, 0.942, 0.735, 
            0.774, 0.829, 0.909, 0.829, 0.813, 
            0.86, 0.897, 0.817, 0.831, 0.818, 
            0.861, 0.787, 0.789),
        cutip = cms.vdouble(0.0246, 0.076, 0.0966, 0.0885, 0.441, 
            0.205, 0.0292, 0.0293, 0.0619, 0.0251, 
            0.159, 0.0815, 7.29, 0.0106, 5.76, 
            6.89, 1.27, 5.89),
        cutisotk = cms.vdouble(24.3, 8.45, 14.4, 27.8, 6.02, 
            10.5, 14.1, 10.2, 14.5, 19.1, 
            6.1, 14.1, 8.59, 8.33, 8.3, 
            8.93, 8.6, 16.0),
        cutsee = cms.vdouble(0.0172, 0.0115, 0.0143, 0.0344, 0.0295, 
            0.0304, 0.0145, 0.0108, 0.0128, 0.0347, 
            0.0307, 0.0316, 0.018, 0.011, 0.0132, 
            0.0349, 0.031, 0.0327),
        cutdphi = cms.vdouble(0.0372, 0.114, 0.118, 0.0488, 0.117, 
            0.119, 0.0606, 0.0548, 0.117, 0.07, 
            0.0355, 0.117, 0.088, 0.045, 0.118, 
            0.0919, 0.0236, 0.0515),
        cutisoecal = cms.vdouble(33.4, 28.1, 7.32, 27.4, 7.33, 
            21.7, 93.8, 102.0, 12.1, 26.0, 
            8.91, 10.0, 16.1, 31.3, 16.9, 
            15.4, 13.3, 37.7)
    ),
    classbasedlooseEleIDCutsV03 = cms.PSet(
        cutdetain = cms.vdouble(0.00989, 0.00484, 0.0146, 0.0146, 0.00902, 
            0.0172, 0.0137, 0.0477, 0.0275, 0.00967, 
            0.00377, 0.00924, 0.013, 0.00666, 0.0123, 
            0.0125, 0.0228, 0.0112, 0.0106, 0.0038, 
            0.00897, 0.0139, 0.00667, 0.0122, 0.0122, 
            0.0193, 0.00239),
        cutiso_sum = cms.vdouble(31.5, 10.3, 8.8, 11.0, 6.13, 
            6.94, 7.52, 9.0, 3.5, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0431, 0.0767, 0.139, 0.101, 0.149, 
            0.154, 0.932, 0.15, 0.124, 0.0238, 
            0.0467, 0.0759, 0.0369, 0.147, 0.0986, 
            0.0626, 0.195, 0.116, 0.0122, 0.0125, 
            0.0693, 0.0162, 0.089, 0.0673, 0.0467, 
            0.0651, 0.0221),
        cuthoe = cms.vdouble(0.166, 0.0771, 0.144, 0.37, 0.0497, 
            0.139, 0.401, 2.68, 0.516, 0.234, 
            0.0556, 0.144, 0.368, 0.031, 0.12, 
            0.602, 2.01, 1.05, 0.104, 0.063, 
            0.0565, 0.38, 0.0192, 0.0294, 0.537, 
            4.65, 1.87),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 2.5, 2.5, 1.5, 2.5, 
            1.5, 1.5, 1.5, 1.5, 0.5, 
            2.5, 2.5, 0.5, 2.5, 1.5, 
            0.5, 1.5, 1.5, 0.5, 2.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(28.9, 15.3, 12.0, 18.3, 7.17, 
            9.42, 11.0, 9.81, 3.94, 22.7, 
            15.9, 12.3, 17.0, 7.58, 8.89, 
            15.2, 12.7, 6.17, 20.8, 21.2, 
            17.2, 15.5, 9.37, 10.6, 19.8, 
            22.1, 15.6),
        cutdcotdist = cms.vdouble(0.0393, 0.0392, 0.0397, 0.0394, 0.0393, 
            0.039, 0.0378, 0.0388, 0.0382, 0.0385, 
            0.0167, 0.00325, 0.0394, 0.0387, 0.0388, 
            0.0227, 0.0258, 0.0127, 0.0298, 0.03, 
            0.00946, 0.039, 0.0231, 0.0278, 0.00162, 
            0.0367, 0.0199),
        cutsee = cms.vdouble(0.0175, 0.0127, 0.0177, 0.0373, 0.0314, 
            0.0329, 0.0157, 0.0409, 0.14, 0.0169, 
            0.0106, 0.0142, 0.0363, 0.0322, 0.0354, 
            0.0117, 0.0372, 28.2, 0.0171, 0.0113, 
            0.014, 0.0403, 0.0323, 0.0411, 0.0104, 
            0.0436, 0.0114),
        cuteseedopcor = cms.vdouble(0.78, 0.302, 0.483, 0.904, 0.168, 
            0.645, 0.108, 0.284, 0.324, 0.591, 
            0.286, 0.488, 0.813, 0.791, 0.672, 
            0.398, 0.834, 0.878, 0.515, 0.937, 
            0.806, 0.816, 0.85, 0.507, 0.367, 
            0.83, 0.648),
        cutdphiin = cms.vdouble(0.041, 0.275, 0.365, 0.047, 0.273, 
            0.296, 0.329, 0.465, 0.627, 0.0581, 
            0.0954, 0.327, 0.0702, 0.0582, 0.279, 
            0.117, 0.318, 0.246, 0.0821, 0.052, 
            0.292, 0.116, 0.0435, 0.312, 0.118, 
            0.296, 0.0459),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 12.0, 12.0, 
            12.0, 12.0, 12.0, 12.0, 12.0, 
            12.0, 12.5)
    ),
    classbasedlooseEleIDCutsV04 = cms.PSet(
        cutdetain = cms.vdouble(0.00989, 0.00484, 0.0146, 0.0146, 0.00902, 
            0.0172, 0.0137, 0.0477, 0.0275, 0.00967, 
            0.00377, 0.00924, 0.013, 0.00666, 0.0123, 
            0.0125, 0.0228, 0.0112, 0.0106, 0.0038, 
            0.00897, 0.0139, 0.00667, 0.0122, 0.0122, 
            0.0193, 0.00239),
        cutiso_sum = cms.vdouble(31.5, 10.3, 8.8, 11.0, 6.13, 
            6.94, 7.52, 9.0, 3.5, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0431, 0.0767, 0.139, 0.101, 0.149, 
            0.154, 0.932, 0.15, 0.124, 0.0238, 
            0.0467, 0.0759, 0.0369, 0.147, 0.0986, 
            0.0626, 0.195, 0.116, 0.0122, 0.0125, 
            0.0693, 0.0162, 0.089, 0.0673, 0.0467, 
            0.0651, 0.0221),
        cuthoe = cms.vdouble(0.166, 0.0771, 0.144, 0.37, 0.0497, 
            0.139, 0.401, 2.68, 0.516, 0.234, 
            0.0556, 0.144, 0.368, 0.031, 0.12, 
            0.602, 2.01, 1.05, 0.104, 0.063, 
            0.0565, 0.38, 0.0192, 0.0294, 0.537, 
            4.65, 1.87),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 2.5, 2.5, 1.5, 2.5, 
            1.5, 1.5, 1.5, 1.5, 0.5, 
            2.5, 2.5, 0.5, 2.5, 1.5, 
            0.5, 1.5, 1.5, 0.5, 2.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(28.9, 15.3, 12.0, 18.3, 7.17, 
            9.42, 11.0, 9.81, 3.94, 22.7, 
            15.9, 12.3, 17.0, 7.58, 8.89, 
            15.2, 12.7, 6.17, 20.8, 21.2, 
            17.2, 15.5, 9.37, 10.6, 19.8, 
            22.1, 15.6),
        cutdcotdist = cms.vdouble(0.0393, 0.0392, 0.0397, 0.0394, 0.0393, 
            0.039, 0.0378, 0.0388, 0.0382, 0.0385, 
            0.0167, 0.00325, 0.0394, 0.0387, 0.0388, 
            0.0227, 0.0258, 0.0127, 0.0298, 0.03, 
            0.00946, 0.039, 0.0231, 0.0278, 0.00162, 
            0.0367, 0.0199),
        cutsee = cms.vdouble(0.0175, 0.0127, 0.0177, 0.0373, 0.0314, 
            0.0329, 0.0157, 0.0409, 0.14, 0.0169, 
            0.0106, 0.0142, 0.0363, 0.0322, 0.0354, 
            0.0117, 0.0372, 28.2, 0.0171, 0.0113, 
            0.014, 0.0403, 0.0323, 0.0411, 0.0104, 
            0.0436, 0.0114),
        cuteseedopcor = cms.vdouble(0.78, 0.302, 0.483, 0.904, 0.168, 
            0.645, 0.108, 0.284, 0.324, 0.591, 
            0.286, 0.488, 0.813, 0.791, 0.672, 
            0.398, 0.834, 0.878, 0.515, 0.937, 
            0.806, 0.816, 0.85, 0.507, 0.367, 
            0.83, 0.648),
        cutdphiin = cms.vdouble(0.041, 0.275, 0.365, 0.047, 0.273, 
            0.296, 0.329, 0.465, 0.627, 0.0581, 
            0.0954, 0.327, 0.0702, 0.0582, 0.279, 
            0.117, 0.318, 0.246, 0.0821, 0.052, 
            0.292, 0.116, 0.0435, 0.312, 0.118, 
            0.296, 0.0459),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 12.0, 12.0, 
            12.0, 12.0, 12.0, 12.0, 12.0, 
            12.0, 12.5)
    ),
    classbasedlooseEleIDCutsV06 = cms.PSet(
        cutdetain = cms.vdouble(0.0137, 0.00678, 0.0241, 0.0187, 0.0161, 
            0.0224, 0.0252, 0.0308, 0.0273),
        cutiso_sum = cms.vdouble(33.0, 17.0, 17.9, 18.8, 8.55, 
            12.5, 17.6, 18.5, 2.98),
        cutip_gsf = cms.vdouble(0.0551, 0.0765, 0.143, 0.0874, 0.594, 
            0.37, 0.0913, 1.15, 0.231),
        cutip_gsfl = cms.vdouble(0.0186, 0.0759, 0.138, 0.0473, 0.62, 
            0.304, 0.109, 0.775, 0.0479),
        cuthoe = cms.vdouble(0.247, 0.137, 0.147, 0.371, 0.0588, 
            0.147, 0.52, 0.452, 0.404),
        cutiso_sumoetl = cms.vdouble(11.3, 9.05, 9.07, 9.94, 5.25, 
            6.15, 10.7, 10.8, 4.4),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 4.5, 3.5, 3.5),
        cuthoel = cms.vdouble(0.236, 0.126, 0.147, 0.375, 0.0392, 
            0.145, 0.365, 0.383, 0.384),
        cutdphiin = cms.vdouble(0.0897, 0.262, 0.353, 0.116, 0.357, 
            0.319, 0.342, 0.404, 0.336),
        cutseel = cms.vdouble(0.0164, 0.0118, 0.015, 0.0523, 0.0326, 
            0.0456, 0.0185, 0.0589, 0.0544),
        cutiso_sumoet = cms.vdouble(34.5, 12.7, 12.1, 19.9, 6.35, 
            8.85, 14.0, 10.5, 9.74),
        cutsee = cms.vdouble(0.0176, 0.0125, 0.0181, 0.0415, 0.0364, 
            0.0418, 0.0146, 0.0678, 0.133),
        cuteseedopcor = cms.vdouble(0.63, 0.82, 0.401, 0.718, 0.4, 
            0.458, 0.15, 0.664, 0.373),
        cutdphiinl = cms.vdouble(0.0747, 0.25, 0.356, 0.0956, 0.347, 
            0.326, 0.333, 0.647, 0.289),
        cutdetainl = cms.vdouble(0.0124, 0.00503, 0.0257, 0.0228, 0.0118, 
            0.0178, 0.0188, 0.14, 0.024)
    ),
    src = cms.InputTag("gsfElectrons"),
    robusttightEleIDCuts = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    classbasedtightEleIDCuts = cms.PSet(
        cutdetain = cms.vdouble(0.0116, 0.00449, 0.00938, 0.0184, 0.00678, 
            0.0109, 0.0252, 0.0268, 0.0139),
        cutiso_sum = cms.vdouble(15.5, 12.2, 12.2, 11.7, 7.16, 
            9.71, 8.66, 11.9, 2.98),
        cutip_gsf = cms.vdouble(0.0131, 0.0586, 0.0839, 0.0366, 0.452, 
            0.204, 0.0913, 0.0802, 0.0731),
        cutip_gsfl = cms.vdouble(0.0119, 0.0527, 0.0471, 0.0212, 0.233, 
            0.267, 0.109, 0.122, 0.0479),
        cuthoe = cms.vdouble(0.215, 0.0608, 0.147, 0.369, 0.0349, 
            0.102, 0.52, 0.422, 0.404),
        cutiso_sumoetl = cms.vdouble(6.21, 6.81, 5.3, 5.39, 2.73, 
            4.73, 4.84, 3.46, 3.73),
        cutfmishits = cms.vdouble(1.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 1.5, 2.5, 0.5),
        cuthoel = cms.vdouble(0.228, 0.0836, 0.143, 0.37, 0.0392, 
            0.0979, 0.3, 0.381, 0.339),
        cutdphiin = cms.vdouble(0.0897, 0.0993, 0.295, 0.0979, 0.151, 
            0.252, 0.341, 0.308, 0.328),
        cutseel = cms.vdouble(0.0132, 0.0117, 0.0112, 0.0387, 0.0281, 
            0.0287, 0.00987, 0.0296, 0.0544),
        cutiso_sumoet = cms.vdouble(11.9, 7.81, 6.28, 8.92, 4.65, 
            5.49, 9.36, 8.84, 5.94),
        cutsee = cms.vdouble(0.0145, 0.0116, 0.012, 0.039, 0.0297, 
            0.0311, 0.00987, 0.0347, 0.0917),
        cuteseedopcor = cms.vdouble(0.637, 0.943, 0.742, 0.748, 0.763, 
            0.631, 0.214, 0.873, 0.473),
        cutdphiinl = cms.vdouble(0.061, 0.14, 0.286, 0.0921, 0.197, 
            0.24, 0.333, 0.303, 0.258),
        cutdetainl = cms.vdouble(0.00816, 0.00401, 0.0081, 0.019, 0.00588, 
            0.00893, 0.0171, 0.0434, 0.0143)
    ),
    algorithm = cms.string('eIDCB'),
    robusthighenergyEleIDCuts = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCuts = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV02 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV03 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV00 = cms.PSet(
        barrel = cms.vdouble(0.115, 0.014, 0.09, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.15, 0.0275, 0.092, 0.0105, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV01 = cms.PSet(
        barrel = cms.vdouble(0.075, 0.0132, 0.058, 0.0077, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.083, 0.027, 0.042, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    additionalCategories = cms.bool(True),
    etBinning = cms.bool(True)
)


process.eidTight = cms.EDProducer("EleIdCutBasedExtProducer",
    electronQuality = cms.string('tight'),
    classbasedtightEleIDCutsV02 = cms.PSet(
        cutisohcal = cms.vdouble(10.9, 7.01, 8.75, 3.51, 7.75, 
            1.62, 11.6, 9.9, 4.97, 5.33, 
            3.18, 2.32, 0.164, 5.46, 12.0, 
            0.00604, 4.1, 0.000628),
        cutmishits = cms.vdouble(5.5, 1.5, 0.5, 1.5, 2.5, 
            0.5, 3.5, 5.5, 0.5, 0.5, 
            0.5, 0.5, 0.5, 1.5, 0.5, 
            0.5, 0.5, 0.5),
        cuthoe = cms.vdouble(0.0871, 0.0289, 0.0783, 0.0946, 0.0245, 
            0.0363, 0.0671, 0.048, 0.0614, 0.0924, 
            0.0158, 0.049, 0.0382, 0.0915, 0.0451, 
            0.0452, 0.00196, 0.0043),
        cutdeta = cms.vdouble(0.00915, 0.00302, 0.0061, 0.0135, 0.00565, 
            0.00793, 0.0102, 0.00266, 0.0106, 0.00903, 
            0.00766, 0.00723, 0.0116, 0.00203, 0.00659, 
            0.0148, 0.00555, 0.0128),
        cuteopin = cms.vdouble(0.878, 0.859, 0.874, 0.944, 0.737, 
            0.773, 0.86, 0.967, 0.917, 0.812, 
            0.915, 1.01, 0.847, 0.953, 0.979, 
            0.841, 0.771, 1.09),
        cutip = cms.vdouble(0.0239, 0.027, 0.0768, 0.0231, 0.178, 
            0.0957, 0.0102, 0.0168, 0.043, 0.0166, 
            0.0594, 0.0308, 2.1, 0.00527, 3.17, 
            4.91, 0.769, 5.9),
        cutisotk = cms.vdouble(6.53, 4.6, 6.0, 8.63, 3.11, 
            7.77, 5.42, 4.81, 4.06, 6.47, 
            2.8, 3.45, 5.29, 5.18, 15.4, 
            5.38, 4.47, 0.0347),
        cutsee = cms.vdouble(0.0131, 0.0106, 0.0115, 0.0306, 0.028, 
            0.0293, 0.0131, 0.0106, 0.0115, 0.0317, 
            0.029, 0.0289, 0.0142, 0.0106, 0.0103, 
            0.035, 0.0296, 0.0333),
        cutdphi = cms.vdouble(0.0369, 0.0307, 0.117, 0.0475, 0.0216, 
            0.117, 0.0372, 0.0246, 0.0426, 0.0612, 
            0.0142, 0.039, 0.0737, 0.0566, 0.0359, 
            0.0187, 0.012, 0.0358),
        cutisoecal = cms.vdouble(20.0, 27.2, 4.48, 13.5, 4.56, 
            3.19, 12.2, 13.1, 7.42, 7.67, 
            4.12, 4.85, 10.1, 12.4, 11.1, 
            11.0, 10.6, 13.4)
    ),
    classbasedtightEleIDCutsV03 = cms.PSet(
        cutdetain = cms.vdouble(0.00811, 0.00341, 0.00633, 0.0103, 0.00667, 
            0.01, 0.0106, 0.0145, 0.0163, 0.0076, 
            0.00259, 0.00511, 0.00941, 0.0043, 0.00857, 
            0.012, 0.0169, 0.00172, 0.00861, 0.00362, 
            0.00601, 0.00925, 0.00489, 0.00832, 0.0119, 
            0.0169, 0.000996),
        cutiso_sum = cms.vdouble(11.8, 8.31, 6.26, 6.18, 3.28, 
            4.38, 4.17, 5.4, 1.57, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0213, 0.0422, 0.0632, 0.0361, 0.073, 
            0.126, 0.171, 0.119, 0.0372, 0.0131, 
            0.0146, 0.0564, 0.0152, 0.0222, 0.0268, 
            0.0314, 0.0884, 0.00374, 0.00852, 0.00761, 
            0.0143, 0.0106, 0.0127, 0.0119, 0.0123, 
            0.0235, 0.00363),
        cuthoe = cms.vdouble(0.0783, 0.0387, 0.105, 0.118, 0.0227, 
            0.062, 0.13, 2.47, 0.38, 0.0888, 
            0.0503, 0.0955, 0.0741, 0.015, 0.03, 
            0.589, 1.13, 0.612, 0.0494, 0.0461, 
            0.0292, 0.0369, 0.0113, 0.0145, 0.124, 
            2.05, 0.61),
        cutfmishits = cms.vdouble(2.5, 1.5, 1.5, 1.5, 1.5, 
            0.5, 2.5, 0.5, 0.5, 2.5, 
            1.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5, -0.5, 2.5, 1.5, 
            0.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(13.7, 11.6, 7.14, 9.98, 3.52, 
            4.87, 6.24, 7.96, 2.53, 11.2, 
            11.9, 7.88, 8.16, 5.58, 5.03, 
            11.4, 8.15, 5.79, 10.4, 11.1, 
            10.4, 7.47, 5.08, 5.9, 11.8, 
            14.1, 11.7),
        cutdcotdist = cms.vdouble(0.0393, 0.0256, 0.00691, 0.0394, 0.0386, 
            0.039, 0.0325, 0.0384, 0.0382, 0.0245, 
            0.000281, 5.46e-05, 0.0342, 0.0232, 0.00107, 
            0.0178, 0.0193, 0.000758, 0.000108, 0.0248, 
            0.000458, 0.0129, 0.00119, 0.0182, 4.53e-05, 
            0.0189, 0.000928),
        cutsee = cms.vdouble(0.0143, 0.0105, 0.0123, 0.0324, 0.0307, 
            0.0301, 0.0109, 0.027, 0.0292, 0.0133, 
            0.0104, 0.0116, 0.0332, 0.0296, 0.031, 
            0.00981, 0.0307, 0.072, 0.0149, 0.0105, 
            0.011, 0.0342, 0.0307, 0.0303, 0.00954, 
            0.0265, 0.0101),
        cuteseedopcor = cms.vdouble(0.784, 0.366, 0.57, 0.911, 0.298, 
            0.645, 0.51, 0.497, 0.932, 0.835, 
            0.968, 0.969, 0.923, 0.898, 0.98, 
            0.63, 0.971, 1.0, 0.515, 0.963, 
            0.986, 0.823, 0.879, 1.01, 0.931, 
            0.937, 1.05),
        cutdphiin = cms.vdouble(0.0404, 0.0499, 0.263, 0.042, 0.0484, 
            0.241, 0.242, 0.231, 0.286, 0.0552, 
            0.0338, 0.154, 0.0623, 0.0183, 0.0392, 
            0.0547, 0.0588, 0.00654, 0.042, 0.0217, 
            0.0885, 0.0445, 0.0141, 0.0234, 0.065, 
            0.0258, 0.0346),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 13.7, 13.2, 
            13.6, 14.2, 14.1, 13.9, 12.9, 
            14.9, 17.7)
    ),
    classbasedtightEleIDCutsV00 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.032, 0.016, 0.0525, 0.09, 0.025, 
            0.035, 0.065, 0.092),
        hOverE = cms.vdouble(0.05, 0.042, 0.045, 0.0, 0.055, 
            0.037, 0.05, 0.0),
        sigmaEtaEta = cms.vdouble(0.0125, 0.011, 0.01, 0.0, 0.0265, 
            0.0252, 0.026, 0.0),
        deltaEtaIn = cms.vdouble(0.0055, 0.003, 0.0065, 0.0, 0.006, 
            0.0055, 0.0075, 0.0),
        eSeedOverPin = cms.vdouble(0.24, 0.94, 0.11, 0.0, 0.32, 
            0.83, 0.0, 0.0)
    ),
    classbasedtightEleIDCutsV01 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.0225, 0.0114, 0.0234, 0.039, 0.0215, 
            0.0095, 0.0148, 0.0167),
        hOverE = cms.vdouble(0.056, 0.0221, 0.037, 0.0, 0.0268, 
            0.0102, 0.0104, 0.0),
        sigmaEtaEta = cms.vdouble(0.0095, 0.0094, 0.0094, 0.0, 0.026, 
            0.0257, 0.0246, 0.0),
        deltaEtaIn = cms.vdouble(0.0043, 0.00282, 0.0036, 0.0, 0.0066, 
            0.0049, 0.0041, 0.0),
        eSeedOverPin = cms.vdouble(0.32, 0.94, 0.221, 0.0, 0.74, 
            0.89, 0.66, 0.0)
    ),
    classbasedtightEleIDCutsV06 = cms.PSet(
        cutdetain = cms.vdouble(0.0116, 0.00449, 0.00938, 0.0184, 0.00678, 
            0.0109, 0.0252, 0.0268, 0.0139),
        cutiso_sum = cms.vdouble(15.5, 12.2, 12.2, 11.7, 7.16, 
            9.71, 8.66, 11.9, 2.98),
        cutip_gsf = cms.vdouble(0.0131, 0.0586, 0.0839, 0.0366, 0.452, 
            0.204, 0.0913, 0.0802, 0.0731),
        cutip_gsfl = cms.vdouble(0.0119, 0.0527, 0.0471, 0.0212, 0.233, 
            0.267, 0.109, 0.122, 0.0479),
        cuthoe = cms.vdouble(0.215, 0.0608, 0.147, 0.369, 0.0349, 
            0.102, 0.52, 0.422, 0.404),
        cutiso_sumoetl = cms.vdouble(6.21, 6.81, 5.3, 5.39, 2.73, 
            4.73, 4.84, 3.46, 3.73),
        cutfmishits = cms.vdouble(1.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 1.5, 2.5, 0.5),
        cuthoel = cms.vdouble(0.228, 0.0836, 0.143, 0.37, 0.0392, 
            0.0979, 0.3, 0.381, 0.339),
        cutdphiin = cms.vdouble(0.0897, 0.0993, 0.295, 0.0979, 0.151, 
            0.252, 0.341, 0.308, 0.328),
        cutseel = cms.vdouble(0.0132, 0.0117, 0.0112, 0.0387, 0.0281, 
            0.0287, 0.00987, 0.0296, 0.0544),
        cutiso_sumoet = cms.vdouble(11.9, 7.81, 6.28, 8.92, 4.65, 
            5.49, 9.36, 8.84, 5.94),
        cutsee = cms.vdouble(0.0145, 0.0116, 0.012, 0.039, 0.0297, 
            0.0311, 0.00987, 0.0347, 0.0917),
        cuteseedopcor = cms.vdouble(0.637, 0.943, 0.742, 0.748, 0.763, 
            0.631, 0.214, 0.873, 0.473),
        cutdphiinl = cms.vdouble(0.061, 0.14, 0.286, 0.0921, 0.197, 
            0.24, 0.333, 0.303, 0.258),
        cutdetainl = cms.vdouble(0.00816, 0.00401, 0.0081, 0.019, 0.00588, 
            0.00893, 0.0171, 0.0434, 0.0143)
    ),
    classbasedtightEleIDCutsV04 = cms.PSet(
        cutdetain = cms.vdouble(0.00811, 0.00341, 0.00633, 0.0103, 0.00667, 
            0.01, 0.0106, 0.0145, 0.0163, 0.0076, 
            0.00259, 0.00511, 0.00941, 0.0043, 0.00857, 
            0.012, 0.0169, 0.00172, 0.00861, 0.00362, 
            0.00601, 0.00925, 0.00489, 0.00832, 0.0119, 
            0.0169, 0.000996),
        cutiso_sum = cms.vdouble(11.8, 8.31, 6.26, 6.18, 3.28, 
            4.38, 4.17, 5.4, 1.57, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0213, 0.0422, 0.0632, 0.0361, 0.073, 
            0.126, 0.171, 0.119, 0.0372, 0.0131, 
            0.0146, 0.0564, 0.0152, 0.0222, 0.0268, 
            0.0314, 0.0884, 0.00374, 0.00852, 0.00761, 
            0.0143, 0.0106, 0.0127, 0.0119, 0.0123, 
            0.0235, 0.00363),
        cuthoe = cms.vdouble(0.0783, 0.0387, 0.105, 0.118, 0.0227, 
            0.062, 0.13, 2.47, 0.38, 0.0888, 
            0.0503, 0.0955, 0.0741, 0.015, 0.03, 
            0.589, 1.13, 0.612, 0.0494, 0.0461, 
            0.0292, 0.0369, 0.0113, 0.0145, 0.124, 
            2.05, 0.61),
        cutfmishits = cms.vdouble(2.5, 1.5, 1.5, 1.5, 1.5, 
            0.5, 2.5, 0.5, 0.5, 2.5, 
            1.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5, -0.5, 2.5, 1.5, 
            0.5, 0.5, 0.5, 0.5, 0.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(13.7, 11.6, 7.14, 9.98, 3.52, 
            4.87, 6.24, 7.96, 2.53, 11.2, 
            11.9, 7.88, 8.16, 5.58, 5.03, 
            11.4, 8.15, 5.79, 10.4, 11.1, 
            10.4, 7.47, 5.08, 5.9, 11.8, 
            14.1, 11.7),
        cutdcotdist = cms.vdouble(0.0393, 0.0256, 0.00691, 0.0394, 0.0386, 
            0.039, 0.0325, 0.0384, 0.0382, 0.0245, 
            0.000281, 5.46e-05, 0.0342, 0.0232, 0.00107, 
            0.0178, 0.0193, 0.000758, 0.000108, 0.0248, 
            0.000458, 0.0129, 0.00119, 0.0182, 4.53e-05, 
            0.0189, 0.000928),
        cutsee = cms.vdouble(0.0143, 0.0105, 0.0123, 0.0324, 0.0307, 
            0.0301, 0.0109, 0.027, 0.0292, 0.0133, 
            0.0104, 0.0116, 0.0332, 0.0296, 0.031, 
            0.00981, 0.0307, 0.072, 0.0149, 0.0105, 
            0.011, 0.0342, 0.0307, 0.0303, 0.00954, 
            0.0265, 0.0101),
        cuteseedopcor = cms.vdouble(0.784, 0.366, 0.57, 0.911, 0.298, 
            0.645, 0.51, 0.497, 0.932, 0.835, 
            0.968, 0.969, 0.923, 0.898, 0.98, 
            0.63, 0.971, 1.0, 0.515, 0.963, 
            0.986, 0.823, 0.879, 1.01, 0.931, 
            0.937, 1.05),
        cutdphiin = cms.vdouble(0.0404, 0.0499, 0.263, 0.042, 0.0484, 
            0.241, 0.242, 0.231, 0.286, 0.0552, 
            0.0338, 0.154, 0.0623, 0.0183, 0.0392, 
            0.0547, 0.0588, 0.00654, 0.042, 0.0217, 
            0.0885, 0.0445, 0.0141, 0.0234, 0.065, 
            0.0258, 0.0346),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 13.7, 13.2, 
            13.6, 14.2, 14.1, 13.9, 12.9, 
            14.9, 17.7)
    ),
    electronIDType = cms.string('classbased'),
    robusttightEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    electronVersion = cms.string(''),
    robusttightEleIDCutsV00 = cms.PSet(
        barrel = cms.vdouble(0.015, 0.0092, 0.02, 0.0025, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.018, 0.025, 0.02, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusttightEleIDCutsV01 = cms.PSet(
        barrel = cms.vdouble(0.01, 0.0099, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.01, 0.028, 0.02, 0.0066, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    reducedBarrelRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    robusttightEleIDCutsV03 = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    reducedEndcapRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    verticesCollection = cms.InputTag("offlinePrimaryVerticesWithBS"),
    classbasedlooseEleIDCuts = cms.PSet(
        cutdetain = cms.vdouble(0.0137, 0.00678, 0.0241, 0.0187, 0.0161, 
            0.0224, 0.0252, 0.0308, 0.0273),
        cutiso_sum = cms.vdouble(33.0, 17.0, 17.9, 18.8, 8.55, 
            12.5, 17.6, 18.5, 2.98),
        cutip_gsf = cms.vdouble(0.0551, 0.0765, 0.143, 0.0874, 0.594, 
            0.37, 0.0913, 1.15, 0.231),
        cutip_gsfl = cms.vdouble(0.0186, 0.0759, 0.138, 0.0473, 0.62, 
            0.304, 0.109, 0.775, 0.0479),
        cuthoe = cms.vdouble(0.247, 0.137, 0.147, 0.371, 0.0588, 
            0.147, 0.52, 0.452, 0.404),
        cutiso_sumoetl = cms.vdouble(11.3, 9.05, 9.07, 9.94, 5.25, 
            6.15, 10.7, 10.8, 4.4),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 4.5, 3.5, 3.5),
        cuthoel = cms.vdouble(0.236, 0.126, 0.147, 0.375, 0.0392, 
            0.145, 0.365, 0.383, 0.384),
        cutdphiin = cms.vdouble(0.0897, 0.262, 0.353, 0.116, 0.357, 
            0.319, 0.342, 0.404, 0.336),
        cutseel = cms.vdouble(0.0164, 0.0118, 0.015, 0.0523, 0.0326, 
            0.0456, 0.0185, 0.0589, 0.0544),
        cutiso_sumoet = cms.vdouble(34.5, 12.7, 12.1, 19.9, 6.35, 
            8.85, 14.0, 10.5, 9.74),
        cutsee = cms.vdouble(0.0176, 0.0125, 0.0181, 0.0415, 0.0364, 
            0.0418, 0.0146, 0.0678, 0.133),
        cuteseedopcor = cms.vdouble(0.63, 0.82, 0.401, 0.718, 0.4, 
            0.458, 0.15, 0.664, 0.373),
        cutdphiinl = cms.vdouble(0.0747, 0.25, 0.356, 0.0956, 0.347, 
            0.326, 0.333, 0.647, 0.289),
        cutdetainl = cms.vdouble(0.0124, 0.00503, 0.0257, 0.0228, 0.0118, 
            0.0178, 0.0188, 0.14, 0.024)
    ),
    robusttightEleIDCutsV02 = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV01 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.0275, 0.09, 0.007, -1, 
            -1, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV00 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.011, 0.09, 0.005, -1, 
            -1, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.1, 0.0275, 0.09, 0.007, -1, 
            -1, 9999.0, 9999.0, 0, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV03 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robusthighenergyEleIDCutsV02 = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    classbasedlooseEleIDCutsV00 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.05, 0.025, 0.053, 0.09, 0.07, 
            0.03, 0.092, 0.092),
        hOverE = cms.vdouble(0.115, 0.1, 0.055, 0.0, 0.145, 
            0.12, 0.15, 0.0),
        sigmaEtaEta = cms.vdouble(0.014, 0.012, 0.0115, 0.0, 0.0275, 
            0.0265, 0.0265, 0.0),
        deltaEtaIn = cms.vdouble(0.009, 0.0045, 0.0085, 0.0, 0.0105, 
            0.0068, 0.01, 0.0),
        eSeedOverPin = cms.vdouble(0.11, 0.91, 0.11, 0.0, 0.0, 
            0.85, 0.0, 0.0)
    ),
    classbasedlooseEleIDCutsV01 = cms.PSet(
        deltaPhiIn = cms.vdouble(0.053, 0.0189, 0.059, 0.099, 0.0278, 
            0.0157, 0.042, 0.08),
        hOverE = cms.vdouble(0.076, 0.033, 0.07, 0.0, 0.083, 
            0.0148, 0.033, 0.0),
        sigmaEtaEta = cms.vdouble(0.0101, 0.0095, 0.0097, 0.0, 0.0271, 
            0.0267, 0.0259, 0.0),
        deltaEtaIn = cms.vdouble(0.0078, 0.00259, 0.0062, 0.0, 0.0078, 
            0.0061, 0.0061, 0.0),
        eSeedOverPin = cms.vdouble(0.3, 0.92, 0.211, 0.0, 0.42, 
            0.88, 0.68, 0.0)
    ),
    classbasedlooseEleIDCutsV02 = cms.PSet(
        cutisohcal = cms.vdouble(13.5, 9.93, 7.56, 14.8, 8.1, 
            10.8, 42.7, 20.1, 9.11, 10.4, 
            6.89, 5.59, 8.53, 9.59, 24.2, 
            2.78, 8.67, 0.288),
        cutmishits = cms.vdouble(5.5, 1.5, 5.5, 2.5, 2.5, 
            2.5, 3.5, 5.5, 0.5, 1.5, 
            2.5, 0.5, 1.5, 1.5, 0.5, 
            0.5, 0.5, 0.5),
        cuthoe = cms.vdouble(0.0887, 0.0934, 0.0949, 0.0986, 0.0431, 
            0.0878, 0.097, 0.0509, 0.098, 0.0991, 
            0.0321, 0.0928, 0.0663, 0.0717, 0.0966, 
            0.0758, 0.0149, 0.0131),
        cutdeta = cms.vdouble(0.00958, 0.00406, 0.0122, 0.0137, 0.00837, 
            0.0127, 0.011, 0.00336, 0.00977, 0.015, 
            0.00675, 0.0109, 0.014, 0.00508, 0.0109, 
            0.0146, 0.00506, 0.0127),
        cuteopin = cms.vdouble(0.878, 0.802, 0.814, 0.942, 0.735, 
            0.774, 0.829, 0.909, 0.829, 0.813, 
            0.86, 0.897, 0.817, 0.831, 0.818, 
            0.861, 0.787, 0.789),
        cutip = cms.vdouble(0.0246, 0.076, 0.0966, 0.0885, 0.441, 
            0.205, 0.0292, 0.0293, 0.0619, 0.0251, 
            0.159, 0.0815, 7.29, 0.0106, 5.76, 
            6.89, 1.27, 5.89),
        cutisotk = cms.vdouble(24.3, 8.45, 14.4, 27.8, 6.02, 
            10.5, 14.1, 10.2, 14.5, 19.1, 
            6.1, 14.1, 8.59, 8.33, 8.3, 
            8.93, 8.6, 16.0),
        cutsee = cms.vdouble(0.0172, 0.0115, 0.0143, 0.0344, 0.0295, 
            0.0304, 0.0145, 0.0108, 0.0128, 0.0347, 
            0.0307, 0.0316, 0.018, 0.011, 0.0132, 
            0.0349, 0.031, 0.0327),
        cutdphi = cms.vdouble(0.0372, 0.114, 0.118, 0.0488, 0.117, 
            0.119, 0.0606, 0.0548, 0.117, 0.07, 
            0.0355, 0.117, 0.088, 0.045, 0.118, 
            0.0919, 0.0236, 0.0515),
        cutisoecal = cms.vdouble(33.4, 28.1, 7.32, 27.4, 7.33, 
            21.7, 93.8, 102.0, 12.1, 26.0, 
            8.91, 10.0, 16.1, 31.3, 16.9, 
            15.4, 13.3, 37.7)
    ),
    classbasedlooseEleIDCutsV03 = cms.PSet(
        cutdetain = cms.vdouble(0.00989, 0.00484, 0.0146, 0.0146, 0.00902, 
            0.0172, 0.0137, 0.0477, 0.0275, 0.00967, 
            0.00377, 0.00924, 0.013, 0.00666, 0.0123, 
            0.0125, 0.0228, 0.0112, 0.0106, 0.0038, 
            0.00897, 0.0139, 0.00667, 0.0122, 0.0122, 
            0.0193, 0.00239),
        cutiso_sum = cms.vdouble(31.5, 10.3, 8.8, 11.0, 6.13, 
            6.94, 7.52, 9.0, 3.5, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0431, 0.0767, 0.139, 0.101, 0.149, 
            0.154, 0.932, 0.15, 0.124, 0.0238, 
            0.0467, 0.0759, 0.0369, 0.147, 0.0986, 
            0.0626, 0.195, 0.116, 0.0122, 0.0125, 
            0.0693, 0.0162, 0.089, 0.0673, 0.0467, 
            0.0651, 0.0221),
        cuthoe = cms.vdouble(0.166, 0.0771, 0.144, 0.37, 0.0497, 
            0.139, 0.401, 2.68, 0.516, 0.234, 
            0.0556, 0.144, 0.368, 0.031, 0.12, 
            0.602, 2.01, 1.05, 0.104, 0.063, 
            0.0565, 0.38, 0.0192, 0.0294, 0.537, 
            4.65, 1.87),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 2.5, 2.5, 1.5, 2.5, 
            1.5, 1.5, 1.5, 1.5, 0.5, 
            2.5, 2.5, 0.5, 2.5, 1.5, 
            0.5, 1.5, 1.5, 0.5, 2.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(28.9, 15.3, 12.0, 18.3, 7.17, 
            9.42, 11.0, 9.81, 3.94, 22.7, 
            15.9, 12.3, 17.0, 7.58, 8.89, 
            15.2, 12.7, 6.17, 20.8, 21.2, 
            17.2, 15.5, 9.37, 10.6, 19.8, 
            22.1, 15.6),
        cutdcotdist = cms.vdouble(0.0393, 0.0392, 0.0397, 0.0394, 0.0393, 
            0.039, 0.0378, 0.0388, 0.0382, 0.0385, 
            0.0167, 0.00325, 0.0394, 0.0387, 0.0388, 
            0.0227, 0.0258, 0.0127, 0.0298, 0.03, 
            0.00946, 0.039, 0.0231, 0.0278, 0.00162, 
            0.0367, 0.0199),
        cutsee = cms.vdouble(0.0175, 0.0127, 0.0177, 0.0373, 0.0314, 
            0.0329, 0.0157, 0.0409, 0.14, 0.0169, 
            0.0106, 0.0142, 0.0363, 0.0322, 0.0354, 
            0.0117, 0.0372, 28.2, 0.0171, 0.0113, 
            0.014, 0.0403, 0.0323, 0.0411, 0.0104, 
            0.0436, 0.0114),
        cuteseedopcor = cms.vdouble(0.78, 0.302, 0.483, 0.904, 0.168, 
            0.645, 0.108, 0.284, 0.324, 0.591, 
            0.286, 0.488, 0.813, 0.791, 0.672, 
            0.398, 0.834, 0.878, 0.515, 0.937, 
            0.806, 0.816, 0.85, 0.507, 0.367, 
            0.83, 0.648),
        cutdphiin = cms.vdouble(0.041, 0.275, 0.365, 0.047, 0.273, 
            0.296, 0.329, 0.465, 0.627, 0.0581, 
            0.0954, 0.327, 0.0702, 0.0582, 0.279, 
            0.117, 0.318, 0.246, 0.0821, 0.052, 
            0.292, 0.116, 0.0435, 0.312, 0.118, 
            0.296, 0.0459),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 12.0, 12.0, 
            12.0, 12.0, 12.0, 12.0, 12.0, 
            12.0, 12.5)
    ),
    classbasedlooseEleIDCutsV04 = cms.PSet(
        cutdetain = cms.vdouble(0.00989, 0.00484, 0.0146, 0.0146, 0.00902, 
            0.0172, 0.0137, 0.0477, 0.0275, 0.00967, 
            0.00377, 0.00924, 0.013, 0.00666, 0.0123, 
            0.0125, 0.0228, 0.0112, 0.0106, 0.0038, 
            0.00897, 0.0139, 0.00667, 0.0122, 0.0122, 
            0.0193, 0.00239),
        cutiso_sum = cms.vdouble(31.5, 10.3, 8.8, 11.0, 6.13, 
            6.94, 7.52, 9.0, 3.5, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 
            100000.0, 100000.0),
        cutip_gsf = cms.vdouble(0.0431, 0.0767, 0.139, 0.101, 0.149, 
            0.154, 0.932, 0.15, 0.124, 0.0238, 
            0.0467, 0.0759, 0.0369, 0.147, 0.0986, 
            0.0626, 0.195, 0.116, 0.0122, 0.0125, 
            0.0693, 0.0162, 0.089, 0.0673, 0.0467, 
            0.0651, 0.0221),
        cuthoe = cms.vdouble(0.166, 0.0771, 0.144, 0.37, 0.0497, 
            0.139, 0.401, 2.68, 0.516, 0.234, 
            0.0556, 0.144, 0.368, 0.031, 0.12, 
            0.602, 2.01, 1.05, 0.104, 0.063, 
            0.0565, 0.38, 0.0192, 0.0294, 0.537, 
            4.65, 1.87),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 2.5, 2.5, 1.5, 2.5, 
            1.5, 1.5, 1.5, 1.5, 0.5, 
            2.5, 2.5, 0.5, 2.5, 1.5, 
            0.5, 1.5, 1.5, 0.5, 2.5, 
            0.5, 0.5),
        cutiso_sumoet = cms.vdouble(28.9, 15.3, 12.0, 18.3, 7.17, 
            9.42, 11.0, 9.81, 3.94, 22.7, 
            15.9, 12.3, 17.0, 7.58, 8.89, 
            15.2, 12.7, 6.17, 20.8, 21.2, 
            17.2, 15.5, 9.37, 10.6, 19.8, 
            22.1, 15.6),
        cutdcotdist = cms.vdouble(0.0393, 0.0392, 0.0397, 0.0394, 0.0393, 
            0.039, 0.0378, 0.0388, 0.0382, 0.0385, 
            0.0167, 0.00325, 0.0394, 0.0387, 0.0388, 
            0.0227, 0.0258, 0.0127, 0.0298, 0.03, 
            0.00946, 0.039, 0.0231, 0.0278, 0.00162, 
            0.0367, 0.0199),
        cutsee = cms.vdouble(0.0175, 0.0127, 0.0177, 0.0373, 0.0314, 
            0.0329, 0.0157, 0.0409, 0.14, 0.0169, 
            0.0106, 0.0142, 0.0363, 0.0322, 0.0354, 
            0.0117, 0.0372, 28.2, 0.0171, 0.0113, 
            0.014, 0.0403, 0.0323, 0.0411, 0.0104, 
            0.0436, 0.0114),
        cuteseedopcor = cms.vdouble(0.78, 0.302, 0.483, 0.904, 0.168, 
            0.645, 0.108, 0.284, 0.324, 0.591, 
            0.286, 0.488, 0.813, 0.791, 0.672, 
            0.398, 0.834, 0.878, 0.515, 0.937, 
            0.806, 0.816, 0.85, 0.507, 0.367, 
            0.83, 0.648),
        cutdphiin = cms.vdouble(0.041, 0.275, 0.365, 0.047, 0.273, 
            0.296, 0.329, 0.465, 0.627, 0.0581, 
            0.0954, 0.327, 0.0702, 0.0582, 0.279, 
            0.117, 0.318, 0.246, 0.0821, 0.052, 
            0.292, 0.116, 0.0435, 0.312, 0.118, 
            0.296, 0.0459),
        cutet = cms.vdouble(-100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, -100000.0, -100000.0, 
            -100000.0, -100000.0, -100000.0, 12.0, 12.0, 
            12.0, 12.0, 12.0, 12.0, 12.0, 
            12.0, 12.5)
    ),
    classbasedlooseEleIDCutsV06 = cms.PSet(
        cutdetain = cms.vdouble(0.0137, 0.00678, 0.0241, 0.0187, 0.0161, 
            0.0224, 0.0252, 0.0308, 0.0273),
        cutiso_sum = cms.vdouble(33.0, 17.0, 17.9, 18.8, 8.55, 
            12.5, 17.6, 18.5, 2.98),
        cutip_gsf = cms.vdouble(0.0551, 0.0765, 0.143, 0.0874, 0.594, 
            0.37, 0.0913, 1.15, 0.231),
        cutip_gsfl = cms.vdouble(0.0186, 0.0759, 0.138, 0.0473, 0.62, 
            0.304, 0.109, 0.775, 0.0479),
        cuthoe = cms.vdouble(0.247, 0.137, 0.147, 0.371, 0.0588, 
            0.147, 0.52, 0.452, 0.404),
        cutiso_sumoetl = cms.vdouble(11.3, 9.05, 9.07, 9.94, 5.25, 
            6.15, 10.7, 10.8, 4.4),
        cutfmishits = cms.vdouble(4.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 4.5, 3.5, 3.5),
        cuthoel = cms.vdouble(0.236, 0.126, 0.147, 0.375, 0.0392, 
            0.145, 0.365, 0.383, 0.384),
        cutdphiin = cms.vdouble(0.0897, 0.262, 0.353, 0.116, 0.357, 
            0.319, 0.342, 0.404, 0.336),
        cutseel = cms.vdouble(0.0164, 0.0118, 0.015, 0.0523, 0.0326, 
            0.0456, 0.0185, 0.0589, 0.0544),
        cutiso_sumoet = cms.vdouble(34.5, 12.7, 12.1, 19.9, 6.35, 
            8.85, 14.0, 10.5, 9.74),
        cutsee = cms.vdouble(0.0176, 0.0125, 0.0181, 0.0415, 0.0364, 
            0.0418, 0.0146, 0.0678, 0.133),
        cuteseedopcor = cms.vdouble(0.63, 0.82, 0.401, 0.718, 0.4, 
            0.458, 0.15, 0.664, 0.373),
        cutdphiinl = cms.vdouble(0.0747, 0.25, 0.356, 0.0956, 0.347, 
            0.326, 0.333, 0.647, 0.289),
        cutdetainl = cms.vdouble(0.0124, 0.00503, 0.0257, 0.0228, 0.0118, 
            0.0178, 0.0188, 0.14, 0.024)
    ),
    src = cms.InputTag("gsfElectrons"),
    robusttightEleIDCuts = cms.PSet(
        barrel = cms.vdouble(0.0201, 0.0102, 0.0211, 0.00606, -1, 
            -1, 2.34, 3.24, 4.51, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.00253, 0.0291, 0.022, 0.0032, -1, 
            -1, 0.826, 2.7, 0.255, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    classbasedtightEleIDCuts = cms.PSet(
        cutdetain = cms.vdouble(0.0116, 0.00449, 0.00938, 0.0184, 0.00678, 
            0.0109, 0.0252, 0.0268, 0.0139),
        cutiso_sum = cms.vdouble(15.5, 12.2, 12.2, 11.7, 7.16, 
            9.71, 8.66, 11.9, 2.98),
        cutip_gsf = cms.vdouble(0.0131, 0.0586, 0.0839, 0.0366, 0.452, 
            0.204, 0.0913, 0.0802, 0.0731),
        cutip_gsfl = cms.vdouble(0.0119, 0.0527, 0.0471, 0.0212, 0.233, 
            0.267, 0.109, 0.122, 0.0479),
        cuthoe = cms.vdouble(0.215, 0.0608, 0.147, 0.369, 0.0349, 
            0.102, 0.52, 0.422, 0.404),
        cutiso_sumoetl = cms.vdouble(6.21, 6.81, 5.3, 5.39, 2.73, 
            4.73, 4.84, 3.46, 3.73),
        cutfmishits = cms.vdouble(1.5, 1.5, 1.5, 2.5, 2.5, 
            1.5, 1.5, 2.5, 0.5),
        cuthoel = cms.vdouble(0.228, 0.0836, 0.143, 0.37, 0.0392, 
            0.0979, 0.3, 0.381, 0.339),
        cutdphiin = cms.vdouble(0.0897, 0.0993, 0.295, 0.0979, 0.151, 
            0.252, 0.341, 0.308, 0.328),
        cutseel = cms.vdouble(0.0132, 0.0117, 0.0112, 0.0387, 0.0281, 
            0.0287, 0.00987, 0.0296, 0.0544),
        cutiso_sumoet = cms.vdouble(11.9, 7.81, 6.28, 8.92, 4.65, 
            5.49, 9.36, 8.84, 5.94),
        cutsee = cms.vdouble(0.0145, 0.0116, 0.012, 0.039, 0.0297, 
            0.0311, 0.00987, 0.0347, 0.0917),
        cuteseedopcor = cms.vdouble(0.637, 0.943, 0.742, 0.748, 0.763, 
            0.631, 0.214, 0.873, 0.473),
        cutdphiinl = cms.vdouble(0.061, 0.14, 0.286, 0.0921, 0.197, 
            0.24, 0.333, 0.303, 0.258),
        cutdetainl = cms.vdouble(0.00816, 0.00401, 0.0081, 0.019, 0.00588, 
            0.00893, 0.0171, 0.0434, 0.0143)
    ),
    algorithm = cms.string('eIDCB'),
    robusthighenergyEleIDCuts = cms.PSet(
        barrel = cms.vdouble(0.05, 9999, 0.09, 0.005, 0.94, 
            0.83, 7.5, 2, 0.03, 9999.0, 
            0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.05, 0.03, 0.09, 0.007, -1, 
            -1, 15, 2.5, 0.03, 2.5, 
            0, 0.5, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCuts = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV02 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV03 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV00 = cms.PSet(
        barrel = cms.vdouble(0.115, 0.014, 0.09, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.15, 0.0275, 0.092, 0.0105, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV01 = cms.PSet(
        barrel = cms.vdouble(0.075, 0.0132, 0.058, 0.0077, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.083, 0.027, 0.042, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    robustlooseEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.05, 0.0103, 0.8, 0.00688, -1, 
            -1, 7.33, 4.68, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0),
        endcap = cms.vdouble(0.0389, 0.0307, 0.7, 0.00944, -1, 
            -1, 7.76, 3.09, 2.23, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 9999, -1, 0, 
            0)
    ),
    additionalCategories = cms.bool(True),
    etBinning = cms.bool(True)
)


process.eleIsoDepositEcalFromHits = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("gsfElectrons"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        isolationVariable = cms.string('et'),
        spikeIdThreshold = cms.double(0.95),
        tryBoth = cms.bool(True),
        recHitFlagsToBeExcluded = cms.vint32(3, 4, 8, 9),
        ComponentName = cms.string('EgammaRecHitExtractor'),
        endcapEcalHits = cms.InputTag("reducedEcalRecHitsEE"),
        intStrip = cms.double(0.0),
        intRadius = cms.double(0.0),
        severityRecHitThreshold = cms.double(5.0),
        severityLevelCut = cms.int32(3),
        energyMin = cms.double(0.08),
        extRadius = cms.double(0.6),
        subtractSuperClusterEnergy = cms.bool(False),
        spikeIdString = cms.string('kSwissCross'),
        vetoClustered = cms.bool(False),
        etMin = cms.double(0.0),
        DepositLabel = cms.untracked.string(''),
        barrelEcalHits = cms.InputTag("reducedEcalRecHitsEB")
    )
)


process.eleIsoDepositHcalFromTowers = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("gsfElectrons"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        caloTowers = cms.InputTag("towerMaker"),
        ComponentName = cms.string('EgammaTowerExtractor'),
        hcalDepth = cms.int32(-1),
        intRadius = cms.double(0.0),
        extRadius = cms.double(0.6),
        DepositLabel = cms.untracked.string(''),
        etMin = cms.double(-999.0)
    )
)


process.eleIsoDepositTk = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("gsfElectrons"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(0.2),
        dzOption = cms.string('vz'),
        BeamSpotLabel = cms.InputTag("offlineBeamSpot"),
        ComponentName = cms.string('EgammaTrackExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(9999.0),
        Chi2Prob_Min = cms.double(-1.0),
        DR_Veto = cms.double(0.0),
        NHits_Min = cms.uint32(0),
        Chi2Ndof_Max = cms.double(1e+64),
        Pt_Min = cms.double(-1.0),
        DepositLabel = cms.untracked.string(''),
        BeamlineOption = cms.string('BeamSpotFromEvent'),
        inputTrackCollection = cms.InputTag("generalTracks")
    )
)


process.eleIsoFromDepsEcalFromHitsByCrystal = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        src = cms.InputTag("eleIsoDepositEcalFromHits"),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring('NumCrystalVeto(3.0)', 
            'NumCrystalEtaPhiVeto(1.5,9999.0)', 
            'EcalBarrel:AbsThresholdFromTransverse(0.08)', 
            'EcalEndcaps:AbsThreshold(0.100)'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
    ))
)


process.eleIsoFromDepsHcalFromTowers = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        src = cms.InputTag("eleIsoDepositHcalFromTowers"),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring('0.15'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
    ))
)


process.eleIsoFromDepsTk = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        src = cms.InputTag("eleIsoDepositTk"),
        deltaR = cms.double(0.3),
        weight = cms.string('1'),
        vetos = cms.vstring('RectangularEtaPhiVeto(-0.015,0.015,-0.5,0.5)', 
            'Threshold(0.7)'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
    ))
)


process.electronMatch = cms.EDProducer("MCMatcher",
    src = cms.InputTag("gsfElectrons"),
    maxDPtRel = cms.double(999999.0),
    mcPdgId = cms.vint32(11),
    mcStatus = cms.vint32(1),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = cms.double(0.2),
    checkCharge = cms.bool(False),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("genParticles")
)


process.electronMatchPF = cms.EDProducer("MCMatcher",
    src = cms.InputTag("gsfElectrons"),
    maxDPtRel = cms.double(999999.0),
    mcPdgId = cms.vint32(11),
    mcStatus = cms.vint32(1),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = cms.double(0.2),
    checkCharge = cms.bool(False),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("genParticles")
)


process.eventCountProducer = cms.EDProducer("EventCountProducer")


process.fixedConePFTauDecayModeIndexProducer = cms.EDProducer("PFRecoTauDecayModeIndexProducer",
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    PFTauProducer = cms.InputTag("fixedConePFTauProducer"),
    PFTauDecayModeProducer = cms.InputTag("fixedConePFTauDecayModeProducer")
)


process.fixedConePFTauDecayModeProducer = cms.EDProducer("PFRecoTauDecayModeDeterminator",
    mergeByBestMatch = cms.bool(True),
    refitTracks = cms.bool(False),
    maxPiZeroMass = cms.double(0.2),
    mergeLowPtPhotonsFirst = cms.bool(True),
    setMergedPi0Mass = cms.bool(True),
    setChargedPionMass = cms.bool(True),
    filterPhotons = cms.bool(True),
    minPtFractionSinglePhotons = cms.double(0.1),
    minPtFractionPiZeroes = cms.double(0.15),
    maxNbrOfIterations = cms.int32(10),
    filterTwoProngs = cms.bool(True),
    minPtFractionForSecondProng = cms.double(0.1),
    maxDistance = cms.double(0.01),
    setPi0Mass = cms.bool(True),
    maxPhotonsToMerge = cms.uint32(2),
    PFTauProducer = cms.InputTag("fixedConePFTauProducer")
)


process.fixedConePFTauDiscriminationAgainstElectron = cms.EDProducer("PFRecoTauDiscriminationAgainstElectron",
    ApplyCut_ElectronPreID_2D = cms.bool(False),
    ElecPreID0_HOverPLead_minValue = cms.double(0.05),
    PFTauProducer = cms.InputTag("fixedConePFTauProducer"),
    ApplyCut_ElectronPreID = cms.bool(False),
    ApplyCut_HcalTotOverPLead = cms.bool(False),
    EOverPLead_minValue = cms.double(0.8),
    ElecPreID1_EOverPLead_maxValue = cms.double(0.8),
    HcalMaxOverPLead_minValue = cms.double(0.1),
    ApplyCut_EmFraction = cms.bool(False),
    EmFraction_maxValue = cms.double(0.9),
    ApplyCut_PFElectronMVA = cms.bool(True),
    PFElectronMVA_maxValue = cms.double(-0.1),
    ApplyCut_HcalMaxOverPLead = cms.bool(False),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("fixedConePFTauDiscriminationByLeadingTrackFinding")
        )
    ),
    Hcal3x3OverPLead_minValue = cms.double(0.1),
    ElecPreID1_HOverPLead_minValue = cms.double(0.15),
    ElecPreID0_EOverPLead_maxValue = cms.double(0.95),
    BremsRecoveryEOverPLead_minValue = cms.double(0.8),
    ApplyCut_EcalCrackCut = cms.bool(False),
    EOverPLead_maxValue = cms.double(1.8),
    HcalTotOverPLead_minValue = cms.double(0.1),
    ApplyCut_BremsRecoveryEOverPLead = cms.bool(False),
    ApplyCut_Hcal3x3OverPLead = cms.bool(False),
    ApplyCut_EOverPLead = cms.bool(False),
    BremsRecoveryEOverPLead_maxValue = cms.double(1.8)
)


process.fixedConePFTauDiscriminationAgainstMuon = cms.EDProducer("PFRecoTauDiscriminationAgainstMuon",
    a = cms.double(0.5),
    c = cms.double(0.0),
    b = cms.double(0.5),
    PFTauProducer = cms.InputTag("fixedConePFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("fixedConePFTauDiscriminationByLeadingTrackFinding")
        )
    ),
    discriminatorOption = cms.string('noSegMatch')
)


process.fixedConePFTauDiscriminationByECALIsolation = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("fixedConePFTauProducer"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("fixedConePFTauDiscriminationByLeadingTrackFinding")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(False),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.fixedConePFTauDiscriminationByECALIsolationUsingLeadingPion = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("fixedConePFTauProducer"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadPion = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("fixedConePFTauDiscriminationByLeadingTrackFinding")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(False),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.fixedConePFTauDiscriminationByIsolation = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("fixedConePFTauProducer"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("fixedConePFTauDiscriminationByLeadingTrackFinding")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.fixedConePFTauDiscriminationByIsolationUsingLeadingPion = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("fixedConePFTauProducer"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadPion = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("fixedConePFTauDiscriminationByLeadingTrackFinding")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.fixedConePFTauDiscriminationByLeadingPionPtCut = cms.EDProducer("PFRecoTauDiscriminationByLeadingObjectPtCut",
    MinPtLeadingObject = cms.double(5.0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    PFTauProducer = cms.InputTag("fixedConePFTauProducer"),
    UseOnlyChargedHadrons = cms.bool(False)
)


process.fixedConePFTauDiscriminationByLeadingTrackFinding = cms.EDProducer("PFRecoTauDiscriminationByLeadingObjectPtCut",
    MinPtLeadingObject = cms.double(0.0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    PFTauProducer = cms.InputTag("fixedConePFTauProducer"),
    UseOnlyChargedHadrons = cms.bool(True)
)


process.fixedConePFTauDiscriminationByLeadingTrackPtCut = cms.EDProducer("PFRecoTauDiscriminationByLeadingObjectPtCut",
    MinPtLeadingObject = cms.double(5.0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    PFTauProducer = cms.InputTag("fixedConePFTauProducer"),
    UseOnlyChargedHadrons = cms.bool(True)
)


process.fixedConePFTauDiscriminationByTrackIsolation = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("fixedConePFTauProducer"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("fixedConePFTauDiscriminationByLeadingTrackFinding")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.fixedConePFTauDiscriminationByTrackIsolationUsingLeadingPion = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("fixedConePFTauProducer"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadPion = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("fixedConePFTauDiscriminationByLeadingTrackFinding")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.fixedConePFTauProducer = cms.EDProducer("PFRecoTauProducer",
    Rphi = cms.double(2.0),
    LeadTrack_minPt = cms.double(0.0),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    ECALSignalConeSizeFormula = cms.string('0.15'),
    TrackerIsolConeMetric = cms.string('DR'),
    TrackerSignalConeMetric = cms.string('DR'),
    EcalStripSumE_deltaPhiOverQ_minValue = cms.double(-0.1),
    MaxEtInEllipse = cms.double(2.0),
    MatchingConeMetric = cms.string('DR'),
    TrackerSignalConeSizeFormula = cms.string('0.07'),
    MatchingConeSizeFormula = cms.string('0.1'),
    TrackerIsolConeSize_min = cms.double(0.0),
    MatchingConeSize_min = cms.double(0.0),
    ElectronPreIDProducer = cms.InputTag("elecpreid"),
    ChargedHadrCandLeadChargedHadrCand_tksmaxDZ = cms.double(0.2),
    TrackerIsolConeSize_max = cms.double(0.6),
    TrackerSignalConeSize_max = cms.double(0.07),
    HCALIsolConeMetric = cms.string('DR'),
    AddEllipseGammas = cms.bool(False),
    maximumForElectrionPreIDOutput = cms.double(-0.1),
    TrackerSignalConeSize_min = cms.double(0.0),
    ECALIsolConeSize_max = cms.double(0.6),
    HCALIsolConeSizeFormula = cms.string('0.50'),
    Track_IsolAnnulus_minNhits = cms.uint32(3),
    smearedPVsigmaZ = cms.double(0.005),
    AreaMetric_recoElements_maxabsEta = cms.double(2.5),
    HCALSignalConeMetric = cms.string('DR'),
    ElecPreIDLeadTkMatch_maxDR = cms.double(0.01),
    ChargedHadrCand_IsolAnnulus_minNhits = cms.uint32(0),
    PFTauTagInfoProducer = cms.InputTag("pfRecoTauTagInfoProducer"),
    ECALIsolConeMetric = cms.string('DR'),
    ECALIsolConeSizeFormula = cms.string('0.50'),
    UseChargedHadrCandLeadChargedHadrCand_tksDZconstraint = cms.bool(True),
    Algorithm = cms.string('ConeBased'),
    JetPtMin = cms.double(0.0),
    ECALSignalConeMetric = cms.string('DR'),
    EcalStripSumE_deltaPhiOverQ_maxValue = cms.double(0.5),
    HCALSignalConeSize_max = cms.double(0.6),
    ECALSignalConeSize_min = cms.double(0.0),
    EcalStripSumE_minClusEnergy = cms.double(0.1),
    EcalStripSumE_deltaEta = cms.double(0.03),
    TrackerIsolConeSizeFormula = cms.string('0.50'),
    LeadPFCand_minPt = cms.double(5.0),
    HCALSignalConeSize_min = cms.double(0.0),
    ECALSignalConeSize_max = cms.double(0.6),
    HCALSignalConeSizeFormula = cms.string('0.10'),
    TrackLeadTrack_maxDZ = cms.double(0.2),
    DataType = cms.string('AOD'),
    ECALIsolConeSize_min = cms.double(0.0),
    UseTrackLeadTrackDZconstraint = cms.bool(True),
    smearedPVsigmaY = cms.double(0.0015),
    HCALIsolConeSize_max = cms.double(0.6),
    smearedPVsigmaX = cms.double(0.0015),
    MatchingConeSize_max = cms.double(0.6),
    HCALIsolConeSize_min = cms.double(0.0)
)


process.gamIsoDepositEcalFromHits = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("photons"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        isolationVariable = cms.string('et'),
        spikeIdThreshold = cms.double(0.95),
        tryBoth = cms.bool(True),
        recHitFlagsToBeExcluded = cms.vint32(3, 4, 8, 9),
        ComponentName = cms.string('EgammaRecHitExtractor'),
        endcapEcalHits = cms.InputTag("reducedEcalRecHitsEE"),
        intStrip = cms.double(0.0),
        intRadius = cms.double(0.0),
        severityRecHitThreshold = cms.double(5.0),
        severityLevelCut = cms.int32(3),
        energyMin = cms.double(0.08),
        extRadius = cms.double(0.6),
        subtractSuperClusterEnergy = cms.bool(False),
        spikeIdString = cms.string('kSwissCross'),
        vetoClustered = cms.bool(False),
        detector = cms.string('Ecal'),
        etMin = cms.double(0.0),
        DepositLabel = cms.untracked.string(''),
        barrelEcalHits = cms.InputTag("reducedEcalRecHitsEB")
    )
)


process.gamIsoDepositHcalFromTowers = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("photons"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        caloTowers = cms.InputTag("towerMaker"),
        ComponentName = cms.string('EgammaTowerExtractor'),
        hcalDepth = cms.int32(-1),
        intRadius = cms.double(0.0),
        extRadius = cms.double(0.6),
        DepositLabel = cms.untracked.string(''),
        etMin = cms.double(-999.0)
    )
)


process.gamIsoDepositTk = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("photons"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(0.2),
        dzOption = cms.string('vz'),
        BeamSpotLabel = cms.InputTag("offlineBeamSpot"),
        ComponentName = cms.string('EgammaTrackExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(9999.0),
        Chi2Prob_Min = cms.double(-1.0),
        DR_Veto = cms.double(0.0),
        NHits_Min = cms.uint32(0),
        Chi2Ndof_Max = cms.double(1e+64),
        Pt_Min = cms.double(-1.0),
        DepositLabel = cms.untracked.string(''),
        BeamlineOption = cms.string('BeamSpotFromEvent'),
        inputTrackCollection = cms.InputTag("generalTracks")
    )
)


process.gamIsoFromDepsEcalFromHits = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        src = cms.InputTag("gamIsoDepositEcalFromHits"),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring('EcalBarrel:0.045', 
            'EcalBarrel:RectangularEtaPhiVeto(-0.02,0.02,-0.5,0.5)', 
            'EcalBarrel:AbsThresholdFromTransverse(0.080)', 
            'EcalEndcaps:0.070', 
            'EcalEndcaps:RectangularEtaPhiVeto(-0.02,0.02,-0.5,0.5)', 
            'EcalEndcaps:AbsThreshold(0.100)'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
    ))
)


process.gamIsoFromDepsHcalFromTowers = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        src = cms.InputTag("gamIsoDepositHcalFromTowers"),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring('0.15'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
    ))
)


process.gamIsoFromDepsTk = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        src = cms.InputTag("gamIsoDepositTk"),
        deltaR = cms.double(0.3),
        weight = cms.string('1'),
        vetos = cms.vstring('RectangularEtaPhiVeto(-0.015,0.015,-0.5,0.5)', 
            'Threshold(1.0)'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
    ))
)


process.genParticlesForJets = cms.EDProducer("InputGenJetsParticleSelector",
    src = cms.InputTag("genParticles"),
    ignoreParticleIDs = cms.vuint32(1000022, 1000012, 1000014, 1000016, 2000012, 
        2000014, 2000016, 1000039, 5100039, 4000012, 
        4000014, 4000016, 9900012, 9900014, 9900016, 
        39),
    partonicFinalState = cms.bool(False),
    excludeResonances = cms.bool(True),
    excludeFromResonancePids = cms.vuint32(12, 13, 14, 16),
    tausAsJets = cms.bool(False)
)


process.genParticlesForJetsNoMuNoNu = cms.EDProducer("InputGenJetsParticleSelector",
    src = cms.InputTag("genParticles"),
    ignoreParticleIDs = cms.vuint32(1000022, 1000012, 1000014, 1000016, 2000012, 
        2000014, 2000016, 1000039, 5100039, 4000012, 
        4000014, 4000016, 9900012, 9900014, 9900016, 
        39, 12, 13, 14, 16),
    partonicFinalState = cms.bool(False),
    excludeResonances = cms.bool(True),
    excludeFromResonancePids = cms.vuint32(12, 13, 14, 16),
    tausAsJets = cms.bool(False)
)


process.genParticlesForJetsNoNu = cms.EDProducer("InputGenJetsParticleSelector",
    src = cms.InputTag("genParticles"),
    ignoreParticleIDs = cms.vuint32(1000022, 1000012, 1000014, 1000016, 2000012, 
        2000014, 2000016, 1000039, 5100039, 4000012, 
        4000014, 4000016, 9900012, 9900014, 9900016, 
        39, 12, 14, 16),
    partonicFinalState = cms.bool(False),
    excludeResonances = cms.bool(True),
    excludeFromResonancePids = cms.vuint32(12, 13, 14, 16),
    tausAsJets = cms.bool(False)
)


process.ghostTrackBJetTags = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('ghostTrack'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfos"), cms.InputTag("ghostTrackVertexTagInfos"))
)


process.ghostTrackVertexTagInfos = cms.EDProducer("SecondaryVertexProducer",
    beamSpotTag = cms.InputTag("offlineBeamSpot"),
    vertexReco = cms.PSet(
        primcut = cms.double(2.0),
        seccut = cms.double(4.0),
        maxFitChi2 = cms.double(10.0),
        fitType = cms.string('RefitGhostTrackWithVertices'),
        mergeThreshold = cms.double(3.0),
        finder = cms.string('gtvr')
    ),
    vertexSelection = cms.PSet(
        sortCriterium = cms.string('dist3dError')
    ),
    constraint = cms.string('BeamSpot'),
    vertexCuts = cms.PSet(
        distSig3dMax = cms.double(99999.9),
        fracPV = cms.double(0.65),
        distVal2dMax = cms.double(2.5),
        useTrackWeights = cms.bool(True),
        maxDeltaRToJetAxis = cms.double(0.5),
        v0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        ),
        distSig2dMin = cms.double(3.0),
        multiplicityMin = cms.uint32(1),
        massMax = cms.double(6.5),
        distSig2dMax = cms.double(99999.9),
        distVal3dMax = cms.double(99999.9),
        minimumTrackWeight = cms.double(0.5),
        distVal3dMin = cms.double(-99999.9),
        distVal2dMin = cms.double(0.01),
        distSig3dMin = cms.double(-99999.9)
    ),
    trackIPTagInfos = cms.InputTag("impactParameterTagInfos"),
    minimumTrackWeight = cms.double(0.5),
    usePVError = cms.bool(True),
    trackSelection = cms.PSet(
        totalHitsMin = cms.uint32(8),
        jetDeltaRMax = cms.double(0.3),
        qualityClass = cms.string('highPurity'),
        pixelHitsMin = cms.uint32(2),
        sip3dSigMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        maxDistToAxis = cms.double(0.2),
        sip2dValMax = cms.double(99999.9),
        maxDecayLen = cms.double(99999.9),
        ptMin = cms.double(1.0),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        sip2dValMin = cms.double(-99999.9),
        normChi2Max = cms.double(99999.9)
    ),
    trackSort = cms.string('sip3dSig')
)


process.gk5GenJets = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    src = cms.InputTag("genParticlesForJets"),
    doAreaFastjet = cms.bool(False),
    doPVCorrection = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    doRhoFastjet = cms.bool(False),
    srcPVs = cms.InputTag(""),
    inputEtMin = cms.double(0.0),
    nSigmaPU = cms.double(1.0),
    radiusPU = cms.double(0.5),
    Rho_EtaMax = cms.double(4.5),
    GhostArea = cms.double(0.01),
    jetPtMin = cms.double(3.0),
    jetType = cms.string('GenJet'),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    jetAlgorithm = cms.string('GeneralizedKt'),
    rParam = cms.double(0.5)
)


process.gk5GenJetsNoMuNoNu = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('GeneralizedKt'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.5),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJetsNoMuNoNu"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.gk5GenJetsNoNu = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('GeneralizedKt'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.5),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJetsNoNu"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.gk7GenJets = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('GeneralizedKt'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.7),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJets"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.gk7GenJetsNoMuNoNu = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('GeneralizedKt'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.7),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJetsNoMuNoNu"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.gk7GenJetsNoNu = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('GeneralizedKt'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.7),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJetsNoNu"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.hiGenParticlesForJets = cms.EDProducer("InputGenJetsParticleSelector",
    src = cms.InputTag("hiGenParticles"),
    ignoreParticleIDs = cms.vuint32(1000022, 1000012, 1000014, 1000016, 2000012, 
        2000014, 2000016, 1000039, 5100039, 4000012, 
        4000014, 4000016, 9900012, 9900014, 9900016, 
        39),
    partonicFinalState = cms.bool(False),
    excludeResonances = cms.bool(True),
    excludeFromResonancePids = cms.vuint32(12, 13, 14, 16),
    tausAsJets = cms.bool(False)
)


process.hpsPFRecoTauProducer = cms.EDProducer("PFRecoTauProducer",
    doOneProngTwoStrips = cms.bool(True),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    minimumSignalCone = cms.double(0.05),
    leadPionThreshold = cms.double(1.0),
    gammaIsolationConeSize = cms.double(0.5),
    stripPtThreshold = cms.double(1.0),
    candOverlapCriterion = cms.string('Isolation'),
    stripEtaAssociationDistance = cms.double(0.05),
    ElectronPreIDProducer = cms.InputTag("elecpreid"),
    oneProngTwoStripsPi0MassWindow = cms.vdouble(0.05, 0.2),
    doThreeProng = cms.bool(True),
    doOneProngStrip = cms.bool(True),
    smearedPVsigmaZ = cms.double(0.005),
    PFTauTagInfoProducer = cms.InputTag("pfRecoTauTagInfoProducer"),
    oneProngStripMassWindow = cms.vdouble(0.3, 1.3),
    maximumSignalCone = cms.double(0.1),
    coneMetric = cms.string('DR'),
    emMergingAlgorithm = cms.string('StripBased'),
    chargeHadrIsolationConeSize = cms.double(0.5),
    Algorithm = cms.string('HPS'),
    JetPtMin = cms.double(0.0),
    doOneProng = cms.bool(True),
    useIsolationAnnulus = cms.bool(False),
    threeProngMassWindow = cms.vdouble(0.8, 1.5),
    tauPtThreshold = cms.double(15.0),
    stripCandidatesPdgIds = cms.vint32(22, 11),
    stripPhiAssociationDistance = cms.double(0.2),
    neutrHadrIsolationConeSize = cms.double(0.5),
    smearedPVsigmaY = cms.double(0.0015),
    smearedPVsigmaX = cms.double(0.0015),
    coneSizeFormula = cms.string('2.8/ET'),
    matchingCone = cms.double(0.1),
    oneProngTwoStripsMassWindow = cms.vdouble(0.4, 1.2)
)


process.hpsPFTauDiscriminationAgainstElectron = cms.EDProducer("PFRecoTauDiscriminationAgainstElectron",
    ApplyCut_ElectronPreID_2D = cms.bool(False),
    ElecPreID0_HOverPLead_minValue = cms.double(0.05),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    ApplyCut_ElectronPreID = cms.bool(False),
    ApplyCut_HcalTotOverPLead = cms.bool(False),
    EOverPLead_minValue = cms.double(0.8),
    ElecPreID1_EOverPLead_maxValue = cms.double(0.8),
    HcalMaxOverPLead_minValue = cms.double(0.1),
    ApplyCut_EmFraction = cms.bool(False),
    EmFraction_maxValue = cms.double(0.9),
    ApplyCut_PFElectronMVA = cms.bool(True),
    PFElectronMVA_maxValue = cms.double(-0.1),
    ApplyCut_HcalMaxOverPLead = cms.bool(False),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    Hcal3x3OverPLead_minValue = cms.double(0.1),
    ElecPreID1_HOverPLead_minValue = cms.double(0.15),
    ElecPreID0_EOverPLead_maxValue = cms.double(0.95),
    BremsRecoveryEOverPLead_minValue = cms.double(0.8),
    ApplyCut_EcalCrackCut = cms.bool(False),
    EOverPLead_maxValue = cms.double(1.8),
    HcalTotOverPLead_minValue = cms.double(0.1),
    ApplyCut_BremsRecoveryEOverPLead = cms.bool(False),
    ApplyCut_Hcal3x3OverPLead = cms.bool(False),
    ApplyCut_EOverPLead = cms.bool(False),
    BremsRecoveryEOverPLead_maxValue = cms.double(1.8)
)


process.hpsPFTauDiscriminationAgainstMuon = cms.EDProducer("PFRecoTauDiscriminationAgainstMuon",
    a = cms.double(0.5),
    c = cms.double(0.0),
    b = cms.double(0.5),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    discriminatorOption = cms.string('noSegMatch')
)


process.hpsPFTauDiscriminationByDecayModeFinding = cms.EDProducer("PFRecoTauDiscriminationByLeadingObjectPtCut",
    MinPtLeadingObject = cms.double(0.0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    UseOnlyChargedHadrons = cms.bool(True)
)


process.hpsPFTauDiscriminationByLooseIsolation = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.hpsPFTauDiscriminationByMediumIsolation = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.8),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.8),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.8),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.hpsPFTauDiscriminationByTightIsolation = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.hpsPFTauProducer = cms.EDProducer("PFRecoTauProducer",
    doOneProngTwoStrips = cms.bool(True),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    minimumSignalCone = cms.double(0.05),
    leadPionThreshold = cms.double(1.0),
    gammaIsolationConeSize = cms.double(0.5),
    stripPtThreshold = cms.double(1.0),
    candOverlapCriterion = cms.string('Isolation'),
    stripEtaAssociationDistance = cms.double(0.05),
    ElectronPreIDProducer = cms.InputTag("elecpreid"),
    oneProngTwoStripsPi0MassWindow = cms.vdouble(0.05, 0.2),
    doThreeProng = cms.bool(True),
    doOneProngStrip = cms.bool(True),
    smearedPVsigmaZ = cms.double(0.005),
    PFTauTagInfoProducer = cms.InputTag("pfRecoTauTagInfoProducer"),
    oneProngStripMassWindow = cms.vdouble(0.3, 1.3),
    maximumSignalCone = cms.double(0.1),
    coneMetric = cms.string('DR'),
    emMergingAlgorithm = cms.string('StripBased'),
    chargeHadrIsolationConeSize = cms.double(0.5),
    Algorithm = cms.string('HPS'),
    JetPtMin = cms.double(0.0),
    doOneProng = cms.bool(True),
    useIsolationAnnulus = cms.bool(False),
    threeProngMassWindow = cms.vdouble(0.8, 1.5),
    tauPtThreshold = cms.double(15.0),
    stripCandidatesPdgIds = cms.vint32(22, 11),
    stripPhiAssociationDistance = cms.double(0.2),
    neutrHadrIsolationConeSize = cms.double(0.5),
    smearedPVsigmaY = cms.double(0.0015),
    smearedPVsigmaX = cms.double(0.0015),
    coneSizeFormula = cms.string('2.8/ET'),
    matchingCone = cms.double(0.1),
    oneProngTwoStripsMassWindow = cms.vdouble(0.4, 1.2)
)


process.ic5PFJetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
    tracks = cms.InputTag("generalTracks"),
    coneSize = cms.double(0.5),
    jets = cms.InputTag("iterativeCone5PFJets")
)


process.impactParameterMVABJetTags = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('impactParameterMVAComputer'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfos"))
)


process.impactParameterTagInfos = cms.EDProducer("TrackIPProducer",
    maximumTransverseImpactParameter = cms.double(0.2),
    minimumNumberOfHits = cms.int32(8),
    minimumTransverseMomentum = cms.double(1.0),
    primaryVertex = cms.InputTag("offlinePrimaryVertices"),
    maximumLongitudinalImpactParameter = cms.double(17.0),
    computeProbabilities = cms.bool(True),
    ghostTrackPriorDeltaR = cms.double(0.03),
    jetTracks = cms.InputTag("ak5JetTracksAssociatorAtVertex"),
    jetDirectionUsingGhostTrack = cms.bool(False),
    minimumNumberOfPixelHits = cms.int32(2),
    jetDirectionUsingTracks = cms.bool(False),
    computeGhostTrack = cms.bool(True),
    useTrackQuality = cms.bool(False),
    maximumChiSquared = cms.double(5.0)
)


process.impactParameterTagInfosAK5JPT = cms.EDProducer("TrackIPProducer",
    maximumTransverseImpactParameter = cms.double(0.2),
    minimumNumberOfHits = cms.int32(8),
    minimumTransverseMomentum = cms.double(1.0),
    primaryVertex = cms.InputTag("offlinePrimaryVertices"),
    maximumLongitudinalImpactParameter = cms.double(17.0),
    computeGhostTrack = cms.bool(True),
    ghostTrackPriorDeltaR = cms.double(0.03),
    jetTracks = cms.InputTag("jetTracksAssociatorAtVertexAK5JPT"),
    jetDirectionUsingGhostTrack = cms.bool(False),
    minimumNumberOfPixelHits = cms.int32(2),
    jetDirectionUsingTracks = cms.bool(False),
    computeProbabilities = cms.bool(True),
    useTrackQuality = cms.bool(False),
    maximumChiSquared = cms.double(5.0)
)


process.impactParameterTagInfosAK5PF = cms.EDProducer("TrackIPProducer",
    maximumTransverseImpactParameter = cms.double(0.2),
    minimumNumberOfHits = cms.int32(8),
    minimumTransverseMomentum = cms.double(1.0),
    primaryVertex = cms.InputTag("offlinePrimaryVertices"),
    maximumLongitudinalImpactParameter = cms.double(17.0),
    computeGhostTrack = cms.bool(True),
    ghostTrackPriorDeltaR = cms.double(0.03),
    jetTracks = cms.InputTag("jetTracksAssociatorAtVertexAK5PF"),
    jetDirectionUsingGhostTrack = cms.bool(False),
    minimumNumberOfPixelHits = cms.int32(2),
    jetDirectionUsingTracks = cms.bool(False),
    computeProbabilities = cms.bool(True),
    useTrackQuality = cms.bool(False),
    maximumChiSquared = cms.double(5.0)
)


process.impactParameterTagInfosAOD = cms.EDProducer("TrackIPProducer",
    maximumTransverseImpactParameter = cms.double(0.2),
    minimumNumberOfHits = cms.int32(8),
    minimumTransverseMomentum = cms.double(1.0),
    primaryVertex = cms.InputTag("offlinePrimaryVertices"),
    maximumLongitudinalImpactParameter = cms.double(17.0),
    computeGhostTrack = cms.bool(True),
    ghostTrackPriorDeltaR = cms.double(0.03),
    jetTracks = cms.InputTag("jetTracksAssociatorAtVertex"),
    jetDirectionUsingGhostTrack = cms.bool(False),
    minimumNumberOfPixelHits = cms.int32(2),
    jetDirectionUsingTracks = cms.bool(False),
    computeProbabilities = cms.bool(True),
    useTrackQuality = cms.bool(False),
    maximumChiSquared = cms.double(5.0)
)


process.impactParameterTagInfosAODPF = cms.EDProducer("TrackIPProducer",
    maximumTransverseImpactParameter = cms.double(0.2),
    minimumNumberOfHits = cms.int32(8),
    minimumTransverseMomentum = cms.double(1.0),
    primaryVertex = cms.InputTag("offlinePrimaryVertices"),
    maximumLongitudinalImpactParameter = cms.double(17.0),
    computeGhostTrack = cms.bool(True),
    ghostTrackPriorDeltaR = cms.double(0.03),
    jetTracks = cms.InputTag("jetTracksAssociatorAtVertexPF"),
    jetDirectionUsingGhostTrack = cms.bool(False),
    minimumNumberOfPixelHits = cms.int32(2),
    jetDirectionUsingTracks = cms.bool(False),
    computeProbabilities = cms.bool(True),
    useTrackQuality = cms.bool(False),
    maximumChiSquared = cms.double(5.0)
)


process.isoDepElectronWithCharged = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("pfSelectedElectrons"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedHadrons"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
    )
)


process.isoDepElectronWithChargedPF = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("pfSelectedElectronsPF"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedHadronsPF"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
    )
)


process.isoDepElectronWithNeutral = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("pfSelectedElectrons"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllNeutralHadrons"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
    )
)


process.isoDepElectronWithNeutralPF = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("pfSelectedElectronsPF"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllNeutralHadronsPF"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
    )
)


process.isoDepElectronWithPhotons = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("pfSelectedElectrons"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllPhotons"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
    )
)


process.isoDepElectronWithPhotonsPF = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("pfSelectedElectronsPF"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllPhotonsPF"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
    )
)


process.isoDepMuonWithCharged = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("pfSelectedMuons"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedHadrons"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
    )
)


process.isoDepMuonWithChargedPF = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("pfSelectedMuonsPF"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedHadronsPF"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
    )
)


process.isoDepMuonWithNeutral = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("pfSelectedMuons"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllNeutralHadrons"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
    )
)


process.isoDepMuonWithNeutralPF = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("pfSelectedMuonsPF"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllNeutralHadronsPF"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
    )
)


process.isoDepMuonWithPhotons = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("pfSelectedMuons"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllPhotons"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
    )
)


process.isoDepMuonWithPhotonsPF = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("pfSelectedMuonsPF"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllPhotonsPF"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
    )
)


process.isoDeposits = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag(""),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag(""),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
    )
)


process.isoValElectronWithCharged = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        src = cms.InputTag("isoDepElectronWithCharged"),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring(),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
    ))
)


process.isoValElectronWithChargedPF = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        src = cms.InputTag("isoDepElectronWithChargedPF"),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring(),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
    ))
)


process.isoValElectronWithNeutral = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        src = cms.InputTag("isoDepElectronWithNeutral"),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring('Threshold(0.5)'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
    ))
)


process.isoValElectronWithNeutralPF = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        src = cms.InputTag("isoDepElectronWithNeutralPF"),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring('Threshold(0.5)'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
    ))
)


process.isoValElectronWithPhotons = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        src = cms.InputTag("isoDepElectronWithPhotons"),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring('Threshold(0.5)'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
    ))
)


process.isoValElectronWithPhotonsPF = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        src = cms.InputTag("isoDepElectronWithPhotonsPF"),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring('Threshold(0.5)'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
    ))
)


process.isoValMuonWithCharged = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        src = cms.InputTag("isoDepMuonWithCharged"),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring(),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
    ))
)


process.isoValMuonWithChargedPF = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        src = cms.InputTag("isoDepMuonWithChargedPF"),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring(),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
    ))
)


process.isoValMuonWithNeutral = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        src = cms.InputTag("isoDepMuonWithNeutral"),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring('Threshold(0.5)'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
    ))
)


process.isoValMuonWithNeutralPF = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        src = cms.InputTag("isoDepMuonWithNeutralPF"),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring('Threshold(0.5)'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
    ))
)


process.isoValMuonWithPhotons = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        src = cms.InputTag("isoDepMuonWithPhotons"),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring('Threshold(0.5)'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
    ))
)


process.isoValMuonWithPhotonsPF = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        src = cms.InputTag("isoDepMuonWithPhotonsPF"),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring('Threshold(0.5)'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
    ))
)


process.iterativeCone3HiGenJets = cms.EDProducer("SubEventGenJetProducer",
    Active_Area_Repeats = cms.int32(5),
    src = cms.InputTag("hiGenParticlesForJets"),
    doAreaFastjet = cms.bool(True),
    doPVCorrection = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    doRhoFastjet = cms.bool(False),
    srcPVs = cms.InputTag(""),
    inputEtMin = cms.double(0.0),
    nSigmaPU = cms.double(1.0),
    radiusPU = cms.double(0.5),
    Rho_EtaMax = cms.double(4.5),
    GhostArea = cms.double(0.01),
    jetPtMin = cms.double(3.0),
    jetType = cms.string('GenJet'),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    jetAlgorithm = cms.string('IterativeCone'),
    rParam = cms.double(0.3)
)


process.iterativeCone4HiGenJets = cms.EDProducer("SubEventGenJetProducer",
    Active_Area_Repeats = cms.int32(5),
    src = cms.InputTag("hiGenParticlesForJets"),
    doAreaFastjet = cms.bool(True),
    doPVCorrection = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    doRhoFastjet = cms.bool(False),
    srcPVs = cms.InputTag(""),
    inputEtMin = cms.double(0.0),
    nSigmaPU = cms.double(1.0),
    radiusPU = cms.double(0.5),
    Rho_EtaMax = cms.double(4.5),
    GhostArea = cms.double(0.01),
    jetPtMin = cms.double(3.0),
    jetType = cms.string('GenJet'),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    jetAlgorithm = cms.string('IterativeCone'),
    rParam = cms.double(0.4)
)


process.iterativeCone5GenJets = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    src = cms.InputTag("genParticlesForJets"),
    doAreaFastjet = cms.bool(False),
    doPVCorrection = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    doRhoFastjet = cms.bool(False),
    srcPVs = cms.InputTag(""),
    inputEtMin = cms.double(0.0),
    nSigmaPU = cms.double(1.0),
    radiusPU = cms.double(0.5),
    Rho_EtaMax = cms.double(4.5),
    GhostArea = cms.double(0.01),
    jetPtMin = cms.double(3.0),
    jetType = cms.string('GenJet'),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    jetAlgorithm = cms.string('IterativeCone'),
    rParam = cms.double(0.5)
)


process.iterativeCone5GenJetsNoMuNoNu = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('IterativeCone'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.5),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJetsNoMuNoNu"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.iterativeCone5GenJetsNoNu = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('IterativeCone'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.5),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJetsNoNu"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.iterativeCone5HiGenJets = cms.EDProducer("SubEventGenJetProducer",
    Active_Area_Repeats = cms.int32(5),
    src = cms.InputTag("hiGenParticlesForJets"),
    doAreaFastjet = cms.bool(True),
    doPVCorrection = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    doRhoFastjet = cms.bool(False),
    srcPVs = cms.InputTag(""),
    inputEtMin = cms.double(0.0),
    nSigmaPU = cms.double(1.0),
    radiusPU = cms.double(0.5),
    Rho_EtaMax = cms.double(4.5),
    GhostArea = cms.double(0.01),
    jetPtMin = cms.double(3.0),
    jetType = cms.string('GenJet'),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    jetAlgorithm = cms.string('IterativeCone'),
    rParam = cms.double(0.5)
)


process.iterativeCone7HiGenJets = cms.EDProducer("SubEventGenJetProducer",
    Active_Area_Repeats = cms.int32(5),
    src = cms.InputTag("hiGenParticlesForJets"),
    doAreaFastjet = cms.bool(True),
    doPVCorrection = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    doRhoFastjet = cms.bool(False),
    srcPVs = cms.InputTag(""),
    inputEtMin = cms.double(0.0),
    nSigmaPU = cms.double(1.0),
    radiusPU = cms.double(0.5),
    Rho_EtaMax = cms.double(4.5),
    GhostArea = cms.double(0.01),
    jetPtMin = cms.double(3.0),
    jetType = cms.string('GenJet'),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    jetAlgorithm = cms.string('IterativeCone'),
    rParam = cms.double(0.7)
)


process.jetBProbabilityBJetTags = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('jetBProbability'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfos"))
)


process.jetBProbabilityBJetTagsAK5JPT = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('jetBProbability'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAK5JPT"))
)


process.jetBProbabilityBJetTagsAK5PF = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('jetBProbability'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAK5PF"))
)


process.jetBProbabilityBJetTagsAOD = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('jetBProbability'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAOD"))
)


process.jetBProbabilityBJetTagsAODPF = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('jetBProbability'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAODPF"))
)


process.jetProbabilityBJetTags = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('jetProbability'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfos"))
)


process.jetProbabilityBJetTagsAK5JPT = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('jetProbability'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAK5JPT"))
)


process.jetProbabilityBJetTagsAK5PF = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('jetProbability'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAK5PF"))
)


process.jetProbabilityBJetTagsAOD = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('jetProbability'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAOD"))
)


process.jetProbabilityBJetTagsAODPF = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('jetProbability'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAODPF"))
)


process.jetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
    jets = cms.InputTag("ak5CaloJets"),
    tracks = cms.InputTag("generalTracks"),
    coneSize = cms.double(0.5)
)


process.jetTracksAssociatorAtVertexAK5JPT = cms.EDProducer("JetTracksAssociatorAtVertex",
    jets = cms.InputTag("JetPlusTrackZSPCorJetAntiKt5"),
    tracks = cms.InputTag("generalTracks"),
    coneSize = cms.double(0.5)
)


process.jetTracksAssociatorAtVertexAK5PF = cms.EDProducer("JetTracksAssociatorAtVertex",
    jets = cms.InputTag("ak5PFJets"),
    tracks = cms.InputTag("generalTracks"),
    coneSize = cms.double(0.5)
)


process.jetTracksAssociatorAtVertexPF = cms.EDProducer("JetTracksAssociatorAtVertex",
    jets = cms.InputTag("pfNoTauPF"),
    tracks = cms.InputTag("generalTracks"),
    coneSize = cms.double(0.5)
)


process.kt3HiGenJets = cms.EDProducer("SubEventGenJetProducer",
    Active_Area_Repeats = cms.int32(5),
    src = cms.InputTag("hiGenParticlesForJets"),
    doAreaFastjet = cms.bool(True),
    doPVCorrection = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    doRhoFastjet = cms.bool(False),
    srcPVs = cms.InputTag(""),
    inputEtMin = cms.double(0.0),
    nSigmaPU = cms.double(1.0),
    radiusPU = cms.double(0.5),
    Rho_EtaMax = cms.double(4.5),
    GhostArea = cms.double(0.01),
    jetPtMin = cms.double(3.0),
    jetType = cms.string('GenJet'),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    jetAlgorithm = cms.string('Kt'),
    rParam = cms.double(0.3)
)


process.kt4GenJets = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    src = cms.InputTag("genParticlesForJets"),
    doAreaFastjet = cms.bool(False),
    doPVCorrection = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    doRhoFastjet = cms.bool(False),
    srcPVs = cms.InputTag(""),
    inputEtMin = cms.double(0.0),
    nSigmaPU = cms.double(1.0),
    radiusPU = cms.double(0.5),
    Rho_EtaMax = cms.double(4.5),
    GhostArea = cms.double(0.01),
    jetPtMin = cms.double(3.0),
    jetType = cms.string('GenJet'),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    jetAlgorithm = cms.string('Kt'),
    rParam = cms.double(0.4)
)


process.kt4GenJetsNoMuNoNu = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('Kt'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.4),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJetsNoMuNoNu"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.kt4GenJetsNoNu = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('Kt'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.4),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJetsNoNu"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.kt4HiGenJets = cms.EDProducer("SubEventGenJetProducer",
    Active_Area_Repeats = cms.int32(5),
    src = cms.InputTag("hiGenParticlesForJets"),
    doAreaFastjet = cms.bool(True),
    doPVCorrection = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    doRhoFastjet = cms.bool(False),
    srcPVs = cms.InputTag(""),
    inputEtMin = cms.double(0.0),
    nSigmaPU = cms.double(1.0),
    radiusPU = cms.double(0.5),
    Rho_EtaMax = cms.double(4.5),
    GhostArea = cms.double(0.01),
    jetPtMin = cms.double(3.0),
    jetType = cms.string('GenJet'),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    jetAlgorithm = cms.string('Kt'),
    rParam = cms.double(0.4)
)


process.kt6GenJets = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('Kt'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.6),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJets"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.kt6GenJetsNoMuNoNu = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('Kt'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.6),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJetsNoMuNoNu"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.kt6GenJetsNoNu = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('Kt'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.6),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJetsNoNu"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.kt6HiGenJets = cms.EDProducer("SubEventGenJetProducer",
    Active_Area_Repeats = cms.int32(5),
    src = cms.InputTag("hiGenParticlesForJets"),
    doAreaFastjet = cms.bool(True),
    doPVCorrection = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    doRhoFastjet = cms.bool(False),
    srcPVs = cms.InputTag(""),
    inputEtMin = cms.double(0.0),
    nSigmaPU = cms.double(1.0),
    radiusPU = cms.double(0.5),
    Rho_EtaMax = cms.double(4.5),
    GhostArea = cms.double(0.01),
    jetPtMin = cms.double(3.0),
    jetType = cms.string('GenJet'),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    jetAlgorithm = cms.string('Kt'),
    rParam = cms.double(0.6)
)


process.metJESCorAK5CaloJet = cms.EDProducer("Type1MET",
    inputUncorJetsLabel = cms.string('ak5CaloJets'),
    jetEMfracLimit = cms.double(0.9),
    metType = cms.string('CaloMET'),
    jetPTthreshold = cms.double(20.0),
    inputUncorMetLabel = cms.string('met'),
    hasMuonsCorr = cms.bool(False),
    useTypeII = cms.bool(False),
    corrector = cms.string('ak5CaloL2L3'),
    UscaleA = cms.double(1.5),
    UscaleB = cms.double(1.8),
    UscaleC = cms.double(-0.06)
)


process.metJESCorAK5CaloJetMuons = cms.EDProducer("MuonMET",
    muonMETDepositValueMapInputTag = cms.InputTag("muonMETValueMapProducer","muCorrData"),
    metTypeInputTag = cms.InputTag("CaloMET"),
    muonsInputTag = cms.InputTag("muons"),
    uncorMETInputTag = cms.InputTag("metJESCorAK5CaloJet")
)


process.metJESCorAK5CaloJetMuonsPF = cms.EDProducer("MuonMET",
    muonMETDepositValueMapInputTag = cms.InputTag("muonMETValueMapProducer","muCorrData"),
    metTypeInputTag = cms.InputTag("CaloMET"),
    muonsInputTag = cms.InputTag("muons"),
    uncorMETInputTag = cms.InputTag("metJESCorAK5CaloJetPF")
)


process.metJESCorAK5CaloJetMuonsTypeII = cms.EDProducer("MuonMET",
    muonMETDepositValueMapInputTag = cms.InputTag("muonMETValueMapProducer","muCorrData"),
    metTypeInputTag = cms.InputTag("CaloMET"),
    muonsInputTag = cms.InputTag("muons"),
    uncorMETInputTag = cms.InputTag("metJESCorAK5CaloJetTypeII")
)


process.metJESCorAK5CaloJetPF = cms.EDProducer("Type1MET",
    inputUncorJetsLabel = cms.string('ak5CaloJets'),
    jetEMfracLimit = cms.double(0.9),
    metType = cms.string('CaloMET'),
    jetPTthreshold = cms.double(20.0),
    inputUncorMetLabel = cms.string('met'),
    hasMuonsCorr = cms.bool(False),
    useTypeII = cms.bool(False),
    corrector = cms.string('ak5CaloL2L3'),
    UscaleA = cms.double(1.5),
    UscaleB = cms.double(1.8),
    UscaleC = cms.double(-0.06)
)


process.metJESCorAK5CaloJetTypeII = cms.EDProducer("Type1MET",
    inputUncorJetsLabel = cms.string('ak5CaloJets'),
    jetEMfracLimit = cms.double(0.9),
    metType = cms.string('CaloMET'),
    jetPTthreshold = cms.double(20.0),
    inputUncorMetLabel = cms.string('met'),
    hasMuonsCorr = cms.bool(False),
    useTypeII = cms.bool(True),
    corrector = cms.string('ak5CaloL2L3'),
    UscaleA = cms.double(1.5),
    UscaleB = cms.double(1.8),
    UscaleC = cms.double(-0.06)
)


process.metJESCorAK5PFJet = cms.EDProducer("Type1MET",
    inputUncorJetsLabel = cms.string('ak5PFJets'),
    jetEMfracLimit = cms.double(0.9),
    metType = cms.string('PFMET'),
    jetPTthreshold = cms.double(10.0),
    inputUncorMetLabel = cms.string('pfMet'),
    hasMuonsCorr = cms.bool(False),
    useTypeII = cms.bool(False),
    corrector = cms.string('ak5PFL2L3'),
    UscaleA = cms.double(1.5),
    UscaleB = cms.double(1.8),
    UscaleC = cms.double(-0.06)
)


process.metJESCorAK5PFTypeI = cms.EDProducer("Type1MET",
    inputUncorJetsLabel = cms.string('patJetsPF'),
    jetEMfracLimit = cms.double(0.9),
    metType = cms.string('pat'),
    jetPTthreshold = cms.double(6.0),
    inputUncorMetLabel = cms.string('pfMet'),
    hasMuonsCorr = cms.bool(False),
    useTypeII = cms.bool(False),
    corrector = cms.string('ak5PFL2L3'),
    UscaleA = cms.double(1.5),
    UscaleB = cms.double(1.8),
    UscaleC = cms.double(-0.06)
)


process.metJESCorAK7CaloJet = cms.EDProducer("Type1MET",
    inputUncorJetsLabel = cms.string('ak7CaloJets'),
    jetEMfracLimit = cms.double(0.9),
    metType = cms.string('CaloMET'),
    jetPTthreshold = cms.double(20.0),
    inputUncorMetLabel = cms.string('met'),
    hasMuonsCorr = cms.bool(False),
    useTypeII = cms.bool(False),
    corrector = cms.string('ak7CaloL2L3'),
    UscaleA = cms.double(1.5),
    UscaleB = cms.double(1.8),
    UscaleC = cms.double(-0.06)
)


process.metJESCorIC5CaloJet = cms.EDProducer("Type1MET",
    inputUncorJetsLabel = cms.string('iterativeCone5CaloJets'),
    jetEMfracLimit = cms.double(0.9),
    metType = cms.string('CaloMET'),
    jetPTthreshold = cms.double(20.0),
    inputUncorMetLabel = cms.string('met'),
    hasMuonsCorr = cms.bool(False),
    useTypeII = cms.bool(False),
    corrector = cms.string('ic5CaloL2L3'),
    UscaleA = cms.double(1.5),
    UscaleB = cms.double(1.8),
    UscaleC = cms.double(-0.06)
)


process.metJESCorKT4CaloJet = cms.EDProducer("Type1MET",
    inputUncorJetsLabel = cms.string('kt4CaloJets'),
    jetEMfracLimit = cms.double(0.9),
    metType = cms.string('CaloMET'),
    jetPTthreshold = cms.double(20.0),
    inputUncorMetLabel = cms.string('met'),
    hasMuonsCorr = cms.bool(False),
    useTypeII = cms.bool(False),
    corrector = cms.string('kt4CaloL2L3'),
    UscaleA = cms.double(1.5),
    UscaleB = cms.double(1.8),
    UscaleC = cms.double(-0.06)
)


process.metJESCorKT6CaloJet = cms.EDProducer("Type1MET",
    inputUncorJetsLabel = cms.string('kt6CaloJets'),
    jetEMfracLimit = cms.double(0.9),
    metType = cms.string('CaloMET'),
    jetPTthreshold = cms.double(20.0),
    inputUncorMetLabel = cms.string('met'),
    hasMuonsCorr = cms.bool(False),
    useTypeII = cms.bool(False),
    corrector = cms.string('kt6CaloL2L3'),
    UscaleA = cms.double(1.5),
    UscaleB = cms.double(1.8),
    UscaleC = cms.double(-0.06)
)


process.metJESCorSC5CaloJet = cms.EDProducer("Type1MET",
    inputUncorJetsLabel = cms.string('sisCone5CaloJets'),
    jetEMfracLimit = cms.double(0.9),
    metType = cms.string('CaloMET'),
    jetPTthreshold = cms.double(20.0),
    inputUncorMetLabel = cms.string('met'),
    hasMuonsCorr = cms.bool(False),
    useTypeII = cms.bool(False),
    corrector = cms.string('sisCone5CaloL2L3'),
    UscaleA = cms.double(1.5),
    UscaleB = cms.double(1.8),
    UscaleC = cms.double(-0.06)
)


process.metJESCorSC7CaloJet = cms.EDProducer("Type1MET",
    inputUncorJetsLabel = cms.string('sisCone7CaloJets'),
    jetEMfracLimit = cms.double(0.9),
    metType = cms.string('CaloMET'),
    jetPTthreshold = cms.double(20.0),
    inputUncorMetLabel = cms.string('met'),
    hasMuonsCorr = cms.bool(False),
    useTypeII = cms.bool(False),
    corrector = cms.string('sisCone7CaloL2L3'),
    UscaleA = cms.double(1.5),
    UscaleB = cms.double(1.8),
    UscaleC = cms.double(-0.06)
)


process.muonMETValueMapProducer = cms.EDProducer("MuonMETValueMapProducer",
    TrackAssociatorParameters = cms.PSet(
        muonMaxDistanceSigmaX = cms.double(0.0),
        muonMaxDistanceSigmaY = cms.double(0.0),
        CSCSegmentCollectionLabel = cms.InputTag("cscSegments"),
        dRHcal = cms.double(9999.0),
        dRPreshowerPreselection = cms.double(0.2),
        CaloTowerCollectionLabel = cms.InputTag("towerMaker"),
        useEcal = cms.bool(False),
        dREcal = cms.double(9999.0),
        dREcalPreselection = cms.double(0.05),
        HORecHitCollectionLabel = cms.InputTag("horeco"),
        dRMuon = cms.double(9999.0),
        propagateAllDirections = cms.bool(True),
        muonMaxDistanceX = cms.double(5.0),
        muonMaxDistanceY = cms.double(5.0),
        useHO = cms.bool(False),
        trajectoryUncertaintyTolerance = cms.double(-1.0),
        usePreshower = cms.bool(False),
        DTRecSegment4DCollectionLabel = cms.InputTag("dt4DSegments"),
        EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
        dRHcalPreselection = cms.double(0.2),
        useMuon = cms.bool(False),
        useCalo = cms.bool(True),
        accountForTrajectoryChangeCalo = cms.bool(False),
        EBRecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
        dRMuonPreselection = cms.double(0.2),
        truthMatch = cms.bool(False),
        HBHERecHitCollectionLabel = cms.InputTag("hbhereco"),
        useHcal = cms.bool(False)
    ),
    beamSpotInputTag = cms.InputTag("offlineBeamSpot"),
    minPt = cms.double(10.0),
    maxNormChi2 = cms.double(10.0),
    minnValidStaHits = cms.int32(1),
    useHO = cms.bool(False),
    minnHits = cms.int32(11),
    useTrackAssociatorPositions = cms.bool(True),
    useRecHits = cms.bool(False),
    maxEta = cms.double(2.5),
    maxd0 = cms.double(0.2),
    towerEtThreshold = cms.double(0.3),
    isAlsoTkMu = cms.bool(True),
    muonInputTag = cms.InputTag("muons")
)


process.muonMatch = cms.EDProducer("MCMatcher",
    src = cms.InputTag("muons"),
    maxDPtRel = cms.double(999999.0),
    mcPdgId = cms.vint32(13),
    mcStatus = cms.vint32(1),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = cms.double(0.2),
    checkCharge = cms.bool(False),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("genParticles")
)


process.muonMatchPF = cms.EDProducer("MCMatcher",
    src = cms.InputTag("pfIsolatedMuonsPF"),
    maxDPtRel = cms.double(999999.0),
    mcPdgId = cms.vint32(13),
    mcStatus = cms.vint32(1),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = cms.double(0.2),
    checkCharge = cms.bool(False),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("genParticles")
)


process.patElectronMatch = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
    matchedCuts = cms.string('path("HLT_Ele15_LW_L1R") || path("HLT_Ele15_SW_L1R") || path("HLT_Ele15_SW_CaloEleId_L1R") || path("HLT_Ele15_SW_EleId_L1R") || path("HLT_Ele17_SW_TightEleId_L1R") || path("HLT_Ele17_SW_TighterEleId_L1R_v1") || path("HLT_Ele22_SW_TighterCaloIdIsol_L1R_v1") || path("HLT_Ele22_SW_TighterCaloIdIsol_L1R_v2") || path("HLT_Ele22_SW_TighterEleId_L1R_v1") || path("HLT_Ele22_SW_TighterEleId_L1R_v2") || path("HLT_Ele22_SW_TighterEleId_L1R_v3")'),
    src = cms.InputTag("cleanPatElectrons"),
    maxDPtRel = cms.double(0.5),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = cms.double(0.5),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("patTrigger")
)


process.patElectronMatchPF = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
    matchedCuts = cms.string('path("HLT_Ele15_LW_L1R") || path("HLT_Ele15_SW_L1R") || path("HLT_Ele15_SW_CaloEleId_L1R") || path("HLT_Ele15_SW_EleId_L1R") || path("HLT_Ele17_SW_TightEleId_L1R") || path("HLT_Ele17_SW_TighterEleId_L1R_v1") || path("HLT_Ele22_SW_TighterCaloIdIsol_L1R_v1") || path("HLT_Ele22_SW_TighterCaloIdIsol_L1R_v2") || path("HLT_Ele22_SW_TighterEleId_L1R_v1") || path("HLT_Ele22_SW_TighterEleId_L1R_v2") || path("HLT_Ele22_SW_TighterEleId_L1R_v3")'),
    src = cms.InputTag("selectedPatElectronsPF"),
    maxDPtRel = cms.double(0.5),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = cms.double(0.5),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("patTrigger")
)


process.patElectrons = cms.EDProducer("PATElectronProducer",
    embedHighLevelSelection = cms.bool(True),
    electronSource = cms.InputTag("gsfElectrons"),
    resolutions = cms.PSet(

    ),
    userIsolation = cms.PSet(

    ),
    embedSuperCluster = cms.bool(False),
    embedPFCandidate = cms.bool(True),
    pfElectronSource = cms.InputTag("pfIsolatedElectrons"),
    addElectronID = cms.bool(True),
    efficiencies = cms.PSet(

    ),
    embedGsfTrack = cms.bool(True),
    useParticleFlow = cms.bool(False),
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
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring()
    ),
    embedTrack = cms.bool(False),
    addEfficiencies = cms.bool(False),
    usePV = cms.bool(True),
    pvSrc = cms.InputTag("offlinePrimaryVertices"),
    electronIDSources = cms.PSet(
        simpleEleId80cIso = cms.InputTag("simpleEleId80cIso"),
        simpleEleId60cIso = cms.InputTag("simpleEleId60cIso"),
        eidRobustTight = cms.InputTag("eidRobustTight"),
        eidRobustHighEnergy = cms.InputTag("eidRobustHighEnergy"),
        simpleEleId85cIso = cms.InputTag("simpleEleId85cIso"),
        simpleEleId95cIso = cms.InputTag("simpleEleId95cIso"),
        simpleEleId90relIso = cms.InputTag("simpleEleId90relIso"),
        simpleEleId70relIso = cms.InputTag("simpleEleId70relIso"),
        simpleEleId60relIso = cms.InputTag("simpleEleId60relIso"),
        simpleEleId80relIso = cms.InputTag("simpleEleId80relIso"),
        simpleEleId70cIso = cms.InputTag("simpleEleId70cIso"),
        simpleEleId95relIso = cms.InputTag("simpleEleId95relIso"),
        eidTight = cms.InputTag("eidTight"),
        eidLoose = cms.InputTag("eidLoose"),
        eidRobustLoose = cms.InputTag("eidRobustLoose"),
        simpleEleId85relIso = cms.InputTag("simpleEleId85relIso"),
        simpleEleId90cIso = cms.InputTag("simpleEleId90cIso")
    ),
    genParticleMatch = cms.InputTag(""),
    beamLineSrc = cms.InputTag("offlineBeamSpot"),
    addGenMatch = cms.bool(False),
    addResolutions = cms.bool(False),
    isoDeposits = cms.PSet(

    ),
    embedGenMatch = cms.bool(False)
)


process.patElectronsPF = cms.EDProducer("PATElectronProducer",
    embedHighLevelSelection = cms.bool(True),
    electronSource = cms.InputTag("gsfElectrons"),
    resolutions = cms.PSet(

    ),
    userIsolation = cms.PSet(

    ),
    embedSuperCluster = cms.bool(True),
    embedPFCandidate = cms.bool(True),
    pfElectronSource = cms.InputTag("pfIsolatedElectronsPF"),
    addElectronID = cms.bool(True),
    efficiencies = cms.PSet(

    ),
    embedGsfTrack = cms.bool(True),
    useParticleFlow = cms.bool(True),
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
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring()
    ),
    embedTrack = cms.bool(False),
    addEfficiencies = cms.bool(False),
    usePV = cms.bool(True),
    pvSrc = cms.InputTag("offlinePrimaryVertices"),
    addResolutions = cms.bool(False),
    genParticleMatch = cms.InputTag(""),
    beamLineSrc = cms.InputTag("offlineBeamSpot"),
    addGenMatch = cms.bool(False),
    electronIDSources = cms.PSet(
        eidTight = cms.InputTag("eidTight"),
        eidLoose = cms.InputTag("eidLoose"),
        eidRobustTight = cms.InputTag("eidRobustTight"),
        eidRobustHighEnergy = cms.InputTag("eidRobustHighEnergy"),
        eidRobustLoose = cms.InputTag("eidRobustLoose")
    ),
    isoDeposits = cms.PSet(
        pfNeutralHadrons = cms.InputTag("isoDepElectronWithNeutralPF"),
        pfChargedHadrons = cms.InputTag("isoDepElectronWithChargedPF"),
        pfPhotons = cms.InputTag("isoDepElectronWithPhotonsPF")
    ),
    embedGenMatch = cms.bool(False),
    isolationValues = cms.PSet(
        pfNeutralHadrons = cms.InputTag("isoValElectronWithNeutralPF"),
        pfChargedHadrons = cms.InputTag("isoValElectronWithChargedPF"),
        pfPhotons = cms.InputTag("isoValElectronWithPhotonsPF")
    )
)


process.patHemispheres = cms.EDProducer("PATHemisphereProducer",
    patJets = cms.InputTag("cleanPatJetsAK5Calo"),
    maxTauEta = cms.double(-1),
    maxPhotonEta = cms.double(5),
    minMuonEt = cms.double(7),
    patMuons = cms.InputTag("cleanLayer1Muons"),
    seedMethod = cms.int32(3),
    patElectrons = cms.InputTag("cleanLayer1Electrons"),
    patMets = cms.InputTag("layer1METs"),
    maxMuonEta = cms.double(5),
    minTauEt = cms.double(1000000),
    minPhotonEt = cms.double(200000),
    minElectronEt = cms.double(7),
    patPhotons = cms.InputTag("cleanLayer1Photons"),
    combinationMethod = cms.int32(3),
    maxJetEta = cms.double(5),
    maxElectronEta = cms.double(5),
    minJetEt = cms.double(30),
    patTaus = cms.InputTag("cleanLayer1Taus")
)


process.patJetCharge = cms.EDProducer("JetChargeProducer",
    var = cms.string('Pt'),
    src = cms.InputTag("jetTracksAssociatorAtVertex"),
    exp = cms.double(1.0)
)


process.patJetChargeAK5JPT = cms.EDProducer("JetChargeProducer",
    var = cms.string('Pt'),
    src = cms.InputTag("jetTracksAssociatorAtVertexAK5JPT"),
    exp = cms.double(1.0)
)


process.patJetChargeAK5PF = cms.EDProducer("JetChargeProducer",
    var = cms.string('Pt'),
    src = cms.InputTag("jetTracksAssociatorAtVertexAK5PF"),
    exp = cms.double(1.0)
)


process.patJetChargePF = cms.EDProducer("JetChargeProducer",
    var = cms.string('Pt'),
    src = cms.InputTag("jetTracksAssociatorAtVertexPF"),
    exp = cms.double(1.0)
)


process.patJetCorrFactors = cms.EDProducer("JetCorrFactorsProducer",
    src = cms.InputTag("ak5CaloJets"),
    emf = cms.bool(False),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    levels = cms.vstring('L1Offset', 
        'L2Relative', 
        'L3Absolute',
        'L2L3Residual'
    ),
    payload = cms.string('AK5Calo'),
    flavorType = cms.string('J')
)


process.patJetCorrFactorsAK5JPT = cms.EDProducer("JetCorrFactorsProducer",
    src = cms.InputTag("JetPlusTrackZSPCorJetAntiKt5"),
    emf = cms.bool(False),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    levels = cms.vstring('L1Offset', 
        'L2Relative', 
        'L3Absolute',
        'L2L3Residual'
    ),
    payload = cms.string('AK5JPT'),
    flavorType = cms.string('J')
)


process.patJetCorrFactorsAK5PF = cms.EDProducer("JetCorrFactorsProducer",
    src = cms.InputTag("ak5PFJets"),
    emf = cms.bool(False),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    levels = cms.vstring('L1Offset', 
        'L2Relative', 
        'L3Absolute',
        'L2L3Residual'
    ),
    payload = cms.string('AK5PF'),
    flavorType = cms.string('J')
)


process.patJetCorrFactorsPF = cms.EDProducer("JetCorrFactorsProducer",
    src = cms.InputTag("pfNoTauPF"),
    emf = cms.bool(False),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    levels = cms.vstring('L1Offset',
        'L2Relative', 
        'L3Absolute', 
        'L2L3Residual'
    ),
    payload = cms.string('AK5PF'),
    flavorType = cms.string('J')
)


process.patJetFlavourAssociation = cms.EDProducer("JetFlavourIdentifier",
    srcByReference = cms.InputTag("patJetPartonAssociation"),
    physicsDefinition = cms.bool(False)
)


process.patJetFlavourAssociationAK5JPT = cms.EDProducer("JetFlavourIdentifier",
    srcByReference = cms.InputTag("patJetPartonAssociationAK5JPT"),
    physicsDefinition = cms.bool(False)
)


process.patJetFlavourAssociationAK5PF = cms.EDProducer("JetFlavourIdentifier",
    srcByReference = cms.InputTag("patJetPartonAssociationAK5PF"),
    physicsDefinition = cms.bool(False)
)


process.patJetFlavourAssociationPF = cms.EDProducer("JetFlavourIdentifier",
    srcByReference = cms.InputTag("patJetPartonAssociationPF"),
    physicsDefinition = cms.bool(False)
)


process.patJetGenJetMatch = cms.EDProducer("GenJetMatcher",
    src = cms.InputTag("ak5CaloJets"),
    maxDPtRel = cms.double(999999.0),
    mcPdgId = cms.vint32(),
    mcStatus = cms.vint32(),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = cms.double(0.25),
    checkCharge = cms.bool(False),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("ak5GenJets")
)


process.patJetGenJetMatchAK5JPT = cms.EDProducer("GenJetMatcher",
    src = cms.InputTag("JetPlusTrackZSPCorJetAntiKt5"),
    maxDPtRel = cms.double(999999.0),
    mcPdgId = cms.vint32(),
    mcStatus = cms.vint32(),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = cms.double(0.25),
    checkCharge = cms.bool(False),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("ak5GenJets")
)


process.patJetGenJetMatchAK5PF = cms.EDProducer("GenJetMatcher",
    src = cms.InputTag("ak5PFJets"),
    maxDPtRel = cms.double(999999.0),
    mcPdgId = cms.vint32(),
    mcStatus = cms.vint32(),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = cms.double(0.25),
    checkCharge = cms.bool(False),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("ak5GenJets")
)


process.patJetGenJetMatchPF = cms.EDProducer("GenJetMatcher",
    src = cms.InputTag("pfNoTauPF"),
    maxDPtRel = cms.double(3.0),
    mcPdgId = cms.vint32(),
    mcStatus = cms.vint32(),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = cms.double(0.4),
    checkCharge = cms.bool(False),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("ak5GenJetsNoNu")
)


process.patJetMatchAK5Calo = cms.EDProducer("PATTriggerMatcherDRLessByR",
    matchedCuts = cms.string('path("HLT_Jet15U") || path("HLT_Jet30U") || path("HLT_Jet50U")'),
    src = cms.InputTag("cleanPatJetsAK5Calo"),
    maxDPtRel = cms.double(3.0),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = cms.double(0.4),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("patTrigger")
)


process.patJetMatchAK5JPT = cms.EDProducer("PATTriggerMatcherDRLessByR",
    matchedCuts = cms.string('path("HLT_Jet15U") || path("HLT_Jet30U") || path("HLT_Jet50U")'),
    src = cms.InputTag("cleanPatJetsAK5JPT"),
    maxDPtRel = cms.double(3.0),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = cms.double(0.4),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("patTrigger")
)


process.patJetMatchAK5PF = cms.EDProducer("PATTriggerMatcherDRLessByR",
    matchedCuts = cms.string('path("HLT_Jet15U") || path("HLT_Jet30U") || path("HLT_Jet50U")'),
    src = cms.InputTag("cleanPatJetsAK5PF"),
    maxDPtRel = cms.double(3.0),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = cms.double(0.4),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("patTrigger")
)


process.patJetMatchPF = cms.EDProducer("PATTriggerMatcherDRLessByR",
    matchedCuts = cms.string('path("HLT_Jet15U") || path("HLT_Jet30U") || path("HLT_Jet50U")'),
    src = cms.InputTag("selectedPatJetsPF"),
    maxDPtRel = cms.double(3.0),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = cms.double(0.4),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("patTrigger")
)


process.patJetPartonAssociation = cms.EDProducer("JetPartonMatcher",
    jets = cms.InputTag("ak5CaloJets"),
    coneSizeToAssociate = cms.double(0.3),
    partons = cms.InputTag("patJetPartons")
)


process.patJetPartonAssociationAK5JPT = cms.EDProducer("JetPartonMatcher",
    jets = cms.InputTag("JetPlusTrackZSPCorJetAntiKt5"),
    coneSizeToAssociate = cms.double(0.3),
    partons = cms.InputTag("patJetPartons")
)


process.patJetPartonAssociationAK5PF = cms.EDProducer("JetPartonMatcher",
    jets = cms.InputTag("ak5PFJets"),
    coneSizeToAssociate = cms.double(0.3),
    partons = cms.InputTag("patJetPartons")
)


process.patJetPartonAssociationPF = cms.EDProducer("JetPartonMatcher",
    jets = cms.InputTag("pfNoTauPF"),
    coneSizeToAssociate = cms.double(0.3),
    partons = cms.InputTag("patJetPartonsPF")
)


process.patJetPartonMatch = cms.EDProducer("MCMatcher",
    src = cms.InputTag("ak5CaloJets"),
    maxDPtRel = cms.double(999999.0),
    mcPdgId = cms.vint32(1, 2, 3, 4, 5, 
        21),
    mcStatus = cms.vint32(3),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = cms.double(0.25),
    checkCharge = cms.bool(False),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("genParticles")
)


process.patJetPartonMatchAK5JPT = cms.EDProducer("MCMatcher",
    src = cms.InputTag("JetPlusTrackZSPCorJetAntiKt5"),
    maxDPtRel = cms.double(999999.0),
    mcPdgId = cms.vint32(1, 2, 3, 4, 5, 
        21),
    mcStatus = cms.vint32(3),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = cms.double(0.25),
    checkCharge = cms.bool(False),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("genParticles")
)


process.patJetPartonMatchAK5PF = cms.EDProducer("MCMatcher",
    src = cms.InputTag("ak5PFJets"),
    maxDPtRel = cms.double(999999.0),
    mcPdgId = cms.vint32(1, 2, 3, 4, 5, 
        21),
    mcStatus = cms.vint32(3),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = cms.double(0.25),
    checkCharge = cms.bool(False),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("genParticles")
)


process.patJetPartonMatchPF = cms.EDProducer("MCMatcher",
    src = cms.InputTag("pfNoTauPF"),
    maxDPtRel = cms.double(3.0),
    mcPdgId = cms.vint32(1, 2, 3, 4, 5, 
        21),
    mcStatus = cms.vint32(3),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = cms.double(0.4),
    checkCharge = cms.bool(False),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("genParticles")
)


process.patJetPartons = cms.EDProducer("PartonSelector",
    src = cms.InputTag("genParticles"),
    withLeptons = cms.bool(False)
)


process.patJetPartonsPF = cms.EDProducer("PartonSelector",
    src = cms.InputTag("genParticles"),
    withLeptons = cms.bool(False)
)


process.patJets = cms.EDProducer("PATJetProducer",
    addJetCharge = cms.bool(True),
    addGenJetMatch = cms.bool(False),
    embedPFCandidates = cms.bool(True),
    embedGenJetMatch = cms.bool(False),
    addAssociatedTracks = cms.bool(True),
    partonJetSource = cms.InputTag("NOT_IMPLEMENTED"),
    addGenPartonMatch = cms.bool(False),
    JetPartonMapSource = cms.InputTag(""),
    resolutions = cms.PSet(

    ),
    genPartonMatch = cms.InputTag(""),
    addTagInfos = cms.bool(False),
    addPartonJetMatch = cms.bool(False),
    embedGenPartonMatch = cms.bool(False),
    efficiencies = cms.PSet(

    ),
    genJetMatch = cms.InputTag(""),
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
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring()
    ),
    jetSource = cms.InputTag("ak5CaloJets"),
    addEfficiencies = cms.bool(False),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactors")),
    trackAssociationSource = cms.InputTag("jetTracksAssociatorAtVertex"),
    tagInfoSources = cms.VInputTag(cms.InputTag("impactParameterTagInfosAOD"), cms.InputTag("secondaryVertexTagInfosAOD"), cms.InputTag("softMuonTagInfosAOD")),
    discriminatorSources = cms.VInputTag(cms.InputTag("jetBProbabilityBJetTagsAOD"), cms.InputTag("jetProbabilityBJetTagsAOD"), cms.InputTag("trackCountingHighPurBJetTagsAOD"), cms.InputTag("trackCountingHighEffBJetTagsAOD"), cms.InputTag("simpleSecondaryVertexHighEffBJetTagsAOD"), 
        cms.InputTag("simpleSecondaryVertexHighPurBJetTagsAOD"), cms.InputTag("combinedSecondaryVertexBJetTagsAOD"), cms.InputTag("combinedSecondaryVertexMVABJetTagsAOD"), cms.InputTag("softMuonBJetTagsAOD"), cms.InputTag("softMuonByPtBJetTagsAOD"), 
        cms.InputTag("softMuonByIP3dBJetTagsAOD")),
    addBTagInfo = cms.bool(True),
    embedCaloTowers = cms.bool(True),
    addResolutions = cms.bool(False),
    getJetMCFlavour = cms.bool(False),
    addDiscriminators = cms.bool(True),
    jetChargeSource = cms.InputTag("patJetCharge"),
    addJetCorrFactors = cms.bool(True),
    jetIDMap = cms.InputTag("ak5JetID"),
    addJetID = cms.bool(True)
)


process.patJetsAK5JPT = cms.EDProducer("PATJetProducer",
    addJetCharge = cms.bool(True),
    addGenJetMatch = cms.bool(False),
    embedGenJetMatch = cms.bool(False),
    addAssociatedTracks = cms.bool(True),
    addBTagInfo = cms.bool(True),
    partonJetSource = cms.InputTag("NOT_IMPLEMENTED"),
    addGenPartonMatch = cms.bool(False),
    JetPartonMapSource = cms.InputTag("AK5JPT"),
    resolutions = cms.PSet(

    ),
    genPartonMatch = cms.InputTag("AK5JPT"),
    addTagInfos = cms.bool(False),
    addPartonJetMatch = cms.bool(False),
    embedGenPartonMatch = cms.bool(False),
    efficiencies = cms.PSet(

    ),
    genJetMatch = cms.InputTag("AK5JPT"),
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
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring()
    ),
    jetSource = cms.InputTag("JetPlusTrackZSPCorJetAntiKt5"),
    addEfficiencies = cms.bool(False),
    discriminatorSources = cms.VInputTag(cms.InputTag("jetBProbabilityBJetTagsAK5JPT"), cms.InputTag("jetProbabilityBJetTagsAK5JPT"), cms.InputTag("trackCountingHighPurBJetTagsAK5JPT"), cms.InputTag("trackCountingHighEffBJetTagsAK5JPT"), cms.InputTag("simpleSecondaryVertexHighEffBJetTagsAK5JPT"), 
        cms.InputTag("simpleSecondaryVertexHighPurBJetTagsAK5JPT"), cms.InputTag("combinedSecondaryVertexBJetTagsAK5JPT"), cms.InputTag("combinedSecondaryVertexMVABJetTagsAK5JPT"), cms.InputTag("softMuonBJetTagsAK5JPT"), cms.InputTag("softMuonByPtBJetTagsAK5JPT"), 
        cms.InputTag("softMuonByIP3dBJetTagsAK5JPT")),
    trackAssociationSource = cms.InputTag("jetTracksAssociatorAtVertexAK5JPT"),
    tagInfoSources = cms.VInputTag(cms.InputTag("impactParameterTagInfosAK5JPT"), cms.InputTag("secondaryVertexTagInfosAK5JPT"), cms.InputTag("softMuonTagInfosAK5JPT")),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsAK5JPT")),
    embedPFCandidates = cms.bool(True),
    addJetCorrFactors = cms.bool(True),
    addResolutions = cms.bool(False),
    getJetMCFlavour = cms.bool(False),
    addDiscriminators = cms.bool(True),
    jetChargeSource = cms.InputTag("patJetChargeAK5JPT"),
    embedCaloTowers = cms.bool(True),
    jetIDMap = cms.InputTag("ak5JetID"),
    addJetID = cms.bool(True)
)


process.patJetsAK5PF = cms.EDProducer("PATJetProducer",
    addJetCharge = cms.bool(True),
    addGenJetMatch = cms.bool(False),
    embedGenJetMatch = cms.bool(False),
    addAssociatedTracks = cms.bool(True),
    addBTagInfo = cms.bool(True),
    partonJetSource = cms.InputTag("NOT_IMPLEMENTED"),
    addGenPartonMatch = cms.bool(False),
    JetPartonMapSource = cms.InputTag("AK5PF"),
    resolutions = cms.PSet(

    ),
    genPartonMatch = cms.InputTag("AK5PF"),
    addTagInfos = cms.bool(False),
    addPartonJetMatch = cms.bool(False),
    embedGenPartonMatch = cms.bool(False),
    efficiencies = cms.PSet(

    ),
    genJetMatch = cms.InputTag("AK5PF"),
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
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring()
    ),
    jetSource = cms.InputTag("ak5PFJets"),
    addEfficiencies = cms.bool(False),
    discriminatorSources = cms.VInputTag(cms.InputTag("jetBProbabilityBJetTagsAK5PF"), cms.InputTag("jetProbabilityBJetTagsAK5PF"), cms.InputTag("trackCountingHighPurBJetTagsAK5PF"), cms.InputTag("trackCountingHighEffBJetTagsAK5PF"), cms.InputTag("simpleSecondaryVertexHighEffBJetTagsAK5PF"), 
        cms.InputTag("simpleSecondaryVertexHighPurBJetTagsAK5PF"), cms.InputTag("combinedSecondaryVertexBJetTagsAK5PF"), cms.InputTag("combinedSecondaryVertexMVABJetTagsAK5PF"), cms.InputTag("softMuonBJetTagsAK5PF"), cms.InputTag("softMuonByPtBJetTagsAK5PF"), 
        cms.InputTag("softMuonByIP3dBJetTagsAK5PF")),
    trackAssociationSource = cms.InputTag("jetTracksAssociatorAtVertexAK5PF"),
    tagInfoSources = cms.VInputTag(cms.InputTag("impactParameterTagInfosAK5PF"), cms.InputTag("secondaryVertexTagInfosAK5PF"), cms.InputTag("softMuonTagInfosAK5PF")),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsAK5PF")),
    embedPFCandidates = cms.bool(True),
    addJetCorrFactors = cms.bool(True),
    addResolutions = cms.bool(False),
    getJetMCFlavour = cms.bool(False),
    addDiscriminators = cms.bool(True),
    jetChargeSource = cms.InputTag("patJetChargeAK5PF"),
    embedCaloTowers = cms.bool(True),
    jetIDMap = cms.InputTag("ak5JetID"),
    addJetID = cms.bool(False)
)


process.patJetsPF = cms.EDProducer("PATJetProducer",
    addJetCharge = cms.bool(True),
    addGenJetMatch = cms.bool(False),
    embedGenJetMatch = cms.bool(False),
    addAssociatedTracks = cms.bool(True),
    addBTagInfo = cms.bool(True),
    partonJetSource = cms.InputTag("NOT_IMPLEMENTED"),
    addGenPartonMatch = cms.bool(False),
    JetPartonMapSource = cms.InputTag(""),
    resolutions = cms.PSet(

    ),
    genPartonMatch = cms.InputTag(""),
    addTagInfos = cms.bool(False),
    addPartonJetMatch = cms.bool(False),
    embedGenPartonMatch = cms.bool(False),
    efficiencies = cms.PSet(

    ),
    genJetMatch = cms.InputTag(""),
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
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring()
    ),
    jetSource = cms.InputTag("pfNoTauPF"),
    addEfficiencies = cms.bool(False),
    discriminatorSources = cms.VInputTag(cms.InputTag("jetBProbabilityBJetTagsAODPF"), cms.InputTag("jetProbabilityBJetTagsAODPF"), cms.InputTag("trackCountingHighPurBJetTagsAODPF"), cms.InputTag("trackCountingHighEffBJetTagsAODPF"), cms.InputTag("simpleSecondaryVertexHighEffBJetTagsAODPF"), 
        cms.InputTag("simpleSecondaryVertexHighPurBJetTagsAODPF"), cms.InputTag("combinedSecondaryVertexBJetTagsAODPF"), cms.InputTag("combinedSecondaryVertexMVABJetTagsAODPF"), cms.InputTag("softMuonBJetTagsAODPF"), cms.InputTag("softMuonByPtBJetTagsAODPF"), 
        cms.InputTag("softMuonByIP3dBJetTagsAODPF")),
    trackAssociationSource = cms.InputTag("jetTracksAssociatorAtVertexPF"),
    tagInfoSources = cms.VInputTag(cms.InputTag("impactParameterTagInfosAODPF"), cms.InputTag("secondaryVertexTagInfosAODPF"), cms.InputTag("softMuonTagInfosAODPF")),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsPF")),
    embedPFCandidates = cms.bool(False),
    addJetCorrFactors = cms.bool(True),
    addResolutions = cms.bool(False),
    getJetMCFlavour = cms.bool(False),
    addDiscriminators = cms.bool(True),
    jetChargeSource = cms.InputTag("patJetChargePF"),
    embedCaloTowers = cms.bool(False),
    jetIDMap = cms.InputTag("AK5JetID"),
    addJetID = cms.bool(True)
)


process.patMETs = cms.EDProducer("PATMETProducer",
    metSource = cms.InputTag("metJESCorAK5CaloJetMuons"),
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
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring()
    ),
    addResolutions = cms.bool(False),
    addEfficiencies = cms.bool(False),
    genMETSource = cms.InputTag(""),
    efficiencies = cms.PSet(

    ),
    addGenMET = cms.bool(False),
    addMuonCorrections = cms.bool(True),
    muonSource = cms.InputTag("muons"),
    resolutions = cms.PSet(

    )
)


process.patMETsAK5Calo = cms.EDProducer("PATMETProducer",
    metSource = cms.InputTag("metJESCorAK5CaloJetMuons"),
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
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring()
    ),
    addResolutions = cms.bool(False),
    addEfficiencies = cms.bool(False),
    genMETSource = cms.InputTag(""),
    efficiencies = cms.PSet(

    ),
    addGenMET = cms.bool(False),
    addMuonCorrections = cms.bool(True),
    muonSource = cms.InputTag("muons"),
    resolutions = cms.PSet(

    )
)


process.patMETsAK5CaloTypeII = cms.EDProducer("PATMETProducer",
    metSource = cms.InputTag("metJESCorAK5CaloJetMuonsTypeII"),
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
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring()
    ),
    addResolutions = cms.bool(False),
    muonSource = cms.InputTag("muons"),
    addEfficiencies = cms.bool(False),
    genMETSource = cms.InputTag(""),
    efficiencies = cms.PSet(

    ),
    addGenMET = cms.bool(False),
    addMuonCorrections = cms.bool(True),
    resolutions = cms.PSet(

    )
)


process.patMETsPF = cms.EDProducer("PATMETProducer",
    metSource = cms.InputTag("pfMETPF"),
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
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring()
    ),
    addResolutions = cms.bool(False),
    muonSource = cms.InputTag("muons"),
    addEfficiencies = cms.bool(False),
    genMETSource = cms.InputTag(""),
    efficiencies = cms.PSet(

    ),
    addGenMET = cms.bool(False),
    addMuonCorrections = cms.bool(False),
    resolutions = cms.PSet(

    )
)


process.patMETsTC = cms.EDProducer("PATMETProducer",
    metSource = cms.InputTag("tcMet"),
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
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring()
    ),
    addResolutions = cms.bool(False),
    muonSource = cms.InputTag("muons"),
    addEfficiencies = cms.bool(False),
    genMETSource = cms.InputTag(""),
    efficiencies = cms.PSet(

    ),
    addGenMET = cms.bool(False),
    addMuonCorrections = cms.bool(True),
    resolutions = cms.PSet(

    )
)


process.patMETsTypeIPF = cms.EDProducer("PATMETProducer",
    metSource = cms.InputTag("metJESCorAK5PFTypeI"),
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
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring()
    ),
    addResolutions = cms.bool(False),
    addEfficiencies = cms.bool(False),
    genMETSource = cms.InputTag(""),
    efficiencies = cms.PSet(

    ),
    addGenMET = cms.bool(False),
    addMuonCorrections = cms.bool(False),
    muonSource = cms.InputTag("muons"),
    resolutions = cms.PSet(

    )
)


process.patMHTs = cms.EDProducer("PATMHTProducer",
    verbose = cms.double(0.0),
    muonEtaMax = cms.double(2.5),
    jetTag = cms.untracked.InputTag("allLayer1Jets"),
    eleEtaMax = cms.double(3.0),
    noHF = cms.bool(False),
    muonTag = cms.untracked.InputTag("allLayer1Muons"),
    CaloTowerTag = cms.InputTag("towerMaker"),
    elePhiUncertaintyParameter0 = cms.double(0.01),
    uncertaintyScaleFactor = cms.double(1.0),
    muonPtMin = cms.double(10.0),
    eleEtUncertaintyParameter0 = cms.double(0.01),
    useHO = cms.bool(False),
    jetEtUncertaintyParameter2 = cms.double(0.033),
    jetEtUncertaintyParameter1 = cms.double(1.25),
    jetEMfracMax = cms.double(0.9),
    jetPhiUncertaintyParameter2 = cms.double(0.023),
    jetPhiUncertaintyParameter0 = cms.double(4.75),
    jetPhiUncertaintyParameter1 = cms.double(-0.426),
    tauTag = cms.untracked.InputTag("allLayer1Taus"),
    jetEtUncertaintyParameter0 = cms.double(5.6),
    electronTag = cms.untracked.InputTag("allLayer1Electrons"),
    jetEtaMax = cms.double(5.0),
    elePtMin = cms.double(10.0),
    jetPtMin = cms.double(20.0),
    muonEtUncertaintyParameter0 = cms.double(0.01),
    photonTag = cms.untracked.InputTag("allLayer1Photons"),
    muonPhiUncertaintyParameter0 = cms.double(0.01),
    controlledUncertainty = cms.bool(True),
    towerEtThreshold = cms.double(0.5)
)


process.patMuonMatch = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
    matchedCuts = cms.string('path("HLT_Mu15_v1") || path("HLT_Mu_9") || path("HLT_Mu_11")'),
    src = cms.InputTag("cleanPatMuons"),
    maxDPtRel = cms.double(0.5),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = cms.double(0.5),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("patTrigger")
)


process.patMuonMatchPF = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
    matchedCuts = cms.string('path("HLT_Mu15_v1") || path("HLT_Mu_9") || path("HLT_Mu_11")'),
    src = cms.InputTag("selectedPatMuonsPF"),
    maxDPtRel = cms.double(0.5),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = cms.double(0.5),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("patTriggerPF")
)


process.patMuons = cms.EDProducer("PATMuonProducer",
    embedTpfmsMuon = cms.bool(True),
    embedHighLevelSelection = cms.bool(True),
    embedCaloMETMuonCorrs = cms.bool(True),
    caloMETMuonCorrs = cms.InputTag("muonMETValueMapProducer","muCorrData"),
    resolutions = cms.PSet(

    ),
    userIsolation = cms.PSet(

    ),
    embedPFCandidate = cms.bool(True),
    pfMuonSource = cms.InputTag("pfIsolatedMuons"),
    efficiencies = cms.PSet(

    ),
    embedStandAloneMuon = cms.bool(True),
    useParticleFlow = cms.bool(False),
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
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring()
    ),
    embedTrack = cms.bool(False),
    addEfficiencies = cms.bool(False),
    usePV = cms.bool(True),
    embedTcMETMuonCorrs = cms.bool(True),
    pvSrc = cms.InputTag("offlinePrimaryVertices"),
    addTeVRefits = cms.bool(True),
    embedCombinedMuon = cms.bool(True),
    genParticleMatch = cms.InputTag(""),
    beamLineSrc = cms.InputTag("offlineBeamSpot"),
    muonSource = cms.InputTag("muons"),
    addGenMatch = cms.bool(False),
    addResolutions = cms.bool(False),
    tpfmsSrc = cms.InputTag("tevMuons","firstHit"),
    pickySrc = cms.InputTag("tevMuons","picky"),
    isoDeposits = cms.PSet(

    ),
    embedGenMatch = cms.bool(False),
    tcMETMuonCorrs = cms.InputTag("muonTCMETValueMapProducer","muCorrData"),
    embedPickyMuon = cms.bool(True)
)


process.patMuonsPF = cms.EDProducer("PATMuonProducer",
    embedTpfmsMuon = cms.bool(True),
    embedHighLevelSelection = cms.bool(True),
    embedCaloMETMuonCorrs = cms.bool(True),
    caloMETMuonCorrs = cms.InputTag("muonMETValueMapProducer","muCorrData"),
    addResolutions = cms.bool(False),
    resolutions = cms.PSet(

    ),
    userIsolation = cms.PSet(

    ),
    embedPFCandidate = cms.bool(True),
    pfMuonSource = cms.InputTag("pfIsolatedMuonsPF"),
    efficiencies = cms.PSet(

    ),
    embedStandAloneMuon = cms.bool(True),
    useParticleFlow = cms.bool(True),
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
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring()
    ),
    embedTrack = cms.bool(False),
    addEfficiencies = cms.bool(False),
    usePV = cms.bool(True),
    embedTcMETMuonCorrs = cms.bool(True),
    pvSrc = cms.InputTag("offlinePrimaryVertices"),
    addTeVRefits = cms.bool(True),
    embedCombinedMuon = cms.bool(True),
    genParticleMatch = cms.InputTag(""),
    beamLineSrc = cms.InputTag("offlineBeamSpot"),
    addGenMatch = cms.bool(False),
    muonSource = cms.InputTag("muons"),
    tpfmsSrc = cms.InputTag("tevMuons","firstHit"),
    pickySrc = cms.InputTag("tevMuons","picky"),
    isoDeposits = cms.PSet(
        pfNeutralHadrons = cms.InputTag("isoDepMuonWithNeutralPF"),
        pfChargedHadrons = cms.InputTag("isoDepMuonWithChargedPF"),
        pfPhotons = cms.InputTag("isoDepMuonWithPhotonsPF")
    ),
    embedGenMatch = cms.bool(False),
    tcMETMuonCorrs = cms.InputTag("muonTCMETValueMapProducer","muCorrData"),
    embedPickyMuon = cms.bool(True),
    isolationValues = cms.PSet(
        pfNeutralHadrons = cms.InputTag("isoValMuonWithNeutralPF"),
        pfChargedHadrons = cms.InputTag("isoValMuonWithChargedPF"),
        pfPhotons = cms.InputTag("isoValMuonWithPhotonsPF")
    )
)


process.patPFParticlesPF = cms.EDProducer("PATPFParticleProducer",
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
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring()
    ),
    addGenMatch = cms.bool(False),
    addResolutions = cms.bool(False),
    addEfficiencies = cms.bool(False),
    efficiencies = cms.PSet(

    ),
    embedGenMatch = cms.bool(False),
    pfCandidateSource = cms.InputTag("pfNoJetPF"),
    resolutions = cms.PSet(

    ),
    genParticleMatch = cms.InputTag("")
)


process.patPhotonMatch = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
    matchedCuts = cms.string('path("HLT_Photon30")'),
    src = cms.InputTag("cleanPatPhotons"),
    maxDPtRel = cms.double(0.5),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = cms.double(0.5),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("patTrigger")
)


process.patPhotons = cms.EDProducer("PATPhotonProducer",
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
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring()
    ),
    addGenMatch = cms.bool(False),
    addResolutions = cms.bool(False),
    addEfficiencies = cms.bool(False),
    photonIDSources = cms.PSet(
        PhotonCutBasedIDTight = cms.InputTag("PhotonIDProd","PhotonCutBasedIDTight"),
        PhotonCutBasedIDLoose = cms.InputTag("PhotonIDProd","PhotonCutBasedIDLoose")
    ),
    isoDeposits = cms.PSet(

    ),
    efficiencies = cms.PSet(

    ),
    embedSuperCluster = cms.bool(False),
    embedGenMatch = cms.bool(False),
    resolutions = cms.PSet(

    ),
    addPhotonID = cms.bool(True),
    photonSource = cms.InputTag("photons"),
    userIsolation = cms.PSet(

    ),
    genParticleMatch = cms.InputTag("")
)


process.patPhotonsPF = cms.EDProducer("PATPhotonProducer",
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
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring()
    ),
    addGenMatch = cms.bool(False),
    addResolutions = cms.bool(False),
    addEfficiencies = cms.bool(False),
    photonIDSources = cms.PSet(
        PhotonCutBasedIDTight = cms.InputTag("PhotonIDProd","PhotonCutBasedIDTight"),
        PhotonCutBasedIDLoose = cms.InputTag("PhotonIDProd","PhotonCutBasedIDLoose")
    ),
    isoDeposits = cms.PSet(

    ),
    efficiencies = cms.PSet(

    ),
    embedSuperCluster = cms.bool(True),
    embedGenMatch = cms.bool(False),
    resolutions = cms.PSet(

    ),
    genParticleMatch = cms.InputTag(""),
    photonSource = cms.InputTag("photons"),
    userIsolation = cms.PSet(

    ),
    addPhotonID = cms.bool(True)
)


process.patTauMatch = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
    matchedCuts = cms.string('path("HLT_Jet15U") || path("HLT_Jet30U") || path("HLT_Jet50U")'),
    src = cms.InputTag("cleanPatTaus"),
    maxDPtRel = cms.double(0.5),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = cms.double(0.5),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("patTrigger")
)


process.patTauMatchPF = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
    matchedCuts = cms.string('path("HLT_Jet15U") || path("HLT_Jet30U") || path("HLT_Jet50U")'),
    src = cms.InputTag("selectedPatTausPF"),
    maxDPtRel = cms.double(0.5),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = cms.double(0.5),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("patTrigger")
)


process.patTaus = cms.EDProducer("PATTauProducer",
    tauIDSources = cms.PSet(
        trackIsolation = cms.InputTag("shrinkingConePFTauDiscriminationByTrackIsolation"),
        byTaNCfrHalfPercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrHalfPercent"),
        byTaNC = cms.InputTag("shrinkingConePFTauDiscriminationByTaNC"),
        byTaNCfrOnePercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrOnePercent"),
        byTaNCfrQuarterPercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrQuarterPercent"),
        leadingPionPtCut = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingPionPtCut"),
        leadingTrackFinding = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFinding"),
        byIsolation = cms.InputTag("shrinkingConePFTauDiscriminationByIsolation"),
        byTaNCfrTenthPercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrTenthPercent"),
        ecalIsolationUsingLeadingPion = cms.InputTag("shrinkingConePFTauDiscriminationByECALIsolationUsingLeadingPion"),
        byIsolationUsingLeadingPion = cms.InputTag("shrinkingConePFTauDiscriminationByIsolationUsingLeadingPion"),
        leadingTrackPtCut = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackPtCut"),
        trackIsolationUsingLeadingPion = cms.InputTag("shrinkingConePFTauDiscriminationByTrackIsolationUsingLeadingPion"),
        againstMuon = cms.InputTag("shrinkingConePFTauDiscriminationAgainstMuon"),
        againstElectron = cms.InputTag("shrinkingConePFTauDiscriminationAgainstElectron"),
        ecalIsolation = cms.InputTag("shrinkingConePFTauDiscriminationByECALIsolation")
    ),
    addGenJetMatch = cms.bool(False),
    embedGenJetMatch = cms.bool(False),
    embedLeadTrack = cms.bool(False),
    embedLeadPFCand = cms.bool(False),
    decayModeSrc = cms.InputTag("shrinkingConePFTauDecayModeProducer"),
    embedSignalPFChargedHadrCands = cms.bool(False),
    resolutions = cms.PSet(

    ),
    userIsolation = cms.PSet(
        pfAllParticles = cms.PSet(
            threshold = cms.double(0.0),
            src = cms.InputTag("tauIsoDepositPFCandidates"),
            deltaR = cms.double(0.5)
        ),
        pfNeutralHadron = cms.PSet(
            threshold = cms.double(0.0),
            src = cms.InputTag("tauIsoDepositPFNeutralHadrons"),
            deltaR = cms.double(0.5)
        ),
        pfChargedHadron = cms.PSet(
            threshold = cms.double(0.0),
            src = cms.InputTag("tauIsoDepositPFChargedHadrons"),
            deltaR = cms.double(0.5)
        ),
        pfGamma = cms.PSet(
            threshold = cms.double(0.0),
            src = cms.InputTag("tauIsoDepositPFGammas"),
            deltaR = cms.double(0.5)
        )
    ),
    addDecayMode = cms.bool(True),
    embedIsolationPFGammaCands = cms.bool(False),
    embedSignalPFGammaCands = cms.bool(False),
    efficiencies = cms.PSet(

    ),
    genJetMatch = cms.InputTag(""),
    embedIsolationPFCands = cms.bool(False),
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
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring()
    ),
    embedSignalPFCands = cms.bool(False),
    addEfficiencies = cms.bool(False),
    embedSignalTracks = cms.bool(False),
    tauSource = cms.InputTag("shrinkingConePFTauProducer"),
    embedIsolationPFNeutralHadrCands = cms.bool(False),
    addTauID = cms.bool(True),
    genParticleMatch = cms.InputTag(""),
    addGenMatch = cms.bool(False),
    addResolutions = cms.bool(False),
    embedIsolationPFChargedHadrCands = cms.bool(False),
    embedIsolationTracks = cms.bool(False),
    embedSignalPFNeutralHadrCands = cms.bool(False),
    isoDeposits = cms.PSet(
        pfAllParticles = cms.InputTag("tauIsoDepositPFCandidates"),
        pfNeutralHadron = cms.InputTag("tauIsoDepositPFNeutralHadrons"),
        pfChargedHadron = cms.InputTag("tauIsoDepositPFChargedHadrons"),
        pfGamma = cms.InputTag("tauIsoDepositPFGammas")
    ),
    embedLeadPFChargedHadrCand = cms.bool(False),
    embedGenMatch = cms.bool(False),
    embedLeadPFNeutralCand = cms.bool(False)
)


process.patTausPF = cms.EDProducer("PATTauProducer",
    tauIDSources = cms.PSet(
        leadingTrackFinding = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFindingPF"),
        leadingTrackPtCut = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackPtCutPF"),
        trackIsolation = cms.InputTag("shrinkingConePFTauDiscriminationByTrackIsolationPF"),
        ecalIsolation = cms.InputTag("shrinkingConePFTauDiscriminationByECALIsolationPF"),
        byIsolation = cms.InputTag("shrinkingConePFTauDiscriminationByIsolationPF"),
        againstElectron = cms.InputTag("shrinkingConePFTauDiscriminationAgainstElectronPF"),
        againstMuon = cms.InputTag("shrinkingConePFTauDiscriminationAgainstMuonPF"),
        leadingPionPtCut = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingPionPtCutPF"),
        trackIsolationUsingLeadingPion = cms.InputTag("shrinkingConePFTauDiscriminationByTrackIsolationUsingLeadingPionPF"),
        ecalIsolationUsingLeadingPion = cms.InputTag("shrinkingConePFTauDiscriminationByECALIsolationUsingLeadingPionPF"),
        byIsolationUsingLeadingPion = cms.InputTag("shrinkingConePFTauDiscriminationByIsolationUsingLeadingPionPF"),
        byTaNC = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCPF"),
        byTaNCfrOnePercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrOnePercentPF"),
        byTaNCfrHalfPercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrHalfPercentPF"),
        byTaNCfrQuarterPercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrQuarterPercentPF"),
        byTaNCfrTenthPercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrTenthPercentPF")
    ),
    addGenJetMatch = cms.bool(False),
    embedGenJetMatch = cms.bool(False),
    embedLeadTrack = cms.bool(False),
    embedLeadPFCand = cms.bool(False),
    decayModeSrc = cms.InputTag("shrinkingConePFTauDecayModeProducerPF"),
    embedSignalPFChargedHadrCands = cms.bool(False),
    resolutions = cms.PSet(

    ),
    userIsolation = cms.PSet(
        pfAllParticles = cms.PSet(
            threshold = cms.double(0.0),
            src = cms.InputTag("tauIsoDepositPFCandidatesPF"),
            deltaR = cms.double(0.5)
        ),
        pfNeutralHadron = cms.PSet(
            threshold = cms.double(0.0),
            src = cms.InputTag("tauIsoDepositPFNeutralHadronsPF"),
            deltaR = cms.double(0.5)
        ),
        pfChargedHadron = cms.PSet(
            threshold = cms.double(0.0),
            src = cms.InputTag("tauIsoDepositPFChargedHadronsPF"),
            deltaR = cms.double(0.5)
        ),
        pfGamma = cms.PSet(
            threshold = cms.double(0.0),
            src = cms.InputTag("tauIsoDepositPFGammasPF"),
            deltaR = cms.double(0.5)
        )
    ),
    addDecayMode = cms.bool(True),
    embedIsolationPFGammaCands = cms.bool(False),
    embedSignalPFGammaCands = cms.bool(False),
    efficiencies = cms.PSet(

    ),
    genJetMatch = cms.InputTag(""),
    addEfficiencies = cms.bool(False),
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
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring()
    ),
    embedSignalPFCands = cms.bool(False),
    embedIsolationPFCands = cms.bool(False),
    embedSignalTracks = cms.bool(False),
    tauSource = cms.InputTag("pfTausPF"),
    addResolutions = cms.bool(False),
    addTauID = cms.bool(True),
    genParticleMatch = cms.InputTag(""),
    addGenMatch = cms.bool(False),
    embedIsolationPFNeutralHadrCands = cms.bool(False),
    embedIsolationPFChargedHadrCands = cms.bool(False),
    embedIsolationTracks = cms.bool(False),
    embedSignalPFNeutralHadrCands = cms.bool(False),
    isoDeposits = cms.PSet(
        pfAllParticles = cms.InputTag("tauIsoDepositPFCandidatesPF"),
        pfNeutralHadron = cms.InputTag("tauIsoDepositPFNeutralHadronsPF"),
        pfChargedHadron = cms.InputTag("tauIsoDepositPFChargedHadronsPF"),
        pfGamma = cms.InputTag("tauIsoDepositPFGammasPF")
    ),
    embedLeadPFChargedHadrCand = cms.bool(False),
    embedGenMatch = cms.bool(False),
    embedLeadPFNeutralCand = cms.bool(False)
)


process.patTrigger = cms.EDProducer("PATTriggerProducer",
    processName = cms.string(options.hltname),
    onlyStandAlone = cms.bool(False)
)


process.patTriggerEvent = cms.EDProducer("PATTriggerEventProducer",
    patTriggerMatches = cms.VInputTag(),
    processName = cms.string(options.hltname),
    patTriggerProducer = cms.InputTag("patTrigger")
)


process.patTriggerEventPF = cms.EDProducer("PATTriggerEventProducer",
    patTriggerMatches = cms.VInputTag(),
    patTriggerProducer = cms.InputTag("patTriggerPF"),
    processName = cms.string(options.hltname)
)


process.patTriggerPF = cms.EDProducer("PATTriggerProducer",
    processName = cms.string(options.hltname),
    onlyStandAlone = cms.bool(False)
)


process.pfJetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
    tracks = cms.InputTag("generalTracks"),
    coneSize = cms.double(0.5),
    jets = cms.InputTag("pfJets")
)


process.pfJetTracksAssociatorAtVertexPF = cms.EDProducer("JetTracksAssociatorAtVertex",
    jets = cms.InputTag("pfJetsPF"),
    tracks = cms.InputTag("generalTracks"),
    coneSize = cms.double(0.5)
)


process.pfJets = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(1),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(5.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('PFJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('AntiKt'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.5),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("pfNoElectron"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.pfJetsPF = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(1),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(5.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('PFJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('AntiKt'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.5),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("pfNoElectronPF"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.pfMET = cms.EDProducer("METProducer",
    HB_EtResPar = cms.vdouble(0.0, 1.22, 0.05),
    EE_EtResPar = cms.vdouble(0.2, 0.03, 0.005),
    HF_PhiResPar = cms.vdouble(0.05022),
    PF_PhiResType7 = cms.vdouble(0.02511),
    HE_PhiResPar = cms.vdouble(0.02511),
    PF_PhiResType2 = cms.vdouble(0.002),
    PF_PhiResType3 = cms.vdouble(0.002),
    HF_EtResPar = cms.vdouble(0.0, 1.82, 0.09),
    PF_PhiResType1 = cms.vdouble(0.002),
    PF_PhiResType6 = cms.vdouble(0.02511),
    HO_EtResPar = cms.vdouble(0.0, 1.3, 0.005),
    PF_PhiResType4 = cms.vdouble(0.005),
    PF_PhiResType5 = cms.vdouble(0.02511),
    EB_EtResPar = cms.vdouble(0.2, 0.03, 0.005),
    EB_PhiResPar = cms.vdouble(0.00502),
    HE_EtResPar = cms.vdouble(0.0, 1.3, 0.05),
    HO_PhiResPar = cms.vdouble(0.02511),
    EE_PhiResPar = cms.vdouble(0.02511),
    HB_PhiResPar = cms.vdouble(0.02511),
    PF_EtResType5 = cms.vdouble(0.0, 1.22, 0.05),
    PF_EtResType4 = cms.vdouble(0.2, 0.03, 0.005),
    PF_EtResType7 = cms.vdouble(0.0, 1.22, 0.05),
    PF_EtResType6 = cms.vdouble(0.0, 1.22, 0.05),
    PF_EtResType1 = cms.vdouble(0.05, 0, 0),
    PF_EtResType3 = cms.vdouble(0.05, 0, 0),
    PF_EtResType2 = cms.vdouble(0.05, 0, 0),
    src = cms.InputTag("particleFlow"),
    METType = cms.string('PFMET'),
    calculateSignificance = cms.bool(True),
    alias = cms.string('PFMET'),
    noHF = cms.bool(False),
    globalThreshold = cms.double(0.0),
    InputType = cms.string('PFCandidateCollection')
)


process.pfMETPF = cms.EDProducer("METProducer",
    HB_EtResPar = cms.vdouble(0.0, 1.22, 0.05),
    EE_EtResPar = cms.vdouble(0.2, 0.03, 0.005),
    HF_PhiResPar = cms.vdouble(0.05022),
    InputType = cms.string('PFCandidateCollection'),
    HE_PhiResPar = cms.vdouble(0.02511),
    HB_PhiResPar = cms.vdouble(0.02511),
    noHF = cms.bool(False),
    PF_PhiResType2 = cms.vdouble(0.002),
    PF_PhiResType3 = cms.vdouble(0.002),
    HF_EtResPar = cms.vdouble(0.0, 1.82, 0.09),
    PF_PhiResType1 = cms.vdouble(0.002),
    PF_PhiResType6 = cms.vdouble(0.02511),
    HO_EtResPar = cms.vdouble(0.0, 1.3, 0.005),
    PF_PhiResType4 = cms.vdouble(0.005),
    PF_PhiResType5 = cms.vdouble(0.02511),
    METType = cms.string('PFMET'),
    EB_EtResPar = cms.vdouble(0.2, 0.03, 0.005),
    globalThreshold = cms.double(0.0),
    EB_PhiResPar = cms.vdouble(0.00502),
    src = cms.InputTag("particleFlow"),
    PF_EtResType6 = cms.vdouble(0.0, 1.22, 0.05),
    HO_PhiResPar = cms.vdouble(0.02511),
    calculateSignificance = cms.bool(True),
    EE_PhiResPar = cms.vdouble(0.02511),
    PF_PhiResType7 = cms.vdouble(0.02511),
    alias = cms.string('PFMET'),
    PF_EtResType5 = cms.vdouble(0.0, 1.22, 0.05),
    PF_EtResType4 = cms.vdouble(0.2, 0.03, 0.005),
    PF_EtResType7 = cms.vdouble(0.0, 1.22, 0.05),
    HE_EtResPar = cms.vdouble(0.0, 1.3, 0.05),
    PF_EtResType1 = cms.vdouble(0.05, 0, 0),
    PF_EtResType3 = cms.vdouble(0.05, 0, 0),
    PF_EtResType2 = cms.vdouble(0.05, 0, 0)
)


process.pfNoElectron = cms.EDProducer("TPPFCandidatesOnPFCandidates",
    bottomCollection = cms.InputTag("pfNoMuon"),
    enable = cms.bool(True),
    topCollection = cms.InputTag("pfIsolatedElectrons"),
    name = cms.untracked.string('noElectron'),
    verbose = cms.untracked.bool(False)
)


process.pfNoElectronPF = cms.EDProducer("TPPFCandidatesOnPFCandidates",
    bottomCollection = cms.InputTag("pfNoMuonPF"),
    enable = cms.bool(True),
    topCollection = cms.InputTag("pfIsolatedElectronsPF"),
    name = cms.untracked.string('noElectron'),
    verbose = cms.untracked.bool(False)
)


process.pfNoJet = cms.EDProducer("TPPFJetsOnPFCandidates",
    bottomCollection = cms.InputTag("pfNoElectron"),
    enable = cms.bool(True),
    topCollection = cms.InputTag("pfJets"),
    name = cms.untracked.string('noJet'),
    verbose = cms.untracked.bool(False)
)


process.pfNoJetPF = cms.EDProducer("TPPFJetsOnPFCandidates",
    bottomCollection = cms.InputTag("pfNoElectronPF"),
    enable = cms.bool(True),
    topCollection = cms.InputTag("pfJetsPF"),
    name = cms.untracked.string('noJet'),
    verbose = cms.untracked.bool(False)
)


process.pfNoMuon = cms.EDProducer("TPPFCandidatesOnPFCandidates",
    bottomCollection = cms.InputTag("pfNoPileUp"),
    enable = cms.bool(True),
    topCollection = cms.InputTag("pfIsolatedMuons"),
    name = cms.untracked.string('noMuon'),
    verbose = cms.untracked.bool(False)
)


process.pfNoMuonPF = cms.EDProducer("TPPFCandidatesOnPFCandidates",
    bottomCollection = cms.InputTag("pfNoPileUpPF"),
    enable = cms.bool(True),
    topCollection = cms.InputTag("pfIsolatedMuonsPF"),
    name = cms.untracked.string('noMuon'),
    verbose = cms.untracked.bool(False)
)


process.pfNoPileUp = cms.EDProducer("TPPileUpPFCandidatesOnPFCandidates",
    bottomCollection = cms.InputTag("particleFlow"),
    enable = cms.bool(True),
    topCollection = cms.InputTag("pfPileUp"),
    name = cms.untracked.string('pileUpOnPFCandidates'),
    verbose = cms.untracked.bool(False)
)


process.pfNoPileUpPF = cms.EDProducer("TPPileUpPFCandidatesOnPFCandidates",
    bottomCollection = cms.InputTag("particleFlow"),
    enable = cms.bool(True),
    topCollection = cms.InputTag("pfPileUpPF"),
    name = cms.untracked.string('pileUpOnPFCandidates'),
    verbose = cms.untracked.bool(False)
)


process.pfNoTau = cms.EDProducer("TPPFTausOnPFJets",
    bottomCollection = cms.InputTag("pfJets"),
    enable = cms.bool(True),
    topCollection = cms.InputTag("pfTaus"),
    name = cms.untracked.string('noTau'),
    verbose = cms.untracked.bool(False)
)


process.pfNoTauPF = cms.EDProducer("TPPFTausOnPFJets",
    bottomCollection = cms.InputTag("pfJetsPF"),
    enable = cms.bool(True),
    topCollection = cms.InputTag("pfTausPF"),
    name = cms.untracked.string('noTau'),
    verbose = cms.untracked.bool(False)
)


process.pfPileUp = cms.EDProducer("PFPileUp",
    Enable = cms.bool(True),
    PFCandidates = cms.InputTag("particleFlow"),
    verbose = cms.untracked.bool(False),
    Vertices = cms.InputTag("offlinePrimaryVertices")
)


process.pfPileUpPF = cms.EDProducer("PFPileUp",
    Enable = cms.bool(True),
    PFCandidates = cms.InputTag("particleFlow"),
    verbose = cms.untracked.bool(False),
    Vertices = cms.InputTag("offlinePrimaryVertices")
)


process.pfRecoTauDiscriminationAgainstElectron = cms.EDProducer("PFRecoTauDiscriminationAgainstElectron",
    ApplyCut_ElectronPreID_2D = cms.bool(False),
    ElecPreID0_HOverPLead_minValue = cms.double(0.05),
    PFTauProducer = cms.InputTag("pfRecoTauProducer"),
    ApplyCut_ElectronPreID = cms.bool(False),
    ApplyCut_HcalTotOverPLead = cms.bool(False),
    EOverPLead_minValue = cms.double(0.8),
    ElecPreID1_EOverPLead_maxValue = cms.double(0.8),
    HcalMaxOverPLead_minValue = cms.double(0.1),
    ApplyCut_EmFraction = cms.bool(False),
    EmFraction_maxValue = cms.double(0.9),
    ApplyCut_PFElectronMVA = cms.bool(True),
    PFElectronMVA_maxValue = cms.double(-0.1),
    ApplyCut_HcalMaxOverPLead = cms.bool(False),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("pfRecoTauDiscriminationByLeadingTrackFinding")
        )
    ),
    Hcal3x3OverPLead_minValue = cms.double(0.1),
    ElecPreID1_HOverPLead_minValue = cms.double(0.15),
    ElecPreID0_EOverPLead_maxValue = cms.double(0.95),
    BremsRecoveryEOverPLead_minValue = cms.double(0.8),
    ApplyCut_EcalCrackCut = cms.bool(False),
    EOverPLead_maxValue = cms.double(1.8),
    HcalTotOverPLead_minValue = cms.double(0.1),
    ApplyCut_BremsRecoveryEOverPLead = cms.bool(False),
    ApplyCut_Hcal3x3OverPLead = cms.bool(False),
    ApplyCut_EOverPLead = cms.bool(False),
    BremsRecoveryEOverPLead_maxValue = cms.double(1.8)
)


process.pfRecoTauDiscriminationAgainstMuon = cms.EDProducer("PFRecoTauDiscriminationAgainstMuon",
    a = cms.double(0.5),
    c = cms.double(0.0),
    b = cms.double(0.5),
    PFTauProducer = cms.InputTag("pfRecoTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("pfRecoTauDiscriminationByLeadingTrackFinding")
        )
    ),
    discriminatorOption = cms.string('noSegMatch')
)


process.pfRecoTauDiscriminationByECALIsolation = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("pfRecoTauProducer"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("pfRecoTauDiscriminationByLeadingTrackFinding")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(False),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.pfRecoTauDiscriminationByECALIsolationUsingLeadingPion = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("pfRecoTauProducer"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadPion = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("pfRecoTauDiscriminationByLeadingTrackFinding")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(False),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.pfRecoTauDiscriminationByIsolation = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("pfRecoTauProducer"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("pfRecoTauDiscriminationByLeadingTrackFinding")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.pfRecoTauDiscriminationByIsolationUsingLeadingPion = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("pfRecoTauProducer"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadPion = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("pfRecoTauDiscriminationByLeadingTrackFinding")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.pfRecoTauDiscriminationByLeadingPionPtCut = cms.EDProducer("PFRecoTauDiscriminationByLeadingObjectPtCut",
    MinPtLeadingObject = cms.double(5.0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    PFTauProducer = cms.InputTag("pfRecoTauProducer"),
    UseOnlyChargedHadrons = cms.bool(False)
)


process.pfRecoTauDiscriminationByLeadingTrackFinding = cms.EDProducer("PFRecoTauDiscriminationByLeadingObjectPtCut",
    MinPtLeadingObject = cms.double(0.0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    PFTauProducer = cms.InputTag("pfRecoTauProducer"),
    UseOnlyChargedHadrons = cms.bool(True)
)


process.pfRecoTauDiscriminationByLeadingTrackPtCut = cms.EDProducer("PFRecoTauDiscriminationByLeadingObjectPtCut",
    MinPtLeadingObject = cms.double(5.0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    PFTauProducer = cms.InputTag("pfRecoTauProducer"),
    UseOnlyChargedHadrons = cms.bool(True)
)


process.pfRecoTauDiscriminationByTrackIsolation = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("pfRecoTauProducer"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("pfRecoTauDiscriminationByLeadingTrackFinding")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.pfRecoTauDiscriminationByTrackIsolationUsingLeadingPion = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("pfRecoTauProducer"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadPion = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("pfRecoTauDiscriminationByLeadingTrackFinding")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.pfRecoTauProducer = cms.EDProducer("PFRecoTauProducer",
    Rphi = cms.double(2.0),
    LeadTrack_minPt = cms.double(0.0),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    ECALSignalConeSizeFormula = cms.string('0.15'),
    TrackerIsolConeMetric = cms.string('DR'),
    TrackerSignalConeMetric = cms.string('DR'),
    EcalStripSumE_deltaPhiOverQ_minValue = cms.double(-0.1),
    MaxEtInEllipse = cms.double(2.0),
    MatchingConeMetric = cms.string('DR'),
    TrackerSignalConeSizeFormula = cms.string('0.07'),
    MatchingConeSizeFormula = cms.string('0.1'),
    TrackerIsolConeSize_min = cms.double(0.0),
    MatchingConeSize_min = cms.double(0.0),
    ElectronPreIDProducer = cms.InputTag("elecpreid"),
    ChargedHadrCandLeadChargedHadrCand_tksmaxDZ = cms.double(0.2),
    TrackerIsolConeSize_max = cms.double(0.6),
    TrackerSignalConeSize_max = cms.double(0.07),
    HCALIsolConeMetric = cms.string('DR'),
    AddEllipseGammas = cms.bool(False),
    maximumForElectrionPreIDOutput = cms.double(-0.1),
    TrackerSignalConeSize_min = cms.double(0.0),
    ECALIsolConeSize_max = cms.double(0.6),
    HCALIsolConeSizeFormula = cms.string('0.50'),
    Track_IsolAnnulus_minNhits = cms.uint32(3),
    smearedPVsigmaZ = cms.double(0.005),
    AreaMetric_recoElements_maxabsEta = cms.double(2.5),
    HCALSignalConeMetric = cms.string('DR'),
    ElecPreIDLeadTkMatch_maxDR = cms.double(0.01),
    ChargedHadrCand_IsolAnnulus_minNhits = cms.uint32(0),
    PFTauTagInfoProducer = cms.InputTag("pfRecoTauTagInfoProducer"),
    ECALIsolConeMetric = cms.string('DR'),
    ECALIsolConeSizeFormula = cms.string('0.50'),
    UseChargedHadrCandLeadChargedHadrCand_tksDZconstraint = cms.bool(True),
    Algorithm = cms.string('ConeBased'),
    JetPtMin = cms.double(0.0),
    ECALSignalConeMetric = cms.string('DR'),
    EcalStripSumE_deltaPhiOverQ_maxValue = cms.double(0.5),
    HCALSignalConeSize_max = cms.double(0.6),
    ECALSignalConeSize_min = cms.double(0.0),
    EcalStripSumE_minClusEnergy = cms.double(0.1),
    EcalStripSumE_deltaEta = cms.double(0.03),
    TrackerIsolConeSizeFormula = cms.string('0.50'),
    LeadPFCand_minPt = cms.double(5.0),
    HCALSignalConeSize_min = cms.double(0.0),
    ECALSignalConeSize_max = cms.double(0.6),
    HCALSignalConeSizeFormula = cms.string('0.10'),
    TrackLeadTrack_maxDZ = cms.double(0.2),
    DataType = cms.string('AOD'),
    ECALIsolConeSize_min = cms.double(0.0),
    UseTrackLeadTrackDZconstraint = cms.bool(True),
    smearedPVsigmaY = cms.double(0.0015),
    HCALIsolConeSize_max = cms.double(0.6),
    smearedPVsigmaX = cms.double(0.0015),
    MatchingConeSize_max = cms.double(0.6),
    HCALIsolConeSize_min = cms.double(0.0)
)


process.pfRecoTauTagInfoProducer = cms.EDProducer("PFRecoTauTagInfoProducer",
    tkminTrackerHitsn = cms.int32(3),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    tkmaxChi2 = cms.double(100.0),
    ChargedHadrCand_AssociationCone = cms.double(0.8),
    ChargedHadrCand_tkminTrackerHitsn = cms.int32(3),
    ChargedHadrCand_tkmaxChi2 = cms.double(100.0),
    tkPVmaxDZ = cms.double(0.2),
    GammaCand_EcalclusMinEt = cms.double(0.5),
    tkminPixelHitsn = cms.int32(0),
    tkminPt = cms.double(0.5),
    PFCandidateProducer = cms.InputTag("pfNoElectron"),
    ChargedHadrCand_tkminPt = cms.double(0.5),
    ChargedHadrCand_tkmaxipt = cms.double(0.03),
    ChargedHadrCand_tkminPixelHitsn = cms.int32(0),
    UsePVconstraint = cms.bool(True),
    NeutrHadrCand_HcalclusMinEt = cms.double(1.0),
    PFJetTracksAssociatorProducer = cms.InputTag("pfJetTracksAssociatorAtVertex"),
    smearedPVsigmaY = cms.double(0.0015),
    smearedPVsigmaX = cms.double(0.0015),
    smearedPVsigmaZ = cms.double(0.005),
    ChargedHadrCand_tkPVmaxDZ = cms.double(0.2),
    tkmaxipt = cms.double(0.03)
)


process.pfRecoTauTagInfoProducerInsideOut = cms.EDProducer("PFRecoTauTagInfoProducer",
    tkminTrackerHitsn = cms.int32(3),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    tkmaxChi2 = cms.double(100.0),
    ChargedHadrCand_AssociationCone = cms.double(1.0),
    ChargedHadrCand_tkminTrackerHitsn = cms.int32(3),
    ChargedHadrCand_tkmaxChi2 = cms.double(100.0),
    tkPVmaxDZ = cms.double(0.2),
    GammaCand_EcalclusMinEt = cms.double(0.5),
    tkminPixelHitsn = cms.int32(0),
    tkminPt = cms.double(0.5),
    PFCandidateProducer = cms.InputTag("particleFlow"),
    ChargedHadrCand_tkminPt = cms.double(0.5),
    ChargedHadrCand_tkmaxipt = cms.double(0.03),
    ChargedHadrCand_tkminPixelHitsn = cms.int32(0),
    UsePVconstraint = cms.bool(True),
    NeutrHadrCand_HcalclusMinEt = cms.double(1.0),
    PFJetTracksAssociatorProducer = cms.InputTag("insideOutJetTracksAssociatorAtVertex"),
    smearedPVsigmaY = cms.double(0.0015),
    smearedPVsigmaX = cms.double(0.0015),
    smearedPVsigmaZ = cms.double(0.005),
    ChargedHadrCand_tkPVmaxDZ = cms.double(0.2),
    tkmaxipt = cms.double(0.03)
)


process.pfRecoTauTagInfoProducerPF = cms.EDProducer("PFRecoTauTagInfoProducer",
    tkminTrackerHitsn = cms.int32(3),
    tkminPt = cms.double(0.5),
    tkmaxChi2 = cms.double(100.0),
    ChargedHadrCand_AssociationCone = cms.double(0.8),
    ChargedHadrCand_tkminTrackerHitsn = cms.int32(3),
    ChargedHadrCand_tkmaxChi2 = cms.double(100.0),
    tkPVmaxDZ = cms.double(0.2),
    GammaCand_EcalclusMinEt = cms.double(0.5),
    tkminPixelHitsn = cms.int32(0),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFCandidateProducer = cms.InputTag("pfNoElectronPF"),
    ChargedHadrCand_tkminPt = cms.double(0.5),
    ChargedHadrCand_tkmaxipt = cms.double(0.03),
    ChargedHadrCand_tkminPixelHitsn = cms.int32(0),
    UsePVconstraint = cms.bool(True),
    NeutrHadrCand_HcalclusMinEt = cms.double(1.0),
    PFJetTracksAssociatorProducer = cms.InputTag("pfJetTracksAssociatorAtVertexPF"),
    smearedPVsigmaY = cms.double(0.0015),
    smearedPVsigmaX = cms.double(0.0015),
    smearedPVsigmaZ = cms.double(0.005),
    ChargedHadrCand_tkPVmaxDZ = cms.double(0.2),
    tkmaxipt = cms.double(0.03)
)


process.pfTauDecayMode = cms.EDProducer("PFRecoTauDecayModeDeterminator",
    mergeByBestMatch = cms.bool(True),
    refitTracks = cms.bool(False),
    maxPiZeroMass = cms.double(0.2),
    mergeLowPtPhotonsFirst = cms.bool(True),
    setMergedPi0Mass = cms.bool(True),
    setChargedPionMass = cms.bool(True),
    filterPhotons = cms.bool(True),
    minPtFractionSinglePhotons = cms.double(0.1),
    minPtFractionPiZeroes = cms.double(0.15),
    maxNbrOfIterations = cms.int32(10),
    filterTwoProngs = cms.bool(True),
    minPtFractionForSecondProng = cms.double(0.1),
    maxDistance = cms.double(0.01),
    setPi0Mass = cms.bool(True),
    maxPhotonsToMerge = cms.uint32(2),
    PFTauProducer = cms.InputTag("pfRecoTauProducer")
)


process.pfTauDecayModeHighEfficiency = cms.EDProducer("PFRecoTauDecayModeDeterminator",
    mergeByBestMatch = cms.bool(True),
    refitTracks = cms.bool(False),
    maxPiZeroMass = cms.double(0.2),
    mergeLowPtPhotonsFirst = cms.bool(True),
    setMergedPi0Mass = cms.bool(True),
    setChargedPionMass = cms.bool(True),
    filterPhotons = cms.bool(True),
    minPtFractionSinglePhotons = cms.double(0.1),
    minPtFractionPiZeroes = cms.double(0.15),
    maxNbrOfIterations = cms.int32(10),
    filterTwoProngs = cms.bool(True),
    minPtFractionForSecondProng = cms.double(0.1),
    maxDistance = cms.double(0.01),
    setPi0Mass = cms.bool(True),
    maxPhotonsToMerge = cms.uint32(2),
    PFTauProducer = cms.InputTag("pfRecoTauProducerHighEfficiency")
)


process.pfTauDecayModeIndexProducer = cms.EDProducer("PFRecoTauDecayModeIndexProducer",
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    PFTauProducer = cms.InputTag("pfRecoTauProducer"),
    PFTauDecayModeProducer = cms.InputTag("pfRecoTauDecayModeProducer")
)


process.pfTauDecayModeInsideOut = cms.EDProducer("PFRecoTauDecayModeDeterminator",
    mergeByBestMatch = cms.bool(True),
    refitTracks = cms.bool(False),
    maxPiZeroMass = cms.double(0.2),
    mergeLowPtPhotonsFirst = cms.bool(True),
    setMergedPi0Mass = cms.bool(True),
    setChargedPionMass = cms.bool(True),
    filterPhotons = cms.bool(True),
    minPtFractionSinglePhotons = cms.double(0.1),
    minPtFractionPiZeroes = cms.double(0.15),
    maxNbrOfIterations = cms.int32(10),
    filterTwoProngs = cms.bool(True),
    minPtFractionForSecondProng = cms.double(0.1),
    maxDistance = cms.double(0.01),
    setPi0Mass = cms.bool(True),
    maxPhotonsToMerge = cms.uint32(2),
    PFTauProducer = cms.InputTag("pfRecoTauProducerInsideOut")
)


process.pfTausBaseDiscriminationByIsolationPF = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducerPF"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("pfTausBaseDiscriminationByLeadingTrackFindingPF")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.pfTausBaseDiscriminationByLeadingPionPtCutPF = cms.EDProducer("PFRecoTauDiscriminationByLeadingObjectPtCut",
    MinPtLeadingObject = cms.double(5.0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    UseOnlyChargedHadrons = cms.bool(False),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducerPF")
)


process.pfTausBaseDiscriminationByLeadingTrackFindingPF = cms.EDProducer("PFRecoTauDiscriminationByLeadingObjectPtCut",
    MinPtLeadingObject = cms.double(0.0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    UseOnlyChargedHadrons = cms.bool(True),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducerPF")
)


process.pfTausDiscriminationByIsolation = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducer"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("pfTausDiscriminationByLeadingTrackFinding")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.pfTausDiscriminationByIsolationPF = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducerPF"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("pfTausDiscriminationByLeadingTrackFindingPF")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.pfTausDiscriminationByLeadingPionPtCut = cms.EDProducer("PFRecoTauDiscriminationByLeadingObjectPtCut",
    MinPtLeadingObject = cms.double(5.0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    UseOnlyChargedHadrons = cms.bool(False),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducer")
)


process.pfTausDiscriminationByLeadingPionPtCutPF = cms.EDProducer("PFRecoTauDiscriminationByLeadingObjectPtCut",
    MinPtLeadingObject = cms.double(5.0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducerPF"),
    UseOnlyChargedHadrons = cms.bool(False)
)


process.pfTausDiscriminationByLeadingTrackFinding = cms.EDProducer("PFRecoTauDiscriminationByLeadingObjectPtCut",
    MinPtLeadingObject = cms.double(0.0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    UseOnlyChargedHadrons = cms.bool(True),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducer")
)


process.pfTausDiscriminationByLeadingTrackFindingPF = cms.EDProducer("PFRecoTauDiscriminationByLeadingObjectPtCut",
    MinPtLeadingObject = cms.double(0.0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducerPF"),
    UseOnlyChargedHadrons = cms.bool(True)
)


process.photonMatch = cms.EDProducer("MCMatcher",
    src = cms.InputTag("photons"),
    maxDPtRel = cms.double(1.0),
    mcPdgId = cms.vint32(22),
    mcStatus = cms.vint32(1),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = cms.double(0.2),
    checkCharge = cms.bool(True),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("genParticles")
)


process.photonMatchPF = cms.EDProducer("MCMatcher",
    src = cms.InputTag("photons"),
    maxDPtRel = cms.double(1.0),
    mcPdgId = cms.vint32(22),
    mcStatus = cms.vint32(1),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = cms.double(0.2),
    checkCharge = cms.bool(True),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("genParticles")
)


process.secondaryVertexNegativeTagInfos = cms.EDProducer("SecondaryVertexProducer",
    trackSelection = cms.PSet(
        totalHitsMin = cms.uint32(8),
        jetDeltaRMax = cms.double(0.3),
        qualityClass = cms.string('highPurity'),
        pixelHitsMin = cms.uint32(2),
        sip3dSigMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        maxDistToAxis = cms.double(0.2),
        sip2dValMax = cms.double(99999.9),
        maxDecayLen = cms.double(99999.9),
        ptMin = cms.double(1.0),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        sip2dValMin = cms.double(-99999.9),
        normChi2Max = cms.double(99999.9)
    ),
    vertexSelection = cms.PSet(
        sortCriterium = cms.string('dist3dError')
    ),
    vertexCuts = cms.PSet(
        distSig3dMax = cms.double(99999.9),
        fracPV = cms.double(0.65),
        distVal2dMax = cms.double(-0.01),
        useTrackWeights = cms.bool(True),
        maxDeltaRToJetAxis = cms.double(-0.5),
        v0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        ),
        distSig2dMin = cms.double(-99999.9),
        multiplicityMin = cms.uint32(2),
        massMax = cms.double(6.5),
        distSig2dMax = cms.double(-3.0),
        distVal3dMax = cms.double(99999.9),
        minimumTrackWeight = cms.double(0.5),
        distVal3dMin = cms.double(-99999.9),
        distVal2dMin = cms.double(-2.5),
        distSig3dMin = cms.double(-99999.9)
    ),
    vertexReco = cms.PSet(
        seccut = cms.double(6.0),
        primcut = cms.double(1.8),
        smoothing = cms.bool(False),
        weightthreshold = cms.double(0.001),
        minweight = cms.double(0.5),
        finder = cms.string('avr')
    ),
    beamSpotTag = cms.InputTag("offlineBeamSpot"),
    constraint = cms.string('BeamSpot'),
    trackIPTagInfos = cms.InputTag("impactParameterTagInfos"),
    minimumTrackWeight = cms.double(0.5),
    usePVError = cms.bool(True),
    trackSort = cms.string('sip3dSig')
)


process.secondaryVertexTagInfos = cms.EDProducer("SecondaryVertexProducer",
    trackSelection = cms.PSet(
        totalHitsMin = cms.uint32(8),
        jetDeltaRMax = cms.double(0.3),
        qualityClass = cms.string('highPurity'),
        pixelHitsMin = cms.uint32(2),
        sip3dSigMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        maxDistToAxis = cms.double(0.2),
        sip2dValMax = cms.double(99999.9),
        maxDecayLen = cms.double(99999.9),
        ptMin = cms.double(1.0),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        sip2dValMin = cms.double(-99999.9),
        normChi2Max = cms.double(99999.9)
    ),
    vertexSelection = cms.PSet(
        sortCriterium = cms.string('dist3dError')
    ),
    vertexCuts = cms.PSet(
        distSig3dMax = cms.double(99999.9),
        fracPV = cms.double(0.65),
        distVal2dMax = cms.double(2.5),
        useTrackWeights = cms.bool(True),
        maxDeltaRToJetAxis = cms.double(0.5),
        v0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        ),
        distSig2dMin = cms.double(3.0),
        multiplicityMin = cms.uint32(2),
        massMax = cms.double(6.5),
        distSig2dMax = cms.double(99999.9),
        distVal3dMax = cms.double(99999.9),
        minimumTrackWeight = cms.double(0.5),
        distVal3dMin = cms.double(-99999.9),
        distVal2dMin = cms.double(0.01),
        distSig3dMin = cms.double(-99999.9)
    ),
    vertexReco = cms.PSet(
        seccut = cms.double(6.0),
        primcut = cms.double(1.8),
        smoothing = cms.bool(False),
        weightthreshold = cms.double(0.001),
        minweight = cms.double(0.5),
        finder = cms.string('avr')
    ),
    beamSpotTag = cms.InputTag("offlineBeamSpot"),
    constraint = cms.string('BeamSpot'),
    trackIPTagInfos = cms.InputTag("impactParameterTagInfos"),
    minimumTrackWeight = cms.double(0.5),
    usePVError = cms.bool(True),
    trackSort = cms.string('sip3dSig')
)


process.secondaryVertexTagInfosAK5JPT = cms.EDProducer("SecondaryVertexProducer",
    beamSpotTag = cms.InputTag("offlineBeamSpot"),
    vertexReco = cms.PSet(
        seccut = cms.double(6.0),
        primcut = cms.double(1.8),
        smoothing = cms.bool(False),
        weightthreshold = cms.double(0.001),
        minweight = cms.double(0.5),
        finder = cms.string('avr')
    ),
    vertexSelection = cms.PSet(
        sortCriterium = cms.string('dist3dError')
    ),
    constraint = cms.string('BeamSpot'),
    vertexCuts = cms.PSet(
        distSig3dMax = cms.double(99999.9),
        fracPV = cms.double(0.65),
        distVal2dMax = cms.double(2.5),
        useTrackWeights = cms.bool(True),
        maxDeltaRToJetAxis = cms.double(0.5),
        v0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        ),
        distSig2dMin = cms.double(3.0),
        multiplicityMin = cms.uint32(2),
        massMax = cms.double(6.5),
        distSig2dMax = cms.double(99999.9),
        distVal3dMax = cms.double(99999.9),
        minimumTrackWeight = cms.double(0.5),
        distVal3dMin = cms.double(-99999.9),
        distVal2dMin = cms.double(0.01),
        distSig3dMin = cms.double(-99999.9)
    ),
    trackIPTagInfos = cms.InputTag("impactParameterTagInfosAK5JPT"),
    minimumTrackWeight = cms.double(0.5),
    usePVError = cms.bool(True),
    trackSelection = cms.PSet(
        totalHitsMin = cms.uint32(8),
        jetDeltaRMax = cms.double(0.3),
        qualityClass = cms.string('highPurity'),
        pixelHitsMin = cms.uint32(2),
        sip3dSigMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        maxDistToAxis = cms.double(0.2),
        sip2dValMax = cms.double(99999.9),
        maxDecayLen = cms.double(99999.9),
        ptMin = cms.double(1.0),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        sip2dValMin = cms.double(-99999.9),
        normChi2Max = cms.double(99999.9)
    ),
    trackSort = cms.string('sip3dSig')
)


process.secondaryVertexTagInfosAK5PF = cms.EDProducer("SecondaryVertexProducer",
    beamSpotTag = cms.InputTag("offlineBeamSpot"),
    vertexReco = cms.PSet(
        seccut = cms.double(6.0),
        primcut = cms.double(1.8),
        smoothing = cms.bool(False),
        weightthreshold = cms.double(0.001),
        minweight = cms.double(0.5),
        finder = cms.string('avr')
    ),
    vertexSelection = cms.PSet(
        sortCriterium = cms.string('dist3dError')
    ),
    constraint = cms.string('BeamSpot'),
    vertexCuts = cms.PSet(
        distSig3dMax = cms.double(99999.9),
        fracPV = cms.double(0.65),
        distVal2dMax = cms.double(2.5),
        useTrackWeights = cms.bool(True),
        maxDeltaRToJetAxis = cms.double(0.5),
        v0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        ),
        distSig2dMin = cms.double(3.0),
        multiplicityMin = cms.uint32(2),
        massMax = cms.double(6.5),
        distSig2dMax = cms.double(99999.9),
        distVal3dMax = cms.double(99999.9),
        minimumTrackWeight = cms.double(0.5),
        distVal3dMin = cms.double(-99999.9),
        distVal2dMin = cms.double(0.01),
        distSig3dMin = cms.double(-99999.9)
    ),
    trackIPTagInfos = cms.InputTag("impactParameterTagInfosAK5PF"),
    minimumTrackWeight = cms.double(0.5),
    usePVError = cms.bool(True),
    trackSelection = cms.PSet(
        totalHitsMin = cms.uint32(8),
        jetDeltaRMax = cms.double(0.3),
        qualityClass = cms.string('highPurity'),
        pixelHitsMin = cms.uint32(2),
        sip3dSigMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        maxDistToAxis = cms.double(0.2),
        sip2dValMax = cms.double(99999.9),
        maxDecayLen = cms.double(99999.9),
        ptMin = cms.double(1.0),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        sip2dValMin = cms.double(-99999.9),
        normChi2Max = cms.double(99999.9)
    ),
    trackSort = cms.string('sip3dSig')
)


process.secondaryVertexTagInfosAOD = cms.EDProducer("SecondaryVertexProducer",
    beamSpotTag = cms.InputTag("offlineBeamSpot"),
    vertexReco = cms.PSet(
        seccut = cms.double(6.0),
        primcut = cms.double(1.8),
        smoothing = cms.bool(False),
        weightthreshold = cms.double(0.001),
        minweight = cms.double(0.5),
        finder = cms.string('avr')
    ),
    vertexSelection = cms.PSet(
        sortCriterium = cms.string('dist3dError')
    ),
    constraint = cms.string('BeamSpot'),
    vertexCuts = cms.PSet(
        distSig3dMax = cms.double(99999.9),
        fracPV = cms.double(0.65),
        distVal2dMax = cms.double(2.5),
        useTrackWeights = cms.bool(True),
        maxDeltaRToJetAxis = cms.double(0.5),
        v0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        ),
        distSig2dMin = cms.double(3.0),
        multiplicityMin = cms.uint32(2),
        massMax = cms.double(6.5),
        distSig2dMax = cms.double(99999.9),
        distVal3dMax = cms.double(99999.9),
        minimumTrackWeight = cms.double(0.5),
        distVal3dMin = cms.double(-99999.9),
        distVal2dMin = cms.double(0.01),
        distSig3dMin = cms.double(-99999.9)
    ),
    trackIPTagInfos = cms.InputTag("impactParameterTagInfosAOD"),
    minimumTrackWeight = cms.double(0.5),
    usePVError = cms.bool(True),
    trackSelection = cms.PSet(
        totalHitsMin = cms.uint32(8),
        jetDeltaRMax = cms.double(0.3),
        qualityClass = cms.string('highPurity'),
        pixelHitsMin = cms.uint32(2),
        sip3dSigMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        maxDistToAxis = cms.double(0.2),
        sip2dValMax = cms.double(99999.9),
        maxDecayLen = cms.double(99999.9),
        ptMin = cms.double(1.0),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        sip2dValMin = cms.double(-99999.9),
        normChi2Max = cms.double(99999.9)
    ),
    trackSort = cms.string('sip3dSig')
)


process.secondaryVertexTagInfosAODPF = cms.EDProducer("SecondaryVertexProducer",
    beamSpotTag = cms.InputTag("offlineBeamSpot"),
    vertexReco = cms.PSet(
        seccut = cms.double(6.0),
        primcut = cms.double(1.8),
        smoothing = cms.bool(False),
        weightthreshold = cms.double(0.001),
        minweight = cms.double(0.5),
        finder = cms.string('avr')
    ),
    vertexSelection = cms.PSet(
        sortCriterium = cms.string('dist3dError')
    ),
    constraint = cms.string('BeamSpot'),
    vertexCuts = cms.PSet(
        distSig3dMax = cms.double(99999.9),
        fracPV = cms.double(0.65),
        distVal2dMax = cms.double(2.5),
        useTrackWeights = cms.bool(True),
        maxDeltaRToJetAxis = cms.double(0.5),
        v0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        ),
        distSig2dMin = cms.double(3.0),
        multiplicityMin = cms.uint32(2),
        massMax = cms.double(6.5),
        distSig2dMax = cms.double(99999.9),
        distVal3dMax = cms.double(99999.9),
        minimumTrackWeight = cms.double(0.5),
        distVal3dMin = cms.double(-99999.9),
        distVal2dMin = cms.double(0.01),
        distSig3dMin = cms.double(-99999.9)
    ),
    trackIPTagInfos = cms.InputTag("impactParameterTagInfosAODPF"),
    minimumTrackWeight = cms.double(0.5),
    usePVError = cms.bool(True),
    trackSelection = cms.PSet(
        totalHitsMin = cms.uint32(8),
        jetDeltaRMax = cms.double(0.3),
        qualityClass = cms.string('highPurity'),
        pixelHitsMin = cms.uint32(2),
        sip3dSigMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        maxDistToAxis = cms.double(0.2),
        sip2dValMax = cms.double(99999.9),
        maxDecayLen = cms.double(99999.9),
        ptMin = cms.double(1.0),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        sip2dValMin = cms.double(-99999.9),
        normChi2Max = cms.double(99999.9)
    ),
    trackSort = cms.string('sip3dSig')
)


process.selectedPatElectronsTriggerMatchPF = cms.EDProducer("PATTriggerMatchElectronEmbedder",
    matches = cms.VInputTag("patElectronMatchPF"),
    src = cms.InputTag("selectedPatElectronsPF")
)


process.selectedPatJetsTriggerMatchPF = cms.EDProducer("PATTriggerMatchJetEmbedder",
    matches = cms.VInputTag("patJetMatchPF"),
    src = cms.InputTag("selectedPatJetsPF")
)


process.selectedPatMuonsTriggerMatchPF = cms.EDProducer("PATTriggerMatchMuonEmbedder",
    matches = cms.VInputTag("patMuonMatchPF"),
    src = cms.InputTag("selectedPatMuonsPF")
)


process.selectedPatTausTriggerMatchPF = cms.EDProducer("PATTriggerMatchTauEmbedder",
    matches = cms.VInputTag("patTauMatchPF"),
    src = cms.InputTag("selectedPatTausPF")
)


process.shrinkingConePFTauDecayModeIndexProducer = cms.EDProducer("PFRecoTauDecayModeIndexProducer",
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducer"),
    PFTauDecayModeProducer = cms.InputTag("shrinkingConePFTauDecayModeProducer")
)


process.shrinkingConePFTauDecayModeIndexProducerPF = cms.EDProducer("PFRecoTauDecayModeIndexProducer",
    PFTauDecayModeProducer = cms.InputTag("shrinkingConePFTauDecayModeProducerPF"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    PFTauProducer = cms.InputTag("pfTausPF")
)


process.shrinkingConePFTauDecayModeProducer = cms.EDProducer("PFRecoTauDecayModeDeterminator",
    mergeByBestMatch = cms.bool(True),
    refitTracks = cms.bool(False),
    maxPiZeroMass = cms.double(0.2),
    mergeLowPtPhotonsFirst = cms.bool(True),
    setMergedPi0Mass = cms.bool(True),
    setChargedPionMass = cms.bool(True),
    filterPhotons = cms.bool(True),
    minPtFractionSinglePhotons = cms.double(0.1),
    minPtFractionPiZeroes = cms.double(0.15),
    maxNbrOfIterations = cms.int32(10),
    filterTwoProngs = cms.bool(True),
    minPtFractionForSecondProng = cms.double(0.1),
    maxDistance = cms.double(0.01),
    setPi0Mass = cms.bool(True),
    maxPhotonsToMerge = cms.uint32(2),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducer")
)


process.shrinkingConePFTauDecayModeProducerPF = cms.EDProducer("PFRecoTauDecayModeDeterminator",
    maxPiZeroMass = cms.double(0.2),
    refitTracks = cms.bool(False),
    maxNbrOfIterations = cms.int32(10),
    PFTauProducer = cms.InputTag("pfTausPF"),
    mergeLowPtPhotonsFirst = cms.bool(True),
    setMergedPi0Mass = cms.bool(True),
    setChargedPionMass = cms.bool(True),
    filterPhotons = cms.bool(True),
    minPtFractionPiZeroes = cms.double(0.15),
    maxPhotonsToMerge = cms.uint32(2),
    filterTwoProngs = cms.bool(True),
    mergeByBestMatch = cms.bool(True),
    minPtFractionForSecondProng = cms.double(0.1),
    maxDistance = cms.double(0.01),
    setPi0Mass = cms.bool(True),
    minPtFractionSinglePhotons = cms.double(0.1)
)


process.shrinkingConePFTauDiscriminationAgainstElectron = cms.EDProducer("PFRecoTauDiscriminationAgainstElectron",
    ApplyCut_ElectronPreID_2D = cms.bool(False),
    ElecPreID0_HOverPLead_minValue = cms.double(0.05),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducer"),
    ApplyCut_ElectronPreID = cms.bool(False),
    ApplyCut_HcalTotOverPLead = cms.bool(False),
    EOverPLead_minValue = cms.double(0.8),
    ElecPreID1_EOverPLead_maxValue = cms.double(0.8),
    HcalMaxOverPLead_minValue = cms.double(0.1),
    ApplyCut_EmFraction = cms.bool(False),
    EmFraction_maxValue = cms.double(0.9),
    ApplyCut_PFElectronMVA = cms.bool(True),
    PFElectronMVA_maxValue = cms.double(-0.1),
    ApplyCut_HcalMaxOverPLead = cms.bool(False),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFinding")
        )
    ),
    Hcal3x3OverPLead_minValue = cms.double(0.1),
    ElecPreID1_HOverPLead_minValue = cms.double(0.15),
    ElecPreID0_EOverPLead_maxValue = cms.double(0.95),
    BremsRecoveryEOverPLead_minValue = cms.double(0.8),
    ApplyCut_EcalCrackCut = cms.bool(False),
    EOverPLead_maxValue = cms.double(1.8),
    HcalTotOverPLead_minValue = cms.double(0.1),
    ApplyCut_BremsRecoveryEOverPLead = cms.bool(False),
    ApplyCut_Hcal3x3OverPLead = cms.bool(False),
    ApplyCut_EOverPLead = cms.bool(False),
    BremsRecoveryEOverPLead_maxValue = cms.double(1.8)
)


process.shrinkingConePFTauDiscriminationAgainstElectronPF = cms.EDProducer("PFRecoTauDiscriminationAgainstElectron",
    ApplyCut_ElectronPreID_2D = cms.bool(False),
    ElecPreID0_HOverPLead_minValue = cms.double(0.05),
    PFTauProducer = cms.InputTag("pfTausPF"),
    ApplyCut_ElectronPreID = cms.bool(False),
    ApplyCut_HcalTotOverPLead = cms.bool(False),
    EOverPLead_minValue = cms.double(0.8),
    ElecPreID1_EOverPLead_maxValue = cms.double(0.8),
    HcalMaxOverPLead_minValue = cms.double(0.1),
    ApplyCut_EmFraction = cms.bool(False),
    EmFraction_maxValue = cms.double(0.9),
    ApplyCut_PFElectronMVA = cms.bool(True),
    PFElectronMVA_maxValue = cms.double(-0.1),
    ApplyCut_HcalMaxOverPLead = cms.bool(False),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFindingPF")
        )
    ),
    Hcal3x3OverPLead_minValue = cms.double(0.1),
    ElecPreID1_HOverPLead_minValue = cms.double(0.15),
    ElecPreID0_EOverPLead_maxValue = cms.double(0.95),
    BremsRecoveryEOverPLead_minValue = cms.double(0.8),
    ApplyCut_EcalCrackCut = cms.bool(False),
    EOverPLead_maxValue = cms.double(1.8),
    HcalTotOverPLead_minValue = cms.double(0.1),
    ApplyCut_BremsRecoveryEOverPLead = cms.bool(False),
    ApplyCut_Hcal3x3OverPLead = cms.bool(False),
    ApplyCut_EOverPLead = cms.bool(False),
    BremsRecoveryEOverPLead_maxValue = cms.double(1.8)
)


process.shrinkingConePFTauDiscriminationAgainstMuon = cms.EDProducer("PFRecoTauDiscriminationAgainstMuon",
    a = cms.double(0.5),
    c = cms.double(0.0),
    b = cms.double(0.5),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFinding")
        )
    ),
    discriminatorOption = cms.string('noSegMatch')
)


process.shrinkingConePFTauDiscriminationAgainstMuonPF = cms.EDProducer("PFRecoTauDiscriminationAgainstMuon",
    a = cms.double(0.5),
    c = cms.double(0.0),
    b = cms.double(0.5),
    PFTauProducer = cms.InputTag("pfTausPF"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFindingPF")
        )
    ),
    discriminatorOption = cms.string('noSegMatch')
)


process.shrinkingConePFTauDiscriminationByECALIsolation = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducer"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFinding")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(False),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.shrinkingConePFTauDiscriminationByECALIsolationPF = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("pfTausPF"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFindingPF")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(False),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.shrinkingConePFTauDiscriminationByECALIsolationUsingLeadingPion = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducer"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadPion = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFinding")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(False),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.shrinkingConePFTauDiscriminationByECALIsolationUsingLeadingPionPF = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("pfTausPF"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadPion = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFindingPF")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(False),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.shrinkingConePFTauDiscriminationByIsolation = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducer"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFinding")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.shrinkingConePFTauDiscriminationByIsolationPF = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("pfTausPF"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFindingPF")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.shrinkingConePFTauDiscriminationByIsolationUsingLeadingPion = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducer"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadPion = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFinding")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.shrinkingConePFTauDiscriminationByIsolationUsingLeadingPionPF = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("pfTausPF"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadPion = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFindingPF")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.shrinkingConePFTauDiscriminationByLeadingPionPtCut = cms.EDProducer("PFRecoTauDiscriminationByLeadingObjectPtCut",
    MinPtLeadingObject = cms.double(5.0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducer"),
    UseOnlyChargedHadrons = cms.bool(False)
)


process.shrinkingConePFTauDiscriminationByLeadingPionPtCutPF = cms.EDProducer("PFRecoTauDiscriminationByLeadingObjectPtCut",
    MinPtLeadingObject = cms.double(5.0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    UseOnlyChargedHadrons = cms.bool(False),
    PFTauProducer = cms.InputTag("pfTausPF")
)


process.shrinkingConePFTauDiscriminationByLeadingTrackFinding = cms.EDProducer("PFRecoTauDiscriminationByLeadingObjectPtCut",
    MinPtLeadingObject = cms.double(0.0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducer"),
    UseOnlyChargedHadrons = cms.bool(True)
)


process.shrinkingConePFTauDiscriminationByLeadingTrackFindingPF = cms.EDProducer("PFRecoTauDiscriminationByLeadingObjectPtCut",
    MinPtLeadingObject = cms.double(0.0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    UseOnlyChargedHadrons = cms.bool(True),
    PFTauProducer = cms.InputTag("pfTausPF")
)


process.shrinkingConePFTauDiscriminationByLeadingTrackPtCut = cms.EDProducer("PFRecoTauDiscriminationByLeadingObjectPtCut",
    MinPtLeadingObject = cms.double(5.0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducer"),
    UseOnlyChargedHadrons = cms.bool(True)
)


process.shrinkingConePFTauDiscriminationByLeadingTrackPtCutPF = cms.EDProducer("PFRecoTauDiscriminationByLeadingObjectPtCut",
    MinPtLeadingObject = cms.double(5.0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    UseOnlyChargedHadrons = cms.bool(True),
    PFTauProducer = cms.InputTag("pfTausPF")
)


process.shrinkingConePFTauDiscriminationByTaNC = cms.EDProducer("PFTauMVADiscriminator",
    RemapOutput = cms.bool(True),
    pfTauDecayModeSrc = cms.InputTag("shrinkingConePFTauDecayModeProducer"),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducer"),
    computers = cms.VPSet(cms.PSet(
        applyIsolation = cms.bool(False),
        cut = cms.double(-10.0),
        computerName = cms.string('OneProngNoPiZero'),
        decayModeIndices = cms.vint32(0)
    ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(-10.0),
            computerName = cms.string('OneProngOnePiZero'),
            decayModeIndices = cms.vint32(1)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(-10.0),
            computerName = cms.string('OneProngTwoPiZero'),
            decayModeIndices = cms.vint32(2)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(-10.0),
            computerName = cms.string('ThreeProngNoPiZero'),
            decayModeIndices = cms.vint32(10)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(-10.0),
            computerName = cms.string('ThreeProngOnePiZero'),
            decayModeIndices = cms.vint32(11)
        )),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFinding")
        )
    ),
    prefailValue = cms.double(-2.0),
    dbLabel = cms.string(''),
    MakeBinaryDecision = cms.bool(False)
)


process.shrinkingConePFTauDiscriminationByTaNCPF = cms.EDProducer("PFTauMVADiscriminator",
    RemapOutput = cms.bool(True),
    pfTauDecayModeSrc = cms.InputTag("shrinkingConePFTauDecayModeProducerPF"),
    PFTauProducer = cms.InputTag("pfTausPF"),
    computers = cms.VPSet(cms.PSet(
        applyIsolation = cms.bool(False),
        cut = cms.double(-10.0),
        computerName = cms.string('OneProngNoPiZero'),
        decayModeIndices = cms.vint32(0)
    ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(-10.0),
            computerName = cms.string('OneProngOnePiZero'),
            decayModeIndices = cms.vint32(1)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(-10.0),
            computerName = cms.string('OneProngTwoPiZero'),
            decayModeIndices = cms.vint32(2)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(-10.0),
            computerName = cms.string('ThreeProngNoPiZero'),
            decayModeIndices = cms.vint32(10)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(-10.0),
            computerName = cms.string('ThreeProngOnePiZero'),
            decayModeIndices = cms.vint32(11)
        )),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFindingPF")
        )
    ),
    prefailValue = cms.double(-2.0),
    dbLabel = cms.string(''),
    MakeBinaryDecision = cms.bool(False)
)


process.shrinkingConePFTauDiscriminationByTaNCfrHalfPercent = cms.EDProducer("PFTauDecayModeCutMultiplexer",
    RemapOutput = cms.bool(True),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducer"),
    computers = cms.VPSet(cms.PSet(
        applyIsolation = cms.bool(False),
        cut = cms.double(0.9087875),
        computerName = cms.string('OneProngNoPiZero'),
        decayModeIndices = cms.vint32(0)
    ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.8707145),
            computerName = cms.string('OneProngOnePiZero'),
            decayModeIndices = cms.vint32(1)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.921793),
            computerName = cms.string('OneProngTwoPiZero'),
            decayModeIndices = cms.vint32(2)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.9451915),
            computerName = cms.string('ThreeProngNoPiZero'),
            decayModeIndices = cms.vint32(10)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.9488565),
            computerName = cms.string('ThreeProngOnePiZero'),
            decayModeIndices = cms.vint32(11)
        )),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFinding")
        )
    ),
    PFTauDiscriminantToMultiplex = cms.InputTag("shrinkingConePFTauDiscriminationByTaNC"),
    PFTauDecayModeSrc = cms.InputTag("shrinkingConePFTauDecayModeIndexProducer"),
    MakeBinaryDecision = cms.bool(True)
)


process.shrinkingConePFTauDiscriminationByTaNCfrHalfPercentPF = cms.EDProducer("PFTauDecayModeCutMultiplexer",
    RemapOutput = cms.bool(True),
    PFTauProducer = cms.InputTag("pfTausPF"),
    computers = cms.VPSet(cms.PSet(
        applyIsolation = cms.bool(False),
        cut = cms.double(0.9087875),
        computerName = cms.string('OneProngNoPiZero'),
        decayModeIndices = cms.vint32(0)
    ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.8707145),
            computerName = cms.string('OneProngOnePiZero'),
            decayModeIndices = cms.vint32(1)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.921793),
            computerName = cms.string('OneProngTwoPiZero'),
            decayModeIndices = cms.vint32(2)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.9451915),
            computerName = cms.string('ThreeProngNoPiZero'),
            decayModeIndices = cms.vint32(10)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.9488565),
            computerName = cms.string('ThreeProngOnePiZero'),
            decayModeIndices = cms.vint32(11)
        )),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFindingPF")
        )
    ),
    PFTauDiscriminantToMultiplex = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCPF"),
    PFTauDecayModeSrc = cms.InputTag("shrinkingConePFTauDecayModeIndexProducerPF"),
    MakeBinaryDecision = cms.bool(True)
)


process.shrinkingConePFTauDiscriminationByTaNCfrOnePercent = cms.EDProducer("PFTauDecayModeCutMultiplexer",
    RemapOutput = cms.bool(True),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducer"),
    computers = cms.VPSet(cms.PSet(
        applyIsolation = cms.bool(False),
        cut = cms.double(0.7649845),
        computerName = cms.string('OneProngNoPiZero'),
        decayModeIndices = cms.vint32(0)
    ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.699697),
            computerName = cms.string('OneProngOnePiZero'),
            decayModeIndices = cms.vint32(1)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.772231),
            computerName = cms.string('OneProngTwoPiZero'),
            decayModeIndices = cms.vint32(2)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.905071),
            computerName = cms.string('ThreeProngNoPiZero'),
            decayModeIndices = cms.vint32(10)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.9238995),
            computerName = cms.string('ThreeProngOnePiZero'),
            decayModeIndices = cms.vint32(11)
        )),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFinding")
        )
    ),
    PFTauDiscriminantToMultiplex = cms.InputTag("shrinkingConePFTauDiscriminationByTaNC"),
    PFTauDecayModeSrc = cms.InputTag("shrinkingConePFTauDecayModeIndexProducer"),
    MakeBinaryDecision = cms.bool(True)
)


process.shrinkingConePFTauDiscriminationByTaNCfrOnePercentPF = cms.EDProducer("PFTauDecayModeCutMultiplexer",
    RemapOutput = cms.bool(True),
    PFTauProducer = cms.InputTag("pfTausPF"),
    computers = cms.VPSet(cms.PSet(
        applyIsolation = cms.bool(False),
        cut = cms.double(0.7649845),
        computerName = cms.string('OneProngNoPiZero'),
        decayModeIndices = cms.vint32(0)
    ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.699697),
            computerName = cms.string('OneProngOnePiZero'),
            decayModeIndices = cms.vint32(1)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.772231),
            computerName = cms.string('OneProngTwoPiZero'),
            decayModeIndices = cms.vint32(2)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.905071),
            computerName = cms.string('ThreeProngNoPiZero'),
            decayModeIndices = cms.vint32(10)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.9238995),
            computerName = cms.string('ThreeProngOnePiZero'),
            decayModeIndices = cms.vint32(11)
        )),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFindingPF")
        )
    ),
    PFTauDiscriminantToMultiplex = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCPF"),
    PFTauDecayModeSrc = cms.InputTag("shrinkingConePFTauDecayModeIndexProducerPF"),
    MakeBinaryDecision = cms.bool(True)
)


process.shrinkingConePFTauDiscriminationByTaNCfrQuarterPercent = cms.EDProducer("PFTauDecayModeCutMultiplexer",
    RemapOutput = cms.bool(True),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducer"),
    computers = cms.VPSet(cms.PSet(
        applyIsolation = cms.bool(False),
        cut = cms.double(0.9539685),
        computerName = cms.string('OneProngNoPiZero'),
        decayModeIndices = cms.vint32(0)
    ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.940843),
            computerName = cms.string('OneProngOnePiZero'),
            decayModeIndices = cms.vint32(1)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.9645195),
            computerName = cms.string('OneProngTwoPiZero'),
            decayModeIndices = cms.vint32(2)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.960407),
            computerName = cms.string('ThreeProngNoPiZero'),
            decayModeIndices = cms.vint32(10)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.994065),
            computerName = cms.string('ThreeProngOnePiZero'),
            decayModeIndices = cms.vint32(11)
        )),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFinding")
        )
    ),
    PFTauDiscriminantToMultiplex = cms.InputTag("shrinkingConePFTauDiscriminationByTaNC"),
    PFTauDecayModeSrc = cms.InputTag("shrinkingConePFTauDecayModeIndexProducer"),
    MakeBinaryDecision = cms.bool(True)
)


process.shrinkingConePFTauDiscriminationByTaNCfrQuarterPercentPF = cms.EDProducer("PFTauDecayModeCutMultiplexer",
    RemapOutput = cms.bool(True),
    PFTauProducer = cms.InputTag("pfTausPF"),
    computers = cms.VPSet(cms.PSet(
        applyIsolation = cms.bool(False),
        cut = cms.double(0.9539685),
        computerName = cms.string('OneProngNoPiZero'),
        decayModeIndices = cms.vint32(0)
    ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.940843),
            computerName = cms.string('OneProngOnePiZero'),
            decayModeIndices = cms.vint32(1)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.9645195),
            computerName = cms.string('OneProngTwoPiZero'),
            decayModeIndices = cms.vint32(2)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.960407),
            computerName = cms.string('ThreeProngNoPiZero'),
            decayModeIndices = cms.vint32(10)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.994065),
            computerName = cms.string('ThreeProngOnePiZero'),
            decayModeIndices = cms.vint32(11)
        )),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFindingPF")
        )
    ),
    PFTauDiscriminantToMultiplex = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCPF"),
    PFTauDecayModeSrc = cms.InputTag("shrinkingConePFTauDecayModeIndexProducerPF"),
    MakeBinaryDecision = cms.bool(True)
)


process.shrinkingConePFTauDiscriminationByTaNCfrTenthPercent = cms.EDProducer("PFTauDecayModeCutMultiplexer",
    RemapOutput = cms.bool(True),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducer"),
    computers = cms.VPSet(cms.PSet(
        applyIsolation = cms.bool(False),
        cut = cms.double(0.959384),
        computerName = cms.string('OneProngNoPiZero'),
        decayModeIndices = cms.vint32(0)
    ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.991353),
            computerName = cms.string('OneProngOnePiZero'),
            decayModeIndices = cms.vint32(1)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.9712775),
            computerName = cms.string('OneProngTwoPiZero'),
            decayModeIndices = cms.vint32(2)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.9875015),
            computerName = cms.string('ThreeProngNoPiZero'),
            decayModeIndices = cms.vint32(10)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(1.0234655),
            computerName = cms.string('ThreeProngOnePiZero'),
            decayModeIndices = cms.vint32(11)
        )),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFinding")
        )
    ),
    PFTauDiscriminantToMultiplex = cms.InputTag("shrinkingConePFTauDiscriminationByTaNC"),
    PFTauDecayModeSrc = cms.InputTag("shrinkingConePFTauDecayModeIndexProducer"),
    MakeBinaryDecision = cms.bool(True)
)


process.shrinkingConePFTauDiscriminationByTaNCfrTenthPercentPF = cms.EDProducer("PFTauDecayModeCutMultiplexer",
    RemapOutput = cms.bool(True),
    PFTauProducer = cms.InputTag("pfTausPF"),
    computers = cms.VPSet(cms.PSet(
        applyIsolation = cms.bool(False),
        cut = cms.double(0.959384),
        computerName = cms.string('OneProngNoPiZero'),
        decayModeIndices = cms.vint32(0)
    ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.991353),
            computerName = cms.string('OneProngOnePiZero'),
            decayModeIndices = cms.vint32(1)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.9712775),
            computerName = cms.string('OneProngTwoPiZero'),
            decayModeIndices = cms.vint32(2)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(0.9875015),
            computerName = cms.string('ThreeProngNoPiZero'),
            decayModeIndices = cms.vint32(10)
        ), 
        cms.PSet(
            applyIsolation = cms.bool(False),
            cut = cms.double(1.0234655),
            computerName = cms.string('ThreeProngOnePiZero'),
            decayModeIndices = cms.vint32(11)
        )),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFindingPF")
        )
    ),
    PFTauDiscriminantToMultiplex = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCPF"),
    PFTauDecayModeSrc = cms.InputTag("shrinkingConePFTauDecayModeIndexProducerPF"),
    MakeBinaryDecision = cms.bool(True)
)


process.shrinkingConePFTauDiscriminationByTrackIsolation = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducer"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFinding")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.shrinkingConePFTauDiscriminationByTrackIsolationPF = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("pfTausPF"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFindingPF")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.shrinkingConePFTauDiscriminationByTrackIsolationUsingLeadingPion = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("shrinkingConePFTauProducer"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadPion = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFinding")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.shrinkingConePFTauDiscriminationByTrackIsolationUsingLeadingPionPF = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    PFTauProducer = cms.InputTag("pfTausPF"),
    maximumSumPtCut = cms.double(6.0),
    relativeSumPtCut = cms.double(0.0),
    maximumOccupancy = cms.uint32(0),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadPion = cms.PSet(
            cut = cms.double(0.5),
            Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFindingPF")
        )
    ),
    applyOccupancyCut = cms.bool(True),
    applySumPtCut = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(8),
            minTrackPt = cms.double(1.0),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(1.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        ),
        signalQualityCuts = cms.PSet(
            minTrackHits = cms.uint32(3),
            minTrackPt = cms.double(0.5),
            maxTrackChi2 = cms.double(100.0),
            minTrackPixelHits = cms.uint32(0),
            minGammaEt = cms.double(0.5),
            useTracksInsteadOfPFHadrons = cms.bool(False),
            maxDeltaZ = cms.double(0.2),
            maxTransverseImpactParameter = cms.double(0.03)
        )
    ),
    applyRelativeSumPtCut = cms.bool(False)
)


process.shrinkingConePFTauProducer = cms.EDProducer("PFRecoTauProducer",
    Rphi = cms.double(2.0),
    LeadTrack_minPt = cms.double(0.0),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    ECALSignalConeSizeFormula = cms.string('0.15'),
    TrackerIsolConeMetric = cms.string('DR'),
    TrackerSignalConeMetric = cms.string('DR'),
    EcalStripSumE_deltaPhiOverQ_minValue = cms.double(-0.1),
    MaxEtInEllipse = cms.double(2.0),
    MatchingConeMetric = cms.string('DR'),
    TrackerSignalConeSizeFormula = cms.string('5/ET'),
    MatchingConeSizeFormula = cms.string('0.1'),
    TrackerIsolConeSize_min = cms.double(0.0),
    MatchingConeSize_min = cms.double(0.0),
    ElectronPreIDProducer = cms.InputTag("elecpreid"),
    ChargedHadrCandLeadChargedHadrCand_tksmaxDZ = cms.double(0.2),
    TrackerIsolConeSize_max = cms.double(0.6),
    TrackerSignalConeSize_max = cms.double(0.15),
    HCALIsolConeMetric = cms.string('DR'),
    AddEllipseGammas = cms.bool(False),
    maximumForElectrionPreIDOutput = cms.double(-0.1),
    TrackerSignalConeSize_min = cms.double(0.07),
    ECALIsolConeSize_max = cms.double(0.6),
    HCALIsolConeSizeFormula = cms.string('0.50'),
    Track_IsolAnnulus_minNhits = cms.uint32(3),
    smearedPVsigmaZ = cms.double(0.005),
    AreaMetric_recoElements_maxabsEta = cms.double(2.5),
    HCALSignalConeMetric = cms.string('DR'),
    ElecPreIDLeadTkMatch_maxDR = cms.double(0.01),
    ChargedHadrCand_IsolAnnulus_minNhits = cms.uint32(0),
    PFTauTagInfoProducer = cms.InputTag("pfRecoTauTagInfoProducer"),
    ECALIsolConeMetric = cms.string('DR'),
    ECALIsolConeSizeFormula = cms.string('0.50'),
    UseChargedHadrCandLeadChargedHadrCand_tksDZconstraint = cms.bool(True),
    Algorithm = cms.string('ConeBased'),
    JetPtMin = cms.double(0.0),
    ECALSignalConeMetric = cms.string('DR'),
    EcalStripSumE_deltaPhiOverQ_maxValue = cms.double(0.5),
    HCALSignalConeSize_max = cms.double(0.5),
    ECALSignalConeSize_min = cms.double(0.0),
    EcalStripSumE_minClusEnergy = cms.double(0.1),
    EcalStripSumE_deltaEta = cms.double(0.03),
    TrackerIsolConeSizeFormula = cms.string('0.50'),
    LeadPFCand_minPt = cms.double(5.0),
    HCALSignalConeSize_min = cms.double(0.0),
    ECALSignalConeSize_max = cms.double(0.5),
    HCALSignalConeSizeFormula = cms.string('0.15'),
    TrackLeadTrack_maxDZ = cms.double(0.2),
    DataType = cms.string('AOD'),
    ECALIsolConeSize_min = cms.double(0.0),
    UseTrackLeadTrackDZconstraint = cms.bool(True),
    smearedPVsigmaY = cms.double(0.0015),
    HCALIsolConeSize_max = cms.double(0.6),
    smearedPVsigmaX = cms.double(0.0015),
    MatchingConeSize_max = cms.double(0.6),
    HCALIsolConeSize_min = cms.double(0.0)
)


process.shrinkingConePFTauProducerPF = cms.EDProducer("PFRecoTauProducer",
    Rphi = cms.double(2.0),
    LeadTrack_minPt = cms.double(0.0),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    ECALSignalConeSizeFormula = cms.string('0.15'),
    TrackerIsolConeMetric = cms.string('DR'),
    TrackerSignalConeMetric = cms.string('DR'),
    EcalStripSumE_deltaPhiOverQ_minValue = cms.double(-0.1),
    smearedPVsigmaY = cms.double(0.0015),
    DataType = cms.string('AOD'),
    MatchingConeMetric = cms.string('DR'),
    TrackerSignalConeSizeFormula = cms.string('5/ET'),
    MatchingConeSizeFormula = cms.string('0.1'),
    TrackerIsolConeSize_min = cms.double(0.0),
    MatchingConeSize_min = cms.double(0.0),
    ElectronPreIDProducer = cms.InputTag("elecpreid"),
    ChargedHadrCandLeadChargedHadrCand_tksmaxDZ = cms.double(0.2),
    TrackerIsolConeSize_max = cms.double(0.6),
    TrackerSignalConeSize_max = cms.double(0.15),
    HCALIsolConeMetric = cms.string('DR'),
    AddEllipseGammas = cms.bool(False),
    maximumForElectrionPreIDOutput = cms.double(-0.1),
    TrackerSignalConeSize_min = cms.double(0.07),
    JetPtMin = cms.double(0.0),
    HCALIsolConeSizeFormula = cms.string('0.50'),
    AreaMetric_recoElements_maxabsEta = cms.double(2.5),
    HCALIsolConeSize_max = cms.double(0.6),
    Track_IsolAnnulus_minNhits = cms.uint32(3),
    HCALSignalConeMetric = cms.string('DR'),
    ElecPreIDLeadTkMatch_maxDR = cms.double(0.01),
    PFTauTagInfoProducer = cms.InputTag("pfRecoTauTagInfoProducerPF"),
    ECALIsolConeMetric = cms.string('DR'),
    ECALIsolConeSizeFormula = cms.string('0.50'),
    UseChargedHadrCandLeadChargedHadrCand_tksDZconstraint = cms.bool(True),
    Algorithm = cms.string('ConeBased'),
    ECALIsolConeSize_max = cms.double(0.6),
    ECALSignalConeMetric = cms.string('DR'),
    EcalStripSumE_deltaPhiOverQ_maxValue = cms.double(0.5),
    HCALSignalConeSize_max = cms.double(0.5),
    ECALSignalConeSize_min = cms.double(0.0),
    EcalStripSumE_minClusEnergy = cms.double(0.1),
    EcalStripSumE_deltaEta = cms.double(0.03),
    TrackerIsolConeSizeFormula = cms.string('0.50'),
    LeadPFCand_minPt = cms.double(5.0),
    HCALSignalConeSize_min = cms.double(0.0),
    ECALSignalConeSize_max = cms.double(0.5),
    HCALSignalConeSizeFormula = cms.string('0.15'),
    TrackLeadTrack_maxDZ = cms.double(0.2),
    ChargedHadrCand_IsolAnnulus_minNhits = cms.uint32(0),
    ECALIsolConeSize_min = cms.double(0.0),
    UseTrackLeadTrackDZconstraint = cms.bool(True),
    MaxEtInEllipse = cms.double(2.0),
    smearedPVsigmaX = cms.double(0.0015),
    smearedPVsigmaZ = cms.double(0.005),
    MatchingConeSize_max = cms.double(0.6),
    HCALIsolConeSize_min = cms.double(0.0)
)


process.simpleCutBasedElectronID = cms.EDProducer("EleIdCutBasedExtProducer",
    robust60relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.04, 0.04, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.025, 0.02, 0.02, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust80cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.07, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.03, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.06, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    electronQuality = cms.string('test'),
    algorithm = cms.string('eIDCB'),
    verticesCollection = cms.InputTag("offlineBeamSpot"),
    robust95relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.15, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.15, 2.0, 0.12, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0),
        endcap = cms.vdouble(0.07, 0.03, 0.7, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.08, 0.06, 0.05, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0)
    ),
    robust70relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.03, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.06, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.025, 0.025, 0.02, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust85relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.006, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.09, 0.08, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.04, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.05, 0.025, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust85cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.006, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.09, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.04, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.06, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    robust90cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.12, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.1, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.05, 0.03, 0.7, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.07, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    robust80relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.09, 0.07, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.03, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.04, 0.05, 0.025, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    electronIDType = cms.string('robust'),
    electronVersion = cms.string('V04'),
    src = cms.InputTag("gsfElectrons"),
    robust90relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.12, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.12, 0.09, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.05, 0.03, 0.7, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.06, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust60cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.03, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.02, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust70cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.03, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.04, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.03, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust95cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.15, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.15, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0),
        endcap = cms.vdouble(0.07, 0.03, 0.7, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.1, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0)
    )
)


process.simpleEleId60cIso = cms.EDProducer("EleIdCutBasedExtProducer",
    robust60relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.04, 0.04, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.025, 0.02, 0.02, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust80cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.07, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.03, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.06, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    algorithm = cms.string('eIDCB'),
    verticesCollection = cms.InputTag("offlineBeamSpot"),
    robust95relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.15, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.15, 2.0, 0.12, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0),
        endcap = cms.vdouble(0.07, 0.03, 0.7, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.08, 0.06, 0.05, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0)
    ),
    robust70relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.03, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.06, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.025, 0.025, 0.02, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    robust70cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.03, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.04, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.03, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust85cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.006, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.09, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.04, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.06, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust80relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.09, 0.07, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.03, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.04, 0.05, 0.025, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust90cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.12, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.1, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.05, 0.03, 0.7, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.07, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    electronQuality = cms.string('60cIso'),
    electronIDType = cms.string('robust'),
    electronVersion = cms.string('V04'),
    src = cms.InputTag("gsfElectrons"),
    robust90relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.12, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.12, 0.09, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.05, 0.03, 0.7, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.06, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust60cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.03, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.02, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust85relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.006, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.09, 0.08, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.04, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.05, 0.025, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust95cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.15, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.15, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0),
        endcap = cms.vdouble(0.07, 0.03, 0.7, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.1, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0)
    )
)


process.simpleEleId60relIso = cms.EDProducer("EleIdCutBasedExtProducer",
    robust60relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.04, 0.04, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.025, 0.02, 0.02, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust80cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.07, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.03, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.06, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    algorithm = cms.string('eIDCB'),
    verticesCollection = cms.InputTag("offlineBeamSpot"),
    robust95relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.15, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.15, 2.0, 0.12, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0),
        endcap = cms.vdouble(0.07, 0.03, 0.7, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.08, 0.06, 0.05, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0)
    ),
    robust70relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.03, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.06, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.025, 0.025, 0.02, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    robust70cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.03, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.04, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.03, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust85cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.006, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.09, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.04, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.06, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust80relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.09, 0.07, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.03, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.04, 0.05, 0.025, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust90cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.12, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.1, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.05, 0.03, 0.7, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.07, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    electronQuality = cms.string('60relIso'),
    electronIDType = cms.string('robust'),
    electronVersion = cms.string('V04'),
    src = cms.InputTag("gsfElectrons"),
    robust90relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.12, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.12, 0.09, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.05, 0.03, 0.7, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.06, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust60cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.03, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.02, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust85relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.006, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.09, 0.08, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.04, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.05, 0.025, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust95cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.15, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.15, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0),
        endcap = cms.vdouble(0.07, 0.03, 0.7, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.1, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0)
    )
)


process.simpleEleId70cIso = cms.EDProducer("EleIdCutBasedExtProducer",
    robust60relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.04, 0.04, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.025, 0.02, 0.02, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust80cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.07, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.03, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.06, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    algorithm = cms.string('eIDCB'),
    verticesCollection = cms.InputTag("offlineBeamSpot"),
    robust95relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.15, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.15, 2.0, 0.12, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0),
        endcap = cms.vdouble(0.07, 0.03, 0.7, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.08, 0.06, 0.05, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0)
    ),
    robust70relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.03, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.06, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.025, 0.025, 0.02, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    robust70cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.03, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.04, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.03, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust85cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.006, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.09, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.04, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.06, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust80relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.09, 0.07, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.03, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.04, 0.05, 0.025, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust90cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.12, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.1, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.05, 0.03, 0.7, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.07, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    electronQuality = cms.string('70cIso'),
    electronIDType = cms.string('robust'),
    electronVersion = cms.string('V04'),
    src = cms.InputTag("gsfElectrons"),
    robust90relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.12, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.12, 0.09, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.05, 0.03, 0.7, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.06, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust60cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.03, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.02, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust85relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.006, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.09, 0.08, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.04, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.05, 0.025, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust95cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.15, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.15, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0),
        endcap = cms.vdouble(0.07, 0.03, 0.7, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.1, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0)
    )
)


process.simpleEleId70relIso = cms.EDProducer("EleIdCutBasedExtProducer",
    robust60relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.04, 0.04, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.025, 0.02, 0.02, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust80cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.07, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.03, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.06, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    algorithm = cms.string('eIDCB'),
    verticesCollection = cms.InputTag("offlineBeamSpot"),
    robust95relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.15, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.15, 2.0, 0.12, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0),
        endcap = cms.vdouble(0.07, 0.03, 0.7, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.08, 0.06, 0.05, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0)
    ),
    robust70relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.03, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.06, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.025, 0.025, 0.02, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    robust70cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.03, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.04, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.03, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust85cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.006, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.09, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.04, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.06, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust80relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.09, 0.07, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.03, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.04, 0.05, 0.025, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust90cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.12, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.1, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.05, 0.03, 0.7, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.07, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    electronQuality = cms.string('70relIso'),
    electronIDType = cms.string('robust'),
    electronVersion = cms.string('V04'),
    src = cms.InputTag("gsfElectrons"),
    robust90relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.12, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.12, 0.09, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.05, 0.03, 0.7, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.06, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust60cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.03, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.02, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust85relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.006, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.09, 0.08, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.04, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.05, 0.025, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust95cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.15, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.15, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0),
        endcap = cms.vdouble(0.07, 0.03, 0.7, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.1, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0)
    )
)


process.simpleEleId80cIso = cms.EDProducer("EleIdCutBasedExtProducer",
    robust60relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.04, 0.04, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.025, 0.02, 0.02, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust80cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.07, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.03, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.06, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    algorithm = cms.string('eIDCB'),
    verticesCollection = cms.InputTag("offlineBeamSpot"),
    robust95relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.15, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.15, 2.0, 0.12, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0),
        endcap = cms.vdouble(0.07, 0.03, 0.7, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.08, 0.06, 0.05, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0)
    ),
    robust70relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.03, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.06, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.025, 0.025, 0.02, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    robust70cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.03, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.04, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.03, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust85cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.006, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.09, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.04, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.06, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust80relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.09, 0.07, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.03, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.04, 0.05, 0.025, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust90cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.12, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.1, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.05, 0.03, 0.7, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.07, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    electronQuality = cms.string('80cIso'),
    electronIDType = cms.string('robust'),
    electronVersion = cms.string('V04'),
    src = cms.InputTag("gsfElectrons"),
    robust90relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.12, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.12, 0.09, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.05, 0.03, 0.7, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.06, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust60cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.03, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.02, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust85relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.006, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.09, 0.08, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.04, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.05, 0.025, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust95cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.15, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.15, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0),
        endcap = cms.vdouble(0.07, 0.03, 0.7, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.1, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0)
    )
)


process.simpleEleId80relIso = cms.EDProducer("EleIdCutBasedExtProducer",
    robust60relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.04, 0.04, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.025, 0.02, 0.02, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust80cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.07, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.03, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.06, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    algorithm = cms.string('eIDCB'),
    verticesCollection = cms.InputTag("offlineBeamSpot"),
    robust95relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.15, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.15, 2.0, 0.12, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0),
        endcap = cms.vdouble(0.07, 0.03, 0.7, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.08, 0.06, 0.05, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0)
    ),
    robust70relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.03, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.06, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.025, 0.025, 0.02, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    robust70cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.03, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.04, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.03, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust85cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.006, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.09, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.04, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.06, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust80relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.09, 0.07, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.03, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.04, 0.05, 0.025, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust90cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.12, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.1, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.05, 0.03, 0.7, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.07, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    electronQuality = cms.string('80relIso'),
    electronIDType = cms.string('robust'),
    electronVersion = cms.string('V04'),
    src = cms.InputTag("gsfElectrons"),
    robust90relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.12, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.12, 0.09, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.05, 0.03, 0.7, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.06, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust60cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.03, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.02, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust85relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.006, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.09, 0.08, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.04, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.05, 0.025, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust95cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.15, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.15, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0),
        endcap = cms.vdouble(0.07, 0.03, 0.7, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.1, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0)
    )
)


process.simpleEleId85cIso = cms.EDProducer("EleIdCutBasedExtProducer",
    robust60relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.04, 0.04, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.025, 0.02, 0.02, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust80cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.07, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.03, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.06, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    algorithm = cms.string('eIDCB'),
    verticesCollection = cms.InputTag("offlineBeamSpot"),
    robust95relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.15, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.15, 2.0, 0.12, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0),
        endcap = cms.vdouble(0.07, 0.03, 0.7, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.08, 0.06, 0.05, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0)
    ),
    robust70relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.03, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.06, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.025, 0.025, 0.02, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    robust70cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.03, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.04, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.03, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust85cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.006, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.09, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.04, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.06, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust80relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.09, 0.07, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.03, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.04, 0.05, 0.025, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust90cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.12, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.1, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.05, 0.03, 0.7, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.07, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    electronQuality = cms.string('85cIso'),
    electronIDType = cms.string('robust'),
    electronVersion = cms.string('V04'),
    src = cms.InputTag("gsfElectrons"),
    robust90relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.12, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.12, 0.09, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.05, 0.03, 0.7, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.06, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust60cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.03, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.02, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust85relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.006, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.09, 0.08, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.04, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.05, 0.025, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust95cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.15, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.15, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0),
        endcap = cms.vdouble(0.07, 0.03, 0.7, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.1, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0)
    )
)


process.simpleEleId85relIso = cms.EDProducer("EleIdCutBasedExtProducer",
    robust60relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.04, 0.04, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.025, 0.02, 0.02, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust80cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.07, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.03, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.06, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    algorithm = cms.string('eIDCB'),
    verticesCollection = cms.InputTag("offlineBeamSpot"),
    robust95relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.15, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.15, 2.0, 0.12, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0),
        endcap = cms.vdouble(0.07, 0.03, 0.7, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.08, 0.06, 0.05, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0)
    ),
    robust70relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.03, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.06, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.025, 0.025, 0.02, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    robust70cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.03, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.04, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.03, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust85cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.006, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.09, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.04, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.06, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust80relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.09, 0.07, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.03, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.04, 0.05, 0.025, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust90cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.12, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.1, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.05, 0.03, 0.7, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.07, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    electronQuality = cms.string('85relIso'),
    electronIDType = cms.string('robust'),
    electronVersion = cms.string('V04'),
    src = cms.InputTag("gsfElectrons"),
    robust90relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.12, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.12, 0.09, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.05, 0.03, 0.7, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.06, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust60cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.03, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.02, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust85relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.006, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.09, 0.08, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.04, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.05, 0.025, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust95cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.15, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.15, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0),
        endcap = cms.vdouble(0.07, 0.03, 0.7, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.1, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0)
    )
)


process.simpleEleId90cIso = cms.EDProducer("EleIdCutBasedExtProducer",
    robust60relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.04, 0.04, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.025, 0.02, 0.02, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust80cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.07, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.03, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.06, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    algorithm = cms.string('eIDCB'),
    verticesCollection = cms.InputTag("offlineBeamSpot"),
    robust95relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.15, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.15, 2.0, 0.12, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0),
        endcap = cms.vdouble(0.07, 0.03, 0.7, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.08, 0.06, 0.05, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0)
    ),
    robust70relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.03, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.06, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.025, 0.025, 0.02, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    robust70cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.03, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.04, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.03, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust85cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.006, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.09, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.04, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.06, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust80relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.09, 0.07, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.03, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.04, 0.05, 0.025, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust90cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.12, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.1, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.05, 0.03, 0.7, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.07, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    electronQuality = cms.string('90cIso'),
    electronIDType = cms.string('robust'),
    electronVersion = cms.string('V04'),
    src = cms.InputTag("gsfElectrons"),
    robust90relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.12, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.12, 0.09, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.05, 0.03, 0.7, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.06, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust60cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.03, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.02, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust85relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.006, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.09, 0.08, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.04, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.05, 0.025, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust95cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.15, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.15, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0),
        endcap = cms.vdouble(0.07, 0.03, 0.7, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.1, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0)
    )
)


process.simpleEleId90relIso = cms.EDProducer("EleIdCutBasedExtProducer",
    robust60relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.04, 0.04, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.025, 0.02, 0.02, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust80cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.07, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.03, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.06, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    algorithm = cms.string('eIDCB'),
    verticesCollection = cms.InputTag("offlineBeamSpot"),
    robust95relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.15, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.15, 2.0, 0.12, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0),
        endcap = cms.vdouble(0.07, 0.03, 0.7, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.08, 0.06, 0.05, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0)
    ),
    robust70relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.03, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.06, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.025, 0.025, 0.02, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    robust70cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.03, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.04, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.03, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust85cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.006, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.09, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.04, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.06, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust80relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.09, 0.07, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.03, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.04, 0.05, 0.025, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust90cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.12, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.1, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.05, 0.03, 0.7, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.07, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    electronQuality = cms.string('90relIso'),
    electronIDType = cms.string('robust'),
    electronVersion = cms.string('V04'),
    src = cms.InputTag("gsfElectrons"),
    robust90relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.12, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.12, 0.09, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.05, 0.03, 0.7, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.06, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust60cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.03, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.02, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust85relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.006, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.09, 0.08, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.04, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.05, 0.025, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust95cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.15, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.15, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0),
        endcap = cms.vdouble(0.07, 0.03, 0.7, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.1, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0)
    )
)


process.simpleEleId95cIso = cms.EDProducer("EleIdCutBasedExtProducer",
    robust60relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.04, 0.04, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.025, 0.02, 0.02, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust80cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.07, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.03, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.06, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    algorithm = cms.string('eIDCB'),
    verticesCollection = cms.InputTag("offlineBeamSpot"),
    robust95relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.15, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.15, 2.0, 0.12, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0),
        endcap = cms.vdouble(0.07, 0.03, 0.7, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.08, 0.06, 0.05, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0)
    ),
    robust70relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.03, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.06, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.025, 0.025, 0.02, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    robust70cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.03, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.04, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.03, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust85cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.006, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.09, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.04, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.06, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust80relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.09, 0.07, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.03, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.04, 0.05, 0.025, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust90cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.12, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.1, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.05, 0.03, 0.7, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.07, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    electronQuality = cms.string('95cIso'),
    electronIDType = cms.string('robust'),
    electronVersion = cms.string('V04'),
    src = cms.InputTag("gsfElectrons"),
    robust90relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.12, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.12, 0.09, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.05, 0.03, 0.7, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.06, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust60cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.03, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.02, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust85relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.006, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.09, 0.08, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.04, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.05, 0.025, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust95cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.15, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.15, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0),
        endcap = cms.vdouble(0.07, 0.03, 0.7, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.1, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0)
    )
)


process.simpleEleId95relIso = cms.EDProducer("EleIdCutBasedExtProducer",
    robust60relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.04, 0.04, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.025, 0.02, 0.02, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust80cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.07, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.03, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.06, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    algorithm = cms.string('eIDCB'),
    verticesCollection = cms.InputTag("offlineBeamSpot"),
    robust95relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.15, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.15, 2.0, 0.12, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0),
        endcap = cms.vdouble(0.07, 0.03, 0.7, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.08, 0.06, 0.05, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0)
    ),
    robust70relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.03, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.06, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.025, 0.025, 0.02, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    robust70cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.03, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.04, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.03, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust85cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.006, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.09, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.04, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.06, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust80relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.09, 0.07, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.03, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.04, 0.05, 0.025, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust90cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.12, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.1, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.05, 0.03, 0.7, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.07, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    electronQuality = cms.string('95relIso'),
    electronIDType = cms.string('robust'),
    electronVersion = cms.string('V04'),
    src = cms.InputTag("gsfElectrons"),
    robust90relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.12, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.12, 0.09, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.05, 0.03, 0.7, 0.009, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.06, 0.03, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust60cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.025, 0.01, 0.025, 0.004, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.03, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.02, 0.005, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.02, 0.0, -9999.0, 
            9999.0, 9999.0, 0, -1, 0.02, 
            0.02)
    ),
    robust85relIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.04, 0.01, 0.06, 0.006, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.09, 0.08, 0.1, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02),
        endcap = cms.vdouble(0.025, 0.03, 0.04, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 0.05, 0.05, 0.025, 9999.0, 
            9999.0, 9999.0, 9999.0, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.02, 
            0.02)
    ),
    robust95cIsoEleIDCutsV04 = cms.PSet(
        barrel = cms.vdouble(0.15, 0.01, 0.8, 0.007, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.15, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0),
        endcap = cms.vdouble(0.07, 0.03, 0.7, 0.01, -1, 
            -1, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 
            9999.0, 9999.0, 0.1, 0.0, -9999.0, 
            9999.0, 9999.0, 1, -1, 0.0, 
            0.0)
    )
)


process.simpleSecondaryVertexBJetTags = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('simpleSecondaryVertex2Trk'),
    tagInfos = cms.VInputTag(cms.InputTag("secondaryVertexTagInfos"))
)


process.simpleSecondaryVertexHighEffBJetTags = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('simpleSecondaryVertex2Trk'),
    tagInfos = cms.VInputTag(cms.InputTag("secondaryVertexTagInfos"))
)


process.simpleSecondaryVertexHighEffBJetTagsAK5JPT = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('simpleSecondaryVertex2Trk'),
    tagInfos = cms.VInputTag(cms.InputTag("secondaryVertexTagInfosAK5JPT"))
)


process.simpleSecondaryVertexHighEffBJetTagsAK5PF = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('simpleSecondaryVertex2Trk'),
    tagInfos = cms.VInputTag(cms.InputTag("secondaryVertexTagInfosAK5PF"))
)


process.simpleSecondaryVertexHighEffBJetTagsAOD = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('simpleSecondaryVertex2Trk'),
    tagInfos = cms.VInputTag(cms.InputTag("secondaryVertexTagInfosAOD"))
)


process.simpleSecondaryVertexHighEffBJetTagsAODPF = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('simpleSecondaryVertex2Trk'),
    tagInfos = cms.VInputTag(cms.InputTag("secondaryVertexTagInfosAODPF"))
)


process.simpleSecondaryVertexHighPurBJetTags = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('simpleSecondaryVertex3Trk'),
    tagInfos = cms.VInputTag(cms.InputTag("secondaryVertexTagInfos"))
)


process.simpleSecondaryVertexHighPurBJetTagsAK5JPT = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('simpleSecondaryVertex3Trk'),
    tagInfos = cms.VInputTag(cms.InputTag("secondaryVertexTagInfosAK5JPT"))
)


process.simpleSecondaryVertexHighPurBJetTagsAK5PF = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('simpleSecondaryVertex3Trk'),
    tagInfos = cms.VInputTag(cms.InputTag("secondaryVertexTagInfosAK5PF"))
)


process.simpleSecondaryVertexHighPurBJetTagsAOD = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('simpleSecondaryVertex3Trk'),
    tagInfos = cms.VInputTag(cms.InputTag("secondaryVertexTagInfosAOD"))
)


process.simpleSecondaryVertexHighPurBJetTagsAODPF = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('simpleSecondaryVertex3Trk'),
    tagInfos = cms.VInputTag(cms.InputTag("secondaryVertexTagInfosAODPF"))
)


process.simpleSecondaryVertexNegativeBJetTags = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('simpleSecondaryVertex'),
    tagInfos = cms.VInputTag(cms.InputTag("secondaryVertexNegativeTagInfos"))
)


process.sisCone5GenJets = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    src = cms.InputTag("genParticlesForJets"),
    doAreaFastjet = cms.bool(False),
    doPVCorrection = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    doRhoFastjet = cms.bool(False),
    srcPVs = cms.InputTag(""),
    inputEtMin = cms.double(0.0),
    nSigmaPU = cms.double(1.0),
    radiusPU = cms.double(0.5),
    Rho_EtaMax = cms.double(4.5),
    GhostArea = cms.double(0.01),
    jetPtMin = cms.double(3.0),
    jetType = cms.string('GenJet'),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    jetAlgorithm = cms.string('SISCone'),
    rParam = cms.double(0.5)
)


process.sisCone5GenJetsNoMuNoNu = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('SISCone'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.5),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJetsNoMuNoNu"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.sisCone5GenJetsNoNu = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('SISCone'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.5),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJetsNoNu"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.sisCone7GenJets = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('SISCone'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.7),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJets"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.sisCone7GenJetsNoMuNoNu = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('SISCone'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.7),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJetsNoMuNoNu"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.sisCone7GenJetsNoNu = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    doAreaFastjet = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('GenJet'),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('SISCone'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.5),
    maxBadEcalCells = cms.uint32(9999999),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.7),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("genParticlesForJetsNoNu"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)


process.softElectronBJetTags = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('softElectron'),
    tagInfos = cms.VInputTag(cms.InputTag("softElectronTagInfos"))
)


process.softElectronByIP3dBJetTags = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('softLeptonByIP3d'),
    tagInfos = cms.VInputTag(cms.InputTag("softElectronTagInfos"))
)


process.softElectronByPtBJetTags = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('softLeptonByPt'),
    tagInfos = cms.VInputTag(cms.InputTag("softElectronTagInfos"))
)


process.softElectronCands = cms.EDProducer("SoftElectronCandProducer",
    BarreldRGsfTrackElectronCuts = cms.vdouble(0.0, 0.017),
    BarrelEemPinRatioCuts = cms.vdouble(-0.9, 0.39),
    BarrelMVACuts = cms.vdouble(-0.1, 1.0),
    BarrelPtCuts = cms.vdouble(2.0, 9999.0),
    ForwarddRGsfTrackElectronCuts = cms.vdouble(0.0, 0.006),
    ForwardPtCuts = cms.vdouble(2.0, 9999.0),
    ForwardMVACuts = cms.vdouble(-0.24, 1.0),
    ForwardInverseFBremCuts = cms.vdouble(1.0, 7.01),
    electrons = cms.InputTag("gsfElectrons")
)


process.softElectronSelector = cms.EDProducer("BtagGsfElectronSelector",
    input = cms.InputTag("gsfElectrons"),
    selection = cms.InputTag("eidLoose"),
    cut = cms.double(0.5)
)


process.softElectronTagInfos = cms.EDProducer("SoftLepton",
    muonSelection = cms.uint32(0),
    leptons = cms.InputTag("gsfElectrons"),
    primaryVertex = cms.InputTag("offlinePrimaryVertices"),
    leptonCands = cms.InputTag("softElectronCands"),
    leptonId = cms.InputTag(""),
    refineJetAxis = cms.uint32(0),
    jets = cms.InputTag("ak5CaloJets"),
    leptonDeltaRCut = cms.double(0.4),
    leptonChi2Cut = cms.double(10.0)
)


process.softMuonBJetTags = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('softMuon'),
    tagInfos = cms.VInputTag(cms.InputTag("softMuonTagInfos"))
)


process.softMuonBJetTagsAK5JPT = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('softMuon'),
    tagInfos = cms.VInputTag(cms.InputTag("softMuonTagInfosAK5JPT"))
)


process.softMuonBJetTagsAK5PF = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('softMuon'),
    tagInfos = cms.VInputTag(cms.InputTag("softMuonTagInfosAK5PF"))
)


process.softMuonBJetTagsAOD = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('softMuon'),
    tagInfos = cms.VInputTag(cms.InputTag("softMuonTagInfosAOD"))
)


process.softMuonBJetTagsAODPF = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('softMuon'),
    tagInfos = cms.VInputTag(cms.InputTag("softMuonTagInfosAODPF"))
)


process.softMuonByIP3dBJetTags = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('softLeptonByIP3d'),
    tagInfos = cms.VInputTag(cms.InputTag("softMuonTagInfos"))
)


process.softMuonByIP3dBJetTagsAK5JPT = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('softLeptonByIP3d'),
    tagInfos = cms.VInputTag(cms.InputTag("softMuonTagInfosAK5JPT"))
)


process.softMuonByIP3dBJetTagsAK5PF = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('softLeptonByIP3d'),
    tagInfos = cms.VInputTag(cms.InputTag("softMuonTagInfosAK5PF"))
)


process.softMuonByIP3dBJetTagsAOD = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('softLeptonByIP3d'),
    tagInfos = cms.VInputTag(cms.InputTag("softMuonTagInfosAOD"))
)


process.softMuonByIP3dBJetTagsAODPF = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('softLeptonByIP3d'),
    tagInfos = cms.VInputTag(cms.InputTag("softMuonTagInfosAODPF"))
)


process.softMuonByPtBJetTags = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('softLeptonByPt'),
    tagInfos = cms.VInputTag(cms.InputTag("softMuonTagInfos"))
)


process.softMuonByPtBJetTagsAK5JPT = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('softLeptonByPt'),
    tagInfos = cms.VInputTag(cms.InputTag("softMuonTagInfosAK5JPT"))
)


process.softMuonByPtBJetTagsAK5PF = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('softLeptonByPt'),
    tagInfos = cms.VInputTag(cms.InputTag("softMuonTagInfosAK5PF"))
)


process.softMuonByPtBJetTagsAOD = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('softLeptonByPt'),
    tagInfos = cms.VInputTag(cms.InputTag("softMuonTagInfosAOD"))
)


process.softMuonByPtBJetTagsAODPF = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('softLeptonByPt'),
    tagInfos = cms.VInputTag(cms.InputTag("softMuonTagInfosAODPF"))
)


process.softMuonNoIPBJetTags = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('softMuonNoIP'),
    tagInfos = cms.VInputTag(cms.InputTag("softMuonTagInfos"))
)


process.softMuonTagInfos = cms.EDProducer("SoftLepton",
    muonSelection = cms.uint32(1),
    leptons = cms.InputTag("muons"),
    primaryVertex = cms.InputTag("offlinePrimaryVertices"),
    leptonCands = cms.InputTag(""),
    leptonId = cms.InputTag(""),
    refineJetAxis = cms.uint32(0),
    jets = cms.InputTag("ak5CaloJets"),
    leptonDeltaRCut = cms.double(0.4),
    leptonChi2Cut = cms.double(9999.0)
)


process.softMuonTagInfosAK5JPT = cms.EDProducer("SoftLepton",
    muonSelection = cms.uint32(1),
    leptons = cms.InputTag("muons"),
    primaryVertex = cms.InputTag("offlinePrimaryVertices"),
    leptonCands = cms.InputTag(""),
    leptonId = cms.InputTag(""),
    refineJetAxis = cms.uint32(0),
    jets = cms.InputTag("JetPlusTrackZSPCorJetAntiKt5"),
    leptonDeltaRCut = cms.double(0.4),
    leptonChi2Cut = cms.double(9999.0)
)


process.softMuonTagInfosAK5PF = cms.EDProducer("SoftLepton",
    muonSelection = cms.uint32(1),
    leptons = cms.InputTag("muons"),
    primaryVertex = cms.InputTag("offlinePrimaryVertices"),
    leptonCands = cms.InputTag(""),
    leptonId = cms.InputTag(""),
    refineJetAxis = cms.uint32(0),
    jets = cms.InputTag("ak5PFJets"),
    leptonDeltaRCut = cms.double(0.4),
    leptonChi2Cut = cms.double(9999.0)
)


process.softMuonTagInfosAOD = cms.EDProducer("SoftLepton",
    muonSelection = cms.uint32(1),
    leptons = cms.InputTag("muons"),
    primaryVertex = cms.InputTag("offlinePrimaryVertices"),
    leptonCands = cms.InputTag(""),
    leptonId = cms.InputTag(""),
    refineJetAxis = cms.uint32(0),
    jets = cms.InputTag("ak5CaloJets"),
    leptonDeltaRCut = cms.double(0.4),
    leptonChi2Cut = cms.double(9999.0)
)


process.softMuonTagInfosAODPF = cms.EDProducer("SoftLepton",
    muonSelection = cms.uint32(1),
    leptons = cms.InputTag("muons"),
    primaryVertex = cms.InputTag("offlinePrimaryVertices"),
    leptonCands = cms.InputTag(""),
    leptonId = cms.InputTag(""),
    refineJetAxis = cms.uint32(0),
    jets = cms.InputTag("pfNoTauPF"),
    leptonDeltaRCut = cms.double(0.4),
    leptonChi2Cut = cms.double(9999.0)
)


process.tauGenJetMatch = cms.EDProducer("GenJetMatcher",
    src = cms.InputTag("shrinkingConePFTauProducer"),
    maxDPtRel = cms.double(3.0),
    mcPdgId = cms.vint32(),
    mcStatus = cms.vint32(),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = cms.double(0.1),
    checkCharge = cms.bool(False),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("tauGenJetsSelectorAllHadrons")
)


process.tauGenJetMatchPF = cms.EDProducer("GenJetMatcher",
    src = cms.InputTag("pfTausPF"),
    maxDPtRel = cms.double(3.0),
    mcPdgId = cms.vint32(),
    mcStatus = cms.vint32(),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = cms.double(0.1),
    checkCharge = cms.bool(False),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("tauGenJetsSelectorAllHadronsPF")
)


process.tauGenJets = cms.EDProducer("TauGenJetProducer",
    includeNeutrinos = cms.bool(False),
    GenParticles = cms.InputTag("genParticles"),
    verbose = cms.untracked.bool(False)
)


process.tauGenJetsPF = cms.EDProducer("TauGenJetProducer",
    includeNeutrinos = cms.bool(False),
    GenParticles = cms.InputTag("genParticles"),
    verbose = cms.untracked.bool(False)
)


process.tauIsoDepositPFCandidates = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("shrinkingConePFTauProducer"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(10000.0),
        ComponentName = cms.string('PFTauExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(10000.0),
        dRvetoPFTauSignalConeConstituents = cms.double(0.01),
        tauSource = cms.InputTag("shrinkingConePFTauProducer"),
        DR_Veto = cms.double(0.0),
        DepositLabel = cms.untracked.string(''),
        candidateSource = cms.InputTag("particleFlow"),
        dRmatchPFTau = cms.double(0.1)
    )
)


process.tauIsoDepositPFCandidatesPF = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("pfTausPF"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(10000.0),
        ComponentName = cms.string('PFTauExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(10000.0),
        dRvetoPFTauSignalConeConstituents = cms.double(0.01),
        tauSource = cms.InputTag("pfTausPF"),
        DR_Veto = cms.double(0.0),
        DepositLabel = cms.untracked.string(''),
        candidateSource = cms.InputTag("particleFlow"),
        dRmatchPFTau = cms.double(0.1)
    )
)


process.tauIsoDepositPFChargedHadrons = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("shrinkingConePFTauProducer"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(0.2),
        ComponentName = cms.string('PFTauExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(0.1),
        dRvetoPFTauSignalConeConstituents = cms.double(0.01),
        tauSource = cms.InputTag("shrinkingConePFTauProducer"),
        DR_Veto = cms.double(0.0),
        DepositLabel = cms.untracked.string(''),
        candidateSource = cms.InputTag("pfAllChargedHadrons"),
        dRmatchPFTau = cms.double(0.1)
    )
)


process.tauIsoDepositPFChargedHadronsPF = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("pfTausPF"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(0.2),
        ComponentName = cms.string('PFTauExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(0.1),
        dRvetoPFTauSignalConeConstituents = cms.double(0.01),
        tauSource = cms.InputTag("pfTausPF"),
        DR_Veto = cms.double(0.0),
        DepositLabel = cms.untracked.string(''),
        candidateSource = cms.InputTag("pfAllChargedHadronsPF"),
        dRmatchPFTau = cms.double(0.1)
    )
)


process.tauIsoDepositPFGammas = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("shrinkingConePFTauProducer"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(10000.0),
        ComponentName = cms.string('PFTauExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(10000.0),
        dRvetoPFTauSignalConeConstituents = cms.double(0.01),
        tauSource = cms.InputTag("shrinkingConePFTauProducer"),
        DR_Veto = cms.double(0.0),
        DepositLabel = cms.untracked.string(''),
        candidateSource = cms.InputTag("pfAllPhotons"),
        dRmatchPFTau = cms.double(0.1)
    )
)


process.tauIsoDepositPFGammasPF = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("pfTausPF"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(10000.0),
        ComponentName = cms.string('PFTauExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(10000.0),
        dRvetoPFTauSignalConeConstituents = cms.double(0.01),
        tauSource = cms.InputTag("pfTausPF"),
        DR_Veto = cms.double(0.0),
        DepositLabel = cms.untracked.string(''),
        candidateSource = cms.InputTag("pfAllPhotonsPF"),
        dRmatchPFTau = cms.double(0.1)
    )
)


process.tauIsoDepositPFNeutralHadrons = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("shrinkingConePFTauProducer"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(10000.0),
        ComponentName = cms.string('PFTauExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(10000.0),
        dRvetoPFTauSignalConeConstituents = cms.double(0.01),
        tauSource = cms.InputTag("shrinkingConePFTauProducer"),
        DR_Veto = cms.double(0.0),
        DepositLabel = cms.untracked.string(''),
        candidateSource = cms.InputTag("pfAllNeutralHadrons"),
        dRmatchPFTau = cms.double(0.1)
    )
)


process.tauIsoDepositPFNeutralHadronsPF = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("pfTausPF"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(10000.0),
        ComponentName = cms.string('PFTauExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(10000.0),
        dRvetoPFTauSignalConeConstituents = cms.double(0.01),
        tauSource = cms.InputTag("pfTausPF"),
        DR_Veto = cms.double(0.0),
        DepositLabel = cms.untracked.string(''),
        candidateSource = cms.InputTag("pfAllNeutralHadronsPF"),
        dRmatchPFTau = cms.double(0.1)
    )
)


process.tauMatch = cms.EDProducer("MCMatcher",
    src = cms.InputTag("shrinkingConePFTauProducer"),
    maxDPtRel = cms.double(999.9),
    mcPdgId = cms.vint32(15),
    mcStatus = cms.vint32(2),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = cms.double(0.3),
    checkCharge = cms.bool(False),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("genParticles")
)


process.tauMatchPF = cms.EDProducer("MCMatcher",
    src = cms.InputTag("pfTausPF"),
    maxDPtRel = cms.double(999.9),
    mcPdgId = cms.vint32(15),
    mcStatus = cms.vint32(2),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = cms.double(999.9),
    checkCharge = cms.bool(True),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("genParticles")
)


process.tcRecoTauProducer = cms.EDProducer("TCRecoTauProducer",
    tkminTrackerHitsn = cms.int32(5),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    EtCaloOverTrackMin = cms.double(-0.9),
    EERecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    tkmaxChi2 = cms.double(100.0),
    EtHcalOverTrackMin = cms.double(-0.3),
    CaloRecoTauProducer = cms.InputTag("JPTCaloRecoTauProducer"),
    EBRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    tkminPixelHitsn = cms.int32(0),
    MatchingConeSize = cms.double(0.1),
    TrackCollection = cms.InputTag("generalTracks"),
    HBHERecHitCollection = cms.InputTag("hbhereco"),
    EtCaloOverTrackMax = cms.double(0.0),
    TrackAssociatorParameters = cms.PSet(
        muonMaxDistanceSigmaX = cms.double(0.0),
        muonMaxDistanceSigmaY = cms.double(0.0),
        CSCSegmentCollectionLabel = cms.InputTag("cscSegments"),
        dRHcal = cms.double(9999.0),
        dRPreshowerPreselection = cms.double(0.2),
        CaloTowerCollectionLabel = cms.InputTag("towerMaker"),
        useEcal = cms.bool(True),
        dREcal = cms.double(9999.0),
        dREcalPreselection = cms.double(0.05),
        HORecHitCollectionLabel = cms.InputTag("horeco"),
        dRMuon = cms.double(9999.0),
        propagateAllDirections = cms.bool(True),
        muonMaxDistanceX = cms.double(5.0),
        muonMaxDistanceY = cms.double(5.0),
        useHO = cms.bool(True),
        trajectoryUncertaintyTolerance = cms.double(-1.0),
        usePreshower = cms.bool(False),
        DTRecSegment4DCollectionLabel = cms.InputTag("dt4DSegments"),
        EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
        dRHcalPreselection = cms.double(0.2),
        useMuon = cms.bool(True),
        useCalo = cms.bool(False),
        accountForTrajectoryChangeCalo = cms.bool(False),
        EBRecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
        dRMuonPreselection = cms.double(0.2),
        truthMatch = cms.bool(False),
        HBHERecHitCollectionLabel = cms.InputTag("hbhereco"),
        useHcal = cms.bool(True)
    ),
    SignalConeSize = cms.double(0.2),
    EcalConeSize = cms.double(0.5),
    Track_minPt = cms.double(1.0),
    EtHcalOverTrackMax = cms.double(1.0),
    HFRecHitCollection = cms.InputTag("hfreco"),
    DropRejectedJets = cms.untracked.bool(False),
    HORecHitCollection = cms.InputTag("horeco"),
    DropCaloJets = cms.untracked.bool(False),
    tkmaxipt = cms.double(0.1)
)


process.trackCountingHighEffBJetTags = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('trackCounting3D2nd'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfos"))
)


process.trackCountingHighEffBJetTagsAK5JPT = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('trackCounting3D2nd'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAK5JPT"))
)


process.trackCountingHighEffBJetTagsAK5PF = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('trackCounting3D2nd'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAK5PF"))
)


process.trackCountingHighEffBJetTagsAOD = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('trackCounting3D2nd'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAOD"))
)


process.trackCountingHighEffBJetTagsAODPF = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('trackCounting3D2nd'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAODPF"))
)


process.trackCountingHighPurBJetTags = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('trackCounting3D3rd'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfos"))
)


process.trackCountingHighPurBJetTagsAK5JPT = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('trackCounting3D3rd'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAK5JPT"))
)


process.trackCountingHighPurBJetTagsAK5PF = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('trackCounting3D3rd'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAK5PF"))
)


process.trackCountingHighPurBJetTagsAOD = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('trackCounting3D3rd'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAOD"))
)


process.trackCountingHighPurBJetTagsAODPF = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('trackCounting3D3rd'),
    tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosAODPF"))
)


process.HBHENoiseFilter = cms.EDFilter("HBHENoiseFilter",
    minHPDHits = cms.int32(17),
    minHighEHitTime = cms.double(-9999.0),
    minHPDNoOtherHits = cms.int32(10),
    maxRatio = cms.double(0.96),
    minZeros = cms.int32(10),
    maxHighEHitTime = cms.double(9999.0),
    maxRBXEMF = cms.double(-9999.0),
    minRBXHits = cms.int32(999),
    minRatio = cms.double(0.7)
)


process.countPatElectrons = cms.EDFilter("PATCandViewCountFilter",
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("cleanPatElectrons"),
    minNumber = cms.uint32(0)
)


process.countPatElectronsPF = cms.EDFilter("PATCandViewCountFilter",
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("selectedPatElectronsPF"),
    minNumber = cms.uint32(0)
)


process.countPatJets = cms.EDFilter("PATCandViewCountFilter",
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("cleanPatJetsAK5Calo"),
    minNumber = cms.uint32(0)
)


process.countPatJetsAK5JPT = cms.EDFilter("PATCandViewCountFilter",
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("cleanPatJetsAK5JPT"),
    minNumber = cms.uint32(0)
)


process.countPatJetsAK5PF = cms.EDFilter("PATCandViewCountFilter",
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("cleanPatJetsAK5PF"),
    minNumber = cms.uint32(0)
)


process.countPatJetsPF = cms.EDFilter("PATCandViewCountFilter",
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("selectedPatJetsPF"),
    minNumber = cms.uint32(0)
)


process.countPatLeptons = cms.EDFilter("PATLeptonCountFilter",
    maxNumber = cms.uint32(999999),
    countElectrons = cms.bool(True),
    muonSource = cms.InputTag("cleanPatMuons"),
    minNumber = cms.uint32(0),
    electronSource = cms.InputTag("cleanPatElectrons"),
    tauSource = cms.InputTag("cleanPatTaus"),
    countTaus = cms.bool(False),
    countMuons = cms.bool(True)
)


process.countPatLeptonsPF = cms.EDFilter("PATLeptonCountFilter",
    maxNumber = cms.uint32(999999),
    countElectrons = cms.bool(True),
    muonSource = cms.InputTag("selectedPatMuonsPF"),
    minNumber = cms.uint32(0),
    electronSource = cms.InputTag("selectedPatElectronsPF"),
    tauSource = cms.InputTag("selectedPatTausPF"),
    countTaus = cms.bool(False),
    countMuons = cms.bool(True)
)


process.countPatMuons = cms.EDFilter("PATCandViewCountFilter",
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("cleanPatMuons"),
    minNumber = cms.uint32(0)
)


process.countPatMuonsPF = cms.EDFilter("PATCandViewCountFilter",
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("selectedPatMuonsPF"),
    minNumber = cms.uint32(0)
)


process.countPatPFParticlesPF = cms.EDFilter("PATCandViewCountFilter",
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("patPFParticlesPF"),
    minNumber = cms.uint32(0)
)


process.countPatPhotons = cms.EDFilter("PATCandViewCountFilter",
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("cleanPatPhotons"),
    minNumber = cms.uint32(0)
)


process.countPatPhotonsPF = cms.EDFilter("PATCandViewCountFilter",
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("selectedPatPhotonsPF"),
    minNumber = cms.uint32(0)
)


process.countPatTaus = cms.EDFilter("PATCandViewCountFilter",
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("cleanPatTaus"),
    minNumber = cms.uint32(0)
)


process.countPatTausPF = cms.EDFilter("PATCandViewCountFilter",
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("selectedPatTausPF"),
    minNumber = cms.uint32(0)
)


process.hbheNoise = cms.EDFilter("HBHENoiseFilter",
    minRBXHits = cms.int32(999),
    minHighEHitTime = cms.double(-9999.0),
    minHPDNoOtherHits = cms.int32(10),
    maxRatio = cms.double(0.96),
    minZeros = cms.int32(10),
    maxHighEHitTime = cms.double(9999.0),
    maxRBXEMF = cms.double(-9999.0),
    minHPDHits = cms.int32(17),
    minRatio = cms.double(0.7)
)


process.hltPhysicsDeclared = cms.EDFilter("HLTPhysicsDeclared",
    L1GtReadoutRecordTag = cms.InputTag("hltGtDigis"),
    invert = cms.bool(False)
)


process.pfAllChargedHadrons = cms.EDFilter("PdgIdPFCandidateSelector",
    pdgId = cms.vint32(211, -211, 321, -321, 999211, 
        2212, -2212),
    src = cms.InputTag("pfNoPileUp")
)


process.pfAllChargedHadronsPF = cms.EDFilter("PdgIdPFCandidateSelector",
    pdgId = cms.vint32(211, -211, 321, -321, 999211, 
        2212, -2212),
    src = cms.InputTag("pfNoPileUpPF")
)


process.pfAllElectrons = cms.EDFilter("PdgIdPFCandidateSelector",
    pdgId = cms.vint32(11, -11),
    src = cms.InputTag("pfNoMuon")
)


process.pfAllElectronsPF = cms.EDFilter("PdgIdPFCandidateSelector",
    pdgId = cms.vint32(11, -11),
    src = cms.InputTag("pfNoMuonPF")
)


process.pfAllMuons = cms.EDFilter("PdgIdPFCandidateSelector",
    pdgId = cms.vint32(-13, 13),
    src = cms.InputTag("pfNoPileUp")
)


process.pfAllMuonsPF = cms.EDFilter("PdgIdPFCandidateSelector",
    pdgId = cms.vint32(-13, 13),
    src = cms.InputTag("pfNoPileUpPF")
)


process.pfAllNeutralHadrons = cms.EDFilter("PdgIdPFCandidateSelector",
    pdgId = cms.vint32(111, 130, 310, 2112),
    src = cms.InputTag("pfNoPileUp")
)


process.pfAllNeutralHadronsPF = cms.EDFilter("PdgIdPFCandidateSelector",
    pdgId = cms.vint32(111, 130, 310, 2112),
    src = cms.InputTag("pfNoPileUpPF")
)


process.pfAllPhotons = cms.EDFilter("PdgIdPFCandidateSelector",
    pdgId = cms.vint32(22),
    src = cms.InputTag("pfNoPileUp")
)


process.pfAllPhotonsPF = cms.EDFilter("PdgIdPFCandidateSelector",
    pdgId = cms.vint32(22),
    src = cms.InputTag("pfNoPileUpPF")
)


process.pfElectronsFromVertex = cms.EDFilter("IPCutPFCandidateSelector",
    d0Cut = cms.double(0.2),
    src = cms.InputTag("pfAllElectrons"),
    dzSigCut = cms.double(99.0),
    d0SigCut = cms.double(99.0),
    vertices = cms.InputTag("offlinePrimaryVertices"),
    dzCut = cms.double(0.5)
)


process.pfElectronsFromVertexPF = cms.EDFilter("IPCutPFCandidateSelector",
    d0Cut = cms.double(0.2),
    src = cms.InputTag("pfAllElectronsPF"),
    dzSigCut = cms.double(99.0),
    vertices = cms.InputTag("offlinePrimaryVertices"),
    d0SigCut = cms.double(99.0),
    dzCut = cms.double(0.5)
)


process.pfIsolatedElectrons = cms.EDFilter("IsolatedPFCandidateSelector",
    src = cms.InputTag("pfSelectedElectrons"),
    isRelative = cms.bool(True),
    combinedIsolationCut = cms.double(0.2),
    isCombined = cms.bool(True),
    isolationValueMaps = cms.VInputTag(cms.InputTag("isoValElectronWithCharged"), cms.InputTag("isoValElectronWithNeutral"), cms.InputTag("isoValElectronWithPhotons")),
    isolationCuts = cms.vdouble(10, 10, 10)
)


process.pfIsolatedElectronsPF = cms.EDFilter("IsolatedPFCandidateSelector",
    src = cms.InputTag("pfSelectedElectronsPF"),
    isRelative = cms.bool(True),
    combinedIsolationCut = cms.double(0.2),
    isCombined = cms.bool(True),
    isolationValueMaps = cms.VInputTag(cms.InputTag("isoValElectronWithChargedPF"), cms.InputTag("isoValElectronWithNeutralPF"), cms.InputTag("isoValElectronWithPhotonsPF")),
    isolationCuts = cms.vdouble(10, 10, 10)
)


process.pfIsolatedMuons = cms.EDFilter("IsolatedPFCandidateSelector",
    src = cms.InputTag("pfSelectedMuons"),
    isRelative = cms.bool(True),
    combinedIsolationCut = cms.double(0.15),
    isCombined = cms.bool(True),
    isolationValueMaps = cms.VInputTag(cms.InputTag("isoValMuonWithCharged"), cms.InputTag("isoValMuonWithNeutral"), cms.InputTag("isoValMuonWithPhotons")),
    isolationCuts = cms.vdouble(10, 10, 10)
)


process.pfIsolatedMuonsPF = cms.EDFilter("IsolatedPFCandidateSelector",
    src = cms.InputTag("pfSelectedMuonsPF"),
    isRelative = cms.bool(True),
    combinedIsolationCut = cms.double(0.15),
    isCombined = cms.bool(True),
    isolationValueMaps = cms.VInputTag(cms.InputTag("isoValMuonWithChargedPF"), cms.InputTag("isoValMuonWithNeutralPF"), cms.InputTag("isoValMuonWithPhotonsPF")),
    isolationCuts = cms.vdouble(10, 10, 10)
)


process.pfMuonsFromVertex = cms.EDFilter("IPCutPFCandidateSelector",
    d0Cut = cms.double(0.2),
    src = cms.InputTag("pfAllMuons"),
    dzSigCut = cms.double(99.0),
    d0SigCut = cms.double(99.0),
    vertices = cms.InputTag("offlinePrimaryVertices"),
    dzCut = cms.double(0.5)
)


process.pfMuonsFromVertexPF = cms.EDFilter("IPCutPFCandidateSelector",
    d0Cut = cms.double(0.2),
    src = cms.InputTag("pfAllMuonsPF"),
    dzSigCut = cms.double(99.0),
    vertices = cms.InputTag("offlinePrimaryVertices"),
    d0SigCut = cms.double(99.0),
    dzCut = cms.double(0.5)
)


process.pfSelectedElectrons = cms.EDFilter("GenericPFCandidateSelector",
    src = cms.InputTag("pfElectronsFromVertex"),
    cut = cms.string('pt>5 && gsfTrackRef.isNonnull && gsfTrackRef.trackerExpectedHitsInner.numberOfHits<2')
)


process.pfSelectedElectronsPF = cms.EDFilter("GenericPFCandidateSelector",
    src = cms.InputTag("pfElectronsFromVertexPF"),
    cut = cms.string('pt>5 && gsfTrackRef.isNonnull && gsfTrackRef.trackerExpectedHitsInner.numberOfHits<2')
)


process.pfSelectedMuons = cms.EDFilter("GenericPFCandidateSelector",
    src = cms.InputTag("pfMuonsFromVertex"),
    cut = cms.string('pt>5')
)


process.pfSelectedMuonsPF = cms.EDFilter("GenericPFCandidateSelector",
    src = cms.InputTag("pfMuonsFromVertexPF"),
    cut = cms.string('pt>5')
)


process.pfTauSelector = cms.EDFilter("PFTauSelector",
    discriminators = cms.VPSet(cms.PSet(
        discriminator = cms.InputTag("fixedConePFTauDiscriminationByIsolation"),
        selectionCut = cms.double(0.5)
    )),
    src = cms.InputTag("fixedConePFTauProducer")
)


process.pfTaus = cms.EDFilter("PFTauSelector",
    discriminators = cms.VPSet(cms.PSet(
        discriminator = cms.InputTag("pfTausDiscriminationByLeadingPionPtCut"),
        selectionCut = cms.double(0.5)
    ), 
        cms.PSet(
            discriminator = cms.InputTag("pfTausDiscriminationByIsolation"),
            selectionCut = cms.double(0.5)
        )),
    src = cms.InputTag("shrinkingConePFTauProducer")
)


process.pfTausPF = cms.EDFilter("PFTauSelector",
    discriminators = cms.VPSet(cms.PSet(
        discriminator = cms.InputTag("pfTausBaseDiscriminationByIsolationPF"),
        selectionCut = cms.double(0.5)
    ), 
        cms.PSet(
            discriminator = cms.InputTag("pfTausBaseDiscriminationByLeadingPionPtCutPF"),
            selectionCut = cms.double(0.5)
        )),
    src = cms.InputTag("shrinkingConePFTauProducerPF")
)


process.physicsDeclared = cms.EDFilter("HLTPhysicsDeclared",
    L1GtReadoutRecordTag = cms.InputTag("gtDigis"),
    invert = cms.bool(False)
)


process.primaryVertexFilter = cms.EDFilter("VertexSelector",
    filter = cms.bool(True),
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string('!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2')
)


process.primaryVertexFilter2 = cms.EDFilter("GoodVertexFilter",
    vertexCollection = cms.InputTag("offlinePrimaryVertices"),
    maxd0 = cms.double(2),
    minimumNDOF = cms.uint32(4),
    maxAbsZ = cms.double(24)
)


process.removePKAM = cms.EDFilter("FilterOutScraping",
    debugOn = cms.untracked.bool(False),
    thresh = cms.untracked.double(0.25),
    numtrack = cms.untracked.uint32(10),
    applyfilter = cms.untracked.bool(True)
)


process.selectedPatElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("patElectrons"),
    cut = cms.string('')
)


process.selectedPatElectronsPF = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("patElectronsPF"),
    cut = cms.string('')
)


process.selectedPatJets = cms.EDFilter("PATJetSelector",
    src = cms.InputTag("patJets"),
    cut = cms.string('')
)


process.selectedPatJetsAK5JPT = cms.EDFilter("PATJetSelector",
    src = cms.InputTag("patJetsAK5JPT"),
    cut = cms.string('')
)


process.selectedPatJetsAK5PF = cms.EDFilter("PATJetSelector",
    src = cms.InputTag("patJetsAK5PF"),
    cut = cms.string('')
)


process.selectedPatJetsPF = cms.EDFilter("PATJetSelector",
    src = cms.InputTag("patJetsPF"),
    cut = cms.string('')
)


process.selectedPatMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string('')
)


process.selectedPatMuonsPF = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsPF"),
    cut = cms.string('')
)


process.selectedPatPFParticlesPF = cms.EDFilter("PATPFParticleSelector",
    src = cms.InputTag("patPFParticlesPF"),
    cut = cms.string('')
)


process.selectedPatPhotons = cms.EDFilter("PATPhotonSelector",
    src = cms.InputTag("patPhotons"),
    cut = cms.string('')
)


process.selectedPatPhotonsPF = cms.EDFilter("PATPhotonSelector",
    src = cms.InputTag("patPhotonsPF"),
    cut = cms.string('')
)


process.selectedPatTaus = cms.EDFilter("PATTauSelector",
    src = cms.InputTag("patTaus"),
    cut = cms.string('')
)


process.selectedPatTausPF = cms.EDFilter("PATTauSelector",
    src = cms.InputTag("patTausPF"),
    cut = cms.string('')
)


process.tauGenJetsSelectorAllHadrons = cms.EDFilter("TauGenJetDecayModeSelector",
    filter = cms.bool(False),
    src = cms.InputTag("tauGenJets"),
    select = cms.vstring('oneProng0Pi0', 
        'oneProng1Pi0', 
        'oneProng2Pi0', 
        'oneProngOther', 
        'threeProng0Pi0', 
        'threeProng1Pi0', 
        'threeProngOther', 
        'rare')
)


process.tauGenJetsSelectorAllHadronsPF = cms.EDFilter("TauGenJetDecayModeSelector",
    filter = cms.bool(False),
    src = cms.InputTag("tauGenJetsPF"),
    select = cms.vstring('oneProng0Pi0', 
        'oneProng1Pi0', 
        'oneProng2Pi0', 
        'oneProngOther', 
        'threeProng0Pi0', 
        'threeProng1Pi0', 
        'threeProngOther', 
        'rare')
)


process.analysisNtuplePAT = cms.EDAnalyzer("AnalysisNtuplePAT",
    trackJetParameters = cms.untracked.PSet(
        photonIso = cms.untracked.double(1.5),
        photonPt = cms.untracked.double(10.0),
        usePFJets = cms.untracked.bool(False),
        muonIso = cms.untracked.double(0.05),
        jetTag = cms.untracked.InputTag("cleanPatJetsAK5Track"),
        jetMaxEMF = cms.untracked.double(0.99),
        tauIso = cms.untracked.double(0.1),
        jetMinPt = cms.untracked.double(8.0),
        useCaloJets = cms.untracked.bool(False),
        useJPTJets = cms.untracked.bool(False),
        genJetTag = cms.untracked.InputTag("ak5GenJets"),
        electronIso = cms.untracked.double(0.1),
        prefixJets = cms.untracked.string('Track'),
        electronPt = cms.untracked.double(30.0),
        useTrackJets = cms.untracked.bool(True),
        debugJets = cms.untracked.int32(0),
        muonPt = cms.untracked.double(10.0),
        tauPt = cms.untracked.double(10.0),
        jetMinEMF = cms.untracked.double(0.01),
        doMCJets = cms.untracked.bool(True),
        jetCorTag = cms.untracked.string('AK5Track'),
        jetMaxEta = cms.untracked.double(5.0)
    ),
    mcTruthParameters = cms.untracked.PSet(
        genParticleTag = cms.untracked.InputTag("genParticles"),
        debugTruth = cms.untracked.int32(0)
    ),
    leptonParameters = cms.untracked.PSet(
        muonMaxEt = cms.untracked.double(10.0),
        tauTag = cms.untracked.InputTag("cleanPatTaus"),
        elecTag = cms.untracked.InputTag("cleanPatElectrons"),
        elecMinEt = cms.untracked.double(2.5),
        elecRelIso = cms.untracked.double(0.5),
        muonRelIso = cms.untracked.double(0.1),
        elecMaxEta = cms.untracked.double(5.0),
        muonTag = cms.untracked.InputTag("cleanPatMuons"),
        tauMaxEta = cms.untracked.double(5.0),
        tauRelIso = cms.untracked.double(0.1),
        prefixLeps = cms.untracked.string(''),
        tauMinEt = cms.untracked.double(5.0),
        muonMinEt = cms.untracked.double(2.5),
        tauMaxEt = cms.untracked.double(9999.0),
        muonMaxEta = cms.untracked.double(5.0),
        debugLeps = cms.untracked.int32(0),
        elecMaxEt = cms.untracked.double(15.0)
    ),
    debugDiJets = cms.untracked.int32(0),
    pfmetTypeIParameters = cms.untracked.PSet(
        metTag = cms.untracked.InputTag("patMETsTypeIPF"),
        prefixMET = cms.untracked.string('PFTypeI'),
        doMCMET = cms.untracked.bool(True),
        debugMET = cms.untracked.int32(0),
        genMETTag = cms.untracked.InputTag("genMetTrue")
    ),
    vertexParameters = cms.untracked.PSet(
        minVtxNdof = cms.untracked.int32(4),
        beamspotTag = cms.untracked.InputTag("offlineBeamSpot"),
        maxVtxRho = cms.untracked.double(2.0),
        maxVtxZ = cms.untracked.double(24.0),
        vtxTag = cms.untracked.InputTag("offlinePrimaryVertices"),
        minVtxTrks = cms.untracked.int32(3),
        maxVtxChi2 = cms.untracked.double(999),
        minNVtx = cms.untracked.int32(1),
        debugVtx = cms.untracked.int32(0)
    ),
    doMCTruth = cms.untracked.bool(False),
    pfphotonParameters = cms.untracked.PSet(
        prefixPhots = cms.untracked.string('PF'),
        photRelIso = cms.untracked.double(1.5),
        photTag = cms.untracked.InputTag("selectedPatPhotonsPF"),
        debugPhots = cms.untracked.int32(0),
        photMinEt = cms.untracked.double(2.5),
        photMaxEta = cms.untracked.double(5.0),
        photMaxEt = cms.untracked.double(10000.0)
    ),
    caloJetParameters = cms.untracked.PSet(
        photonIso = cms.untracked.double(1.5),
        photonPt = cms.untracked.double(10.0),
        usePFJets = cms.untracked.bool(False),
        muonIso = cms.untracked.double(0.05),
        jetTag = cms.untracked.InputTag("cleanPatJetsAK5Calo"),
        jetMaxEMF = cms.untracked.double(0.99),
        tauIso = cms.untracked.double(0.1),
        jetMinPt = cms.untracked.double(10.0),
        useCaloJets = cms.untracked.bool(True),
        useJPTJets = cms.untracked.bool(False),
        genJetTag = cms.untracked.InputTag("ak5GenJets"),
        electronIso = cms.untracked.double(0.1),
        prefixJets = cms.untracked.string('Calo'),
        electronPt = cms.untracked.double(30.0),
        useTrackJets = cms.untracked.bool(False),
        debugJets = cms.untracked.int32(0),
        muonPt = cms.untracked.double(10.0),
        tauPt = cms.untracked.double(10.0),
        jetMinEMF = cms.untracked.double(0.01),
        doMCJets = cms.untracked.bool(False),
        jetCorTag = cms.untracked.string('AK5Calo'),
        jetMaxEta = cms.untracked.double(5.0)
    ),
    pfJetParameters = cms.untracked.PSet(
        photonIso = cms.untracked.double(1.5),
        photonPt = cms.untracked.double(10.0),
        usePFJets = cms.untracked.bool(True),
        muonIso = cms.untracked.double(0.05),
        jetTag = cms.untracked.InputTag("cleanPatJetsAK5PF"),
        jetMaxEMF = cms.untracked.double(0.99),
        tauIso = cms.untracked.double(0.1),
        jetMinPt = cms.untracked.double(8.0),
        useCaloJets = cms.untracked.bool(False),
        useJPTJets = cms.untracked.bool(False),
        genJetTag = cms.untracked.InputTag("ak5GenJets"),
        electronIso = cms.untracked.double(0.1),
        prefixJets = cms.untracked.string('PF'),
        electronPt = cms.untracked.double(30.0),
        useTrackJets = cms.untracked.bool(False),
        debugJets = cms.untracked.int32(0),
        muonPt = cms.untracked.double(10.0),
        tauPt = cms.untracked.double(10.0),
        jetMinEMF = cms.untracked.double(0.01),
        doMCJets = cms.untracked.bool(False),
        jetCorTag = cms.untracked.string('AK5PF'),
        jetMaxEta = cms.untracked.double(5.0)
    ),
    triggerParameters = cms.untracked.PSet(
        hlTriggerResults = cms.untracked.InputTag("TriggerResults","","REDIGI"),
        l1TriggerResults = cms.untracked.InputTag("gtDigis"),
        getL1Info = cms.untracked.bool(False),
        getHLTfromConfig = cms.untracked.bool(False),
        debugTriggers = cms.untracked.int32(0)
    ),
    pfmetParameters = cms.untracked.PSet(
        metTag = cms.untracked.InputTag("patMETsPF"),
        prefixMET = cms.untracked.string('PF'),
        doMCMET = cms.untracked.bool(False),
        debugMET = cms.untracked.int32(0),
        genMETTag = cms.untracked.InputTag("genMetTrue")
    ),
    jptJetParameters = cms.untracked.PSet(
        photonIso = cms.untracked.double(1.5),
        photonPt = cms.untracked.double(10.0),
        usePFJets = cms.untracked.bool(False),
        muonIso = cms.untracked.double(0.05),
        jetTag = cms.untracked.InputTag("cleanPatJetsAK5JPT"),
        jetMaxEMF = cms.untracked.double(0.99),
        tauIso = cms.untracked.double(0.1),
        jetMinPt = cms.untracked.double(8.0),
        useCaloJets = cms.untracked.bool(False),
        useJPTJets = cms.untracked.bool(True),
        genJetTag = cms.untracked.InputTag("ak5GenJets"),
        electronIso = cms.untracked.double(0.1),
        prefixJets = cms.untracked.string('JPT'),
        electronPt = cms.untracked.double(30.0),
        useTrackJets = cms.untracked.bool(False),
        debugJets = cms.untracked.int32(0),
        muonPt = cms.untracked.double(10.0),
        tauPt = cms.untracked.double(10.0),
        jetMinEMF = cms.untracked.double(0.01),
        doMCJets = cms.untracked.bool(False),
        jetCorTag = cms.untracked.string('AK5JPT'),
        jetMaxEta = cms.untracked.double(5.0)
    ),
    calometTypeIIParameters = cms.untracked.PSet(
        metTag = cms.untracked.InputTag("patMETsAK5CaloTypeII"),
        prefixMET = cms.untracked.string('CaloTypeII'),
        doMCMET = cms.untracked.bool(False),
        debugMET = cms.untracked.int32(0),
        genMETTag = cms.untracked.InputTag("genMetCalo")
    ),
    tcmetParameters = cms.untracked.PSet(
        metTag = cms.untracked.InputTag("patMETsTC"),
        prefixMET = cms.untracked.string('TC'),
        doMCMET = cms.untracked.bool(False),
        debugMET = cms.untracked.int32(0),
        genMETTag = cms.untracked.InputTag("genMetCalo")
    ),
    trackParameters = cms.untracked.PSet(
        debugTrack = cms.untracked.int32(0),
        doMCTracks = cms.untracked.bool(False),
        trackTag = cms.untracked.InputTag("generalTracks")
    ),
    pf2patJetParameters = cms.untracked.PSet(
        photonIso = cms.untracked.double(1.5),
        photonPt = cms.untracked.double(10.0),
        usePFJets = cms.untracked.bool(True),
        muonIso = cms.untracked.double(0.05),
        jetTag = cms.untracked.InputTag("selectedPatJetsPF"),
        jetMaxEMF = cms.untracked.double(0.99),
        tauIso = cms.untracked.double(0.1),
        jetMinPt = cms.untracked.double(8.0),
        useCaloJets = cms.untracked.bool(False),
        useJPTJets = cms.untracked.bool(False),
        genJetTag = cms.untracked.InputTag("ak5GenJets"),
        electronIso = cms.untracked.double(0.1),
        prefixJets = cms.untracked.string('PF2PAT'),
        electronPt = cms.untracked.double(30.0),
        useTrackJets = cms.untracked.bool(False),
        debugJets = cms.untracked.int32(0),
        muonPt = cms.untracked.double(10.0),
        tauPt = cms.untracked.double(10.0),
        jetMinEMF = cms.untracked.double(0.01),
        doMCJets = cms.untracked.bool(False),
        jetCorTag = cms.untracked.string('AK5PF'),
        jetMaxEta = cms.untracked.double(5.0)
    ),
    photonParameters = cms.untracked.PSet(
        prefixPhots = cms.untracked.string(''),
        photRelIso = cms.untracked.double(1.5),
        photTag = cms.untracked.InputTag("cleanPatPhotons"),
        debugPhots = cms.untracked.int32(0),
        photMinEt = cms.untracked.double(2.5),
        photMaxEta = cms.untracked.double(5.0),
        photMaxEt = cms.untracked.double(10000.0)
    ),
    pfleptonParameters = cms.untracked.PSet(
        muonMaxEt = cms.untracked.double(10.0),
        tauTag = cms.untracked.InputTag("selectedPatTausPF"),
        elecTag = cms.untracked.InputTag("selectedPatElectronsPF"),
        elecMinEt = cms.untracked.double(2.5),
        elecRelIso = cms.untracked.double(0.5),
        muonRelIso = cms.untracked.double(0.1),
        elecMaxEta = cms.untracked.double(5.0),
        muonTag = cms.untracked.InputTag("selectedPatMuonsPF"),
        tauMaxEta = cms.untracked.double(5.0),
        tauRelIso = cms.untracked.double(0.1),
        prefixLeps = cms.untracked.string('PF'),
        tauMinEt = cms.untracked.double(5.0),
        muonMinEt = cms.untracked.double(2.5),
        tauMaxEt = cms.untracked.double(9999.0),
        muonMaxEta = cms.untracked.double(5.0),
        debugLeps = cms.untracked.int32(0),
        elecMaxEt = cms.untracked.double(15.0)
    ),
    calometParameters = cms.untracked.PSet(
        metTag = cms.untracked.InputTag("patMETsAK5Calo"),
        prefixMET = cms.untracked.string('CaloTypeI'),
        doMCMET = cms.untracked.bool(False),
        debugMET = cms.untracked.int32(0),
        genMETTag = cms.untracked.InputTag("genMetCalo")
    )
)


process.cleanPatCandidateSummary = cms.EDAnalyzer("CandidateSummaryTable",
    logName = cms.untracked.string('cleanPatCandidates|PATSummaryTables'),
    candidates = cms.VInputTag(cms.InputTag("cleanPatElectrons"), cms.InputTag("cleanPatMuons"), cms.InputTag("cleanPatTaus"), cms.InputTag("cleanPatPhotons"), cms.InputTag("cleanPatJetsAK5Calo"), 
        cms.InputTag("cleanPatJetsAK5JPT"), cms.InputTag("cleanPatJetsAK5PF"))
)


process.cleanPatCandidateSummaryPF = cms.EDAnalyzer("CandidateSummaryTable",
    candidates = cms.VInputTag(cms.InputTag("cleanPatElectronsPF"), cms.InputTag("cleanPatMuonsPF"), cms.InputTag("cleanPatTausPF"), cms.InputTag("cleanPatPhotonsPF"), cms.InputTag("cleanPatJetsPF")),
    logName = cms.untracked.string('cleanPatCandidates|PATSummaryTables')
)


process.patCandidateSummary = cms.EDAnalyzer("CandidateSummaryTable",
    logName = cms.untracked.string('patCandidates|PATSummaryTables'),
    candidates = cms.VInputTag(cms.InputTag("patElectrons"), cms.InputTag("patMuons"), cms.InputTag("patTaus"), cms.InputTag("patPhotons"), cms.InputTag("patJets"), 
        cms.InputTag("patMETsTC"), cms.InputTag("patMETsAK5Calo"), cms.InputTag("patMHTsAK5Calo"), cms.InputTag("patJetsAK5JPT"), cms.InputTag("patJetsAK5PF"))
)


process.patCandidateSummaryPF = cms.EDAnalyzer("CandidateSummaryTable",
    candidates = cms.VInputTag(cms.InputTag("patElectronsPF"), cms.InputTag("patMuonsPF"), cms.InputTag("patTausPF"), cms.InputTag("patPhotonsPF"), cms.InputTag("patJetsPF"), 
        cms.InputTag("patMETsPF"), cms.InputTag("patPFParticlesPF")),
    logName = cms.untracked.string('patCandidates|PATSummaryTables')
)


process.selectedPatCandidateSummary = cms.EDAnalyzer("CandidateSummaryTable",
    logName = cms.untracked.string('selectedPatCanddiates|PATSummaryTables'),
    candidates = cms.VInputTag(cms.InputTag("selectedPatElectrons"), cms.InputTag("selectedPatMuons"), cms.InputTag("selectedPatTaus"), cms.InputTag("selectedPatPhotons"), cms.InputTag("selectedPatJets"), 
        cms.InputTag("selectedPatJetsAK5JPT"), cms.InputTag("selectedPatJetsAK5PF"))
)


process.selectedPatCandidateSummaryPF = cms.EDAnalyzer("CandidateSummaryTable",
    candidates = cms.VInputTag(cms.InputTag("selectedPatElectronsPF"), cms.InputTag("selectedPatMuonsPF"), cms.InputTag("selectedPatTausPF"), cms.InputTag("selectedPatPhotonsPF"), cms.InputTag("selectedPatJetsPF"), 
        cms.InputTag("selectedPatPFParticlesPF")),
    logName = cms.untracked.string('selectedPatCanddiates|PATSummaryTables')
)


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(options.outputfile),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
    ),
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_selectedPatJets*_*_*', 
        'drop patJets_selectedPatJets*_*_*', 
        'drop *_selectedPatJets_pfCandidates_*', 
        'drop *_*PF_caloTowers_*', 
        'drop *_*JPT_pfCandidates_*', 
        'drop *_*Calo_pfCandidates_*', 
        'keep *_cleanPatPhotons*_*_*', 
        'keep *_cleanPatElectrons*_*_*', 
        'keep *_cleanPatMuons*_*_*', 
        'keep *_cleanPatTaus*_*_*', 
        'keep *_cleanPatJets*_*_*', 
        'keep *_patMETs*_*_*', 
        'keep *_cleanPatHemispheres*_*_*', 
        'keep *_cleanPatPFParticles*_*_*', 
        'keep *_cleanPatTrackCands*_*_*', 
        'keep recoGenParticles_genParticles*_*_*', 
        'keep GenEventInfoProduct_*_*_*', 
        'keep GenRunInfoProduct_*_*_*', 
        'keep recoTracks_generalTracks*_*_*', 
        'keep *_towerMaker_*_*', 
        'keep *_offlineBeamSpot_*_*', 
        'keep *_offlinePrimaryVertices*_*_*', 
        'keep edmTriggerResults_TriggerResults*_*_*', 
        'keep *_hltTriggerSummaryAOD_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep edmConditionsIn*Block_conditionsInEdm_*_*', 
        'keep patTriggerAlgorithms_patTrigger_*_*', 
        'keep patTriggerConditions_patTrigger_*_*', 
        'keep patTriggerObjects_patTrigger_*_*', 
        'keep patTriggerFilters_patTrigger_*_*', 
        'keep patTriggerPaths_patTrigger_*_*', 
        'keep *_patTriggerEvent_*_*', 
        'keep *_*PatPhotons*TriggerMatch_*_*', 
        'keep *_*PatElectrons*TriggerMatch_*_*', 
        'keep *_*PatMuons*TriggerMatch_*_*', 
        'keep *_*PatTaus*TriggerMatch_*_*', 
        'keep *_*PatJets*TriggerMatch_*_*', 
        'keep *_patMETs*TriggerMatch_*_*', 
        'keep *_selectedPatMuonsPF_*_*', 
        'keep *_selectedPatElectronsPF_*_*', 
        'keep *_selectedPatTausPF_*_*', 
        'keep *_selectedPatJetsPF_*_*', 
        'keep *_cleanPatMuonsPF_*_*', 
        'keep *_cleanPatElectronsPF_*_*', 
        'keep *_cleanPatTausPF_*_*', 
        'keep *_cleanPatJetsPF_*_*', 
        'keep L1GlobalTriggerObjectMapRecord_*_*_*', 
        'keep L1GlobalTriggerReadoutRecord_*_*_*', 
        'keep recoGenJets_*GenJets*_*_*', 
        'keep recoGenMETs_*_*_*', 
        'keep edmMergeableCounter_eventCountProducer_*_*', 
        'keep recoRecoChargedRefCandidates_trackRefsForJets_*_*', 
        'keep recoTrackJets_ak5TrackJets_*_*', 
        'keep *_electronMergedSeeds_*_*', 
        'keep *_Conversions_*_*', 
        'keep recoPFCandidates_particleFlow_*_*', 
        'keep recoSuperClusters_corrected*_*_*', 
        'keep recoSuperClusters_pfElectronTranslator_*_*', 
        'keep *_gsfElectronCores_*_*', 
        'keep *_photonCore_*_*', 
        'keep recoConversions_conversions_*_*', 
        'keep recoTracks_*onversions_*_*', 
        'keep HcalNoiseSummary_*_*_*', 
        'keep *BeamHaloSummary_*_*_*', 
        'keep *_MEtoEDMConverter_*_PAT'),
    splitLevel = cms.untracked.int32(99),
    overrideInputFileSplitLevels = cms.untracked.bool(True),
    dropMetaData = cms.untracked.string('DROPPED')
)


process.JetPlusTrackCorrections = cms.Sequence(process.JPTeidTight*process.JetPlusTrackZSPCorJetIcone5)


process.patElectronTrackIsolation = cms.Sequence(process.eleIsoDepositTk*process.eleIsoFromDepsTk)


process.recoJetId = cms.Sequence(process.ak5JetID)


process.pfMuonIsoDepositsSequencePF = cms.Sequence(process.isoDepMuonWithChargedPF+process.isoDepMuonWithNeutralPF+process.isoDepMuonWithPhotonsPF)


process.patShrinkingConePFTauDiscriminationPF = cms.Sequence(process.shrinkingConePFTauDiscriminationByLeadingTrackFindingPF+process.shrinkingConePFTauDiscriminationByLeadingTrackPtCutPF+process.shrinkingConePFTauDiscriminationByLeadingPionPtCutPF+process.shrinkingConePFTauDiscriminationByIsolationPF+process.shrinkingConePFTauDiscriminationByTrackIsolationPF+process.shrinkingConePFTauDiscriminationByECALIsolationPF+process.shrinkingConePFTauDiscriminationByIsolationUsingLeadingPionPF+process.shrinkingConePFTauDiscriminationByTrackIsolationUsingLeadingPionPF+process.shrinkingConePFTauDiscriminationByECALIsolationUsingLeadingPionPF+process.shrinkingConePFTauDiscriminationAgainstElectronPF+process.shrinkingConePFTauDiscriminationAgainstMuonPF+process.shrinkingConePFTauDecayModeProducerPF+process.shrinkingConePFTauDecayModeIndexProducerPF+process.shrinkingConePFTauDiscriminationByTaNCPF+process.shrinkingConePFTauDiscriminationByTaNCfrOnePercentPF+process.shrinkingConePFTauDiscriminationByTaNCfrHalfPercentPF+process.shrinkingConePFTauDiscriminationByTaNCfrQuarterPercentPF+process.shrinkingConePFTauDiscriminationByTaNCfrTenthPercentPF)


process.makePatPhotons = cms.Sequence(process.patPhotons)


process.JetPlusTrackCorrectionsSisCone5 = cms.Sequence(process.JPTeidTight*process.JetPlusTrackZSPCorJetSiscone5)


process.recoAllGenJetsNoMuNoNu = cms.Sequence(process.sisCone5GenJetsNoMuNoNu+process.sisCone7GenJetsNoMuNoNu+process.kt4GenJetsNoMuNoNu+process.kt6GenJetsNoMuNoNu+process.iterativeCone5GenJetsNoMuNoNu+process.ak5GenJetsNoMuNoNu+process.ak7GenJetsNoMuNoNu+process.gk5GenJetsNoMuNoNu+process.gk7GenJetsNoMuNoNu+process.ca4GenJetsNoMuNoNu+process.ca6GenJetsNoMuNoNu)


process.btaggingJetTagsAK5JPT = cms.Sequence(process.jetBProbabilityBJetTagsAK5JPT+process.jetProbabilityBJetTagsAK5JPT+process.trackCountingHighPurBJetTagsAK5JPT+process.trackCountingHighEffBJetTagsAK5JPT+process.simpleSecondaryVertexHighEffBJetTagsAK5JPT+process.simpleSecondaryVertexHighPurBJetTagsAK5JPT+process.combinedSecondaryVertexBJetTagsAK5JPT+process.combinedSecondaryVertexMVABJetTagsAK5JPT+process.softMuonBJetTagsAK5JPT+process.softMuonByPtBJetTagsAK5JPT+process.softMuonByIP3dBJetTagsAK5JPT)


process.cleanPatCandidatesPF = cms.Sequence(process.cleanPatMuonsPF+process.cleanPatElectronsPF+process.cleanPatPhotonsPF+process.cleanPatTausPF+process.cleanPatJetsPF+process.cleanPatCandidateSummaryPF)


process.patJetFlavourId = cms.Sequence(process.patJetPartons*process.patJetPartonAssociation*process.patJetFlavourAssociation)


process.patTriggerSequence = cms.Sequence(process.patTrigger*process.patElectronMatch*process.cleanPatElectronsTriggerMatch*process.patMuonMatch*process.cleanPatMuonsTriggerMatch*process.patPhotonMatch*process.cleanPatPhotonsTriggerMatch*process.patJetMatchAK5Calo*process.cleanPatJetsAK5CaloTriggerMatch*process.patJetMatchAK5JPT*process.cleanPatJetsAK5JPTTriggerMatch*process.patJetMatchAK5PF*process.cleanPatJetsAK5PFTriggerMatch*process.patTauMatch*process.cleanPatTausTriggerMatch)


process.ak5JTA = cms.Sequence(process.ak5JetTracksAssociatorAtVertex*process.ak5JetTracksAssociatorAtCaloFace*process.ak5JetExtender)


process.patJetFlavourIdPF = cms.Sequence(process.patJetPartonsPF+process.patJetPartonAssociationPF+process.patJetFlavourAssociationPF)


process.patElectronHcalIsolation = cms.Sequence(process.eleIsoDepositHcalFromTowers*process.eleIsoFromDepsHcalFromTowers)


process.pfMuonIsolationFromDepositsSequencePF = cms.Sequence(process.isoValMuonWithChargedPF+process.isoValMuonWithNeutralPF+process.isoValMuonWithPhotonsPF)


process.pfTausBaseSequencePF = cms.Sequence(process.shrinkingConePFTauProducerPF+process.pfTausBaseDiscriminationByLeadingTrackFindingPF+process.pfTausBaseDiscriminationByIsolationPF+process.pfTausBaseDiscriminationByLeadingPionPtCutPF)


process.patJetCorrectionsPF = cms.Sequence(process.patJetCorrFactorsPF)


process.recoGenJets = cms.Sequence(process.kt4GenJets+process.kt6GenJets+process.iterativeCone5GenJets+process.ak5GenJets+process.ak7GenJets)


process.countPatCandidatesPF = cms.Sequence(process.countPatElectronsPF+process.countPatMuonsPF+process.countPatTausPF+process.countPatLeptonsPF+process.countPatJetsPF+process.countPatPFParticlesPF)


process.selectedPatCandidatesPF = cms.Sequence(process.selectedPatElectronsPF+process.selectedPatMuonsPF+process.selectedPatTausPF+process.selectedPatJetsPF+process.selectedPatPFParticlesPF+process.selectedPatCandidateSummaryPF)


process.cleanPatCandidates = cms.Sequence(process.cleanPatMuons*process.cleanPatElectrons*process.cleanPatPhotons*process.cleanPatTaus*process.cleanPatJetsAK5Calo*process.cleanPatJetsAK5PF*process.cleanPatJetsAK5JPT*process.cleanPatCandidateSummary)


process.jptRecoTauProducer = cms.Sequence(process.JPTeidTight*process.TCTauJetPlusTrackZSPCorJetAntiKt5*process.JPTAntiKt5JetTracksAssociatorAtVertex*process.caloRecoTauTagInfoProducer*process.JPTCaloRecoTauProducer)


process.simpleEleIdSequence = cms.Sequence(process.simpleEleId95relIso+process.simpleEleId90relIso+process.simpleEleId85relIso+process.simpleEleId80relIso+process.simpleEleId70relIso+process.simpleEleId60relIso+process.simpleEleId95cIso+process.simpleEleId90cIso+process.simpleEleId85cIso+process.simpleEleId80cIso+process.simpleEleId70cIso+process.simpleEleId60cIso)


process.patElectronId = cms.Sequence(process.eidRobustHighEnergy)


process.patHPSPFTauDiscrimination = cms.Sequence(process.hpsPFTauDiscriminationByDecayModeFinding+process.hpsPFTauDiscriminationByLooseIsolation+process.hpsPFTauDiscriminationByMediumIsolation+process.hpsPFTauDiscriminationByTightIsolation+process.hpsPFTauDiscriminationAgainstElectron+process.hpsPFTauDiscriminationAgainstMuon)


process.btaggingJetTagsAOD = cms.Sequence(process.jetBProbabilityBJetTagsAOD+process.jetProbabilityBJetTagsAOD+process.trackCountingHighPurBJetTagsAOD+process.trackCountingHighEffBJetTagsAOD+process.simpleSecondaryVertexHighEffBJetTagsAOD+process.simpleSecondaryVertexHighPurBJetTagsAOD+process.combinedSecondaryVertexBJetTagsAOD+process.combinedSecondaryVertexMVABJetTagsAOD+process.softMuonBJetTagsAOD+process.softMuonByPtBJetTagsAOD+process.softMuonByIP3dBJetTagsAOD)


process.recoAllGenJetsNoNu = cms.Sequence(process.sisCone5GenJetsNoNu+process.sisCone7GenJetsNoNu+process.kt4GenJetsNoNu+process.kt6GenJetsNoNu+process.iterativeCone5GenJetsNoNu+process.ak5GenJetsNoNu+process.ak7GenJetsNoNu+process.gk5GenJetsNoNu+process.gk7GenJetsNoNu+process.ca4GenJetsNoNu+process.ca6GenJetsNoNu)


process.pfElectronIsoDepositsSequencePF = cms.Sequence(process.isoDepElectronWithChargedPF+process.isoDepElectronWithNeutralPF+process.isoDepElectronWithPhotonsPF)


process.doAnalysisNtuplePAT = cms.Sequence(process.analysisNtuplePAT)


process.makePatMuonsPF = cms.Sequence(process.patMuonsPF)


process.pfNoPileUpSequence = cms.Sequence(process.pfPileUp+process.pfNoPileUp)


process.btaggingJetTagsAK5PF = cms.Sequence(process.jetBProbabilityBJetTagsAK5PF+process.jetProbabilityBJetTagsAK5PF+process.trackCountingHighPurBJetTagsAK5PF+process.trackCountingHighEffBJetTagsAK5PF+process.simpleSecondaryVertexHighEffBJetTagsAK5PF+process.simpleSecondaryVertexHighPurBJetTagsAK5PF+process.combinedSecondaryVertexBJetTagsAK5PF+process.combinedSecondaryVertexMVABJetTagsAK5PF+process.softMuonBJetTagsAK5PF+process.softMuonByPtBJetTagsAK5PF+process.softMuonByIP3dBJetTagsAK5PF)


process.RunTanc = cms.Sequence(process.shrinkingConePFTauDiscriminationByTaNCfrOnePercent+process.shrinkingConePFTauDiscriminationByTaNCfrHalfPercent+process.shrinkingConePFTauDiscriminationByTaNCfrQuarterPercent+process.shrinkingConePFTauDiscriminationByTaNCfrTenthPercent)


process.patPFCandidateIsoDepositSelection = cms.Sequence(process.pfNoPileUpSequence*(process.pfAllNeutralHadrons+process.pfAllChargedHadrons+process.pfAllPhotons))


process.JetPlusTrackCorrectionsAntiKt5 = cms.Sequence(process.JPTeidTight*process.JetPlusTrackZSPCorJetAntiKt5)


process.patTriggerEventSequence = cms.Sequence(process.patTriggerEvent)


process.cleanupFilterMC = cms.Sequence(process.removePKAM+process.primaryVertexFilter+process.hbheNoise)


process.pfJetSequence = cms.Sequence(process.pfJets)


process.patCaloTauDiscrimination = cms.Sequence(process.caloRecoTauDiscriminationAgainstElectron+process.caloRecoTauDiscriminationByIsolation+process.caloRecoTauDiscriminationByLeadingTrackFinding+process.caloRecoTauDiscriminationByLeadingTrackPtCut)


process.genForPF2PATSequence = cms.Sequence(process.genParticlesForJetsNoNu+process.ak5GenJetsNoNu+process.ak7GenJetsNoNu)


process.countPatCandidates = cms.Sequence(process.countPatElectrons+process.countPatMuons+process.countPatTaus+process.countPatLeptons+process.countPatPhotons+process.countPatJets*process.countPatJetsAK5PF*process.countPatJetsAK5JPT)


process.genJetParticles = cms.Sequence(process.genParticlesForJets)


process.btaggingTagInfosAK5JPT = cms.Sequence(process.impactParameterTagInfosAK5JPT+process.secondaryVertexTagInfosAK5JPT+process.softMuonTagInfosAK5JPT+process.btaggingJetTagsAK5JPT)


process.produceAndDiscriminateHPSPFTaus = cms.Sequence(process.hpsPFTauProducer*process.hpsPFTauDiscriminationByDecayModeFinding*process.hpsPFTauDiscriminationByLooseIsolation*process.hpsPFTauDiscriminationByMediumIsolation*process.hpsPFTauDiscriminationByTightIsolation*process.hpsPFTauDiscriminationAgainstElectron*process.hpsPFTauDiscriminationAgainstMuon)


process.selectedPatCandidates = cms.Sequence(process.selectedPatElectrons+process.selectedPatMuons+process.selectedPatTaus+process.selectedPatPhotons+process.selectedPatJets*process.selectedPatJetsAK5PF*process.selectedPatJetsAK5JPT+process.selectedPatCandidateSummary)


process.pfTauSequencePF = cms.Sequence(process.pfJetTracksAssociatorAtVertexPF+process.pfRecoTauTagInfoProducerPF+process.pfTausBaseSequencePF+process.pfTausPF)


process.makePatMuons = cms.Sequence(process.patMuons)


process.theNtupler = cms.Sequence(process.doAnalysisNtuplePAT)


process.pfSortByTypeSequence = cms.Sequence(process.pfAllNeutralHadrons+process.pfAllChargedHadrons+process.pfAllPhotons+process.pfAllElectrons+process.pfAllMuons)


process.produceAndDiscriminateShrinkingConePFTaus = cms.Sequence(process.shrinkingConePFTauProducer*process.shrinkingConePFTauDecayModeProducer*process.shrinkingConePFTauDecayModeIndexProducer*process.shrinkingConePFTauDiscriminationByLeadingTrackFinding*process.shrinkingConePFTauDiscriminationByLeadingTrackPtCut*process.shrinkingConePFTauDiscriminationByLeadingPionPtCut*process.shrinkingConePFTauDiscriminationByIsolation*process.shrinkingConePFTauDiscriminationByTrackIsolation*process.shrinkingConePFTauDiscriminationByECALIsolation*process.shrinkingConePFTauDiscriminationByIsolationUsingLeadingPion*process.shrinkingConePFTauDiscriminationByTrackIsolationUsingLeadingPion*process.shrinkingConePFTauDiscriminationByECALIsolationUsingLeadingPion*process.shrinkingConePFTauDiscriminationAgainstElectron*process.shrinkingConePFTauDiscriminationAgainstMuon)


process.pfElectronIsoDepositsSequence = cms.Sequence(process.isoDepElectronWithCharged+process.isoDepElectronWithNeutral+process.isoDepElectronWithPhotons)


process.patShrinkingConePFTauDiscrimination = cms.Sequence(process.shrinkingConePFTauDiscriminationByLeadingTrackFinding+process.shrinkingConePFTauDiscriminationByLeadingTrackPtCut+process.shrinkingConePFTauDiscriminationByLeadingPionPtCut+process.shrinkingConePFTauDiscriminationByIsolation+process.shrinkingConePFTauDiscriminationByTrackIsolation+process.shrinkingConePFTauDiscriminationByECALIsolation+process.shrinkingConePFTauDiscriminationByIsolationUsingLeadingPion+process.shrinkingConePFTauDiscriminationByTrackIsolationUsingLeadingPion+process.shrinkingConePFTauDiscriminationByECALIsolationUsingLeadingPion+process.shrinkingConePFTauDiscriminationAgainstElectron+process.shrinkingConePFTauDiscriminationAgainstMuon+process.shrinkingConePFTauDecayModeProducer+process.shrinkingConePFTauDecayModeIndexProducer+process.shrinkingConePFTauDiscriminationByTaNC+process.shrinkingConePFTauDiscriminationByTaNCfrOnePercent+process.shrinkingConePFTauDiscriminationByTaNCfrHalfPercent+process.shrinkingConePFTauDiscriminationByTaNCfrQuarterPercent+process.shrinkingConePFTauDiscriminationByTaNCfrTenthPercent)


process.patPFTauIsolationPF = cms.Sequence(process.tauIsoDepositPFCandidatesPF+process.tauIsoDepositPFChargedHadronsPF+process.tauIsoDepositPFNeutralHadronsPF+process.tauIsoDepositPFGammasPF)


process.btaggingTagInfosAOD = cms.Sequence(process.impactParameterTagInfosAOD+process.secondaryVertexTagInfosAOD+process.softMuonTagInfosAOD+process.btaggingJetTagsAOD)


process.patPhotonTrackIsolation = cms.Sequence(process.gamIsoDepositTk*process.gamIsoFromDepsTk)


process.recoAllGenJets = cms.Sequence(process.sisCone5GenJets+process.sisCone7GenJets+process.kt4GenJets+process.kt6GenJets+process.iterativeCone5GenJets+process.ak5GenJets+process.ak7GenJets+process.gk5GenJets+process.gk7GenJets+process.ca4GenJets+process.ca6GenJets)


process.patPhotonEcalIsolation = cms.Sequence(process.gamIsoDepositEcalFromHits*process.gamIsoFromDepsEcalFromHits)


process.pfJetSequencePF = cms.Sequence(process.pfJetsPF)


process.pfMuonIsoDepositsSequence = cms.Sequence(process.isoDepMuonWithCharged+process.isoDepMuonWithNeutral+process.isoDepMuonWithPhotons)


process.pfNoPileUpSequencePF = cms.Sequence(process.pfPileUpPF+process.pfNoPileUpPF)


process.patFixedConePFTauDiscrimination = cms.Sequence(process.fixedConePFTauDiscriminationByLeadingTrackFinding+process.fixedConePFTauDiscriminationByLeadingTrackPtCut+process.fixedConePFTauDiscriminationByLeadingPionPtCut+process.fixedConePFTauDiscriminationByIsolation+process.fixedConePFTauDiscriminationByTrackIsolation+process.fixedConePFTauDiscriminationByECALIsolation+process.fixedConePFTauDiscriminationByIsolationUsingLeadingPion+process.fixedConePFTauDiscriminationByTrackIsolationUsingLeadingPion+process.fixedConePFTauDiscriminationByECALIsolationUsingLeadingPion+process.fixedConePFTauDiscriminationAgainstElectron+process.fixedConePFTauDiscriminationAgainstMuon)


process.makePatMETsPF = cms.Sequence(process.patMETsPF+process.metJESCorAK5PFTypeI+process.patMETsTypeIPF)


process.pfMuonIsolationFromDepositsSequence = cms.Sequence(process.isoValMuonWithCharged+process.isoValMuonWithNeutral+process.isoValMuonWithPhotons)


process.produceAndDiscriminateFixedConePFTaus = cms.Sequence(process.fixedConePFTauProducer*process.fixedConePFTauDiscriminationByLeadingTrackFinding*process.fixedConePFTauDiscriminationByLeadingTrackPtCut*process.fixedConePFTauDiscriminationByLeadingPionPtCut*process.fixedConePFTauDiscriminationByIsolation*process.fixedConePFTauDiscriminationByTrackIsolation*process.fixedConePFTauDiscriminationByECALIsolation*process.fixedConePFTauDiscriminationByIsolationUsingLeadingPion*process.fixedConePFTauDiscriminationByTrackIsolationUsingLeadingPion*process.fixedConePFTauDiscriminationByECALIsolationUsingLeadingPion*process.fixedConePFTauDiscriminationAgainstElectron*process.fixedConePFTauDiscriminationAgainstMuon)


process.makePatElectronsPF = cms.Sequence(process.patElectronsPF)


process.JetPlusTrackCorrectionsIcone5 = cms.Sequence(process.JPTeidTight*process.JetPlusTrackZSPCorJetIcone5)


process.btagging = cms.Sequence(process.impactParameterTagInfos*(process.trackCountingHighEffBJetTags+process.trackCountingHighPurBJetTags+process.jetProbabilityBJetTags+process.jetBProbabilityBJetTags+process.secondaryVertexTagInfos*(process.simpleSecondaryVertexHighEffBJetTags+process.simpleSecondaryVertexHighPurBJetTags+process.combinedSecondaryVertexBJetTags+process.combinedSecondaryVertexMVABJetTags)+process.ghostTrackVertexTagInfos*process.ghostTrackBJetTags)+process.softElectronCands*process.softElectronTagInfos*(process.softElectronByIP3dBJetTags+process.softElectronByPtBJetTags)+process.softMuonTagInfos*(process.softMuonBJetTags+process.softMuonByIP3dBJetTags+process.softMuonByPtBJetTags))


process.tautagging = cms.Sequence(process.jptRecoTauProducer*process.caloRecoTauProducer*process.caloRecoTauDiscriminationByLeadingTrackFinding*process.caloRecoTauDiscriminationByLeadingTrackPtCut*process.caloRecoTauDiscriminationByTrackIsolation*process.caloRecoTauDiscriminationByECALIsolation*process.caloRecoTauDiscriminationByIsolation*process.caloRecoTauDiscriminationAgainstElectron*process.caloRecoTauDiscriminationAgainstMuon)


process.pfTausBaseSequence = cms.Sequence(process.shrinkingConePFTauProducer+process.pfTausDiscriminationByLeadingTrackFinding+process.pfTausDiscriminationByLeadingPionPtCut+process.pfTausDiscriminationByIsolation)


process.patPFTauIsolation = cms.Sequence(process.tauIsoDepositPFCandidates*process.tauIsoDepositPFChargedHadrons*process.tauIsoDepositPFNeutralHadrons*process.tauIsoDepositPFGammas)


process.patJetCorrections = cms.Sequence(process.patJetCorrFactors*process.patJetCorrFactorsAK5PF*process.patJetCorrFactorsAK5JPT)


process.hiRecoGenJets = cms.Sequence(process.iterativeCone5HiGenJets+process.iterativeCone7HiGenJets+process.ak5HiGenJets+process.ak7HiGenJets+process.kt4HiGenJets+process.kt6HiGenJets)


process.makePatMHTs = cms.Sequence(process.patMHTs)


process.pfMuonIsolationSequencePF = cms.Sequence(process.pfMuonIsoDepositsSequencePF+process.pfMuonIsolationFromDepositsSequencePF)


process.pfTauSequence = cms.Sequence(process.pfJetTracksAssociatorAtVertex+process.pfRecoTauTagInfoProducer+process.pfTausBaseSequence+process.pfTaus)


process.MetType1Corrections = cms.Sequence(process.metJESCorIC5CaloJet*process.metJESCorKT4CaloJet*process.metJESCorKT6CaloJet*process.metJESCorAK5CaloJet*process.metJESCorAK7CaloJet*process.metJESCorSC5CaloJet*process.metJESCorSC7CaloJet)


process.eIdSequence = cms.Sequence(process.eidRobustLoose+process.eidRobustTight+process.eidRobustHighEnergy+process.eidLoose+process.eidTight)


process.patElectronEcalIsolation = cms.Sequence(process.eleIsoDepositEcalFromHits*process.eleIsoFromDepsEcalFromHitsByCrystal)


process.btaggingJetTagsAODPF = cms.Sequence(process.jetBProbabilityBJetTagsAODPF+process.jetProbabilityBJetTagsAODPF+process.trackCountingHighPurBJetTagsAODPF+process.trackCountingHighEffBJetTagsAODPF+process.simpleSecondaryVertexHighEffBJetTagsAODPF+process.simpleSecondaryVertexHighPurBJetTagsAODPF+process.combinedSecondaryVertexBJetTagsAODPF+process.combinedSecondaryVertexMVABJetTagsAODPF+process.softMuonBJetTagsAODPF+process.softMuonByPtBJetTagsAODPF+process.softMuonByIP3dBJetTagsAODPF)


process.pfMuonIsolationSequence = cms.Sequence(process.pfMuonIsoDepositsSequence+process.pfMuonIsolationFromDepositsSequence)


process.makePatPhotonsPF = cms.Sequence()


process.pfElectronIsolationFromDepositsSequencePF = cms.Sequence(process.isoValElectronWithChargedPF+process.isoValElectronWithNeutralPF+process.isoValElectronWithPhotonsPF)


process.patMETCorrectionsPF = cms.Sequence(process.metJESCorAK5CaloJetPF+process.metJESCorAK5CaloJetMuonsPF)


process.patPhotonHcalIsolation = cms.Sequence(process.gamIsoDepositHcalFromTowers*process.gamIsoFromDepsHcalFromTowers)


process.patMETCorrections = cms.Sequence((process.metJESCorAK5CaloJetTypeII*process.metJESCorAK5CaloJetMuonsTypeII+process.metJESCorAK5CaloJet)*process.metJESCorAK5CaloJetMuons)


process.pfElectronIsolationFromDepositsSequence = cms.Sequence(process.isoValElectronWithCharged+process.isoValElectronWithNeutral+process.isoValElectronWithPhotons)


process.produceShrinkingConeDiscriminationByTauNeuralClassifier = cms.Sequence(process.shrinkingConePFTauDiscriminationByTaNC*process.shrinkingConePFTauDiscriminationByTaNCfrOnePercent*process.shrinkingConePFTauDiscriminationByTaNCfrHalfPercent*process.shrinkingConePFTauDiscriminationByTaNCfrQuarterPercent*process.shrinkingConePFTauDiscriminationByTaNCfrTenthPercent)


process.makePatTausPF = cms.Sequence(process.patPFTauIsolationPF+process.patShrinkingConePFTauDiscriminationPF*process.patTausPF)


process.btaggingAOD = cms.Sequence(process.btaggingTagInfosAOD)


process.pfMuonSequence = cms.Sequence(process.pfAllMuons+process.pfMuonsFromVertex+process.pfSelectedMuons+process.pfMuonIsolationSequence+process.pfIsolatedMuons)


process.cleanupFilterData = cms.Sequence(process.physicsDeclared+process.cleanupFilterMC)


process.btaggingAK5JPT = cms.Sequence(process.btaggingTagInfosAK5JPT)


process.makePatTaus = cms.Sequence(process.patPFCandidateIsoDepositSelection*process.patPFTauIsolation*(process.shrinkingConePFTauDecayModeProducer+process.patTaus))


process.makePatElectrons = cms.Sequence(process.simpleEleIdSequence+process.patElectrons)


process.PFTau = cms.Sequence(process.ak5PFJetTracksAssociatorAtVertex*process.pfRecoTauTagInfoProducer*process.produceAndDiscriminateShrinkingConePFTaus+process.produceShrinkingConeDiscriminationByTauNeuralClassifier+process.produceAndDiscriminateFixedConePFTaus+process.produceAndDiscriminateHPSPFTaus)


process.btaggingTagInfosAK5PF = cms.Sequence(process.impactParameterTagInfosAK5PF+process.secondaryVertexTagInfosAK5PF+process.softMuonTagInfosAK5PF+process.btaggingJetTagsAK5PF)


process.cleanEvents = cms.Sequence(process.cleanupFilterData)


process.pfMuonSequencePF = cms.Sequence(process.pfAllMuonsPF+process.pfMuonsFromVertexPF+process.pfSelectedMuonsPF+process.pfMuonIsolationSequencePF+process.pfIsolatedMuonsPF)


process.patElectronIsolation = cms.Sequence(process.patElectronTrackIsolation+process.patElectronEcalIsolation+process.patElectronHcalIsolation)


process.patJetMETCorrections = cms.Sequence(process.patJetCorrections)


process.patPhotonIsolation = cms.Sequence(process.patPhotonTrackIsolation+process.patPhotonEcalIsolation+process.patPhotonHcalIsolation)


process.btaggingTagInfosAODPF = cms.Sequence(process.impactParameterTagInfosAODPF+process.secondaryVertexTagInfosAODPF+process.softMuonTagInfosAODPF+process.btaggingJetTagsAODPF)


process.selectClean = cms.Sequence(process.cleanEvents)


process.patPFCandidateIsoDepositSelectionPF = cms.Sequence(process.pfNoPileUpSequencePF+process.pfAllNeutralHadronsPF+process.pfAllChargedHadronsPF+process.pfAllPhotonsPF)


process.makePatMETs = cms.Sequence(process.patMETCorrections*(process.patMETsAK5Calo+process.patMETsAK5CaloTypeII)*process.patMETsTC)


process.pfElectronIsolationSequencePF = cms.Sequence(process.pfElectronIsoDepositsSequencePF+process.pfElectronIsolationFromDepositsSequencePF)


process.pfElectronIsolationSequence = cms.Sequence(process.pfElectronIsoDepositsSequence+process.pfElectronIsolationFromDepositsSequence)


process.pfElectronSequence = cms.Sequence(process.pfAllElectrons+process.pfElectronsFromVertex+process.pfSelectedElectrons+process.pfElectronIsolationSequence+process.pfIsolatedElectrons)


process.btaggingAK5PF = cms.Sequence(process.btaggingTagInfosAK5PF)


process.pfElectronSequencePF = cms.Sequence(process.pfAllElectronsPF+process.pfElectronsFromVertexPF+process.pfSelectedElectronsPF+process.pfElectronIsolationSequencePF+process.pfIsolatedElectronsPF)


process.btaggingAODPF = cms.Sequence(process.btaggingTagInfosAODPF)


process.PF2PAT = cms.Sequence(process.pfMET+process.pfNoPileUpSequence+process.pfAllNeutralHadrons+process.pfAllChargedHadrons+process.pfAllPhotons+process.pfMuonSequence+process.pfNoMuon+process.pfElectronSequence+process.pfNoElectron+process.pfJetSequence+process.pfNoJet+process.pfTauSequence+process.pfNoTau)


process.makePatJetsPF = cms.Sequence(process.patJetCorrectionsPF+process.jetTracksAssociatorAtVertexPF+process.btaggingAODPF+process.patJetChargePF+process.patJetsPF)


process.patCandidatesPF = cms.Sequence(process.makePatElectronsPF+process.makePatMuonsPF+process.makePatTausPF+process.makePatJetsPF+process.makePatMETsPF+process.patPFParticlesPF+process.patCandidateSummaryPF)


process.makePatJets = cms.Sequence(process.patJetCorrections*(process.jetTracksAssociatorAtVertex+process.btaggingAOD+process.jetTracksAssociatorAtVertexAK5JPT+process.btaggingAK5JPT+(process.jetTracksAssociatorAtVertexAK5PF+process.btaggingAK5PF+process.patJetCharge*process.patJetChargeAK5PF)*process.patJetChargeAK5JPT)*process.patJets*process.patJetsAK5PF*process.patJetsAK5JPT)


process.PF2PATPF = cms.Sequence(process.pfMETPF+process.pfNoPileUpSequencePF+process.pfAllNeutralHadronsPF+process.pfAllChargedHadronsPF+process.pfAllPhotonsPF+process.pfMuonSequencePF+process.pfNoMuonPF+process.pfElectronSequencePF+process.pfNoElectronPF+process.pfJetSequencePF+process.pfNoJetPF+process.pfTauSequencePF+process.pfNoTauPF)


process.patDefaultSequencePF = cms.Sequence(process.patCandidatesPF+process.selectedPatCandidatesPF+process.countPatCandidatesPF+process.patTriggerPF+process.patTriggerEventPF+process.patElectronMatchPF+process.selectedPatElectronsTriggerMatchPF+process.patMuonMatchPF+process.selectedPatMuonsTriggerMatchPF+process.patJetMatchPF+process.selectedPatJetsTriggerMatchPF+process.patTauMatchPF+process.selectedPatTausTriggerMatchPF)


process.patCandidates = cms.Sequence(process.makePatElectrons+process.makePatMuons+process.makePatTaus+process.makePatPhotons+process.makePatJets+process.makePatMETs+process.patCandidateSummary)


process.patPF2PATSequencePF = cms.Sequence(process.PF2PATPF+process.patDefaultSequencePF)


process.patDefaultSequence = cms.Sequence(process.patCandidates*process.selectedPatCandidates*process.cleanPatCandidates*process.countPatCandidates*process.patTriggerSequence*process.patTriggerEventSequence)


process.patPF2PATSequence = cms.Sequence(process.PF2PAT+process.patDefaultSequence)


process.susyPatDefaultSequence = cms.Sequence(process.eventCountProducer*process.patDefaultSequence*process.patPF2PATSequencePF)


process.susyStep = cms.Sequence(process.susyPatDefaultSequence)


process.p = cms.Path(process.cleanEvents*process.susyStep*process.theNtupler)

#process.p = cms.Path(process.cleanEvents*process.theNtupler)
#process.p = cms.Path(process.cleanEvents*process.susyStep)
#process.outpath = cms.EndPath(process.out)


process.MessageLogger = cms.Service("MessageLogger",
    suppressInfo = cms.untracked.vstring(),
    debugs = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    suppressDebug = cms.untracked.vstring(),
    cout = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    cerr_stats = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING'),
        output = cms.untracked.string('cerr'),
        optionalPSet = cms.untracked.bool(True)
    ),
    warnings = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    default = cms.untracked.PSet(

    ),
    statistics = cms.untracked.vstring('cerr_stats'),
    cerr = cms.untracked.PSet(
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        noTimeStamps = cms.untracked.bool(False),
        FwkReport = cms.untracked.PSet(
            reportEvery = cms.untracked.int32(10000),
            optionalPSet = cms.untracked.bool(True),
            limit = cms.untracked.int32(10000000)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000)
        ),
        Root_NoDictionary = cms.untracked.PSet(
            optionalPSet = cms.untracked.bool(True),
            limit = cms.untracked.int32(0)
        ),
        threshold = cms.untracked.string('INFO'),
        FwkJob = cms.untracked.PSet(
            optionalPSet = cms.untracked.bool(True),
            limit = cms.untracked.int32(0)
        ),
        FwkSummary = cms.untracked.PSet(
            reportEvery = cms.untracked.int32(1),
            optionalPSet = cms.untracked.bool(True),
            limit = cms.untracked.int32(10000000)
        ),
        optionalPSet = cms.untracked.bool(True),
        PATSummaryTables = cms.untracked.PSet(
            reportEvery = cms.untracked.int32(1),
            limit = cms.untracked.int32(-1)
        )
    ),
    FrameworkJobReport = cms.untracked.PSet(
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        optionalPSet = cms.untracked.bool(True),
        FwkJob = cms.untracked.PSet(
            optionalPSet = cms.untracked.bool(True),
            limit = cms.untracked.int32(10000000)
        )
    ),
    suppressWarning = cms.untracked.vstring(),
    errors = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    destinations = cms.untracked.vstring('warnings', 
        'errors', 
        'infos', 
        'debugs', 
        'cout', 
        'cerr'),
    debugModules = cms.untracked.vstring(),
    infos = cms.untracked.PSet(
        optionalPSet = cms.untracked.bool(True),
        Root_NoDictionary = cms.untracked.PSet(
            optionalPSet = cms.untracked.bool(True),
            limit = cms.untracked.int32(0)
        ),
        placeholder = cms.untracked.bool(True)
    ),
    categories = cms.untracked.vstring('FwkJob', 
        'FwkReport', 
        'FwkSummary', 
        'Root_NoDictionary', 
        'PATSummaryTables'),
    fwkJobReports = cms.untracked.vstring('FrameworkJobReport')
)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.patout)
)


process.CSCGeometryESModule = cms.ESProducer("CSCGeometryESModule",
    appendToDataLabel = cms.string(''),
    useDDD = cms.bool(True),
    debugV = cms.untracked.bool(False),
    useGangedStripsInME1a = cms.bool(True),
    alignmentsLabel = cms.string(''),
    useOnlyWiresInME1a = cms.bool(False),
    useRealWireGeometry = cms.bool(True),
    useCentreTIOffsets = cms.bool(False),
    applyAlignment = cms.bool(True)
)


process.CaloGeometryBuilder = cms.ESProducer("CaloGeometryBuilder",
    SelectedCalos = cms.vstring('HCAL', 
        'ZDC', 
        'CASTOR', 
        'EcalBarrel', 
        'EcalEndcap', 
        'EcalPreshower', 
        'TOWER')
)


process.CaloTopologyBuilder = cms.ESProducer("CaloTopologyBuilder")


process.CaloTowerHardcodeGeometryEP = cms.ESProducer("CaloTowerHardcodeGeometryEP")


process.CastorHardcodeGeometryEP = cms.ESProducer("CastorHardcodeGeometryEP")


process.DTGeometryESModule = cms.ESProducer("DTGeometryESModule",
    appendToDataLabel = cms.string(''),
    fromDDD = cms.bool(True),
    applyAlignment = cms.bool(True),
    alignmentsLabel = cms.string('')
)


process.EcalBarrelGeometryEP = cms.ESProducer("EcalBarrelGeometryEP",
    applyAlignment = cms.bool(False)
)


process.EcalElectronicsMappingBuilder = cms.ESProducer("EcalElectronicsMappingBuilder")


process.EcalEndcapGeometryEP = cms.ESProducer("EcalEndcapGeometryEP",
    applyAlignment = cms.bool(False)
)


process.EcalLaserCorrectionService = cms.ESProducer("EcalLaserCorrectionService")


process.EcalPreshowerGeometryEP = cms.ESProducer("EcalPreshowerGeometryEP",
    applyAlignment = cms.bool(False)
)


process.EcalTrigTowerConstituentsMapBuilder = cms.ESProducer("EcalTrigTowerConstituentsMapBuilder",
    MapFile = cms.untracked.string('Geometry/EcalMapping/data/EndCap_TTMap.txt')
)


process.GlobalTrackingGeometryESProducer = cms.ESProducer("GlobalTrackingGeometryESProducer")


process.HcalHardcodeGeometryEP = cms.ESProducer("HcalHardcodeGeometryEP")


process.HcalTopologyIdealEP = cms.ESProducer("HcalTopologyIdealEP")


process.MuonDetLayerGeometryESProducer = cms.ESProducer("MuonDetLayerGeometryESProducer")


process.MuonNumberingInitialization = cms.ESProducer("MuonNumberingInitialization")


process.ParametrizedMagneticFieldProducer = cms.ESProducer("ParametrizedMagneticFieldProducer",
    version = cms.string('OAE_1103l_071212'),
    parameters = cms.PSet(
        BValue = cms.string('3_8T')
    ),
    label = cms.untracked.string('parametrizedField')
)


process.RPCGeometryESModule = cms.ESProducer("RPCGeometryESModule",
    useDDD = cms.untracked.bool(True),
    compatibiltyWith11 = cms.untracked.bool(True)
)


process.SiStripRecHitMatcherESProducer = cms.ESProducer("SiStripRecHitMatcherESProducer",
    ComponentName = cms.string('StandardMatcher'),
    NSigmaInside = cms.double(3.0)
)


process.SteppingHelixPropagatorAlong = cms.ESProducer("SteppingHelixPropagatorESProducer",
    NoErrorPropagation = cms.bool(False),
    PropagationDirection = cms.string('alongMomentum'),
    useMatVolumes = cms.bool(True),
    useTuningForL2Speed = cms.bool(False),
    useIsYokeFlag = cms.bool(True),
    endcapShiftInZNeg = cms.double(0.0),
    SetVBFPointer = cms.bool(False),
    AssumeNoMaterial = cms.bool(False),
    returnTangentPlane = cms.bool(True),
    useInTeslaFromMagField = cms.bool(False),
    VBFName = cms.string('VolumeBasedMagneticField'),
    useEndcapShiftsInZ = cms.bool(False),
    sendLogWarning = cms.bool(False),
    ComponentName = cms.string('SteppingHelixPropagatorAlong'),
    debug = cms.bool(False),
    ApplyRadX0Correction = cms.bool(True),
    useMagVolumes = cms.bool(True),
    endcapShiftInZPos = cms.double(0.0)
)


process.StripCPEfromTrackAngleESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('StripCPEfromTrackAngle')
)


process.TrackerDigiGeometryESModule = cms.ESProducer("TrackerDigiGeometryESModule",
    appendToDataLabel = cms.string(''),
    fromDDD = cms.bool(True),
    applyAlignment = cms.bool(True),
    alignmentsLabel = cms.string('')
)


process.TrackerGeometricDetESModule = cms.ESProducer("TrackerGeometricDetESModule",
    fromDDD = cms.bool(True)
)


process.TrackerRecoGeometryESProducer = cms.ESProducer("TrackerRecoGeometryESProducer")


process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName = cms.string('TransientTrackBuilder')
)


process.VolumeBasedMagneticFieldESProducer = cms.ESProducer("VolumeBasedMagneticFieldESProducer",
    scalingVolumes = cms.vint32(14100, 14200, 17600, 17800, 17900, 
        18100, 18300, 18400, 18600, 23100, 
        23300, 23400, 23600, 23800, 23900, 
        24100, 28600, 28800, 28900, 29100, 
        29300, 29400, 29600, 28609, 28809, 
        28909, 29109, 29309, 29409, 29609, 
        28610, 28810, 28910, 29110, 29310, 
        29410, 29610, 28611, 28811, 28911, 
        29111, 29311, 29411, 29611),
    scalingFactors = cms.vdouble(1, 1, 0.994, 1.004, 1.004, 
        1.005, 1.004, 1.004, 0.994, 0.965, 
        0.958, 0.958, 0.953, 0.958, 0.958, 
        0.965, 0.918, 0.924, 0.924, 0.906, 
        0.924, 0.924, 0.918, 0.991, 0.998, 
        0.998, 0.978, 0.998, 0.998, 0.991, 
        0.991, 0.998, 0.998, 0.978, 0.998, 
        0.998, 0.991, 0.991, 0.998, 0.998, 
        0.978, 0.998, 0.998, 0.991),
    overrideMasterSector = cms.bool(False),
    useParametrizedTrackerField = cms.bool(True),
    label = cms.untracked.string(''),
    version = cms.string('grid_1103l_090322_3_8t'),
    debugBuilder = cms.untracked.bool(False),
    paramLabel = cms.string('parametrizedField'),
    cacheLastVolume = cms.untracked.bool(True)
)


process.ZdcHardcodeGeometryEP = cms.ESProducer("ZdcHardcodeGeometryEP")


process.caloDetIdAssociator = cms.ESProducer("DetIdAssociatorESProducer",
    ComponentName = cms.string('CaloDetIdAssociator'),
    etaBinSize = cms.double(0.087),
    nEta = cms.int32(70),
    nPhi = cms.int32(72)
)


process.combinedMVA = cms.ESProducer("CombinedMVAJetTagESProducer",
    useCategories = cms.bool(False),
    calibrationRecord = cms.string('CombinedMVA'),
    jetTagComputers = cms.VPSet(cms.PSet(
        discriminator = cms.bool(True),
        variables = cms.bool(False),
        jetTagComputer = cms.string('jetProbability')
    ), 
        cms.PSet(
            discriminator = cms.bool(True),
            variables = cms.bool(False),
            jetTagComputer = cms.string('combinedSecondaryVertex')
        ), 
        cms.PSet(
            discriminator = cms.bool(True),
            variables = cms.bool(False),
            jetTagComputer = cms.string('softMuon')
        ), 
        cms.PSet(
            discriminator = cms.bool(True),
            variables = cms.bool(False),
            jetTagComputer = cms.string('softElectron')
        ))
)


process.combinedSecondaryVertex = cms.ESProducer("CombinedSecondaryVertexESProducer",
    useTrackWeights = cms.bool(True),
    pseudoMultiplicityMin = cms.uint32(2),
    correctVertexMass = cms.bool(True),
    trackSelection = cms.PSet(
        totalHitsMin = cms.uint32(0),
        jetDeltaRMax = cms.double(0.3),
        qualityClass = cms.string('highPurity'),
        pixelHitsMin = cms.uint32(0),
        sip3dSigMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        maxDistToAxis = cms.double(0.07),
        sip2dValMax = cms.double(99999.9),
        maxDecayLen = cms.double(5),
        ptMin = cms.double(0.0),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        sip2dValMin = cms.double(-99999.9),
        normChi2Max = cms.double(99999.9)
    ),
    pseudoVertexV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.05)
    ),
    charmCut = cms.double(1.5),
    vertexFlip = cms.bool(False),
    minimumTrackWeight = cms.double(0.5),
    trackPairV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.03)
    ),
    trackMultiplicityMin = cms.uint32(3),
    trackPseudoSelection = cms.PSet(
        totalHitsMin = cms.uint32(0),
        jetDeltaRMax = cms.double(0.3),
        qualityClass = cms.string('highPurity'),
        pixelHitsMin = cms.uint32(0),
        sip3dSigMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        maxDistToAxis = cms.double(0.07),
        sip2dValMax = cms.double(99999.9),
        maxDecayLen = cms.double(5),
        ptMin = cms.double(0.0),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(2.0),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        sip2dValMin = cms.double(-99999.9),
        normChi2Max = cms.double(99999.9)
    ),
    trackSort = cms.string('sip2dSig'),
    trackFlip = cms.bool(False),
    calibrationRecords = cms.vstring('CombinedSVRecoVertex', 
        'CombinedSVPseudoVertex', 
        'CombinedSVNoVertex'),
    useCategories = cms.bool(True),
    categoryVariableName = cms.string('vertexCategory')
)


process.combinedSecondaryVertexMVA = cms.ESProducer("CombinedSecondaryVertexESProducer",
    useTrackWeights = cms.bool(True),
    pseudoMultiplicityMin = cms.uint32(2),
    correctVertexMass = cms.bool(True),
    trackSelection = cms.PSet(
        totalHitsMin = cms.uint32(0),
        jetDeltaRMax = cms.double(0.3),
        qualityClass = cms.string('highPurity'),
        pixelHitsMin = cms.uint32(0),
        sip3dSigMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        maxDistToAxis = cms.double(0.07),
        sip2dValMax = cms.double(99999.9),
        maxDecayLen = cms.double(5),
        ptMin = cms.double(0.0),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        sip2dValMin = cms.double(-99999.9),
        normChi2Max = cms.double(99999.9)
    ),
    pseudoVertexV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.05)
    ),
    charmCut = cms.double(1.5),
    vertexFlip = cms.bool(False),
    minimumTrackWeight = cms.double(0.5),
    trackPairV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.03)
    ),
    trackMultiplicityMin = cms.uint32(3),
    trackPseudoSelection = cms.PSet(
        totalHitsMin = cms.uint32(0),
        jetDeltaRMax = cms.double(0.3),
        qualityClass = cms.string('highPurity'),
        pixelHitsMin = cms.uint32(0),
        sip3dSigMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        maxDistToAxis = cms.double(0.07),
        sip2dValMax = cms.double(99999.9),
        maxDecayLen = cms.double(5),
        ptMin = cms.double(0.0),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(2.0),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        sip2dValMin = cms.double(-99999.9),
        normChi2Max = cms.double(99999.9)
    ),
    trackSort = cms.string('sip2dSig'),
    trackFlip = cms.bool(False),
    calibrationRecords = cms.vstring('CombinedSVMVARecoVertex', 
        'CombinedSVMVAPseudoVertex', 
        'CombinedSVMVANoVertex'),
    useCategories = cms.bool(True),
    categoryVariableName = cms.string('vertexCategory')
)


process.ecalDetIdAssociator = cms.ESProducer("DetIdAssociatorESProducer",
    ComponentName = cms.string('EcalDetIdAssociator'),
    etaBinSize = cms.double(0.02),
    nEta = cms.int32(300),
    nPhi = cms.int32(360)
)


process.fakeForIdealAlignment = cms.ESProducer("FakeAlignmentProducer",
    appendToDataLabel = cms.string('fakeForIdeal')
)


process.ghostTrack = cms.ESProducer("GhostTrackESProducer",
    trackPairV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.03)
    ),
    charmCut = cms.double(1.5),
    trackSort = cms.string('sip2dSig'),
    minimumTrackWeight = cms.double(0.5),
    trackSelection = cms.PSet(
        totalHitsMin = cms.uint32(0),
        jetDeltaRMax = cms.double(0.3),
        qualityClass = cms.string('highPurity'),
        pixelHitsMin = cms.uint32(0),
        sip3dSigMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        maxDistToAxis = cms.double(0.07),
        sip2dValMax = cms.double(99999.9),
        maxDecayLen = cms.double(5),
        ptMin = cms.double(0.0),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        sip2dValMin = cms.double(-99999.9),
        normChi2Max = cms.double(99999.9)
    ),
    calibrationRecords = cms.vstring('GhostTrackRecoVertex', 
        'GhostTrackPseudoVertex', 
        'GhostTrackNoVertex'),
    useCategories = cms.bool(True),
    categoryVariableName = cms.string('vertexCategory')
)


process.hcalDetIdAssociator = cms.ESProducer("DetIdAssociatorESProducer",
    ComponentName = cms.string('HcalDetIdAssociator'),
    etaBinSize = cms.double(0.087),
    nEta = cms.int32(70),
    nPhi = cms.int32(72)
)


process.hcal_db_producer = cms.ESProducer("HcalDbProducer",
    file = cms.untracked.string(''),
    dump = cms.untracked.vstring('')
)


process.hoDetIdAssociator = cms.ESProducer("DetIdAssociatorESProducer",
    ComponentName = cms.string('HODetIdAssociator'),
    etaBinSize = cms.double(0.087),
    nEta = cms.int32(30),
    nPhi = cms.int32(72)
)


process.idealForDigiCSCGeometry = cms.ESProducer("CSCGeometryESModule",
    appendToDataLabel = cms.string('idealForDigi'),
    useDDD = cms.bool(True),
    debugV = cms.untracked.bool(False),
    useGangedStripsInME1a = cms.bool(True),
    alignmentsLabel = cms.string('fakeForIdeal'),
    useOnlyWiresInME1a = cms.bool(False),
    useRealWireGeometry = cms.bool(True),
    useCentreTIOffsets = cms.bool(False),
    applyAlignment = cms.bool(False)
)


process.idealForDigiDTGeometry = cms.ESProducer("DTGeometryESModule",
    appendToDataLabel = cms.string('idealForDigi'),
    fromDDD = cms.bool(True),
    applyAlignment = cms.bool(False),
    alignmentsLabel = cms.string('fakeForIdeal')
)


process.idealForDigiTrackerGeometry = cms.ESProducer("TrackerDigiGeometryESModule",
    appendToDataLabel = cms.string('idealForDigi'),
    fromDDD = cms.bool(True),
    applyAlignment = cms.bool(False),
    alignmentsLabel = cms.string('fakeForIdeal')
)


process.impactParameterMVAComputer = cms.ESProducer("GenericMVAJetTagESProducer",
    useCategories = cms.bool(False),
    calibrationRecord = cms.string('ImpactParameterMVA')
)


process.jetBProbability = cms.ESProducer("JetBProbabilityESProducer",
    deltaR = cms.double(-1.0),
    maximumDistanceToJetAxis = cms.double(0.07),
    impactParameterType = cms.int32(0),
    trackQualityClass = cms.string('any'),
    trackIpSign = cms.int32(1),
    minimumProbability = cms.double(0.005),
    numberOfBTracks = cms.uint32(4),
    maximumDecayLength = cms.double(5.0)
)


process.jetProbability = cms.ESProducer("JetProbabilityESProducer",
    deltaR = cms.double(0.3),
    maximumDistanceToJetAxis = cms.double(0.07),
    impactParameterType = cms.int32(0),
    trackQualityClass = cms.string('any'),
    trackIpSign = cms.int32(1),
    minimumProbability = cms.double(0.005),
    maximumDecayLength = cms.double(5.0)
)


process.muonDetIdAssociator = cms.ESProducer("DetIdAssociatorESProducer",
    ComponentName = cms.string('MuonDetIdAssociator'),
    includeBadChambers = cms.bool(False),
    etaBinSize = cms.double(0.125),
    nEta = cms.int32(48),
    nPhi = cms.int32(48)
)


process.preshowerDetIdAssociator = cms.ESProducer("DetIdAssociatorESProducer",
    ComponentName = cms.string('PreshowerDetIdAssociator'),
    etaBinSize = cms.double(0.1),
    nEta = cms.int32(60),
    nPhi = cms.int32(30)
)


process.siStripGainESProducer = cms.ESProducer("SiStripGainESProducer",
    printDebug = cms.untracked.bool(False),
    appendToDataLabel = cms.string(''),
    APVGain = cms.VPSet(cms.PSet(
        Record = cms.string('SiStripApvGainRcd'),
        NormalizationFactor = cms.untracked.double(1.0),
        Label = cms.untracked.string('')
    ), 
        cms.PSet(
            Record = cms.string('SiStripApvGain2Rcd'),
            NormalizationFactor = cms.untracked.double(1.0),
            Label = cms.untracked.string('')
        )),
    AutomaticNormalization = cms.bool(False)
)


process.siStripQualityESProducer = cms.ESProducer("SiStripQualityESProducer",
    appendToDataLabel = cms.string(''),
    PrintDebugOutput = cms.bool(False),
    ThresholdForReducedGranularity = cms.double(0.3),
    UseEmptyRunInfo = cms.bool(False),
    ReduceGranularity = cms.bool(False),
    ListOfRecordToMerge = cms.VPSet(cms.PSet(
        record = cms.string('SiStripDetVOffRcd'),
        tag = cms.string('')
    ), 
        cms.PSet(
            record = cms.string('SiStripDetCablingRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('RunInfoRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadChannelRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadFiberRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadModuleRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadStripRcd'),
            tag = cms.string('')
        ))
)


process.simpleSecondaryVertex2Trk = cms.ESProducer("SimpleSecondaryVertexESProducer",
    minTracks = cms.uint32(2),
    unBoost = cms.bool(False),
    useSignificance = cms.bool(True),
    use3d = cms.bool(True)
)


process.simpleSecondaryVertex3Trk = cms.ESProducer("SimpleSecondaryVertexESProducer",
    minTracks = cms.uint32(3),
    unBoost = cms.bool(False),
    useSignificance = cms.bool(True),
    use3d = cms.bool(True)
)


process.sistripconn = cms.ESProducer("SiStripConnectivity")


process.softElectron = cms.ESProducer("ElectronTaggerESProducer",
    ipSign = cms.string('any')
)


process.softLeptonByIP3d = cms.ESProducer("LeptonTaggerByIPESProducer",
    use3d = cms.bool(True),
    ipSign = cms.string('any')
)


process.softLeptonByPt = cms.ESProducer("LeptonTaggerByPtESProducer",
    ipSign = cms.string('any')
)


process.softMuon = cms.ESProducer("MuonTaggerESProducer",
    ipSign = cms.string('any')
)


process.softMuonNoIP = cms.ESProducer("MuonTaggerNoIPESProducer",
    ipSign = cms.string('any')
)


process.trackCounting3D2nd = cms.ESProducer("TrackCountingESProducer",
    deltaR = cms.double(-1.0),
    maximumDistanceToJetAxis = cms.double(0.07),
    impactParameterType = cms.int32(0),
    trackQualityClass = cms.string('any'),
    maximumDecayLength = cms.double(5.0),
    nthTrack = cms.int32(2)
)


process.trackCounting3D3rd = cms.ESProducer("TrackCountingESProducer",
    deltaR = cms.double(-1.0),
    maximumDistanceToJetAxis = cms.double(0.07),
    impactParameterType = cms.int32(0),
    trackQualityClass = cms.string('any'),
    maximumDecayLength = cms.double(5.0),
    nthTrack = cms.int32(3)
)


process.BTagRecord = cms.ESSource("EmptyESSource",
    iovIsRunNotTime = cms.bool(True),
    recordName = cms.string('JetTagComputerRecord'),
    firstValid = cms.vuint32(1)
)


process.GlobalTag = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string('.'),
        enableReadOnlySessionOnUpdateConnection = cms.untracked.bool(False),
        idleConnectionCleanupPeriod = cms.untracked.int32(10),
        messageLevel = cms.untracked.int32(0),
        enablePoolAutomaticCleanUp = cms.untracked.bool(False),
        enableConnectionSharing = cms.untracked.bool(True),
        connectionRetrialTimeOut = cms.untracked.int32(60),
        connectionTimeOut = cms.untracked.int32(60),
        connectionRetrialPeriod = cms.untracked.int32(10)
    ),
    BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
    connect = cms.string('frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'),
    globaltag = cms.string(options.globtag+'::All')
)


process.XMLIdealGeometryESSource = cms.ESSource("XMLIdealGeometryESSource",
    geomXMLFiles = cms.vstring('Geometry/CMSCommonData/data/materials.xml', 
        'Geometry/CMSCommonData/data/rotations.xml', 
        'Geometry/CMSCommonData/data/normal/cmsextent.xml', 
        'Geometry/CMSCommonData/data/cms.xml', 
        'Geometry/CMSCommonData/data/cmsMother.xml', 
        'Geometry/CMSCommonData/data/cmsTracker.xml', 
        'Geometry/CMSCommonData/data/caloBase.xml', 
        'Geometry/CMSCommonData/data/cmsCalo.xml', 
        'Geometry/CMSCommonData/data/muonBase.xml', 
        'Geometry/CMSCommonData/data/cmsMuon.xml', 
        'Geometry/CMSCommonData/data/mgnt.xml', 
        'Geometry/CMSCommonData/data/beampipe.xml', 
        'Geometry/CMSCommonData/data/cmsBeam.xml', 
        'Geometry/CMSCommonData/data/muonMB.xml', 
        'Geometry/CMSCommonData/data/muonMagnet.xml', 
        'Geometry/TrackerCommonData/data/pixfwdMaterials.xml', 
        'Geometry/TrackerCommonData/data/pixfwdCommon.xml', 
        'Geometry/TrackerCommonData/data/pixfwdPlaq.xml', 
        'Geometry/TrackerCommonData/data/pixfwdPlaq1x2.xml', 
        'Geometry/TrackerCommonData/data/pixfwdPlaq1x5.xml', 
        'Geometry/TrackerCommonData/data/pixfwdPlaq2x3.xml', 
        'Geometry/TrackerCommonData/data/pixfwdPlaq2x4.xml', 
        'Geometry/TrackerCommonData/data/pixfwdPlaq2x5.xml', 
        'Geometry/TrackerCommonData/data/pixfwdPanelBase.xml', 
        'Geometry/TrackerCommonData/data/pixfwdPanel.xml', 
        'Geometry/TrackerCommonData/data/pixfwdBlade.xml', 
        'Geometry/TrackerCommonData/data/pixfwdNipple.xml', 
        'Geometry/TrackerCommonData/data/pixfwdDisk.xml', 
        'Geometry/TrackerCommonData/data/pixfwdCylinder.xml', 
        'Geometry/TrackerCommonData/data/pixfwd.xml', 
        'Geometry/TrackerCommonData/data/pixbarmaterial.xml', 
        'Geometry/TrackerCommonData/data/pixbarladder.xml', 
        'Geometry/TrackerCommonData/data/pixbarladderfull.xml', 
        'Geometry/TrackerCommonData/data/pixbarladderhalf.xml', 
        'Geometry/TrackerCommonData/data/pixbarlayer.xml', 
        'Geometry/TrackerCommonData/data/pixbarlayer0.xml', 
        'Geometry/TrackerCommonData/data/pixbarlayer1.xml', 
        'Geometry/TrackerCommonData/data/pixbarlayer2.xml', 
        'Geometry/TrackerCommonData/data/pixbar.xml', 
        'Geometry/TrackerCommonData/data/tibtidcommonmaterial.xml', 
        'Geometry/TrackerCommonData/data/tibmaterial.xml', 
        'Geometry/TrackerCommonData/data/tibmodpar.xml', 
        'Geometry/TrackerCommonData/data/tibmodule0.xml', 
        'Geometry/TrackerCommonData/data/tibmodule0a.xml', 
        'Geometry/TrackerCommonData/data/tibmodule0b.xml', 
        'Geometry/TrackerCommonData/data/tibmodule2.xml', 
        'Geometry/TrackerCommonData/data/tibstringpar.xml', 
        'Geometry/TrackerCommonData/data/tibstring0ll.xml', 
        'Geometry/TrackerCommonData/data/tibstring0lr.xml', 
        'Geometry/TrackerCommonData/data/tibstring0ul.xml', 
        'Geometry/TrackerCommonData/data/tibstring0ur.xml', 
        'Geometry/TrackerCommonData/data/tibstring0.xml', 
        'Geometry/TrackerCommonData/data/tibstring1ll.xml', 
        'Geometry/TrackerCommonData/data/tibstring1lr.xml', 
        'Geometry/TrackerCommonData/data/tibstring1ul.xml', 
        'Geometry/TrackerCommonData/data/tibstring1ur.xml', 
        'Geometry/TrackerCommonData/data/tibstring1.xml', 
        'Geometry/TrackerCommonData/data/tibstring2ll.xml', 
        'Geometry/TrackerCommonData/data/tibstring2lr.xml', 
        'Geometry/TrackerCommonData/data/tibstring2ul.xml', 
        'Geometry/TrackerCommonData/data/tibstring2ur.xml', 
        'Geometry/TrackerCommonData/data/tibstring2.xml', 
        'Geometry/TrackerCommonData/data/tibstring3ll.xml', 
        'Geometry/TrackerCommonData/data/tibstring3lr.xml', 
        'Geometry/TrackerCommonData/data/tibstring3ul.xml', 
        'Geometry/TrackerCommonData/data/tibstring3ur.xml', 
        'Geometry/TrackerCommonData/data/tibstring3.xml', 
        'Geometry/TrackerCommonData/data/tiblayerpar.xml', 
        'Geometry/TrackerCommonData/data/tiblayer0.xml', 
        'Geometry/TrackerCommonData/data/tiblayer1.xml', 
        'Geometry/TrackerCommonData/data/tiblayer2.xml', 
        'Geometry/TrackerCommonData/data/tiblayer3.xml', 
        'Geometry/TrackerCommonData/data/tib.xml', 
        'Geometry/TrackerCommonData/data/tidmaterial.xml', 
        'Geometry/TrackerCommonData/data/tidmodpar.xml', 
        'Geometry/TrackerCommonData/data/tidmodule0.xml', 
        'Geometry/TrackerCommonData/data/tidmodule0r.xml', 
        'Geometry/TrackerCommonData/data/tidmodule0l.xml', 
        'Geometry/TrackerCommonData/data/tidmodule1.xml', 
        'Geometry/TrackerCommonData/data/tidmodule1r.xml', 
        'Geometry/TrackerCommonData/data/tidmodule1l.xml', 
        'Geometry/TrackerCommonData/data/tidmodule2.xml', 
        'Geometry/TrackerCommonData/data/tidringpar.xml', 
        'Geometry/TrackerCommonData/data/tidring0.xml', 
        'Geometry/TrackerCommonData/data/tidring0f.xml', 
        'Geometry/TrackerCommonData/data/tidring0b.xml', 
        'Geometry/TrackerCommonData/data/tidring1.xml', 
        'Geometry/TrackerCommonData/data/tidring1f.xml', 
        'Geometry/TrackerCommonData/data/tidring1b.xml', 
        'Geometry/TrackerCommonData/data/tidring2.xml', 
        'Geometry/TrackerCommonData/data/tid.xml', 
        'Geometry/TrackerCommonData/data/tidf.xml', 
        'Geometry/TrackerCommonData/data/tidb.xml', 
        'Geometry/TrackerCommonData/data/tibtidservices.xml', 
        'Geometry/TrackerCommonData/data/tibtidservicesf.xml', 
        'Geometry/TrackerCommonData/data/tibtidservicesb.xml', 
        'Geometry/TrackerCommonData/data/tobmaterial.xml', 
        'Geometry/TrackerCommonData/data/tobmodpar.xml', 
        'Geometry/TrackerCommonData/data/tobmodule0.xml', 
        'Geometry/TrackerCommonData/data/tobmodule2.xml', 
        'Geometry/TrackerCommonData/data/tobmodule4.xml', 
        'Geometry/TrackerCommonData/data/tobrodpar.xml', 
        'Geometry/TrackerCommonData/data/tobrod0c.xml', 
        'Geometry/TrackerCommonData/data/tobrod0l.xml', 
        'Geometry/TrackerCommonData/data/tobrod0h.xml', 
        'Geometry/TrackerCommonData/data/tobrod0.xml', 
        'Geometry/TrackerCommonData/data/tobrod1l.xml', 
        'Geometry/TrackerCommonData/data/tobrod1h.xml', 
        'Geometry/TrackerCommonData/data/tobrod1.xml', 
        'Geometry/TrackerCommonData/data/tobrod2c.xml', 
        'Geometry/TrackerCommonData/data/tobrod2l.xml', 
        'Geometry/TrackerCommonData/data/tobrod2h.xml', 
        'Geometry/TrackerCommonData/data/tobrod2.xml', 
        'Geometry/TrackerCommonData/data/tobrod3l.xml', 
        'Geometry/TrackerCommonData/data/tobrod3h.xml', 
        'Geometry/TrackerCommonData/data/tobrod3.xml', 
        'Geometry/TrackerCommonData/data/tobrod4c.xml', 
        'Geometry/TrackerCommonData/data/tobrod4l.xml', 
        'Geometry/TrackerCommonData/data/tobrod4h.xml', 
        'Geometry/TrackerCommonData/data/tobrod4.xml', 
        'Geometry/TrackerCommonData/data/tobrod5l.xml', 
        'Geometry/TrackerCommonData/data/tobrod5h.xml', 
        'Geometry/TrackerCommonData/data/tobrod5.xml', 
        'Geometry/TrackerCommonData/data/tob.xml', 
        'Geometry/TrackerCommonData/data/tecmaterial.xml', 
        'Geometry/TrackerCommonData/data/tecmodpar.xml', 
        'Geometry/TrackerCommonData/data/tecmodule0.xml', 
        'Geometry/TrackerCommonData/data/tecmodule0r.xml', 
        'Geometry/TrackerCommonData/data/tecmodule0s.xml', 
        'Geometry/TrackerCommonData/data/tecmodule1.xml', 
        'Geometry/TrackerCommonData/data/tecmodule1r.xml', 
        'Geometry/TrackerCommonData/data/tecmodule1s.xml', 
        'Geometry/TrackerCommonData/data/tecmodule2.xml', 
        'Geometry/TrackerCommonData/data/tecmodule3.xml', 
        'Geometry/TrackerCommonData/data/tecmodule4.xml', 
        'Geometry/TrackerCommonData/data/tecmodule4r.xml', 
        'Geometry/TrackerCommonData/data/tecmodule4s.xml', 
        'Geometry/TrackerCommonData/data/tecmodule5.xml', 
        'Geometry/TrackerCommonData/data/tecmodule6.xml', 
        'Geometry/TrackerCommonData/data/tecpetpar.xml', 
        'Geometry/TrackerCommonData/data/tecring0.xml', 
        'Geometry/TrackerCommonData/data/tecring1.xml', 
        'Geometry/TrackerCommonData/data/tecring2.xml', 
        'Geometry/TrackerCommonData/data/tecring3.xml', 
        'Geometry/TrackerCommonData/data/tecring4.xml', 
        'Geometry/TrackerCommonData/data/tecring5.xml', 
        'Geometry/TrackerCommonData/data/tecring6.xml', 
        'Geometry/TrackerCommonData/data/tecring0f.xml', 
        'Geometry/TrackerCommonData/data/tecring1f.xml', 
        'Geometry/TrackerCommonData/data/tecring2f.xml', 
        'Geometry/TrackerCommonData/data/tecring3f.xml', 
        'Geometry/TrackerCommonData/data/tecring4f.xml', 
        'Geometry/TrackerCommonData/data/tecring5f.xml', 
        'Geometry/TrackerCommonData/data/tecring6f.xml', 
        'Geometry/TrackerCommonData/data/tecring0b.xml', 
        'Geometry/TrackerCommonData/data/tecring1b.xml', 
        'Geometry/TrackerCommonData/data/tecring2b.xml', 
        'Geometry/TrackerCommonData/data/tecring3b.xml', 
        'Geometry/TrackerCommonData/data/tecring4b.xml', 
        'Geometry/TrackerCommonData/data/tecring5b.xml', 
        'Geometry/TrackerCommonData/data/tecring6b.xml', 
        'Geometry/TrackerCommonData/data/tecpetalf.xml', 
        'Geometry/TrackerCommonData/data/tecpetalb.xml', 
        'Geometry/TrackerCommonData/data/tecpetal0.xml', 
        'Geometry/TrackerCommonData/data/tecpetal0f.xml', 
        'Geometry/TrackerCommonData/data/tecpetal0b.xml', 
        'Geometry/TrackerCommonData/data/tecpetal3.xml', 
        'Geometry/TrackerCommonData/data/tecpetal3f.xml', 
        'Geometry/TrackerCommonData/data/tecpetal3b.xml', 
        'Geometry/TrackerCommonData/data/tecpetal6f.xml', 
        'Geometry/TrackerCommonData/data/tecpetal6b.xml', 
        'Geometry/TrackerCommonData/data/tecpetal8f.xml', 
        'Geometry/TrackerCommonData/data/tecpetal8b.xml', 
        'Geometry/TrackerCommonData/data/tecwheel.xml', 
        'Geometry/TrackerCommonData/data/tecwheela.xml', 
        'Geometry/TrackerCommonData/data/tecwheelb.xml', 
        'Geometry/TrackerCommonData/data/tecwheelc.xml', 
        'Geometry/TrackerCommonData/data/tecwheeld.xml', 
        'Geometry/TrackerCommonData/data/tecwheel6.xml', 
        'Geometry/TrackerCommonData/data/tecservices.xml', 
        'Geometry/TrackerCommonData/data/tecbackplate.xml', 
        'Geometry/TrackerCommonData/data/tec.xml', 
        'Geometry/TrackerCommonData/data/trackermaterial.xml', 
        'Geometry/TrackerCommonData/data/tracker.xml', 
        'Geometry/TrackerCommonData/data/trackerpixbar.xml', 
        'Geometry/TrackerCommonData/data/trackerpixfwd.xml', 
        'Geometry/TrackerCommonData/data/trackertibtidservices.xml', 
        'Geometry/TrackerCommonData/data/trackertib.xml', 
        'Geometry/TrackerCommonData/data/trackertid.xml', 
        'Geometry/TrackerCommonData/data/trackertob.xml', 
        'Geometry/TrackerCommonData/data/trackertec.xml', 
        'Geometry/TrackerCommonData/data/trackerbulkhead.xml', 
        'Geometry/TrackerCommonData/data/trackerother.xml', 
        'Geometry/EcalCommonData/data/eregalgo.xml', 
        'Geometry/EcalCommonData/data/ebalgo.xml', 
        'Geometry/EcalCommonData/data/ebcon.xml', 
        'Geometry/EcalCommonData/data/ebrot.xml', 
        'Geometry/EcalCommonData/data/eecon.xml', 
        'Geometry/EcalCommonData/data/eefixed.xml', 
        'Geometry/EcalCommonData/data/eehier.xml', 
        'Geometry/EcalCommonData/data/eealgo.xml', 
        'Geometry/EcalCommonData/data/escon.xml', 
        'Geometry/EcalCommonData/data/esalgo.xml', 
        'Geometry/EcalCommonData/data/eeF.xml', 
        'Geometry/EcalCommonData/data/eeB.xml', 
        'Geometry/HcalCommonData/data/hcalrotations.xml', 
        'Geometry/HcalCommonData/data/hcalalgo.xml', 
        'Geometry/HcalCommonData/data/hcalbarrelalgo.xml', 
        'Geometry/HcalCommonData/data/hcalendcapalgo.xml', 
        'Geometry/HcalCommonData/data/hcalouteralgo.xml', 
        'Geometry/HcalCommonData/data/hcalforwardalgo.xml', 
        'Geometry/HcalCommonData/data/hcalforwardfibre.xml', 
        'Geometry/HcalCommonData/data/hcalforwardmaterial.xml', 
        'Geometry/MuonCommonData/data/mbCommon.xml', 
        'Geometry/MuonCommonData/data/mb1.xml', 
        'Geometry/MuonCommonData/data/mb2.xml', 
        'Geometry/MuonCommonData/data/mb3.xml', 
        'Geometry/MuonCommonData/data/mb4.xml', 
        'Geometry/MuonCommonData/data/muonYoke.xml', 
        'Geometry/MuonCommonData/data/mf.xml', 
        'Geometry/ForwardCommonData/data/forward.xml', 
        'Geometry/ForwardCommonData/data/forwardshield.xml', 
        'Geometry/ForwardCommonData/data/brmrotations.xml', 
        'Geometry/ForwardCommonData/data/brm.xml', 
        'Geometry/ForwardCommonData/data/totemMaterials.xml', 
        'Geometry/ForwardCommonData/data/totemRotations.xml', 
        'Geometry/ForwardCommonData/data/totemt2.xml', 
        'Geometry/ForwardCommonData/data/ionpump.xml', 
        'Geometry/MuonCommonData/data/muonNumbering.xml', 
        'Geometry/TrackerCommonData/data/trackerStructureTopology.xml', 
        'Geometry/TrackerSimData/data/trackersens.xml', 
        'Geometry/TrackerRecoData/data/trackerRecoMaterial.xml', 
        'Geometry/EcalSimData/data/ecalsens.xml', 
        'Geometry/HcalCommonData/data/hcalsens.xml', 
        'Geometry/HcalSimData/data/CaloUtil.xml', 
        'Geometry/HcalSimData/data/hf.xml', 
        'Geometry/HcalSimData/data/hffibre.xml', 
        'Geometry/MuonSimData/data/muonSens.xml', 
        'Geometry/DTGeometryBuilder/data/dtSpecsFilter.xml', 
        'Geometry/CSCGeometryBuilder/data/cscSpecsFilter.xml', 
        'Geometry/CSCGeometryBuilder/data/cscSpecs.xml', 
        'Geometry/RPCGeometryBuilder/data/RPCSpecs.xml', 
        'Geometry/ForwardCommonData/data/brmsens.xml', 
        'Geometry/HcalSimData/data/HcalProdCuts.xml', 
        'Geometry/EcalSimData/data/EcalProdCuts.xml', 
        'Geometry/EcalSimData/data/ESProdCuts.xml', 
        'Geometry/TrackerSimData/data/trackerProdCuts.xml', 
        'Geometry/TrackerSimData/data/trackerProdCutsBEAM.xml', 
        'Geometry/MuonSimData/data/muonProdCuts.xml', 
        'Geometry/ForwardSimData/data/ForwardShieldProdCuts.xml', 
        'Geometry/CMSCommonData/data/FieldParameters.xml'),
    rootNodeName = cms.string('cms:OCMS')
)


process.ak5CaloL1FastL2L3 = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ak5CaloL1Fastjet', 
        'ak5CaloL2Relative', 
        'ak5CaloL3Absolute')
)


process.ak5CaloL1FastL2L3L6 = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ak5CaloL1Fastjet', 
        'ak5CaloL2Relative', 
        'ak5CaloL3Absolute', 
        'ak5CaloL6SLB')
)


process.ak5CaloL1FastL2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ak5CaloL1Fastjet', 
        'ak5CaloL2Relative', 
        'ak5CaloL3Absolute', 
        'ak5CaloResidual')
)


process.ak5CaloL1Fastjet = cms.ESSource("L1FastjetCorrectionService",
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1FastJet'),
    section = cms.string(''),
    srcRho = cms.InputTag("kt6PFJets","rho"),
    era = cms.string('Jec10V1')
)


process.ak5CaloL1L2L3 = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ak5CaloL1Offset', 
        'ak5CaloL2Relative', 
        'ak5CaloL3Absolute')
)


process.ak5CaloL1L2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ak5CaloL1Offset', 
        'ak5CaloL2Relative', 
        'ak5CaloL3Absolute', 
        'ak5CaloResidual')
)


process.ak5CaloL1Offset = cms.ESSource("L1OffsetCorrectionService",
    minVtxNdof = cms.int32(4),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1Offset'),
    section = cms.string(''),
    vertexCollection = cms.string('offlinePrimaryVertices'),
    era = cms.string('Fall10')
)


process.ak5CaloL2L3 = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ak5CaloL2Relative', 
        'ak5CaloL3Absolute')
)


process.ak5CaloL2L3L6 = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ak5CaloL2Relative', 
        'ak5CaloL3Absolute', 
        'ak5CaloL6SLB')
)


process.ak5CaloL2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ak5CaloL2Relative', 
        'ak5CaloL3Absolute', 
        'ak5CaloResidual')
)


process.ak5CaloL2Relative = cms.ESSource("LXXXCorrectionService",
    section = cms.string(''),
    level = cms.string('L2Relative'),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5Calo'),
    era = cms.string('Spring10')
)


process.ak5CaloL3Absolute = cms.ESSource("LXXXCorrectionService",
    section = cms.string(''),
    level = cms.string('L3Absolute'),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5Calo'),
    era = cms.string('Spring10')
)


process.ak5CaloL6SLB = cms.ESSource("L6SLBCorrectionService",
    srcBTagInfoElectron = cms.InputTag("ak5CaloJetsSoftElectronTagInfos"),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    section = cms.string(''),
    srcBTagInfoMuon = cms.InputTag("ak5CaloJetsSoftMuonTagInfos"),
    addMuonToJet = cms.bool(True),
    era = cms.string('')
)


process.ak5CaloResidual = cms.ESSource("LXXXCorrectionService",
    section = cms.string(''),
    level = cms.string('L2L3Residual'),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5Calo'),
    era = cms.string('Spring10DataV2')
)


process.ak5JPTL1FastL2L3 = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ak5JPTL1Fastjet', 
        'ak5JPTL1Offset', 
        'ak5JPTL2Relative', 
        'ak5JPTL3Absolute')
)


process.ak5JPTL1FastL2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ak5JPTL1Fastjet', 
        'ak5L1JPTOffset', 
        'ak5JPTL2Relative', 
        'ak5JPTL3Absolute', 
        'ak5JPTResidual')
)


process.ak5JPTL1Fastjet = cms.ESSource("L1FastjetCorrectionService",
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1FastJet'),
    section = cms.string(''),
    srcRho = cms.InputTag("kt6PFJets","rho"),
    era = cms.string('Jec10V1')
)


process.ak5JPTL1L2L3 = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ak5JPTL1Offset', 
        'ak5L1JPTOffset', 
        'ak5JPTL2Relative', 
        'ak5JPTL3Absolute')
)


process.ak5JPTL1L2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ak5JPTL1Offset', 
        'ak5L1JPTOffset', 
        'ak5JPTL2Relative', 
        'ak5JPTL3Absolute', 
        'ak5JPTResidual')
)


process.ak5JPTL1Offset = cms.ESSource("L1OffsetCorrectionService",
    minVtxNdof = cms.int32(4),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1Offset'),
    section = cms.string(''),
    vertexCollection = cms.string('offlinePrimaryVertices'),
    era = cms.string('Fall10')
)


process.ak5JPTL2L3 = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ak5L1JPTOffset', 
        'ak5JPTL2Relative', 
        'ak5JPTL3Absolute')
)


process.ak5JPTL2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ak5L1JPTOffset', 
        'ak5JPTL2Relative', 
        'ak5JPTL3Absolute', 
        'ak5JPTResidual')
)


process.ak5JPTL2Relative = cms.ESSource("LXXXCorrectionService",
    useCondDB = cms.untracked.bool(True),
    era = cms.string('Summer10'),
    section = cms.string(''),
    algorithm = cms.string('AK5JPT'),
    level = cms.string('L2Relative')
)


process.ak5JPTL3Absolute = cms.ESSource("LXXXCorrectionService",
    useCondDB = cms.untracked.bool(True),
    era = cms.string('Summer10'),
    section = cms.string(''),
    algorithm = cms.string('AK5JPT'),
    level = cms.string('L3Absolute')
)


process.ak5JPTResidual = cms.ESSource("LXXXCorrectionService",
    useCondDB = cms.untracked.bool(True),
    era = cms.string('Spring10DataV2'),
    section = cms.string(''),
    algorithm = cms.string('AK5JPT'),
    level = cms.string('L2L3Residual')
)


process.ak5L1JPTOffset = cms.ESSource("LXXXCorrectionService",
    section = cms.string(''),
    era = cms.string('Summer10'),
    algorithm = cms.string('AK5JPT'),
    level = cms.string('L1JPTOffset')
)


process.ak5PFL1FastL2L3 = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ak5PFL1Fastjet', 
        'ak5PFL2Relative', 
        'ak5PFL3Absolute')
)


process.ak5PFL1FastL2L3L6 = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ak5PFL1Fastjet', 
        'ak5PFL2Relative', 
        'ak5PFL3Absolute', 
        'ak5PFL6SLB')
)


process.ak5PFL1FastL2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ak5PFL1Fastjet', 
        'ak5PFL2Relative', 
        'ak5PFL3Absolute', 
        'ak5PFResidual')
)


process.ak5PFL1Fastjet = cms.ESSource("L1FastjetCorrectionService",
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5PF'),
    level = cms.string('L1FastJet'),
    section = cms.string(''),
    srcRho = cms.InputTag("kt6PFJets","rho"),
    era = cms.string('Jec10V1')
)


process.ak5PFL1L2L3 = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ak5PFL1Offset', 
        'ak5PFL2Relative', 
        'ak5PFL3Absolute')
)


process.ak5PFL1L2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ak5PFL1Offset', 
        'ak5PFL2Relative', 
        'ak5PFL3Absolute', 
        'ak5PFResidual')
)


process.ak5PFL1Offset = cms.ESSource("L1OffsetCorrectionService",
    minVtxNdof = cms.int32(4),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5PF'),
    level = cms.string('L1Offset'),
    section = cms.string(''),
    vertexCollection = cms.string('offlinePrimaryVertices'),
    era = cms.string('Fall10')
)


process.ak5PFL2L3 = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ak5PFL2Relative', 
        'ak5PFL3Absolute')
)


process.ak5PFL2L3L6 = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ak5PFL2Relative', 
        'ak5PFL3Absolute', 
        'ak5PFL6SLB')
)


process.ak5PFL2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ak5PFL2Relative', 
        'ak5PFL3Absolute', 
        'ak5PFResidual')
)


process.ak5PFL2Relative = cms.ESSource("LXXXCorrectionService",
    useCondDB = cms.untracked.bool(True),
    era = cms.string('Spring10'),
    section = cms.string(''),
    algorithm = cms.string('AK5PF'),
    level = cms.string('L2Relative')
)


process.ak5PFL3Absolute = cms.ESSource("LXXXCorrectionService",
    useCondDB = cms.untracked.bool(True),
    era = cms.string('Spring10'),
    section = cms.string(''),
    algorithm = cms.string('AK5PF'),
    level = cms.string('L3Absolute')
)


process.ak5PFL6SLB = cms.ESSource("L6SLBCorrectionService",
    srcBTagInfoElectron = cms.InputTag("ak5PFJetsSoftElectronTagInfos"),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    section = cms.string(''),
    srcBTagInfoMuon = cms.InputTag("ak5PFJetsSoftMuonTagInfos"),
    addMuonToJet = cms.bool(False),
    era = cms.string('')
)


process.ak5PFResidual = cms.ESSource("LXXXCorrectionService",
    useCondDB = cms.untracked.bool(True),
    era = cms.string('Spring10DataV2'),
    section = cms.string(''),
    algorithm = cms.string('AK5PF'),
    level = cms.string('L2L3Residual')
)


process.ak5TrackL1FastL2L3 = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ak5CaloL1Fastjet', 
        'ak5TrackL2Relative', 
        'ak5TrackL3Absolute')
)


process.ak5TrackL2L3 = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ak5TrackL2Relative', 
        'ak5TrackL3Absolute')
)


process.ak5TrackL2Relative = cms.ESSource("LXXXCorrectionService",
    useCondDB = cms.untracked.bool(True),
    era = cms.string('Spring10'),
    section = cms.string(''),
    algorithm = cms.string('AK5TRK'),
    level = cms.string('L2Relative')
)


process.ak5TrackL3Absolute = cms.ESSource("LXXXCorrectionService",
    useCondDB = cms.untracked.bool(True),
    era = cms.string('Spring10'),
    section = cms.string(''),
    algorithm = cms.string('AK5TRK'),
    level = cms.string('L3Absolute')
)


process.ak7CaloL1FastL2L3 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak5CaloL1Fastjet', 
        'ak7CaloL2Relative', 
        'ak7CaloL3Absolute')
)


process.ak7CaloL1FastL2L3L6 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak7CaloL1Offset', 
        'ak7CaloL2Relative', 
        'ak7CaloL3Absolute', 
        'ak7CaloL6SLB')
)


process.ak7CaloL1FastL2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak7CaloL1Fastjet', 
        'ak7CaloL2Relative', 
        'ak7CaloL3Absolute', 
        'ak7CaloResidual')
)


process.ak7CaloL1Fastjet = cms.ESSource("L1FastjetCorrectionService",
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1FastJet'),
    section = cms.string(''),
    srcRho = cms.InputTag("kt6PFJets","rho"),
    era = cms.string('Jec10V1')
)


process.ak7CaloL1L2L3 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak7CaloL1Offset', 
        'ak7CaloL2Relative', 
        'ak7CaloL3Absolute')
)


process.ak7CaloL1L2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak7CaloL1Offset', 
        'ak7CaloL2Relative', 
        'ak7CaloL3Absolute', 
        'ak7CaloResidual')
)


process.ak7CaloL1Offset = cms.ESSource("L1OffsetCorrectionService",
    minVtxNdof = cms.int32(4),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1Offset'),
    section = cms.string(''),
    vertexCollection = cms.string('offlinePrimaryVertices'),
    era = cms.string('Fall10')
)


process.ak7CaloL2L3 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak7CaloL2Relative', 
        'ak7CaloL3Absolute')
)


process.ak7CaloL2L3L6 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak7CaloL2Relative', 
        'ak7CaloL3Absolute', 
        'ak7CaloL6SLB')
)


process.ak7CaloL2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ak7CaloL2Relative', 
        'ak7CaloL3Absolute', 
        'ak7CaloResidual')
)


process.ak7CaloL2Relative = cms.ESSource("LXXXCorrectionService",
    useCondDB = cms.untracked.bool(True),
    era = cms.string('Spring10'),
    section = cms.string(''),
    algorithm = cms.string('AK7Calo'),
    level = cms.string('L2Relative')
)


process.ak7CaloL3Absolute = cms.ESSource("LXXXCorrectionService",
    useCondDB = cms.untracked.bool(True),
    era = cms.string('Spring10'),
    section = cms.string(''),
    algorithm = cms.string('AK7Calo'),
    level = cms.string('L3Absolute')
)


process.ak7CaloL6SLB = cms.ESSource("L6SLBCorrectionService",
    srcBTagInfoElectron = cms.InputTag("ak7CaloJetsSoftElectronTagInfos"),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    section = cms.string(''),
    srcBTagInfoMuon = cms.InputTag("ak7CaloJetsSoftMuonTagInfos"),
    addMuonToJet = cms.bool(True),
    era = cms.string('')
)


process.ak7CaloResidual = cms.ESSource("LXXXCorrectionService",
    useCondDB = cms.untracked.bool(True),
    era = cms.string('Spring10DataV2'),
    section = cms.string(''),
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L2L3Residual')
)


process.ak7JPTL1FastL2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak7JPTL1Fastjet', 
        'ak7L1JPTOffset', 
        'ak7JPTL2Relative', 
        'ak7JPTL3Absolute', 
        'ak7JPTResidual')
)


process.ak7JPTL1Fastjet = cms.ESSource("L1FastjetCorrectionService",
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1FastJet'),
    section = cms.string(''),
    srcRho = cms.InputTag("kt6PFJets","rho"),
    era = cms.string('Jec10V1')
)


process.ak7JPTL1L2L3 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak7JPTL1Offset', 
        'ak7L1JPTOffset', 
        'ak7JPTL2Relative', 
        'ak7JPTL3Absolute')
)


process.ak7JPTL1L2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak7JPTL1Offset', 
        'ak7L1JPTOffset', 
        'ak7JPTL2Relative', 
        'ak7JPTL3Absolute', 
        'ak7JPTResidual')
)


process.ak7JPTL1Offset = cms.ESSource("L1OffsetCorrectionService",
    minVtxNdof = cms.int32(4),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1Offset'),
    section = cms.string(''),
    vertexCollection = cms.string('offlinePrimaryVertices'),
    era = cms.string('Fall10')
)


process.ak7JPTL2L3 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak7L1JPTOffset', 
        'ak7JPTL2Relative', 
        'ak7JPTL3Absolute')
)


process.ak7L1JPTOffset = cms.ESSource("LXXXCorrectionService",
    section = cms.string(''),
    era = cms.string('Summer10'),
    algorithm = cms.string('AK5JPT'),
    level = cms.string('L1JPTOffset')
)


process.ak7PFL1FastL2L3 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak5PFL1Fastjet', 
        'ak7PFL2Relative', 
        'ak7PFL3Absolute')
)


process.ak7PFL1FastL2L3L6 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak5PFL1Fastjet', 
        'ak7PFL2Relative', 
        'ak7PFL3Absolute', 
        'ak7PFL6SLB')
)


process.ak7PFL1FastL2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak7PFL1Fastjet', 
        'ak7PFL2Relative', 
        'ak7PFL3Absolute', 
        'ak7PFResidual')
)


process.ak7PFL1Fastjet = cms.ESSource("L1FastjetCorrectionService",
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5PF'),
    level = cms.string('L1FastJet'),
    section = cms.string(''),
    srcRho = cms.InputTag("kt6PFJets","rho"),
    era = cms.string('Jec10V1')
)


process.ak7PFL1L2L3 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak7PFL1Offset', 
        'ak7PFL2Relative', 
        'ak7PFL3Absolute')
)


process.ak7PFL1L2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak7PFL1Offset', 
        'ak7PFL2Relative', 
        'ak7PFL3Absolute', 
        'ak7PFResidual')
)


process.ak7PFL1Offset = cms.ESSource("L1OffsetCorrectionService",
    minVtxNdof = cms.int32(4),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5PF'),
    level = cms.string('L1Offset'),
    section = cms.string(''),
    vertexCollection = cms.string('offlinePrimaryVertices'),
    era = cms.string('Fall10')
)


process.ak7PFL2L3 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak7PFL2Relative', 
        'ak7PFL3Absolute')
)


process.ak7PFL2L3L6 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak7PFL2Relative', 
        'ak7PFL3Absolute', 
        'ak7PFL6SLB')
)


process.ak7PFL2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ak7PFL2Relative', 
        'ak7PFL3Absolute', 
        'ak7PFResidual')
)


process.ak7PFL2Relative = cms.ESSource("LXXXCorrectionService",
    section = cms.string(''),
    era = cms.string('Spring10'),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK7PF'),
    level = cms.string('L2Relative')
)


process.ak7PFL3Absolute = cms.ESSource("LXXXCorrectionService",
    section = cms.string(''),
    era = cms.string('Spring10'),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK7PF'),
    level = cms.string('L3Absolute')
)


process.ak7PFL6SLB = cms.ESSource("L6SLBCorrectionService",
    srcBTagInfoElectron = cms.InputTag("ak7PFJetsSoftElectronTagInfos"),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    section = cms.string(''),
    srcBTagInfoMuon = cms.InputTag("ak7PFJetsSoftMuonTagInfos"),
    addMuonToJet = cms.bool(False),
    era = cms.string('')
)


process.ak7PFResidual = cms.ESSource("LXXXCorrectionService",
    section = cms.string(''),
    era = cms.string('Spring10DataV2'),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5PF'),
    level = cms.string('L2L3Residual')
)


process.eegeom = cms.ESSource("EmptyESSource",
    iovIsRunNotTime = cms.bool(True),
    recordName = cms.string('EcalMappingRcd'),
    firstValid = cms.vuint32(1)
)


process.es_hardcode = cms.ESSource("HcalHardcodeCalibrations",
    toGet = cms.untracked.vstring('GainWidths')
)


process.ic5CaloL1FastL2L3 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak5CaloL1Fastjet', 
        'ic5CaloL2Relative', 
        'ic5CaloL3Absolute')
)


process.ic5CaloL1FastL2L3L6 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ic5CaloL1Offset', 
        'ic5CaloL2Relative', 
        'ic5CaloL3Absolute', 
        'ic5CaloL6SLB')
)


process.ic5CaloL1FastL2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ic5CaloL1Fastjet', 
        'ic5CaloL2Relative', 
        'ic5CaloL3Absolute', 
        'ic5CaloResidual')
)


process.ic5CaloL1Fastjet = cms.ESSource("L1FastjetCorrectionService",
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1FastJet'),
    section = cms.string(''),
    srcRho = cms.InputTag("kt6PFJets","rho"),
    era = cms.string('Jec10V1')
)


process.ic5CaloL1L2L3 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ic5CaloL1Offset', 
        'ic5CaloL2Relative', 
        'ic5CaloL3Absolute')
)


process.ic5CaloL1L2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ic5CaloL1Offset', 
        'ic5CaloL2Relative', 
        'ic5CaloL3Absolute', 
        'ic5CaloResidual')
)


process.ic5CaloL1Offset = cms.ESSource("L1OffsetCorrectionService",
    minVtxNdof = cms.int32(4),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1Offset'),
    section = cms.string(''),
    vertexCollection = cms.string('offlinePrimaryVertices'),
    era = cms.string('Fall10')
)


process.ic5CaloL2L3 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ic5CaloL2Relative', 
        'ic5CaloL3Absolute')
)


process.ic5CaloL2L3L6 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ic5CaloL2Relative', 
        'ic5CaloL3Absolute', 
        'ic5CaloL6SLB')
)


process.ic5CaloL2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ic5CaloL2Relative', 
        'ic5CaloL3Absolute', 
        'ic5CaloResidual')
)


process.ic5CaloL2Relative = cms.ESSource("LXXXCorrectionService",
    useCondDB = cms.untracked.bool(True),
    era = cms.string('Spring10'),
    section = cms.string(''),
    algorithm = cms.string('IC5Calo'),
    level = cms.string('L2Relative')
)


process.ic5CaloL3Absolute = cms.ESSource("LXXXCorrectionService",
    useCondDB = cms.untracked.bool(True),
    era = cms.string('Spring10'),
    section = cms.string(''),
    algorithm = cms.string('IC5Calo'),
    level = cms.string('L3Absolute')
)


process.ic5CaloL6SLB = cms.ESSource("L6SLBCorrectionService",
    srcBTagInfoElectron = cms.InputTag("ic5CaloJetsSoftElectronTagInfos"),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    section = cms.string(''),
    srcBTagInfoMuon = cms.InputTag("ic5CaloJetsSoftMuonTagInfos"),
    addMuonToJet = cms.bool(True),
    era = cms.string('')
)


process.ic5CaloResidual = cms.ESSource("LXXXCorrectionService",
    useCondDB = cms.untracked.bool(True),
    era = cms.string('Spring10DataV2'),
    section = cms.string(''),
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L2L3Residual')
)


process.ic5PFL1FastL2L3 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak5PFL1Fastjet', 
        'ic5PFL2Relative', 
        'ic5PFL3Absolute')
)


process.ic5PFL1FastL2L3L6 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak5PFL1Fastjet', 
        'ic5PFL2Relative', 
        'ic5PFL3Absolute', 
        'ic5PFL6SLB')
)


process.ic5PFL1FastL2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ic5PFL1Fastjet', 
        'ic5PFL2Relative', 
        'ic5PFL3Absolute', 
        'ic5PFResidual')
)


process.ic5PFL1Fastjet = cms.ESSource("L1FastjetCorrectionService",
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5PF'),
    level = cms.string('L1FastJet'),
    section = cms.string(''),
    srcRho = cms.InputTag("kt6PFJets","rho"),
    era = cms.string('Jec10V1')
)


process.ic5PFL1L2L3 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ic5PFL1Offset', 
        'ic5PFL2Relative', 
        'ic5PFL3Absolute')
)


process.ic5PFL1L2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ic5PFL1Offset', 
        'ic5PFL2Relative', 
        'ic5PFL3Absolute', 
        'ic5PFResidual')
)


process.ic5PFL1Offset = cms.ESSource("L1OffsetCorrectionService",
    minVtxNdof = cms.int32(4),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5PF'),
    level = cms.string('L1Offset'),
    section = cms.string(''),
    vertexCollection = cms.string('offlinePrimaryVertices'),
    era = cms.string('Fall10')
)


process.ic5PFL2L3 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ic5PFL2Relative', 
        'ic5PFL3Absolute')
)


process.ic5PFL2L3L6 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ic5PFL2Relative', 
        'ic5PFL3Absolute', 
        'ic5PFL6SLB')
)


process.ic5PFL2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('ic5PFL2Relative', 
        'ic5PFL3Absolute', 
        'ic5PFResidual')
)


process.ic5PFL2Relative = cms.ESSource("LXXXCorrectionService",
    section = cms.string(''),
    era = cms.string('Spring10'),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('IC5PF'),
    level = cms.string('L2Relative')
)


process.ic5PFL3Absolute = cms.ESSource("LXXXCorrectionService",
    section = cms.string(''),
    era = cms.string('Spring10'),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('IC5PF'),
    level = cms.string('L3Absolute')
)


process.ic5PFL6SLB = cms.ESSource("L6SLBCorrectionService",
    srcBTagInfoElectron = cms.InputTag("ic5PFJetsSoftElectronTagInfos"),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    section = cms.string(''),
    srcBTagInfoMuon = cms.InputTag("ic5PFJetsSoftMuonTagInfos"),
    addMuonToJet = cms.bool(False),
    era = cms.string('')
)


process.ic5PFResidual = cms.ESSource("LXXXCorrectionService",
    section = cms.string(''),
    era = cms.string('Spring10DataV2'),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5PF'),
    level = cms.string('L2L3Residual')
)


process.kt4CaloL1FastL2L3 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak5CaloL1Fastjet', 
        'kt4CaloL2Relative', 
        'kt4CaloL3Absolute')
)


process.kt4CaloL1FastL2L3L6 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('kt4CaloL1Offset', 
        'kt4CaloL2Relative', 
        'kt4CaloL3Absolute', 
        'kt4CaloL6SLB')
)


process.kt4CaloL1FastL2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('kt4CaloL1Fastjet', 
        'kt4CaloL2Relative', 
        'kt4CaloL3Absolute', 
        'kt4CaloResidual')
)


process.kt4CaloL1Fastjet = cms.ESSource("L1FastjetCorrectionService",
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1FastJet'),
    section = cms.string(''),
    srcRho = cms.InputTag("kt6PFJets","rho"),
    era = cms.string('Jec10V1')
)


process.kt4CaloL1L2L3 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('kt4CaloL1Offset', 
        'kt4CaloL2Relative', 
        'kt4CaloL3Absolute')
)


process.kt4CaloL1L2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('kt4CaloL1Offset', 
        'kt4CaloL2Relative', 
        'kt4CaloL3Absolute', 
        'kt4CaloResidual')
)


process.kt4CaloL1Offset = cms.ESSource("L1OffsetCorrectionService",
    minVtxNdof = cms.int32(4),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1Offset'),
    section = cms.string(''),
    vertexCollection = cms.string('offlinePrimaryVertices'),
    era = cms.string('Fall10')
)


process.kt4CaloL2L3 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('kt4CaloL2Relative', 
        'kt4CaloL3Absolute')
)


process.kt4CaloL2L3L6 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('kt4CaloL2Relative', 
        'kt4CaloL3Absolute', 
        'kt4CaloL6SLB')
)


process.kt4CaloL2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('kt4CaloL2Relative', 
        'kt4CaloL3Absolute', 
        'kt4CaloResidual')
)


process.kt4CaloL2Relative = cms.ESSource("LXXXCorrectionService",
    useCondDB = cms.untracked.bool(True),
    era = cms.string('Spring10'),
    section = cms.string(''),
    algorithm = cms.string('KT4Calo'),
    level = cms.string('L2Relative')
)


process.kt4CaloL3Absolute = cms.ESSource("LXXXCorrectionService",
    useCondDB = cms.untracked.bool(True),
    era = cms.string('Spring10'),
    section = cms.string(''),
    algorithm = cms.string('KT4Calo'),
    level = cms.string('L3Absolute')
)


process.kt4CaloL6SLB = cms.ESSource("L6SLBCorrectionService",
    srcBTagInfoElectron = cms.InputTag("kt4CaloJetsSoftElectronTagInfos"),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    section = cms.string(''),
    srcBTagInfoMuon = cms.InputTag("kt4CaloJetsSoftMuonTagInfos"),
    addMuonToJet = cms.bool(True),
    era = cms.string('')
)


process.kt4CaloResidual = cms.ESSource("LXXXCorrectionService",
    useCondDB = cms.untracked.bool(True),
    era = cms.string('Spring10DataV2'),
    section = cms.string(''),
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L2L3Residual')
)


process.kt4PFL1FastL2L3 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak5PFL1Fastjet', 
        'kt4PFL2Relative', 
        'kt4PFL3Absolute')
)


process.kt4PFL1FastL2L3L6 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak5PFL1Fastjet', 
        'kt4PFL2Relative', 
        'kt4PFL3Absolute', 
        'kt4PFL6SLB')
)


process.kt4PFL1FastL2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('kt4PFL1Fastjet', 
        'kt4PFL2Relative', 
        'kt4PFL3Absolute', 
        'kt4PFResidual')
)


process.kt4PFL1Fastjet = cms.ESSource("L1FastjetCorrectionService",
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5PF'),
    level = cms.string('L1FastJet'),
    section = cms.string(''),
    srcRho = cms.InputTag("kt6PFJets","rho"),
    era = cms.string('Jec10V1')
)


process.kt4PFL1L2L3 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('kt4PFL1Offset', 
        'kt4PFL2Relative', 
        'kt4PFL3Absolute')
)


process.kt4PFL1L2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('kt4PFL1Offset', 
        'kt4PFL2Relative', 
        'kt4PFL3Absolute', 
        'kt4PFResidual')
)


process.kt4PFL1Offset = cms.ESSource("L1OffsetCorrectionService",
    minVtxNdof = cms.int32(4),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5PF'),
    level = cms.string('L1Offset'),
    section = cms.string(''),
    vertexCollection = cms.string('offlinePrimaryVertices'),
    era = cms.string('Fall10')
)


process.kt4PFL2L3 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('kt4PFL2Relative', 
        'kt4PFL3Absolute')
)


process.kt4PFL2L3L6 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('kt4PFL2Relative', 
        'kt4PFL3Absolute', 
        'kt4PFL6SLB')
)


process.kt4PFL2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('kt4PFL2Relative', 
        'kt4PFL3Absolute', 
        'kt4PFResidual')
)


process.kt4PFL2Relative = cms.ESSource("LXXXCorrectionService",
    section = cms.string(''),
    era = cms.string('Spring10'),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('KT4PF'),
    level = cms.string('L2Relative')
)


process.kt4PFL3Absolute = cms.ESSource("LXXXCorrectionService",
    section = cms.string(''),
    era = cms.string('Spring10'),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('KT4PF'),
    level = cms.string('L3Absolute')
)


process.kt4PFL6SLB = cms.ESSource("L6SLBCorrectionService",
    srcBTagInfoElectron = cms.InputTag("kt4PFJetsSoftElectronTagInfos"),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    section = cms.string(''),
    srcBTagInfoMuon = cms.InputTag("kt4PFJetsSoftMuonTagInfos"),
    addMuonToJet = cms.bool(False),
    era = cms.string('')
)


process.kt4PFResidual = cms.ESSource("LXXXCorrectionService",
    section = cms.string(''),
    era = cms.string('Spring10DataV2'),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5PF'),
    level = cms.string('L2L3Residual')
)


process.kt6CaloL1FastL2L3 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak5CaloL1Fastjet', 
        'kt6CaloL2Relative', 
        'kt6CaloL3Absolute')
)


process.kt6CaloL1FastL2L3L6 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('kt6CaloL1Offset', 
        'kt6CaloL2Relative', 
        'kt6CaloL3Absolute', 
        'kt6CaloL6SLB')
)


process.kt6CaloL1FastL2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('kt6CaloL1Fastjet', 
        'kt6CaloL2Relative', 
        'kt6CaloL3Absolute', 
        'kt6CaloResidual')
)


process.kt6CaloL1Fastjet = cms.ESSource("L1FastjetCorrectionService",
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1FastJet'),
    section = cms.string(''),
    srcRho = cms.InputTag("kt6PFJets","rho"),
    era = cms.string('Jec10V1')
)


process.kt6CaloL1L2L3 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('kt6CaloL1Offset', 
        'kt6CaloL2Relative', 
        'kt6CaloL3Absolute')
)


process.kt6CaloL1L2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('kt6CaloL1Offset', 
        'kt6CaloL2Relative', 
        'kt6CaloL3Absolute', 
        'kt6CaloResidual')
)


process.kt6CaloL1Offset = cms.ESSource("L1OffsetCorrectionService",
    minVtxNdof = cms.int32(4),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1Offset'),
    section = cms.string(''),
    vertexCollection = cms.string('offlinePrimaryVertices'),
    era = cms.string('Fall10')
)


process.kt6CaloL2L3 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('kt6CaloL2Relative', 
        'kt6CaloL3Absolute')
)


process.kt6CaloL2L3L6 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('kt6CaloL2Relative', 
        'kt6CaloL3Absolute', 
        'kt6CaloL6SLB')
)


process.kt6CaloL2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('kt6CaloL2Relative', 
        'kt6CaloL3Absolute', 
        'kt6CaloResidual')
)


process.kt6CaloL2Relative = cms.ESSource("LXXXCorrectionService",
    useCondDB = cms.untracked.bool(True),
    era = cms.string('Spring10'),
    section = cms.string(''),
    algorithm = cms.string('KT6Calo'),
    level = cms.string('L2Relative')
)


process.kt6CaloL3Absolute = cms.ESSource("LXXXCorrectionService",
    useCondDB = cms.untracked.bool(True),
    era = cms.string('Spring10'),
    section = cms.string(''),
    algorithm = cms.string('KT6Calo'),
    level = cms.string('L3Absolute')
)


process.kt6CaloL6SLB = cms.ESSource("L6SLBCorrectionService",
    srcBTagInfoElectron = cms.InputTag("kt6CaloJetsSoftElectronTagInfos"),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    section = cms.string(''),
    srcBTagInfoMuon = cms.InputTag("kt6CaloJetsSoftMuonTagInfos"),
    addMuonToJet = cms.bool(True),
    era = cms.string('')
)


process.kt6CaloResidual = cms.ESSource("LXXXCorrectionService",
    useCondDB = cms.untracked.bool(True),
    era = cms.string('Spring10DataV2'),
    section = cms.string(''),
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L2L3Residual')
)


process.kt6PFL1FastL2L3 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak5PFL1Fastjet', 
        'kt6PFL2Relative', 
        'kt6PFL3Absolute')
)


process.kt6PFL1FastL2L3L6 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('ak5PFL1Fastjet', 
        'kt6PFL2Relative', 
        'kt6PFL3Absolute', 
        'kt6PFL6SLB')
)


process.kt6PFL1FastL2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('kt6PFL1Fastjet', 
        'kt6PFL2Relative', 
        'kt6PFL3Absolute', 
        'kt6PFResidual')
)


process.kt6PFL1Fastjet = cms.ESSource("L1FastjetCorrectionService",
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5PF'),
    level = cms.string('L1FastJet'),
    section = cms.string(''),
    srcRho = cms.InputTag("kt6PFJets","rho"),
    era = cms.string('Jec10V1')
)


process.kt6PFL1L2L3 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('kt6PFL1Offset', 
        'kt6PFL2Relative', 
        'kt6PFL3Absolute')
)


process.kt6PFL1L2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('kt6PFL1Offset', 
        'kt6PFL2Relative', 
        'kt6PFL3Absolute', 
        'kt6PFResidual')
)


process.kt6PFL1Offset = cms.ESSource("L1OffsetCorrectionService",
    minVtxNdof = cms.int32(4),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5PF'),
    level = cms.string('L1Offset'),
    section = cms.string(''),
    vertexCollection = cms.string('offlinePrimaryVertices'),
    era = cms.string('Fall10')
)


process.kt6PFL2L3 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('kt6PFL2Relative', 
        'kt6PFL3Absolute')
)


process.kt6PFL2L3L6 = cms.ESSource("JetCorrectionServiceChain",
    useCondDB = cms.untracked.bool(True),
    correctors = cms.vstring('kt6PFL2Relative', 
        'kt6PFL3Absolute', 
        'kt6PFL6SLB')
)


process.kt6PFL2L3Residual = cms.ESSource("JetCorrectionServiceChain",
    correctors = cms.vstring('kt6PFL2Relative', 
        'kt6PFL3Absolute', 
        'kt6PFResidual')
)


process.kt6PFL2Relative = cms.ESSource("LXXXCorrectionService",
    section = cms.string(''),
    era = cms.string('Spring10'),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('KT6PF'),
    level = cms.string('L2Relative')
)


process.kt6PFL3Absolute = cms.ESSource("LXXXCorrectionService",
    section = cms.string(''),
    era = cms.string('Spring10'),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('KT6PF'),
    level = cms.string('L3Absolute')
)


process.kt6PFL6SLB = cms.ESSource("L6SLBCorrectionService",
    srcBTagInfoElectron = cms.InputTag("kt6PFJetsSoftElectronTagInfos"),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    section = cms.string(''),
    srcBTagInfoMuon = cms.InputTag("kt6PFJetsSoftMuonTagInfos"),
    addMuonToJet = cms.bool(False),
    era = cms.string('')
)


process.kt6PFResidual = cms.ESSource("LXXXCorrectionService",
    section = cms.string(''),
    era = cms.string('Spring10DataV2'),
    useCondDB = cms.untracked.bool(True),
    algorithm = cms.string('AK5PF'),
    level = cms.string('L2L3Residual')
)


process.magfield = cms.ESSource("XMLIdealGeometryESSource",
    geomXMLFiles = cms.vstring('Geometry/CMSCommonData/data/normal/cmsextent.xml', 
        'Geometry/CMSCommonData/data/cms.xml', 
        'Geometry/CMSCommonData/data/cmsMagneticField.xml', 
        'MagneticField/GeomBuilder/data/MagneticFieldVolumes_1103l.xml', 
        'MagneticField/GeomBuilder/data/MagneticFieldParameters_07_2pi.xml', 
        'Geometry/CMSCommonData/data/materials.xml'),
    rootNodeName = cms.string('cmsMagneticField:MAGF')
)


process.prefer("magfield")

process.AnomalousCellParameters = cms.PSet(
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999)
)

process.CondDBSetup = cms.PSet(
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string('.'),
        enableReadOnlySessionOnUpdateConnection = cms.untracked.bool(False),
        idleConnectionCleanupPeriod = cms.untracked.int32(10),
        messageLevel = cms.untracked.int32(0),
        enablePoolAutomaticCleanUp = cms.untracked.bool(False),
        enableConnectionSharing = cms.untracked.bool(True),
        connectionRetrialTimeOut = cms.untracked.int32(60),
        connectionTimeOut = cms.untracked.int32(60),
        connectionRetrialPeriod = cms.untracked.int32(10)
    )
)

process.GenJetParameters = cms.PSet(
    Active_Area_Repeats = cms.int32(5),
    src = cms.InputTag("genParticlesForJets"),
    doAreaFastjet = cms.bool(False),
    doPVCorrection = cms.bool(False),
    Ghost_EtaMax = cms.double(6.0),
    doRhoFastjet = cms.bool(False),
    srcPVs = cms.InputTag(""),
    inputEtMin = cms.double(0.0),
    doPUOffsetCorr = cms.bool(False),
    nSigmaPU = cms.double(1.0),
    radiusPU = cms.double(0.5),
    Rho_EtaMax = cms.double(4.5),
    jetPtMin = cms.double(3.0),
    jetType = cms.string('GenJet'),
    GhostArea = cms.double(0.01),
    inputEMin = cms.double(0.0)
)

process.JPTZSPCorrectorICone5 = cms.PSet(
    VectorialCorrection = cms.bool(True),
    ElectronIds = cms.InputTag("JPTeidTight"),
    UseMuons = cms.bool(True),
    Muons = cms.InputTag("muons"),
    UseTrackQuality = cms.bool(True),
    JetTracksAssociationAtCaloFace = cms.InputTag("iterativeCone5JetTracksAssociatorAtCaloFace"),
    LeakageMap = cms.string('CondFormats/JetMETObjects/data/CMSSW_362_TrackLeakage.txt'),
    UseOutOfConeTracks = cms.bool(True),
    UseInConeTracks = cms.bool(True),
    UseOutOfVertexTracks = cms.bool(True),
    UseResponseInVecCorr = cms.bool(False),
    ResponseMap = cms.string('CondFormats/JetMETObjects/data/CMSSW_31X_resptowers.txt'),
    EfficiencyMap = cms.string('CondFormats/JetMETObjects/data/CMSSW_362_TrackNonEff.txt'),
    UsePions = cms.bool(True),
    Electrons = cms.InputTag("gsfElectrons"),
    JetSplitMerge = cms.int32(0),
    DzVertexCut = cms.double(0.2),
    MaxJetEta = cms.double(3.0),
    Verbose = cms.bool(True),
    UseElectrons = cms.bool(True),
    JetTracksAssociationAtVertex = cms.InputTag("iterativeCone5JetTracksAssociatorAtVertex"),
    PtErrorQuality = cms.double(0.05),
    TrackQuality = cms.string('highPurity'),
    UseEfficiency = cms.bool(False)
)

process.METSignificance_params = cms.PSet(
    HB_EtResPar = cms.vdouble(0.0, 1.22, 0.05),
    EE_EtResPar = cms.vdouble(0.2, 0.03, 0.005),
    HF_PhiResPar = cms.vdouble(0.05022),
    PF_PhiResType7 = cms.vdouble(0.02511),
    HE_PhiResPar = cms.vdouble(0.02511),
    PF_PhiResType2 = cms.vdouble(0.002),
    PF_PhiResType3 = cms.vdouble(0.002),
    HF_EtResPar = cms.vdouble(0.0, 1.82, 0.09),
    PF_PhiResType1 = cms.vdouble(0.002),
    PF_PhiResType6 = cms.vdouble(0.02511),
    HO_EtResPar = cms.vdouble(0.0, 1.3, 0.005),
    PF_PhiResType4 = cms.vdouble(0.005),
    PF_PhiResType5 = cms.vdouble(0.02511),
    EB_EtResPar = cms.vdouble(0.2, 0.03, 0.005),
    EB_PhiResPar = cms.vdouble(0.00502),
    PF_EtResType6 = cms.vdouble(0.0, 1.22, 0.05),
    HO_PhiResPar = cms.vdouble(0.02511),
    EE_PhiResPar = cms.vdouble(0.02511),
    HB_PhiResPar = cms.vdouble(0.02511),
    PF_EtResType5 = cms.vdouble(0.0, 1.22, 0.05),
    PF_EtResType4 = cms.vdouble(0.2, 0.03, 0.005),
    PF_EtResType7 = cms.vdouble(0.0, 1.22, 0.05),
    HE_EtResPar = cms.vdouble(0.0, 1.3, 0.05),
    PF_EtResType1 = cms.vdouble(0.05, 0, 0),
    PF_EtResType3 = cms.vdouble(0.05, 0, 0),
    PF_EtResType2 = cms.vdouble(0.05, 0, 0)
)

process.OneProngNoPiZero = cms.PSet(
    applyIsolation = cms.bool(False),
    cut = cms.double(-10.0),
    computerName = cms.string('OneProngNoPiZero'),
    decayModeIndices = cms.vint32(0)
)

process.OneProngNoPiZeroIso = cms.PSet(
    applyIsolation = cms.bool(True),
    cut = cms.double(-10.0),
    computerName = cms.string('OneProngNoPiZeroIso'),
    decayModeIndices = cms.vint32(0)
)

process.OneProngOnePiZero = cms.PSet(
    applyIsolation = cms.bool(False),
    cut = cms.double(-10.0),
    computerName = cms.string('OneProngOnePiZero'),
    decayModeIndices = cms.vint32(1)
)

process.OneProngOnePiZeroIso = cms.PSet(
    applyIsolation = cms.bool(True),
    cut = cms.double(-10.0),
    computerName = cms.string('OneProngOnePiZeroIso'),
    decayModeIndices = cms.vint32(1)
)

process.OneProngTwoPiZero = cms.PSet(
    applyIsolation = cms.bool(False),
    cut = cms.double(-10.0),
    computerName = cms.string('OneProngTwoPiZero'),
    decayModeIndices = cms.vint32(2)
)

process.OneProngTwoPiZeroIso = cms.PSet(
    applyIsolation = cms.bool(True),
    cut = cms.double(-10.0),
    computerName = cms.string('OneProngTwoPiZeroIso'),
    decayModeIndices = cms.vint32(2)
)

process.PFTauQualityCuts = cms.PSet(
    isolationQualityCuts = cms.PSet(
        minTrackHits = cms.uint32(8),
        minTrackPt = cms.double(1.0),
        maxTrackChi2 = cms.double(100.0),
        minTrackPixelHits = cms.uint32(0),
        minGammaEt = cms.double(1.5),
        useTracksInsteadOfPFHadrons = cms.bool(False),
        maxDeltaZ = cms.double(0.2),
        maxTransverseImpactParameter = cms.double(0.03)
    ),
    signalQualityCuts = cms.PSet(
        minTrackHits = cms.uint32(3),
        minTrackPt = cms.double(0.5),
        maxTrackChi2 = cms.double(100.0),
        minTrackPixelHits = cms.uint32(0),
        minGammaEt = cms.double(0.5),
        useTracksInsteadOfPFHadrons = cms.bool(False),
        maxDeltaZ = cms.double(0.2),
        maxTransverseImpactParameter = cms.double(0.03)
    )
)

process.SingleNet = cms.PSet(
    applyIsolation = cms.bool(False),
    cut = cms.double(-10.0),
    computerName = cms.string('SingleNet'),
    decayModeIndices = cms.vint32(0, 1, 2, 10, 11)
)

process.SingleNetIso = cms.PSet(
    applyIsolation = cms.bool(True),
    cut = cms.double(-10.0),
    computerName = cms.string('SingleNetIso'),
    decayModeIndices = cms.vint32(0, 1, 2, 10, 11)
)

process.ThreeProngNoPiZero = cms.PSet(
    applyIsolation = cms.bool(False),
    cut = cms.double(-10.0),
    computerName = cms.string('ThreeProngNoPiZero'),
    decayModeIndices = cms.vint32(10)
)

process.ThreeProngNoPiZeroIso = cms.PSet(
    applyIsolation = cms.bool(True),
    cut = cms.double(-10.0),
    computerName = cms.string('ThreeProngNoPiZeroIso'),
    decayModeIndices = cms.vint32(10)
)

process.ThreeProngOnePiZero = cms.PSet(
    applyIsolation = cms.bool(False),
    cut = cms.double(-10.0),
    computerName = cms.string('ThreeProngOnePiZero'),
    decayModeIndices = cms.vint32(11)
)

process.ThreeProngOnePiZeroIso = cms.PSet(
    applyIsolation = cms.bool(True),
    cut = cms.double(-10.0),
    computerName = cms.string('ThreeProngOnePiZeroIso'),
    decayModeIndices = cms.vint32(11)
)

process.TrackAssociatorParameterBlock = cms.PSet(
    TrackAssociatorParameters = cms.PSet(
        muonMaxDistanceSigmaX = cms.double(0.0),
        muonMaxDistanceSigmaY = cms.double(0.0),
        CSCSegmentCollectionLabel = cms.InputTag("cscSegments"),
        dRHcal = cms.double(9999.0),
        dRPreshowerPreselection = cms.double(0.2),
        CaloTowerCollectionLabel = cms.InputTag("towerMaker"),
        useEcal = cms.bool(True),
        dREcal = cms.double(9999.0),
        dREcalPreselection = cms.double(0.05),
        HORecHitCollectionLabel = cms.InputTag("horeco"),
        dRMuon = cms.double(9999.0),
        propagateAllDirections = cms.bool(True),
        muonMaxDistanceX = cms.double(5.0),
        muonMaxDistanceY = cms.double(5.0),
        useHO = cms.bool(True),
        trajectoryUncertaintyTolerance = cms.double(-1.0),
        usePreshower = cms.bool(False),
        DTRecSegment4DCollectionLabel = cms.InputTag("dt4DSegments"),
        EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
        dRHcalPreselection = cms.double(0.2),
        useMuon = cms.bool(True),
        useCalo = cms.bool(False),
        accountForTrajectoryChangeCalo = cms.bool(False),
        EBRecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
        dRMuonPreselection = cms.double(0.2),
        truthMatch = cms.bool(False),
        HBHERecHitCollectionLabel = cms.InputTag("hbhereco"),
        useHcal = cms.bool(True)
    )
)

process.TrackAssociatorParameters = cms.PSet(
    muonMaxDistanceSigmaX = cms.double(0.0),
    muonMaxDistanceSigmaY = cms.double(0.0),
    CSCSegmentCollectionLabel = cms.InputTag("cscSegments"),
    dRHcal = cms.double(9999.0),
    dREcal = cms.double(9999.0),
    CaloTowerCollectionLabel = cms.InputTag("towerMaker"),
    useEcal = cms.bool(True),
    dREcalPreselection = cms.double(0.05),
    HORecHitCollectionLabel = cms.InputTag("horeco"),
    dRMuon = cms.double(9999.0),
    propagateAllDirections = cms.bool(True),
    muonMaxDistanceX = cms.double(5.0),
    muonMaxDistanceY = cms.double(5.0),
    useHO = cms.bool(True),
    trajectoryUncertaintyTolerance = cms.double(-1.0),
    usePreshower = cms.bool(False),
    DTRecSegment4DCollectionLabel = cms.InputTag("dt4DSegments"),
    EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    dRHcalPreselection = cms.double(0.2),
    useMuon = cms.bool(True),
    useCalo = cms.bool(False),
    accountForTrajectoryChangeCalo = cms.bool(False),
    EBRecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    dRMuonPreselection = cms.double(0.2),
    truthMatch = cms.bool(False),
    HBHERecHitCollectionLabel = cms.InputTag("hbhereco"),
    useHcal = cms.bool(True)
)

process.combinedSecondaryVertexCommon = cms.PSet(
    trackPseudoSelection = cms.PSet(
        totalHitsMin = cms.uint32(0),
        jetDeltaRMax = cms.double(0.3),
        qualityClass = cms.string('highPurity'),
        pixelHitsMin = cms.uint32(0),
        sip3dSigMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        maxDistToAxis = cms.double(0.07),
        sip2dValMax = cms.double(99999.9),
        maxDecayLen = cms.double(5),
        ptMin = cms.double(0.0),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(2.0),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        sip2dValMin = cms.double(-99999.9),
        normChi2Max = cms.double(99999.9)
    ),
    trackSelection = cms.PSet(
        totalHitsMin = cms.uint32(0),
        jetDeltaRMax = cms.double(0.3),
        qualityClass = cms.string('highPurity'),
        pixelHitsMin = cms.uint32(0),
        sip3dSigMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        maxDistToAxis = cms.double(0.07),
        sip2dValMax = cms.double(99999.9),
        maxDecayLen = cms.double(5),
        ptMin = cms.double(0.0),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        sip2dValMin = cms.double(-99999.9),
        normChi2Max = cms.double(99999.9)
    ),
    useTrackWeights = cms.bool(True),
    pseudoMultiplicityMin = cms.uint32(2),
    correctVertexMass = cms.bool(True),
    trackPairV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.03)
    ),
    charmCut = cms.double(1.5),
    vertexFlip = cms.bool(False),
    minimumTrackWeight = cms.double(0.5),
    pseudoVertexV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.05)
    ),
    trackMultiplicityMin = cms.uint32(3),
    trackSort = cms.string('sip2dSig'),
    trackFlip = cms.bool(False)
)

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    annotation = cms.untracked.string('SUSY pattuple definition'),
    name = cms.untracked.string('$Source: /cvs_server/repositories/CMSSW/UserCode/JSturdy/SUSY/AnalysisNtuplePAT/test/SusyPAT_DATA_new.py,v $')
)

process.fieldScaling = cms.PSet(
    scalingVolumes = cms.vint32(14100, 14200, 17600, 17800, 17900, 
        18100, 18300, 18400, 18600, 23100, 
        23300, 23400, 23600, 23800, 23900, 
        24100, 28600, 28800, 28900, 29100, 
        29300, 29400, 29600, 28609, 28809, 
        28909, 29109, 29309, 29409, 29609, 
        28610, 28810, 28910, 29110, 29310, 
        29410, 29610, 28611, 28811, 28911, 
        29111, 29311, 29411, 29611),
    scalingFactors = cms.vdouble(1, 1, 0.994, 1.004, 1.004, 
        1.005, 1.004, 1.004, 0.994, 0.965, 
        0.958, 0.958, 0.953, 0.958, 0.958, 
        0.965, 0.918, 0.924, 0.924, 0.906, 
        0.924, 0.924, 0.918, 0.991, 0.998, 
        0.998, 0.978, 0.998, 0.998, 0.991, 
        0.991, 0.998, 0.998, 0.978, 0.998, 
        0.998, 0.991, 0.991, 0.998, 0.998, 
        0.978, 0.998, 0.998, 0.991)
)

process.ghostTrackCommon = cms.PSet(
    trackSelection = cms.PSet(
        totalHitsMin = cms.uint32(0),
        jetDeltaRMax = cms.double(0.3),
        qualityClass = cms.string('highPurity'),
        pixelHitsMin = cms.uint32(0),
        sip3dSigMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        maxDistToAxis = cms.double(0.07),
        sip2dValMax = cms.double(99999.9),
        maxDecayLen = cms.double(5),
        ptMin = cms.double(0.0),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        sip2dValMin = cms.double(-99999.9),
        normChi2Max = cms.double(99999.9)
    ),
    trackPairV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.03)
    ),
    charmCut = cms.double(1.5),
    trackSort = cms.string('sip2dSig'),
    minimumTrackWeight = cms.double(0.5)
)

process.ghostTrackVertexRecoBlock = cms.PSet(
    vertexReco = cms.PSet(
        primcut = cms.double(2.0),
        seccut = cms.double(4.0),
        maxFitChi2 = cms.double(10.0),
        fitType = cms.string('RefitGhostTrackWithVertices'),
        mergeThreshold = cms.double(3.0),
        finder = cms.string('gtvr')
    )
)

process.j2tParametersCALO = cms.PSet(
    trackQuality = cms.string('goodIterative'),
    tracks = cms.InputTag("generalTracks"),
    coneSize = cms.double(0.5),
    extrapolations = cms.InputTag("trackExtrapolator")
)

process.j2tParametersVX = cms.PSet(
    tracks = cms.InputTag("generalTracks"),
    coneSize = cms.double(0.5)
)

process.jetAnalyzerPAT = cms.untracked.PSet(
    photonIso = cms.untracked.double(1.5),
    photonPt = cms.untracked.double(10.0),
    usePFJets = cms.untracked.bool(False),
    muonIso = cms.untracked.double(0.05),
    jetTag = cms.untracked.InputTag("cleanPatJetsAK5"),
    jetMaxEMF = cms.untracked.double(0.99),
    tauIso = cms.untracked.double(0.1),
    jetMinPt = cms.untracked.double(0.0),
    useCaloJets = cms.untracked.bool(False),
    useJPTJets = cms.untracked.bool(False),
    genJetTag = cms.untracked.InputTag("ak5GenJets"),
    electronIso = cms.untracked.double(0.1),
    prefixJets = cms.untracked.string(''),
    electronPt = cms.untracked.double(30.0),
    useTrackJets = cms.untracked.bool(False),
    debugJets = cms.untracked.int32(0),
    muonPt = cms.untracked.double(10.0),
    tauPt = cms.untracked.double(10.0),
    jetMinEMF = cms.untracked.double(0.01),
    doMCJets = cms.untracked.bool(False),
    jetCorTag = cms.untracked.string('AK5Calo'),
    jetMaxEta = cms.untracked.double(5.0)
)

process.leadTrackFinding = cms.PSet(
    cut = cms.double(0.5),
    Producer = cms.InputTag("pfRecoTauDiscriminationByLeadingTrackFinding")
)

process.leptonAnalyzerPAT = cms.untracked.PSet(
    tauMaxEt = cms.untracked.double(9999.0),
    elecMaxEta = cms.untracked.double(5.0),
    tauTag = cms.untracked.InputTag("cleanPatTaus"),
    elecTag = cms.untracked.InputTag("cleanPatElectrons"),
    elecMinEt = cms.untracked.double(2.5),
    tauMinEt = cms.untracked.double(5.0),
    muonRelIso = cms.untracked.double(0.1),
    muonTag = cms.untracked.InputTag("cleanPatMuons"),
    tauMaxEta = cms.untracked.double(5.0),
    tauRelIso = cms.untracked.double(0.1),
    prefixLeps = cms.untracked.string(''),
    elecRelIso = cms.untracked.double(0.5),
    muonMinEt = cms.untracked.double(2.5),
    muonMaxEt = cms.untracked.double(10.0),
    muonMaxEta = cms.untracked.double(5.0),
    debugLeps = cms.untracked.int32(0),
    elecMaxEt = cms.untracked.double(15.0)
)

process.loosePFTauQualityCuts = cms.PSet(
    isolationQualityCuts = cms.PSet(
        minTrackHits = cms.uint32(3),
        minTrackPt = cms.double(0.5),
        maxTrackChi2 = cms.double(100.0),
        minTrackPixelHits = cms.uint32(0),
        minGammaEt = cms.double(0.5),
        useTracksInsteadOfPFHadrons = cms.bool(False),
        maxDeltaZ = cms.double(0.2),
        maxTransverseImpactParameter = cms.double(0.03)
    ),
    signalQualityCuts = cms.PSet(
        minTrackHits = cms.uint32(3),
        minTrackPt = cms.double(0.5),
        maxTrackChi2 = cms.double(100.0),
        minTrackPixelHits = cms.uint32(0),
        minGammaEt = cms.double(0.5),
        useTracksInsteadOfPFHadrons = cms.bool(False),
        maxDeltaZ = cms.double(0.2),
        maxTransverseImpactParameter = cms.double(0.03)
    )
)

process.looseSoftPFElectronCleanerBarrelCuts = cms.PSet(
    BarreldRGsfTrackElectronCuts = cms.vdouble(0.0, 0.017),
    BarrelEemPinRatioCuts = cms.vdouble(-0.9, 0.39),
    BarrelMVACuts = cms.vdouble(-0.1, 1.0),
    BarrelPtCuts = cms.vdouble(2.0, 9999.0)
)

process.looseSoftPFElectronCleanerForwardCuts = cms.PSet(
    ForwarddRGsfTrackElectronCuts = cms.vdouble(0.0, 0.006),
    ForwardPtCuts = cms.vdouble(2.0, 9999.0),
    ForwardMVACuts = cms.vdouble(-0.24, 1.0),
    ForwardInverseFBremCuts = cms.vdouble(1.0, 7.01)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5000)
)

process.mcTruthAnalyzerPAT = cms.untracked.PSet(
    genParticleTag = cms.untracked.InputTag("genParticles"),
    debugTruth = cms.untracked.int32(0)
)

process.mediumPFTauQualityCuts = cms.PSet(
    isolationQualityCuts = cms.PSet(
        minTrackHits = cms.uint32(3),
        minTrackPt = cms.double(0.8),
        maxTrackChi2 = cms.double(100.0),
        minTrackPixelHits = cms.uint32(0),
        minGammaEt = cms.double(0.8),
        useTracksInsteadOfPFHadrons = cms.bool(False),
        maxDeltaZ = cms.double(0.2),
        maxTransverseImpactParameter = cms.double(0.03)
    ),
    signalQualityCuts = cms.PSet(
        minTrackHits = cms.uint32(3),
        minTrackPt = cms.double(0.8),
        maxTrackChi2 = cms.double(100.0),
        minTrackPixelHits = cms.uint32(0),
        minGammaEt = cms.double(0.5),
        useTracksInsteadOfPFHadrons = cms.bool(False),
        maxDeltaZ = cms.double(0.2),
        maxTransverseImpactParameter = cms.double(0.03)
    )
)

process.mediumSoftPFElectronCleanerBarrelCuts = cms.PSet(
    BarreldRGsfTrackElectronCuts = cms.vdouble(0.0, 0.0047),
    BarrelEemPinRatioCuts = cms.vdouble(-0.9, 0.54),
    BarrelMVACuts = cms.vdouble(0.6, 1.0),
    BarrelPtCuts = cms.vdouble(2.0, 9999.0)
)

process.mediumSoftPFElectronCleanerForwardCuts = cms.PSet(
    ForwarddRGsfTrackElectronCuts = cms.vdouble(0.0, 0.003),
    ForwardPtCuts = cms.vdouble(2.0, 9999.0),
    ForwardMVACuts = cms.vdouble(0.37, 1.0),
    ForwardInverseFBremCuts = cms.vdouble(1.0, 20.0)
)

process.metAnalyzerPAT = cms.untracked.PSet(
    metTag = cms.untracked.InputTag("patMETsAK5Calo"),
    prefixMET = cms.untracked.string('Calo'),
    doMCMET = cms.untracked.bool(False),
    debugMET = cms.untracked.int32(0),
    genMETTag = cms.untracked.InputTag("genMetCalo")
)

process.noPrediscriminants = cms.PSet(
    BooleanOperator = cms.string('and')
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False)
)

process.photonAnalyzerPAT = cms.untracked.PSet(
    prefixPhots = cms.untracked.string(''),
    photRelIso = cms.untracked.double(1.5),
    photTag = cms.untracked.InputTag("cleanPatPhotons"),
    debugPhots = cms.untracked.int32(0),
    photMinEt = cms.untracked.double(2.5),
    photMaxEta = cms.untracked.double(5.0),
    photMaxEt = cms.untracked.double(10000.0)
)

process.requireDecayMode = cms.PSet(
    BooleanOperator = cms.string('and'),
    decayMode = cms.PSet(
        cut = cms.double(0.5),
        Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding")
    )
)

process.requireLeadPion = cms.PSet(
    BooleanOperator = cms.string('and'),
    leadPion = cms.PSet(
        cut = cms.double(0.5),
        Producer = cms.InputTag("pfRecoTauDiscriminationByLeadingTrackFinding")
    )
)

process.requireLeadTrack = cms.PSet(
    BooleanOperator = cms.string('and'),
    leadTrack = cms.PSet(
        cut = cms.double(0.5),
        Producer = cms.InputTag("pfRecoTauDiscriminationByLeadingTrackFinding")
    )
)

process.requireLeadTrackCalo = cms.PSet(
    BooleanOperator = cms.string('and'),
    leadTrack = cms.PSet(
        cut = cms.double(0.5),
        Producer = cms.InputTag("caloRecoTauDiscriminationByLeadingTrackFinding")
    )
)

process.shrinkingConeLeadTrackFinding = cms.PSet(
    BooleanOperator = cms.string('and'),
    leadTrack = cms.PSet(
        cut = cms.double(0.5),
        Producer = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFinding")
    )
)

process.standardDecayModeParams = cms.PSet(
    mergeByBestMatch = cms.bool(True),
    refitTracks = cms.bool(False),
    maxNbrOfIterations = cms.int32(10),
    mergeLowPtPhotonsFirst = cms.bool(True),
    setMergedPi0Mass = cms.bool(True),
    setChargedPionMass = cms.bool(True),
    filterPhotons = cms.bool(True),
    minPtFractionPiZeroes = cms.double(0.15),
    maxPhotonsToMerge = cms.uint32(2),
    filterTwoProngs = cms.bool(True),
    maxPiZeroMass = cms.double(0.2),
    minPtFractionForSecondProng = cms.double(0.1),
    maxDistance = cms.double(0.01),
    setPi0Mass = cms.bool(True),
    minPtFractionSinglePhotons = cms.double(0.1)
)

process.tcTauAlgoParameters = cms.PSet(
    tkminTrackerHitsn = cms.int32(5),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    EtCaloOverTrackMin = cms.double(-0.9),
    EERecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    tkmaxChi2 = cms.double(100.0),
    EtHcalOverTrackMin = cms.double(-0.3),
    CaloRecoTauProducer = cms.InputTag("JPTCaloRecoTauProducer"),
    EBRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    tkminPixelHitsn = cms.int32(0),
    MatchingConeSize = cms.double(0.1),
    TrackCollection = cms.InputTag("generalTracks"),
    HBHERecHitCollection = cms.InputTag("hbhereco"),
    EtCaloOverTrackMax = cms.double(0.0),
    TrackAssociatorParameters = cms.PSet(
        muonMaxDistanceSigmaX = cms.double(0.0),
        muonMaxDistanceSigmaY = cms.double(0.0),
        CSCSegmentCollectionLabel = cms.InputTag("cscSegments"),
        dRHcal = cms.double(9999.0),
        dRPreshowerPreselection = cms.double(0.2),
        CaloTowerCollectionLabel = cms.InputTag("towerMaker"),
        useEcal = cms.bool(True),
        dREcal = cms.double(9999.0),
        dREcalPreselection = cms.double(0.05),
        HORecHitCollectionLabel = cms.InputTag("horeco"),
        dRMuon = cms.double(9999.0),
        propagateAllDirections = cms.bool(True),
        muonMaxDistanceX = cms.double(5.0),
        muonMaxDistanceY = cms.double(5.0),
        useHO = cms.bool(True),
        trajectoryUncertaintyTolerance = cms.double(-1.0),
        usePreshower = cms.bool(False),
        DTRecSegment4DCollectionLabel = cms.InputTag("dt4DSegments"),
        EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
        dRHcalPreselection = cms.double(0.2),
        useMuon = cms.bool(True),
        useCalo = cms.bool(False),
        accountForTrajectoryChangeCalo = cms.bool(False),
        EBRecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
        dRMuonPreselection = cms.double(0.2),
        truthMatch = cms.bool(False),
        HBHERecHitCollectionLabel = cms.InputTag("hbhereco"),
        useHcal = cms.bool(True)
    ),
    SignalConeSize = cms.double(0.2),
    EcalConeSize = cms.double(0.5),
    Track_minPt = cms.double(1.0),
    EtHcalOverTrackMax = cms.double(1.0),
    HFRecHitCollection = cms.InputTag("hfreco"),
    DropRejectedJets = cms.untracked.bool(False),
    HORecHitCollection = cms.InputTag("horeco"),
    DropCaloJets = cms.untracked.bool(False),
    tkmaxipt = cms.double(0.1)
)

process.tightSoftPFElectronCleanerBarrelCuts = cms.PSet(
    BarreldRGsfTrackElectronCuts = cms.vdouble(0.0, 0.006),
    BarrelEemPinRatioCuts = cms.vdouble(-0.9, 0.065),
    BarrelMVACuts = cms.vdouble(0.58, 1.0),
    BarrelPtCuts = cms.vdouble(2.0, 9999.0)
)

process.tightSoftPFElectronCleanerForwardCuts = cms.PSet(
    ForwarddRGsfTrackElectronCuts = cms.vdouble(0.0, 0.01),
    ForwardPtCuts = cms.vdouble(2.0, 9999.0),
    ForwardMVACuts = cms.vdouble(0.6, 1.0),
    ForwardInverseFBremCuts = cms.vdouble(1.0, 15.0)
)

process.trackAnalyzerPAT = cms.untracked.PSet(
    debugTrack = cms.untracked.int32(0),
    doMCTracks = cms.untracked.bool(False),
    trackTag = cms.untracked.InputTag("generalTracks")
)

process.trackPseudoSelectionBlock = cms.PSet(
    trackPseudoSelection = cms.PSet(
        totalHitsMin = cms.uint32(0),
        jetDeltaRMax = cms.double(0.3),
        qualityClass = cms.string('highPurity'),
        pixelHitsMin = cms.uint32(0),
        sip3dSigMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        maxDistToAxis = cms.double(0.07),
        sip2dValMax = cms.double(99999.9),
        maxDecayLen = cms.double(5),
        ptMin = cms.double(0.0),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(2.0),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        sip2dValMin = cms.double(-99999.9),
        normChi2Max = cms.double(99999.9)
    )
)

process.trackSelectionBlock = cms.PSet(
    trackSelection = cms.PSet(
        totalHitsMin = cms.uint32(0),
        jetDeltaRMax = cms.double(0.3),
        qualityClass = cms.string('highPurity'),
        pixelHitsMin = cms.uint32(0),
        sip3dSigMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        maxDistToAxis = cms.double(0.07),
        sip2dValMax = cms.double(99999.9),
        maxDecayLen = cms.double(5),
        ptMin = cms.double(0.0),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        sip2dValMin = cms.double(-99999.9),
        normChi2Max = cms.double(99999.9)
    )
)

process.triggerAnalyzerPAT = cms.untracked.PSet(
    hlTriggerResults = cms.untracked.InputTag("TriggerResults","","REDIGI"),
    l1TriggerResults = cms.untracked.InputTag("gtDigis"),
    getL1Info = cms.untracked.bool(False),
    getHLTfromConfig = cms.untracked.bool(False),
    debugTriggers = cms.untracked.int32(0)
)

process.vertexAnalyzerPAT = cms.untracked.PSet(
    minVtxNdof = cms.untracked.int32(4),
    beamspotTag = cms.untracked.InputTag("offlineBeamSpot"),
    maxVtxRho = cms.untracked.double(2.0),
    maxVtxZ = cms.untracked.double(24.0),
    vtxTag = cms.untracked.InputTag("offlinePrimaryVertices"),
    minVtxTrks = cms.untracked.int32(3),
    maxVtxChi2 = cms.untracked.double(999),
    minNVtx = cms.untracked.int32(1),
    debugVtx = cms.untracked.int32(0)
)

process.vertexCutsBlock = cms.PSet(
    vertexCuts = cms.PSet(
        distSig3dMax = cms.double(99999.9),
        fracPV = cms.double(0.65),
        distVal2dMax = cms.double(2.5),
        useTrackWeights = cms.bool(True),
        maxDeltaRToJetAxis = cms.double(0.5),
        v0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        ),
        distSig2dMin = cms.double(3.0),
        multiplicityMin = cms.uint32(2),
        massMax = cms.double(6.5),
        distSig2dMax = cms.double(99999.9),
        distVal3dMax = cms.double(99999.9),
        minimumTrackWeight = cms.double(0.5),
        distVal3dMin = cms.double(-99999.9),
        distVal2dMin = cms.double(0.01),
        distSig3dMin = cms.double(-99999.9)
    )
)

process.vertexRecoBlock = cms.PSet(
    vertexReco = cms.PSet(
        seccut = cms.double(6.0),
        primcut = cms.double(1.8),
        smoothing = cms.bool(False),
        weightthreshold = cms.double(0.001),
        minweight = cms.double(0.5),
        finder = cms.string('avr')
    )
)

process.vertexSelectionBlock = cms.PSet(
    vertexSelection = cms.PSet(
        sortCriterium = cms.string('dist3dError')
    )
)

process.vertexTrackSelectionBlock = cms.PSet(
    trackSelection = cms.PSet(
        totalHitsMin = cms.uint32(8),
        jetDeltaRMax = cms.double(0.3),
        qualityClass = cms.string('highPurity'),
        pixelHitsMin = cms.uint32(2),
        sip3dSigMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        maxDistToAxis = cms.double(0.2),
        sip2dValMax = cms.double(99999.9),
        maxDecayLen = cms.double(99999.9),
        ptMin = cms.double(1.0),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        sip2dValMin = cms.double(-99999.9),
        normChi2Max = cms.double(99999.9)
    )
)

process.MultiNetIso = cms.VPSet(cms.PSet(
    applyIsolation = cms.bool(True),
    cut = cms.double(-10.0),
    computerName = cms.string('OneProngNoPiZeroIso'),
    decayModeIndices = cms.vint32(0)
), 
    cms.PSet(
        applyIsolation = cms.bool(True),
        cut = cms.double(-10.0),
        computerName = cms.string('OneProngOnePiZeroIso'),
        decayModeIndices = cms.vint32(1)
    ), 
    cms.PSet(
        applyIsolation = cms.bool(True),
        cut = cms.double(-10.0),
        computerName = cms.string('OneProngTwoPiZeroIso'),
        decayModeIndices = cms.vint32(2)
    ), 
    cms.PSet(
        applyIsolation = cms.bool(True),
        cut = cms.double(-10.0),
        computerName = cms.string('ThreeProngNoPiZeroIso'),
        decayModeIndices = cms.vint32(10)
    ), 
    cms.PSet(
        applyIsolation = cms.bool(True),
        cut = cms.double(-10.0),
        computerName = cms.string('ThreeProngOnePiZeroIso'),
        decayModeIndices = cms.vint32(11)
    ))

process.SingleNetBasedTauID = cms.VPSet(cms.PSet(
    applyIsolation = cms.bool(False),
    cut = cms.double(-10.0),
    computerName = cms.string('SingleNet'),
    decayModeIndices = cms.vint32(0, 1, 2, 10, 11)
))

process.TaNC = cms.VPSet(cms.PSet(
    applyIsolation = cms.bool(False),
    cut = cms.double(-10.0),
    computerName = cms.string('OneProngNoPiZero'),
    decayModeIndices = cms.vint32(0)
), 
    cms.PSet(
        applyIsolation = cms.bool(False),
        cut = cms.double(-10.0),
        computerName = cms.string('OneProngOnePiZero'),
        decayModeIndices = cms.vint32(1)
    ), 
    cms.PSet(
        applyIsolation = cms.bool(False),
        cut = cms.double(-10.0),
        computerName = cms.string('OneProngTwoPiZero'),
        decayModeIndices = cms.vint32(2)
    ), 
    cms.PSet(
        applyIsolation = cms.bool(False),
        cut = cms.double(-10.0),
        computerName = cms.string('ThreeProngNoPiZero'),
        decayModeIndices = cms.vint32(10)
    ), 
    cms.PSet(
        applyIsolation = cms.bool(False),
        cut = cms.double(-10.0),
        computerName = cms.string('ThreeProngOnePiZero'),
        decayModeIndices = cms.vint32(11)
    ))

