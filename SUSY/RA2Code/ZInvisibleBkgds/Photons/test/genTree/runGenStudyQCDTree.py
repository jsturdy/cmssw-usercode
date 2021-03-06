import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

#===================== Message Logger =============================
process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
            limit = cms.untracked.int32(10),
            reportEvery = cms.untracked.int32(100)
            )
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet(
            wantSummary = cms.untracked.bool(True)
            )

#================= configure poolsource module ===================

#process.load('SusyAnalysis.PhotonAnalysis.PhotonRun2011AMay10ReReco_160404to163869_cfi');
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(FILELIST ))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/user/lpcsusyhad/QCD_Pt_15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12/samantha/QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6/Summer12_PU_S7_START52_V9_RA2_NoCuts_09Aug2012V1/f42b6bbcc7ea08f6809cc21efa861dac/susypat_881_1_wZl.root',
        '/store/user/lpcsusyhad/QCD_Pt_15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12/samantha/QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6/Summer12_PU_S7_START52_V9_RA2_NoCuts_09Aug2012V1/f42b6bbcc7ea08f6809cc21efa861dac/susypat_882_1_hxl.root',
        '/store/user/lpcsusyhad/QCD_Pt_15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12/samantha/QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6/Summer12_PU_S7_START52_V9_RA2_NoCuts_09Aug2012V1/f42b6bbcc7ea08f6809cc21efa861dac/susypat_883_1_RBc.root',
        '/store/user/lpcsusyhad/QCD_Pt_15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12/samantha/QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6/Summer12_PU_S7_START52_V9_RA2_NoCuts_09Aug2012V1/f42b6bbcc7ea08f6809cc21efa861dac/susypat_884_1_Umt.root',
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source.skipEvents = cms.untracked.uint32(0)

#========================= analysis module =====================================

scaleF = -1.0
from RA2Classic.WeightProducer.weightProducer_cfi import weightProducer
process.eventWeight = weightProducer.clone(
    weight = cms.double(scaleF),
    Method = cms.string("PtHat"),
    XS = cms.double(2.99913994E10),
    NumberEvts = cms.double(9998151),
    Lumi = cms.double(10000.),
    Exponent = cms.double(-4.5),
    LumiScale = cms.double( 1.0 ),
)
from RA2Classic.WeightProducer.puWeightProducer_cfi import puWeightProducer
process.puWeight = puWeightProducer.clone(
    weight = cms.double(1.0),
)
from ZInvisibleBkgds.Photons.genstudytree_cfi import *
process.directPhotonsID = genstudytree.clone(
    debug            = cms.bool(False),
    genSrc           = cms.InputTag("zinvBkgdDirectPhotons"),
    debugString      = cms.string("direct photons"),
    ScaleFactor      = cms.double(scaleF),
    recoPhotonSrc    = cms.InputTag("patPhotonsID"),
    recoJetSrc       = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30"),
    htJetSrc         = cms.InputTag("patJetsPFNoPhotonIDSpecialPt50Eta25"),
    bJetSrc          = cms.InputTag("patCSVTJetsPFNoPhotonIDSpecialPt30Eta24"),
    htNoBosonSource  = cms.InputTag("htPFchsNoPhotID"),
    mhtNoBosonSource = cms.InputTag("mhtPFchsNoPhotID"),
)
process.directPhotonsIDPFIso = process.directPhotonsID.clone(
    recoPhotonSrc    = cms.InputTag("patPhotonsIDPFIso"),
    recoJetSrc       = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30"),
    htJetSrc         = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt50Eta25"),
    bJetSrc          = cms.InputTag("patCSVTJetsPFNoPhotonIDPFIsoSpecialPt30Eta24"),
    htNoBosonSource  = cms.InputTag("htPFchsNoPhotIDPFIso"),
    mhtNoBosonSource = cms.InputTag("mhtPFchsNoPhotIDPFIso"),
)
process.directPhotonsIDNoVeto = process.directPhotonsID.clone()
process.directPhotonsIDPFIsoNoVeto = process.directPhotonsIDPFIso.clone()
process.secondaryPhotonsID = process.directPhotonsID.clone(
    genSrc           = cms.InputTag("zinvBkgdSecondaryPhotons"),
    debugString      = cms.string("secondary photons"))
process.secondaryPhotonsIDPFIso = process.directPhotonsIDPFIso.clone(
    genSrc           = cms.InputTag("zinvBkgdSecondaryPhotons"),
    debugString      = cms.string("secondary photons"))
process.secondaryPhotonsIDNoVeto = process.secondaryPhotonsID.clone()
process.secondaryPhotonsIDPFIsoNoVeto = process.secondaryPhotonsIDPFIso.clone()
process.fragmentationPhotonsID = process.directPhotonsID.clone(
    genSrc           = cms.InputTag("zinvBkgdFragmentationPhotons"),
    debugString      = cms.string("fragmentation photons"))
process.fragmentationPhotonsIDPFIso = process.directPhotonsIDPFIso.clone(
    genSrc           = cms.InputTag("zinvBkgdFragmentationPhotons"),
    debugString      = cms.string("fragmentation photons"))
process.fragmentationPhotonsIDNoVeto = process.fragmentationPhotonsID.clone()
process.fragmentationPhotonsIDPFIsoNoVeto = process.fragmentationPhotonsIDPFIso.clone()
#================ configure filters and analysis sequence=======================

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
#from SusyAnalysis.MyAnalysis.filterBoolean_cfi import *
#process.load("SusyAnalysis.MyAnalysis.filterBoolean_cfi")

process.load('ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdObjects_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')
from ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff import *
from ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff import *
from ZInvisibleBkgds.Photons.ZinvPhotonJets_cff import *
from ZInvisibleBkgds.Photons.ZinvBkgdJets_cff import *

process.load('ZInvisibleBkgds.Photons.adduserdata_cfi')
from ZInvisibleBkgds.Photons.addphotonuserdata_cfi import addphotonuserbasic
process.patPhotonsAlt = addphotonuserbasic.clone()

from ZInvisibleBkgds.Photons.photonmap_cfi import *
process.rhoToPhotonMap = photonmap.clone()
from ZInvisibleBkgds.Photons.addphotonuserdata_cfi import *
process.patPhotonsUser1 = addphotonuserdata1.clone()
process.patPhotonsUser1.photonLabel = cms.InputTag("patPhotonsAlt")
process.patPhotonsUser1.userData.userFloats = cms.PSet(
    src = cms.VInputTag(
        cms.InputTag("rhoToPhotonMap")
    )
)
process.patPhotonsUserData = addphotonuserdata2.clone()
process.patPhotonsUserData.photonLabel = cms.InputTag("patPhotonsUser1")

process.countDirectPhotons        = countPhotonsIDPFIso.clone(src = cms.InputTag("zinvBkgdDirectPhotons"))
process.countSecondaryPhotons     = countPhotonsIDPFIso.clone(src = cms.InputTag("zinvBkgdSecondaryPhotons"))
process.countFragmentationPhotons = countPhotonsIDPFIso.clone(src = cms.InputTag("zinvBkgdFragmentationPhotons"))

process.analysisSeq = cms.Sequence(  process.ra2PFchsJets
                                   * process.patPhotonsAlt
                                   * process.rhoToPhotonMap
                                   * process.patPhotonsUser1
                                   * process.patPhotonsUserData
                                   * process.zinvBkgdGenPhotons
#                                   * process.zinvBkgdGenZBosons
                                   * process.photonObjectsPF
                                   * process.zinvBJetsPFNoPhotonIDSpecial
                                   * process.zinvBJetsPFNoPhotonIDPFIsoSpecial
                                   * process.zinvBJetsPF
#                                   * process.countGenBosons
#                                   * process.directPhotonsIDNoVeto
#                                   * process.directPhotonsIDPFIsoNoVeto
#                                   * process.ra2MuonVeto
#                                   * process.ra2ElectronVeto
#                                   * process.directPhotonsID
#                                   * process.directPhotonsIDPFIso
)

process.directAnalysisSeq = cms.Sequence(process.countDirectPhotons
                                       * process.directPhotonsIDNoVeto
                                       * process.directPhotonsIDPFIsoNoVeto
                                       * process.ra2MuonVeto
                                       * process.ra2ElectronVeto
                                       * process.directPhotonsID
                                       * process.directPhotonsIDPFIso
)

process.secondaryAnalysisSeq = cms.Sequence(process.countSecondaryPhotons
                                          * process.secondaryPhotonsIDNoVeto
                                          * process.secondaryPhotonsIDPFIsoNoVeto
                                          * process.ra2MuonVeto
                                          * process.ra2ElectronVeto
                                          * process.secondaryPhotonsID
                                          * process.secondaryPhotonsIDPFIso
)

process.fragmentationAnalysisSeq = cms.Sequence(process.countFragmentationPhotons
                                              * process.fragmentationPhotonsIDNoVeto
                                              * process.fragmentationPhotonsIDPFIsoNoVeto
                                              * process.ra2MuonVeto
                                              * process.ra2ElectronVeto
                                              * process.fragmentationPhotonsID
                                              * process.fragmentationPhotonsIDPFIso
)

#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('qcd_flat_pt_15to3000_gen_tree.root')
)
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('qcdflatcontent.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p1')
    ),
    outputCommands = cms.untracked.vstring('keep *')
)
#=================== run range & HLT filters ===============================
#process.load('SusyAnalysis.PhotonAnalysis.Photon_RunRangeHLTSeq_cfi')

#process.load('SandBox.Utilities.puWeightProducer_cfi')
##process.puWeight.DataPileUpHistFile = "SandBox/Utilities/data/May10_Prompt167151_pudist.root"
#process.puWeight.DataPileUpHistFile = "SandBox/Utilities/data/Cert_160404-177515_JSON.pileup.root"

#============================== configure paths ===============================
#process.p1 = cms.Path( process.analysisSeq )
process.p1 = cms.Path(process.puWeight
                    * process.eventWeight
                    * process.analysisSeq )
process.pdirect        = cms.Path(process.directAnalysisSeq)
process.psecondary     = cms.Path(process.secondaryAnalysisSeq)
process.pfragmentation = cms.Path(process.fragmentationAnalysisSeq)

#process.outpath = cms.EndPath(process.out)
##file = open('wtf_gentree.py','w')
##file.write(str(process.dumpPython()))
##file.close()
