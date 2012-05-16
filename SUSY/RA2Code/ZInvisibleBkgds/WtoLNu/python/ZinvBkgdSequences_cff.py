import FWCore.ParameterSet.Config as cms

from ZInvisibleBkgds.WtoLNu.ZinvBkgdJets_cff import *
from ZInvisibleBkgds.WtoLNu.ZinvBkgdLeptons_cff import *
from ZInvisibleBkgds.WtoLNu.ZinvBkgdObjects_cff import *
from ZInvisibleBkgds.WtoLNu.WLNuHT_cff import *
from ZInvisibleBkgds.WtoLNu.WLNuMHT_cff import *
from ZInvisibleBkgds.WtoLNu.WLNuJetMHTDPhi_cff import *

wenuseq_bveto = cms.Sequence(
    wenuobjs *
    zinvBVetoPt30 *
    countJetsAK5PFPt50Eta25NoPFElectronIDIso *
    countPFElectronsIDIsoWLNu *
    htPFNoElecIDIsoFilter *
    jetMHTPFNoElecIDIsoDPhiFilter *
    mhtPFNoElecIDIsoFilter 
)

wenuseq_bsamp = cms.Sequence(
    wenuobjs *
    countSSVHPTBJetsAK5PFPt30 *
    countJetsAK5PFPt50Eta25NoPFElectronIDIso *
    countPFElectronsIDIsoWLNu *
    htPFNoElecIDIsoFilter *
    jetMHTPFNoElecIDIsoDPhiFilter *
    mhtPFNoElecIDIsoFilter 
)

wenuseq = cms.Sequence(
    wenuobjs *
    countJetsAK5PFPt50Eta25NoPFElectronIDIso *
    countPFElectronsIDIsoWLNu *
    htPFNoElecIDIsoFilter *
    jetMHTPFNoElecIDIsoDPhiFilter *
    mhtPFNoElecIDIsoFilter 
)

#### muon sequences
wmunuseq_bveto = cms.Sequence(
    wmunuobjs *
    zinvBVetoPt30 *
    countJetsAK5PFPt50Eta25NoPFMuonIDIso *
    countPFMuonsIDIsoWLNu *
    htPFNoMuonIDIsoFilter *
    jetMHTPFNoMuonIDIsoDPhiFilter *
    mhtPFNoMuonIDIsoFilter 
)

wmunuseq_bsamp = cms.Sequence(
    wmunuobjs *
    countSSVHPTBJetsAK5PFPt30 *
    countJetsAK5PFPt50Eta25NoPFMuonIDIso *
    countPFMuonsIDIsoWLNu *
    htPFNoMuonIDIsoFilter *
    jetMHTPFNoMuonIDIsoDPhiFilter *
    mhtPFNoMuonIDIsoFilter 
)

wmunuseq = cms.Sequence(
    wmunuobjs *
    countJetsAK5PFPt50Eta25NoPFMuonIDIso *
    countPFMuonsIDIsoWLNu *
    htPFNoMuonIDIsoFilter *
    jetMHTPFNoMuonIDIsoDPhiFilter *
    mhtPFNoMuonIDIsoFilter 
)
