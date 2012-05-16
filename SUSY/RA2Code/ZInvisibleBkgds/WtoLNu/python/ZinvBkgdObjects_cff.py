import FWCore.ParameterSet.Config as cms

from ZInvisibleBkgds.WtoLNu.ZinvBkgdJets_cff import *
from ZInvisibleBkgds.WtoLNu.ZinvBkgdLeptons_cff import *
from ZInvisibleBkgds.WtoLNu.WLNuHT_cff import *
from ZInvisibleBkgds.WtoLNu.WLNuMHT_cff import *

wenuobjs = cms.Sequence(
    electronCleanedPFJets *
    htPFNoElecIDIso  *
    mhtPFNoElecIDIso *
    zinvBJets 
)

wmunuobjs = cms.Sequence(
    muonCleanedPFJets *
    htPFNoMuonIDIso   *
    mhtPFNoMuonIDIso  *
    zinvBJets 
)

wlnuobjs = cms.Sequence(
    muonCleanedPFJets     *
    electronCleanedPFJets *
    htPFNoElecIDIso  *
    mhtPFNoElecIDIso *
    htPFNoMuonIDIso  *
    mhtPFNoMuonIDIso *
    zinvBJets     
)
