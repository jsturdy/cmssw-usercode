
# muon filters
from SandBox.Skims.RA2Leptons_cff import *
from PhysicsTools.PatAlgos.selectionLayer1.muonCountFilter_cfi import *

patMuonsPFIDPt20Eta21 = patMuonsPFID.clone()
patMuonsPFIDPt20Eta21.MinMuPt  = 20.0
patMuonsPFIDPt20Eta21.MaxMuEta = 2.1

patMuonsPFIDIsoPt20Eta21 = patMuonsPFIDIso.clone()
patMuonsPFIDIsoPt20Eta21.MinMuPt  = 20.0
patMuonsPFIDIsoPt20Eta21.MaxMuEta = 2.1

muonsPFPt20Eta21 = cms.Sequence(
    patMuonsPFIDPt20Eta21+
    patMuonsPFIDIsoPt20Eta21
)

countPFMuonsIDWLNu = countPatMuons.clone()
countPFMuonsIDWLNu.src = cms.InputTag('patMuonsPFID')
countPFMuonsIDWLNu.minNumber = cms.uint32(1)
countPFMuonsIDWLNu.maxNumber = cms.uint32(1)

countPFMuonsIDIsoWLNu = countPatMuons.clone()
countPFMuonsIDIsoWLNu.src = cms.InputTag('patMuonsPFIDIso')
countPFMuonsIDIsoWLNu.minNumber = cms.uint32(1)
countPFMuonsIDIsoWLNu.maxNumber = cms.uint32(1)

countPFMuonsIDPt20Eta21WLNu = countPatMuons.clone()
countPFMuonsIDPt20Eta21WLNu.src = cms.InputTag('patMuonsPFIDPt20Eta21')
countPFMuonsIDPt20Eta21WLNu.minNumber = cms.uint32(1)
countPFMuonsIDPt20Eta21WLNu.maxNumber = cms.uint32(1)

countPFMuonsIDIsoPt20Eta21WLNu = countPatMuons.clone()
countPFMuonsIDIsoPt20Eta21WLNu.src = cms.InputTag('patMuonsPFIDIsoPt20Eta21')
countPFMuonsIDIsoPt20Eta21WLNu.minNumber = cms.uint32(1)
countPFMuonsIDIsoPt20Eta21WLNu.maxNumber = cms.uint32(1)

countWtoMuNu = cms.Sequence(
    #countPFMuonsIDIsoWLNu+
    countPFMuonsIDIsoPt20Eta21WLNu    
)
# electrons filters
from PhysicsTools.PatAlgos.selectionLayer1.electronCountFilter_cfi import *

patElectronsPFIDPt20Eta21 = patElectronsPFID.clone()
patElectronsPFIDPt20Eta21.MinElePt  = 20.0
patElectronsPFIDPt20Eta21.MaxEleEta = 2.1

patElectronsPFIDIsoPt20Eta21 = patElectronsPFIDIso.clone()
patElectronsPFIDIsoPt20Eta21.MinElePt  = 20.0
patElectronsPFIDIsoPt20Eta21.MaxEleEta = 2.1

electronsPFPt20Eta21 = cms.Sequence(
    patElectronsPFIDPt20Eta21+
    patElectronsPFIDIsoPt20Eta21
)

countPFElectronsIDWLNu = countPatElectrons.clone()
countPFElectronsIDWLNu.src = cms.InputTag('patElectronsPFID')
countPFElectronsIDWLNu.minNumber = cms.uint32(1)
countPFElectronsIDWLNu.maxNumber = cms.uint32(1)

countPFElectronsIDIsoWLNu = countPatElectrons.clone()
countPFElectronsIDIsoWLNu.src = cms.InputTag('patElectronsPFIDIso')
countPFElectronsIDIsoWLNu.minNumber = cms.uint32(1)
countPFElectronsIDIsoWLNu.maxNumber = cms.uint32(1)

countPFElectronsIDPt20Eta21WLNu = countPatElectrons.clone()
countPFElectronsIDPt20Eta21WLNu.src = cms.InputTag('patElectronsPFIDPt20Eta21')
countPFElectronsIDPt20Eta21WLNu.minNumber = cms.uint32(1)
countPFElectronsIDPt20Eta21WLNu.maxNumber = cms.uint32(1)

countPFElectronsIDIsoPt20Eta21WLNu = countPatElectrons.clone()
countPFElectronsIDIsoPt20Eta21WLNu.src = cms.InputTag('patElectronsPFIDIsoPt20Eta21')
countPFElectronsIDIsoPt20Eta21WLNu.minNumber = cms.uint32(1)
countPFElectronsIDIsoPt20Eta21WLNu.maxNumber = cms.uint32(1)

countWtoENu = cms.Sequence(
    countPFElectronsIDIsoWLNu
    #countPFElectronsIDIsoPt20Eta21WLNu
)
