// -*- C++ -*-
//
// Package:    WtoLNu
// Class:      WtoLNu
// 
/**\class WtoLNu WtoLNu.cc ZInvisibleBkgds/WtoLNu/src/WtoLNu.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Jared Sturdy
//         Created:  Wed Apr 18 16:06:24 CDT 2012
// $Id$
//
//


#include "ZInvisibleBkgds/WtoLNu/interface/WtoLNu.h"


//
// constructors and destructor
//
WtoLNu::WtoLNu(const edm::ParameterSet& conf) :
  debug_         ( conf.getParameter< bool >( "debug" ) ),
  vertexLabel_   ( conf.getParameter< edm::InputTag >( "vertexLabel" ) ),
  //puWeightLabel_ ( conf.getParameter< edm::InputTag >( "puWeightLabel" ) ),
  jetLabel_      ( conf.getParameter< edm::InputTag >( "jetLabel" ) ),
  htLabel_       ( conf.getParameter< edm::InputTag >( "htLabel" ) ),
  mhtLabel_      ( conf.getParameter< edm::InputTag >( "mhtLabel" ) ),
  metLabel_      ( conf.getParameter< edm::InputTag >( "metLabel" ) ),
  muonLabel_     ( conf.getParameter< edm::InputTag >( "muonLabel" ) ),
  electronLabel_ ( conf.getParameter< edm::InputTag >( "electronLabel" ) )
{
  
}


WtoLNu::~WtoLNu()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void WtoLNu::produce(edm::Event& ev, const edm::EventSetup& es)
{
  using namespace edm;

  //soon//edm::Handle<double> thPUWeight;
  //soon//ev.getByLabel(puWeightLabel_,thePUWeight);
  //get the jets
  edm::Handle<pat::JetCollection> jets;
  ev.getByLabel(jetLabel_,jets);
  //get the ht
  edm::Handle<double> theHT;
  ev.getByLabel(htLabel_,theHT);
  //get the mht
  edm::Handle<edm::View<reco::MET> > theMHT;
  ev.getByLabel(mhtLabel_,theMHT);
  //get the met
  edm::Handle<edm::View<pat::MET> > theMET;
  //edm::Handle<reco::PFMET > theMET;
  ev.getByLabel(metLabel_,theMET);
  //get the muons
  edm::Handle<pat::MuonCollection> muons;
  ev.getByLabel(muonLabel_,muons);
  //get the electrons
  edm::Handle<pat::ElectronCollection> electrons;
  ev.getByLabel(electronLabel_,electrons);

  //get the vertices
  edm::Handle<reco::VertexCollection> vertices;
  ev.getByLabel(vertexLabel_,vertices);

  ///compute MTW from Electrons
  //TLorentzVector WCandElec;
  reco::Candidate::LorentzVector WCandElec;
  if (electrons->size() != 0) {
    WCandElec = (*electrons)[0].p4() + (*theMET)[0].p4();
    histos1D_[ "ElecEta" ]->Fill((*electrons)[0].eta());
    histos1D_[ "ElecPhi" ]->Fill((*electrons)[0].phi());
    histos1D_[ "ElecPt"  ]->Fill((*electrons)[0].pt() );
    histos1D_[ "eventMTWElec" ]->Fill(WCandElec.mt()  );
  }

  ///compute MTW from Muons
  //TLorentzVector WCandMuon;
  reco::Candidate::LorentzVector WCandMuon;
  if (muons->size() != 0) {
    WCandMuon = (*muons)[0].p4() + (*theMET)[0].p4();
    histos1D_[ "MuonEta" ]->Fill((*muons)[0].eta()  );
    histos1D_[ "MuonPhi" ]->Fill((*muons)[0].phi()  );
    histos1D_[ "MuonPt"  ]->Fill((*muons)[0].pt()   );
    histos1D_[ "eventMTWMuon" ]->Fill(WCandMuon.mt());
  }

  ////Make histograms
  histos1D_[ "nVertices" ]->Fill(vertices->size());
  histos1D_[ "nJetsPt30" ]->Fill(jets->size()    );

  if (jets->size() > 0) {
    histos1D_[ "Jet1Eta" ]->Fill((*jets)[0].eta());
    histos1D_[ "Jet1Phi" ]->Fill((*jets)[0].phi());
    histos1D_[ "Jet1Pt" ]->Fill((*jets)[0].pt());
    
    if (jets->size() > 3) {
      histos1D_[ "Jet2Eta" ]->Fill((*jets)[1].eta());
      histos1D_[ "Jet2Phi" ]->Fill((*jets)[1].phi());
      histos1D_[ "Jet2Pt" ]->Fill((*jets)[1].pt());
      
      if (jets->size() > 2) {
	histos1D_[ "Jet3Eta" ]->Fill((*jets)[2].eta());
	histos1D_[ "Jet3Phi" ]->Fill((*jets)[2].phi());
	histos1D_[ "Jet3Pt" ]->Fill((*jets)[2].pt());
      }
    }
  }
  histos1D_[ "nElecsIDIso" ]->Fill(electrons->size());
  histos1D_[ "nMuonsIDIso" ]->Fill(muons->size()    );

  histos1D_[ "eventHT" ]->Fill(*theHT);

  histos1D_[ "eventMHT"    ]->Fill((*theMHT)[0].pt() );
  histos1D_[ "eventMHTPhi" ]->Fill((*theMHT)[0].phi());

  histos1D_[ "eventMEff" ]->Fill(*theHT+(*theMHT)[0].pt());

  histos1D_[ "eventMET"    ]->Fill((*theMET)[0].pt() );
  histos1D_[ "eventMETPhi" ]->Fill((*theMET)[0].phi());

  ///2D
  histos2D_[ "eventHTvsMHT"  ]->Fill(*theHT,            (*theMHT)[0].pt());
  histos2D_[ "eventHTvsMET"  ]->Fill(*theHT,            (*theMET)[0].pt());
  histos2D_[ "eventMHTvsMET" ]->Fill((*theMHT)[0].pt(), (*theMET)[0].pt());
 
}

// ------------ method called once each job just before starting event loop  ------------
void WtoLNu::beginJob()
{
  edm::Service< TFileService > fs;
  
  //# vertices
  histos1D_[ "nVertices" ] = fs->make< TH1D >( "nVertices", "# Vertices", 75, 0., 75.);
  histos1D_[ "nVertices" ]->SetXTitle( "# Vertices" );
  histos1D_[ "nVertices" ]->SetYTitle( "# of events" );
  histos1D_[ "nVertices" ]->Sumw2();

  ///Jet variables
  //# jets
  histos1D_[ "nJetsPt30" ] = fs->make< TH1D >( "nJetsPt30", "# Jets", 25, 0., 25.);
  histos1D_[ "nJetsPt30" ]->SetXTitle( "# Jets" );
  histos1D_[ "nJetsPt30" ]->SetYTitle( "# of events" );
  histos1D_[ "nJetsPt30" ]->Sumw2();

  //Eta
  histos1D_[ "Jet1Eta" ] = fs->make< TH1D >( "Jet1Eta", "#eta^{J1}", 50, -5., 5.);
  histos1D_[ "Jet1Eta" ]->SetXTitle( "#eta Leading Jet" );
  histos1D_[ "Jet1Eta" ]->SetYTitle( "# of events" );
  histos1D_[ "Jet1Eta" ]->Sumw2();

  histos1D_[ "Jet2Eta" ] = fs->make< TH1D >( "Jet2Eta", "#eta^{J2}", 50, -5., 5.);
  histos1D_[ "Jet2Eta" ]->SetXTitle( "#eta Leading Jet" );
  histos1D_[ "Jet2Eta" ]->SetYTitle( "# of events" );
  histos1D_[ "Jet2Eta" ]->Sumw2();

  histos1D_[ "Jet3Eta" ] = fs->make< TH1D >( "Jet3Eta", "#eta^{J3}", 50, -5., 5.);
  histos1D_[ "Jet3Eta" ]->SetXTitle( "#eta Leading Jet" );
  histos1D_[ "Jet3Eta" ]->SetYTitle( "# of events" );
  histos1D_[ "Jet3Eta" ]->Sumw2();

  //Phi
  histos1D_[ "Jet1Phi" ] = fs->make< TH1D >( "Jet1Phi", "#phi^{J1}", 50, -M_PI, M_PI);
  histos1D_[ "Jet1Phi" ]->SetXTitle( "#phi Leading Jet" );
  histos1D_[ "Jet1Phi" ]->SetYTitle( "# of events" );
  histos1D_[ "Jet1Phi" ]->Sumw2();

  histos1D_[ "Jet2Phi" ] = fs->make< TH1D >( "Jet2Phi", "#phi^{J2}", 50, -M_PI, M_PI);
  histos1D_[ "Jet2Phi" ]->SetXTitle( "#phi Leading Jet" );
  histos1D_[ "Jet2Phi" ]->SetYTitle( "# of events" );
  histos1D_[ "Jet2Phi" ]->Sumw2();

  histos1D_[ "Jet3Phi" ] = fs->make< TH1D >( "Jet3Phi", "#phi^{J3}", 50, -M_PI, M_PI);
  histos1D_[ "Jet3Phi" ]->SetXTitle( "#phi Leading Jet" );
  histos1D_[ "Jet3Phi" ]->SetYTitle( "# of events" );
  histos1D_[ "Jet3Phi" ]->Sumw2();

  //Pt
  histos1D_[ "Jet1Pt" ] = fs->make< TH1D >( "Jet1Pt", "p_{T}^{J1}", 250, 0., 1500.);
  histos1D_[ "Jet1Pt" ]->SetXTitle( "p_{T} Leading Jet" );
  histos1D_[ "Jet1Pt" ]->SetYTitle( "# of events" );
  histos1D_[ "Jet1Pt" ]->Sumw2();

  histos1D_[ "Jet2Pt" ] = fs->make< TH1D >( "Jet2Pt", "p_{T}^{J2}", 250, 0., 1500.);
  histos1D_[ "Jet2Pt" ]->SetXTitle( "p_{T} Leading Jet" );
  histos1D_[ "Jet2Pt" ]->SetYTitle( "# of events" );
  histos1D_[ "Jet2Pt" ]->Sumw2();

  histos1D_[ "Jet3Pt" ] = fs->make< TH1D >( "Jet3Pt", "p_{T}^{J3}", 250, 0., 1500.);
  histos1D_[ "Jet3Pt" ]->SetXTitle( "p_{T} Leading Jet" );
  histos1D_[ "Jet3Pt" ]->SetYTitle( "# of events" );
  histos1D_[ "Jet3Pt" ]->Sumw2();

  ////Lepton variables
  //# elecs
  histos1D_[ "nElecsIDIso" ] = fs->make< TH1D >( "nElecsIDIso", "# Elecs", 25, 0, 25);
  histos1D_[ "nElecsIDIso" ]->SetXTitle( "# Elecs" );
  histos1D_[ "nElecsIDIso" ]->SetYTitle( "# of events" );
  histos1D_[ "nElecsIDIso" ]->Sumw2();

  //Eta
  histos1D_[ "ElecEta" ] = fs->make< TH1D >( "ElecEta", "#eta", 50, -5., 5.);
  histos1D_[ "ElecEta" ]->SetXTitle( "#eta" );
  histos1D_[ "ElecEta" ]->SetYTitle( "# of events" );
  histos1D_[ "ElecEta" ]->Sumw2();

  //Phi
  histos1D_[ "ElecPhi" ] = fs->make< TH1D >( "ElecPhi", "#phi", 50, -M_PI, M_PI);
  histos1D_[ "ElecPhi" ]->SetXTitle( "#phi" );
  histos1D_[ "ElecPhi" ]->SetYTitle( "# of events" );
  histos1D_[ "ElecPhi" ]->Sumw2();


  //Pt
  histos1D_[ "ElecPt" ] = fs->make< TH1D >( "ElecPt", "p_{T}", 250, 0., 1500.);
  histos1D_[ "ElecPt" ]->SetXTitle( "p_{T}" );
  histos1D_[ "ElecPt" ]->SetYTitle( "# of events" );
  histos1D_[ "ElecPt" ]->Sumw2();

  //# muons
  histos1D_[ "nMuonsIDIso" ] = fs->make< TH1D >( "nMuonsIDIso", "# Muons", 25, 0, 25);
  histos1D_[ "nMuonsIDIso" ]->SetXTitle( "# Muons" );
  histos1D_[ "nMuonsIDIso" ]->SetYTitle( "# of events" );
  histos1D_[ "nMuonsIDIso" ]->Sumw2();

  //Eta
  histos1D_[ "MuonEta" ] = fs->make< TH1D >( "MuonEta", "#eta", 50, -5., 5.);
  histos1D_[ "MuonEta" ]->SetXTitle( "#eta" );
  histos1D_[ "MuonEta" ]->SetYTitle( "# of events" );
  histos1D_[ "MuonEta" ]->Sumw2();

  //Phi
  histos1D_[ "MuonPhi" ] = fs->make< TH1D >( "MuonPhi", "#phi", 50, -M_PI, M_PI);
  histos1D_[ "MuonPhi" ]->SetXTitle( "#phi" );
  histos1D_[ "MuonPhi" ]->SetYTitle( "# of events" );
  histos1D_[ "MuonPhi" ]->Sumw2();


  //Pt
  histos1D_[ "MuonPt" ] = fs->make< TH1D >( "MuonPt", "p_{T}", 250, 0., 1500.);
  histos1D_[ "MuonPt" ]->SetXTitle( "p_{T}" );
  histos1D_[ "MuonPt" ]->SetYTitle( "# of events" );
  histos1D_[ "MuonPt" ]->Sumw2();

  ///Event variables
  //HT
  histos1D_[ "eventHT" ] = fs->make< TH1D >( "eventHT", "H_{T}", 500, 0., 5000.);
  histos1D_[ "eventHT" ]->SetXTitle( "H_{T}" );
  histos1D_[ "eventHT" ]->SetYTitle( "# of events" );
  histos1D_[ "eventHT" ]->Sumw2();

  //MHT
  histos1D_[ "eventMHT" ] = fs->make< TH1D >( "eventMHT", "#slash H_{T}", 300, 0., 1500.);
  histos1D_[ "eventMHT" ]->SetXTitle( "#slash H_{T}" );
  histos1D_[ "eventMHT" ]->SetYTitle( "# of events" );
  histos1D_[ "eventMHT" ]->Sumw2();

  histos1D_[ "eventMHTPhi" ] = fs->make< TH1D >( "eventMHTPhi", "#phi(#slash H_{T})", 50, -M_PI, M_PI);
  histos1D_[ "eventMHTPhi" ]->SetXTitle( "#phi(#slash H_{T})" );
  histos1D_[ "eventMHTPhi" ]->SetYTitle( "# of events" );
  histos1D_[ "eventMHTPhi" ]->Sumw2();

  //MEff
  histos1D_[ "eventMEff" ] = fs->make< TH1D >( "eventMEff", "M_{eff}", 500, 0., 5000.);
  histos1D_[ "eventMEff" ]->SetXTitle( "M_{eff}" );
  histos1D_[ "eventMEff" ]->SetYTitle( "# of events" );
  histos1D_[ "eventMEff" ]->Sumw2();


  //MET
  histos1D_[ "eventMET" ] = fs->make< TH1D >( "eventMET", "#slash E_{T}", 300, 0., 1500.);
  histos1D_[ "eventMET" ]->SetXTitle( "#slash E_{T}" );
  histos1D_[ "eventMET" ]->SetYTitle( "# of events" );
  histos1D_[ "eventMET" ]->Sumw2();

  histos1D_[ "eventMETPhi" ] = fs->make< TH1D >( "eventMETPhi", "#phi(#slash E_{T})", 50, -M_PI, M_PI);
  histos1D_[ "eventMETPhi" ]->SetXTitle( "#phi(#slash E_{T})" );
  histos1D_[ "eventMETPhi" ]->SetYTitle( "# of events" );
  histos1D_[ "eventMETPhi" ]->Sumw2();


  //MTW electron
  histos1D_[ "eventMTWElec" ] = fs->make< TH1D >( "eventMTWElec", "M_{T}(W,e)", 200, 0., 200.);
  histos1D_[ "eventMTWElec" ]->SetXTitle( "M_{T}(W,e)" );
  histos1D_[ "eventMTWElec" ]->SetYTitle( "# of events" );
  histos1D_[ "eventMTWElec" ]->Sumw2();

  //MTW muon
  histos1D_[ "eventMTWMuon" ] = fs->make< TH1D >( "eventMTWMuon", "M_{T}(W,#mu)", 200, 0., 200.);
  histos1D_[ "eventMTWMuon" ]->SetXTitle( "M_{T}(W,#mu)" );
  histos1D_[ "eventMTWMuon" ]->SetYTitle( "# of events" );
  histos1D_[ "eventMTWMuon" ]->Sumw2();


  ////2D histograms
  histos2D_[ "eventHTvsMHT" ] = fs->make< TH2D >( "eventHTvsMHT", "H_{T} vs. #slash H_{T}", 500, 0., 5000., 300, 0., 1500.);
  histos2D_[ "eventHTvsMHT" ]->SetXTitle( "H_{T}" );
  histos2D_[ "eventHTvsMHT" ]->SetYTitle( "#slash H_{T}" );
  histos2D_[ "eventHTvsMHT" ]->Sumw2();

  histos2D_[ "eventHTvsMET" ] = fs->make< TH2D >( "eventHTvsMET", "H_{T} vs. #slash E_{T}", 500, 0., 5000., 300, 0., 1500.);
  histos2D_[ "eventHTvsMET" ]->SetXTitle( "H_{T}" );
  histos2D_[ "eventHTvsMET" ]->SetYTitle( "#slash E_{T}" );
  histos2D_[ "eventHTvsMET" ]->Sumw2();

  histos2D_[ "eventMHTvsMET" ] = fs->make< TH2D >( "eventMHTvsMET", "#slash H_{T} vs. #slash E_{T}", 300, 0., 1500., 300, 0., 1500.);
  histos2D_[ "eventMHTvsMET" ]->SetXTitle( "#slash H_{T}" );
  histos2D_[ "eventMHTvsMET" ]->SetYTitle( "#slash E_{T}" );
  histos2D_[ "eventMHTvsMET" ]->Sumw2();

}

// ------------ method called once each job just after ending the event loop  ------------
void WtoLNu::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void WtoLNu::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void WtoLNu::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void WtoLNu::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void WtoLNu::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void WtoLNu::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(WtoLNu);
