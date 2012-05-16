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


// system include files
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <set>
#include <math.h>
#include <utility>

//ROOT includes
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Run.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

//Used data formats
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//
// class declaration
//

class WtoLNu : public edm::EDProducer {
   public:
      explicit WtoLNu(const edm::ParameterSet&);
      ~WtoLNu();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
  bool debug_;
  
  edm::InputTag vertexLabel_;
  edm::InputTag puWeightLabel_;
  edm::InputTag jetLabel_;
  edm::InputTag htLabel_;
  edm::InputTag mhtLabel_;
  edm::InputTag metLabel_;
  edm::InputTag muonLabel_;
  edm::InputTag electronLabel_ ;

  std::map<std::string, TH1D*> histos1D_;
  std::map<std::string, TH2D*> histos2D_;
  
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//
