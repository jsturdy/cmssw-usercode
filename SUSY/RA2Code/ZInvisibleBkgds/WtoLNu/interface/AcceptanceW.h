// -*- C++ -*-
//
// Package:    WtoLNu
// Class:      WtoLNu
// 
/**\class WtoLNu AcceptanceW.h ZInvisibleBkgds/WtoLNu/interface/AcceptanceW.h

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
#include <memory>

//ROOT includes
#include <TH1.h>
#include <TH2.h>

// Framework include files
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
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

//For MC truth information
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

//
// class declaration
//

class AcceptanceW : public edm::EDProducer {
public:
  explicit AcceptanceW(const edm::ParameterSet&);
  ~AcceptanceW();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
public:
  bool debug_;
  std::string debugString_;
  
  edm::InputTag genLabel_;
  edm::InputTag genMuonLabel_;
  edm::InputTag genElectronLabel_;
  int genStatus_;

  edm::InputTag jetLabel_;
  
  edm::InputTag muonLabel_;
  double muonMinPt_       ;
  double muonMaxEta_      ;
  
  edm::InputTag electronLabel_ ;
  double electronMinPt_        ;
  double electronMaxEta_       ;

  std::map<std::string, TH1D*> histos1D_;
  std::map<std::string, TH2D*> histos2D_;
};

