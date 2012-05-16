// -*- C++ -*-
//
// Package:    WtoLNu
// Class:      WtoLNu
// 
/**\class WtoLNu BTagVeto.h ZInvisibleBkgds/WtoLNu/interface/BTagVeto.h

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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//Used data formats
#include "DataFormats/PatCandidates/interface/Jet.h"


//
// class declaration
//

class BTagVeto : public edm::EDProducer {
public:
  explicit BTagVeto(const edm::ParameterSet&);
  ~BTagVeto();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run&, edm::EventSetup const&);
  virtual void endRun(edm::Run&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  
public:
  double eff_b_hp_mc(double discriminator, std::string& alg);
  double eff_c_hp_mc(double discriminator, std::string& alg);
  double eff_lf_hp(double pT, std::string& alg);
  
  double sf_b_hp(double pT, std::string& alg);
  double sf_c_hp(double pT, std::string& alg);
  double sf_lf_hp(double pT, std::string& alg);
  
  ///
  double eff_b_he_mc(double discriminator, std::string& alg);
  double eff_c_he_mc(double discriminator, std::string& alg);
  double eff_lf_he(double pT, std::string& alg);
  
  double sf_b_he(double pT, std::string& alg);
  double sf_c_he(double pT, std::string& alg);
  double sf_lf_he(double pT, std::string& alg);
  
      // ----------member data ---------------------------
public:
  bool debug_;
  std::string debugString_;
  std::string tagger_;

  double vetoDisc_;
  double signalDisc_;
  edm::InputTag jetLabel_;
  edm::InputTag vetoJetLabel_;
  edm::InputTag signalJetLabel_;

  std::map<std::string, TH1D*> histos1D_;
  std::map<std::string, TH2D*> histos2D_;
private:
  double eff_hem_raw   ;
  double eff_hem_scaled;
  double eff_hpt_raw   ;
  double eff_hpt_scaled;

};

//
// constants, enums and typedefs
//


//
// static data member definitions
//
