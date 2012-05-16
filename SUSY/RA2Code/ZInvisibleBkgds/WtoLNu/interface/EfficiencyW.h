// -*- C++ -*-
//
// Package:    WtoLNu
// Class:      WtoLNu
// 
/**\class WtoLNu EfficiencyW.h ZInvisibleBkgds/WtoLNu/interface/EfficiencyW.h

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
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


//
// class declaration
//

class EfficiencyW : public edm::EDProducer {
   public:
      explicit EfficiencyW(const edm::ParameterSet&);
      ~EfficiencyW();

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
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//
