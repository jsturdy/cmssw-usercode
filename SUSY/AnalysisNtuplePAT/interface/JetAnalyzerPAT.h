#ifndef JETANALYZERPAT
#define JETANALYZERPAT

// System include files
#include <memory>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <set>
#include <utility>

// ROOT includes
#include <TNtuple.h>

// Framework include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// SUSY include files
//#include "SusyAnalysis/EventSelector/interface/SusyEventSelector.h"
//#include "SusyAnalysis/EventSelector/interface/SelectorSequence.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
//Included for cross cleaning
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

//#include "PhysicsTools/Utilities/interface/deltaPhi.h"
//#include "PhysicsTools/Utilities/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Hemisphere.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"


#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

//
// Class declaration
//

class JetAnalyzerPAT {
 public:
  JetAnalyzerPAT(const edm::ParameterSet&, TTree*);
  ~JetAnalyzerPAT();
  
  bool filter(const edm::Event&, const edm::EventSetup&);
  
  void beginRun(const edm::Run& run, const edm::EventSetup& es);
  //*** Plotting
  /// Define all plots
  void bookTTree();

  typedef struct {
    double Ht;
    double MHx;
    double MHy;
    double MHz;
    double MHE;
    double MHt;
  } MYMHT;

  typedef struct {
    bool   JetIDMinimal;
    bool   JetIDLoose;
    bool   JetIDTight;
  
    double JetFem;
    double JetFhad;
    double JetCharge;
    int    JetHemi;
    int    JetNConst;
    //calo/jpt jet specific
    double JetfHPD;
    double JetfRBX;
    double JetN90;
  
    //jpt/pf jet specific
    double JetChargedFem;
    double JetNeutralFem;
    double JetChargedFhad;
    double JetNeutralFhad;
    double JetChargedMult;
    double JetNeutralMult;
    double JetElecMult;
    double JetMuonMult;
  
    //pf jet specific
    double JetChargedFmu;
    double JetChargedFele;
    double JetChargedFpho;
    double JetHFFem;
    double JetHFFhad;
    double JetChargedHadMult;
    double JetNeutralHadMult;
    double JetHFHadMult;
    double JetHFEMMult;
    double JetHFMult;
    double JetPhotonMult;
  } MYJETID;

  typedef struct {
    double JetBTag_TCHE;
    double JetBTag_TCHP;
    double JetBTag_jetProb;
    double JetBTag_jetBProb;
    double JetBTag_SSVHE;
    double JetBTag_SSVHP;
    double JetBTag_CSV;
    double JetBTag_CSVMVA;
    double JetBTag_SoftLepton;
    double JetBTag_SoftLeptonByIP;
    double JetBTag_SoftLeptonByPt;
  } BTAGINFO;


 private:

  // Configuration parameters
  TString prefix_;

  int    minNJets_;                  /// for preselection cuts on jets and to calculate HT and MHT
  double jetMaxEta_;
  double jetMinPt_;

  double htMaxEta_;
  double htMinPt_;

  std::vector<double> selJetMaxEta_;
  std::vector<double> selJetMinPt_;
  std::vector<double> selJetMaxEMF_;
  std::vector<double> selJetMinEMF_;

  int  debug_;
  bool doMCData_;
  bool usePFJets_;
  bool useJPTJets_;
  bool useCaloJets_;
  bool useTrackJets_;

  double electronPt_  ;
  double electronIso_ ;
  double tauPt_       ;
  double tauIso_      ;
  double muonPt_      ;
  double muonIso_     ;
  double photonPt_    ;
  double photonIso_   ;
  

  JetCorrectionUncertainty *jecUnc;
  
  // Data tags
  edm::InputTag jetTag_;
  edm::InputTag genJetTag_;

  std::string   jetCorTag_;

  char logmessage[128];

  // Plots
  TTree * mJetData;      /// Will contain the data passing the jet selection

  std::map<std::string, std::vector<float> > map_s_vf_correctionFactor;
  //std::map<std::string, std::vector<int> >   map_s_vi_JetOverlaps;
  //std::map<std::string, std::vector<int> >   map_s_vi_JetNOverlaps;
  std::vector<int>                           vi_JetElectronOverlaps;
  std::vector<int>                           vi_JetElectronNOverlaps;
  std::vector<int>                           vi_JetMuonOverlaps;
  std::vector<int>                           vi_JetMuonNOverlaps;
  std::vector<int>                           vi_JetTauOverlaps;
  std::vector<int>                           vi_JetTauNOverlaps;
  std::vector<int>                           vi_JetPhotonOverlaps;
  std::vector<int>                           vi_JetPhotonNOverlaps;

  std::vector<reco::Candidate::LorentzVector > v_JetP4;
  std::vector<reco::Candidate::LorentzVector > v_RawJetP4;
  std::vector<reco::Candidate::LorentzVector > v_GenJetP4;

  reco::Candidate::LorentzVector MHtP4;
  reco::Candidate::LorentzVector GenMHtP4;

  std::vector<float> vf_JECUncPlus;
  std::vector<float> vf_JECUncMinus;

  int    i_NJets;

  MYMHT JetMHt;
  double d_Ht;

  MYMHT GenMHt;
  double d_GenHt;

  std::vector<double> vd_JetEtaEtaMoment;
  std::vector<double> vd_JetEtaPhiMoment;
  std::vector<double> vd_JetPhiPhiMoment;

  //JetID variables
  pat::strbitset retmin;
  pat::strbitset retloo;
  pat::strbitset rettig;

  std::vector<int>   vb_JetIDMinimal;
  std::vector<int>   vb_JetIDLoose;
  std::vector<int>   vb_JetIDTight;
  
  std::vector<double> vd_JetEFrac_em;
  std::vector<double> vd_JetEFrac_had;
  std::vector<double> vd_JetCharge;
  std::vector<int>    vi_JetHemi;
  std::vector<int>    vi_JetNConst;

  //calo/jpt jet specific
  std::vector<double> vd_JetfHPD;
  std::vector<double> vd_JetfRBX;
  std::vector<double> vd_JetN90;
  
  //jpt/pf jet specific
  std::vector<double> vd_JetChargedFrac_em;
  std::vector<double> vd_JetNeutralFrac_em;
  std::vector<double> vd_JetChargedFrac_had;
  std::vector<double> vd_JetNeutralFrac_had;

  std::vector<double> vd_JetChargedEn_em;
  std::vector<double> vd_JetNeutralEn_em;
  std::vector<double> vd_JetChargedEn_had;
  std::vector<double> vd_JetNeutralEn_had;

  std::vector<double> vd_JetChargedMult;
  std::vector<double> vd_JetElectronMult;
  std::vector<double> vd_JetMuonMult;
  
  //pf jet specific
  std::vector<double> vd_JetHFMult_had;
  std::vector<double> vd_JetHFMult_em;
  std::vector<double> vd_JetNeutralMult_had;
  std::vector<double> vd_JetChargedMult_had;
  std::vector<double> vd_JetPhotonMult;
  std::vector<double> vd_JetNeutralMult;

  std::vector<double> vd_JetHFFrac_em;
  std::vector<double> vd_JetHFFrac_had;
  std::vector<double> vd_JetEFrac_muon;
  std::vector<double> vd_JetChargedFrac_muon;
  std::vector<double> vd_JetEFrac_electron;
  std::vector<double> vd_JetEFrac_photon;

  std::vector<double> vd_JetHFEn_em;
  std::vector<double> vd_JetHFEn_had;
  std::vector<double> vd_JetEn_muon;
  std::vector<double> vd_JetChargedEn_muon;
  std::vector<double> vd_JetEn_electron;
  std::vector<double> vd_JetEn_photon;

  //std::vector<MYJETID> vjid_JetID;
  
  // track info:
  std::vector<int>    vi_JetTrackNo;
  std::vector<double> vd_JetTrackPhi;
  std::vector<double> vd_JetTrackPhiWeighted;
  std::vector<double> vd_JetTrackPt;
  
  //calo jet corrections
  std::vector<double> vd_JetMCCorrFactor;
  std::vector<double> vd_JetJPTCorrFactor;

  //b-tagging
  std::vector<BTAGINFO> vbtag_JetBtag;
  std::vector<double> vd_JetBTag_TCHE;
  std::vector<double> vd_JetBTag_TCHP;
  std::vector<double> vd_JetBTag_jetProb;
  std::vector<double> vd_JetBTag_jetBProb;
  std::vector<double> vd_JetBTag_SSVHE;
  std::vector<double> vd_JetBTag_SSVHP;
  std::vector<double> vd_JetBTag_CSV;
  std::vector<double> vd_JetBTag_CSVMVA;
  std::vector<double> vd_JetBTag_SoftLepton;
  std::vector<double> vd_JetBTag_SoftLeptonByIP;
  std::vector<double> vd_JetBTag_SoftLeptonByPt;
  
  //Parton level id for gen jets?
  std::vector<reco::Candidate::LorentzVector > v_JetPartonP4;
  std::vector<int>    vi_JetPartonId;
  std::vector<int>    vi_JetPartonStatus;
  std::vector<int>    vi_JetPartonMother;
  std::vector<int>    vi_JetPartonMotherStatus;
  std::vector<int>    vi_JetPartonFlavour;

  //Result of some predefined preselection
  bool   bool_JetPreselection;

 public:
  void maintenance() {
    //Setup the vectors
    map_s_vf_correctionFactor.clear();

    //map_s_vi_JetNOverlaps.clear();
    //map_s_vi_JetOverlaps.clear();

    vi_JetElectronNOverlaps.clear();
    vi_JetElectronOverlaps.clear();
    vi_JetMuonNOverlaps.clear();
    vi_JetMuonOverlaps.clear();
    vi_JetTauNOverlaps.clear();
    vi_JetTauOverlaps.clear();
    vi_JetPhotonNOverlaps.clear();
    vi_JetPhotonOverlaps.clear();

    vf_JECUncPlus.clear();
    vf_JECUncMinus.clear();
    
    v_JetP4.clear();
    v_RawJetP4.clear();
    v_GenJetP4.clear();
    
    vd_JetEtaEtaMoment.clear();
    vd_JetEtaPhiMoment.clear();
    vd_JetPhiPhiMoment.clear();
    
    vb_JetIDMinimal.clear();
    vb_JetIDLoose.clear();
    vb_JetIDTight.clear();
    
    vd_JetEFrac_em.clear();
    vd_JetEFrac_had.clear();
    vd_JetCharge.clear();
    vi_JetHemi.clear();
    vi_JetNConst.clear();
    vd_JetfHPD.clear();
    vd_JetfRBX.clear();
    vd_JetN90.clear();
    
    vd_JetChargedFrac_em.clear();
    vd_JetNeutralFrac_em.clear();
    vd_JetChargedFrac_had.clear();
    vd_JetNeutralFrac_had.clear();

    vd_JetChargedEn_em.clear();
    vd_JetNeutralEn_em.clear();
    vd_JetChargedEn_had.clear();
    vd_JetNeutralEn_had.clear();

    vd_JetChargedMult.clear();
    vd_JetElectronMult.clear();
    vd_JetMuonMult.clear();
    
    vd_JetHFMult_had.clear();
    vd_JetHFMult_em.clear();
    vd_JetNeutralMult_had.clear();
    vd_JetChargedMult_had.clear();
    vd_JetPhotonMult.clear();
    vd_JetNeutralMult.clear();
    
    vd_JetHFFrac_em.clear();
    vd_JetHFFrac_had.clear();
    vd_JetEFrac_muon.clear();
    vd_JetChargedFrac_muon.clear();
    vd_JetEFrac_electron.clear();
    vd_JetEFrac_photon.clear();

    vd_JetHFEn_em.clear();
    vd_JetHFEn_had.clear();
    vd_JetEn_muon.clear();
    vd_JetChargedEn_muon.clear();
    vd_JetEn_electron.clear();
    vd_JetEn_photon.clear();

    vi_JetTrackNo.clear();
    vd_JetTrackPhi.clear();
    vd_JetTrackPhiWeighted.clear();
    vd_JetTrackPt.clear();
    
    vd_JetMCCorrFactor.clear();
    vd_JetJPTCorrFactor.clear();
    
    vd_JetBTag_TCHE.clear();
    vd_JetBTag_TCHP.clear();
    vd_JetBTag_jetProb.clear();
    vd_JetBTag_jetBProb.clear();
    vd_JetBTag_SSVHE.clear();
    vd_JetBTag_SSVHP.clear();
    vd_JetBTag_CSV.clear();
    vd_JetBTag_CSVMVA.clear();
    vd_JetBTag_SoftLepton.clear();
    vd_JetBTag_SoftLeptonByIP.clear();
    vd_JetBTag_SoftLeptonByPt.clear();
    
    v_JetPartonP4.clear();
    vi_JetPartonId.clear();
    vi_JetPartonStatus.clear();
    vi_JetPartonMother.clear();
    vi_JetPartonMotherStatus.clear();
    vi_JetPartonFlavour.clear();

  }
};
#endif
