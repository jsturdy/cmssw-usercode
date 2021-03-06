#ifndef JAREDSUSYANALYZER
#define JAREDSUSYANALYZER

// System include files
#include <memory>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

// ROOT includes
#include <TNtuple.h>

// Framework include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// SUSY include files
//#include "SusyAnalysis/EventSelector/interface/SusyEventSelector.h"
//#include "SusyAnalysis/EventSelector/interface/SelectorSequence.h"

//#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

//#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"

//#include "UserCode/AnalysisTools/test/ALPGENParticleId.cc"

//#include "PhysicsTools/Utilities/interface/deltaPhi.h"
//#include "PhysicsTools/Utilities/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Hemisphere.h"


//
// Class declaration
//
class JaredSusyAnalyzer : public edm::EDAnalyzer {
public:
  explicit JaredSusyAnalyzer(const edm::ParameterSet&);
  ~JaredSusyAnalyzer();
  
private:
  //*** CMSSW interface
  /// Called once per job, at start
  virtual void beginJob(const edm::EventSetup&) ;
  /// Called for each event
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  /// Called once per job, at end
  virtual void endJob();

  /// Print a summary of counts for all selectors
  //  virtual void printSummary(void);
  // Print an HLT trigger report
  virtual void printHLTreport(void); // georgia
 
  //*** Plotting
  /// Define all plots
  virtual void initPlots();
  /// Fill all plots for an event
  //  virtual void fillPlots( const edm::Event&, const SelectorDecisions& );

  //  virtual bool filter(const edm::Event& evt,const edm::EventSetup& iSetup );


private:

  bool matchJetsByCaloTowers( const pat::Jet&, const pat::Jet& );

  // Event information
  double eventWeight_;           ///< Event weight from config. file (if <0, get it from event)
  edm::InputTag weightSource_;   ///< Source for CSA07 event weight producer
  double weight_;                ///< Actual event weight (either config. or event)
  int    processId_;             ///< CSA07 generator process ID
  int run_, event_;
  double jetMaxEta_, jetMinPt_;  /// for preselection cuts on jets and to calculate HT and MHT
  bool doMCData_;                 /// switch to turn off generator level information

  // Counters
  unsigned int nrEventTotalRaw_;          ///< Raw number of events (+1 at each event)
  double nrEventTotalWeighted_;           ///< Weighted #(events)
  std::vector<float> nrEventSelected_;    ///< Selected #(events) for each module
  std::vector<float> nrEventAllButOne_;   ///< All-but-one selected #(events) for each module
  std::vector<float> nrEventCumulative_;  ///< Cumulative selected #(events) for each module
    
  // Plots
  TNtuple* ntuple_;      /// Will contain all the selector information we want to keep
  TTree * mAllData;      /// Will contain the additional di-jet specific data
  //TTree * mJetData;
  //TTree * mMETData;
  //TTree * mElectronData;
  //TTree * mMuonData;
  //TTree * mHLTData;
  //TTree * mVertexData;
  //TTree * mMuonData;
  //TTree * mMuonData;

  //TTree * mSelectorData; /// Will contain the information on the selector decisions
  //TTree* decisionTree_;  ///< Will contain all the decisions for ALL processed events
  //TTree* selectionTree_; ///< Will contain all the information we want to keep
  float* variables_;     ///< Container for the tree variables (from selectors)
  bool*  decisions_;     ///< Container for all selector decisions
  bool   globalDec_;     ///< Global decision for event

  int mTempTreeRun;
  int mTempTreeEvent;


  bool mTempTreeHLT1JET;
  bool mTempTreeHLT2JET;
  bool mTempTreeHLT1MET1HT;
  bool mTempTreeHLT1Muon;

  int nHLT; int HLTArray[200];
  int mTempTreenVtx;
  double mTempTreeVtxChi2[5];
  double mTempTreeVtxNdof[5];
  double mTempTreeVtxNormalizedChi2[5];
  double mTempTreeVtxX[5]; double mTempTreeVtxY[5]; double mTempTreeVtxZ[5];
  double mTempTreeVtxdX[5]; double mTempTreeVtxdY[5]; double mTempTreeVtxdZ[5];

  int    mTempTreeNjets;
  double mTempTreeHt;
  double mTempTreeMHt;
  double mTempTreeJetsEt[50];
  double mTempTreeJetsPt[50];
  double mTempTreeJetsPx[50];
  double mTempTreeJetsPy[50];
  double mTempTreeJetsPz[50];
  double mTempTreeJetsE[50];
  double mTempTreeJetsEta[50];
  double mTempTreeJetsPhi[50];
  double mTempTreeJetsFem[50];
  int    mTempTreeJetPartonFlavour[50];
  int    mTempTreeJetsHemi[50];

  // track info:
  int    mTempTreeJetTrackNo[50];
  double mTempTreeJetTrackPhi[50];
  double mTempTreeJetTrackPhiWeighted[50];
  double mTempTreeJetTrackPt[50];

  double mTempTreeJetMCCorrFactor[50];
  double mTempTreeJetJPTCorrFactor[50];

  double mTempTreeccJetMCCorrFactor[50];
  double mTempTreeccJetJPTCorrFactor[50];

  bool    mTempTreeccJetAssoc[50];
  double  mTempTreeccJetAssoc_E[50];
  double  mTempTreeccJetAssoc_px[50];
  double  mTempTreeccJetAssoc_py[50];
  double  mTempTreeccJetAssoc_pz[50];

  int    mTempTreeJetPartonId[50];
  int    mTempTreeJetPartonMother[50];
  double mTempTreeJetPartonPx[50];
  double mTempTreeJetPartonPy[50];
  double mTempTreeJetPartonPz[50];
  double mTempTreeJetPartonEt[50];
  double mTempTreeJetPartonEnergy[50];
  double mTempTreeJetPartonPhi[50];
  double mTempTreeJetPartonEta[50];

  float mTempTreeJetsBTag_TkCountHighEff[50];
  float mTempTreeJetsBTag_SimpleSecVtx[50];
  float mTempTreeJetsBTag_CombSecVtx[50];

  double mTempTreeGenJetsEt[50];
  double mTempTreeGenJetsPt[50];
  double mTempTreeGenJetsE[50];
  double mTempTreeGenJetsPx[50];
  double mTempTreeGenJetsPy[50];
  double mTempTreeGenJetsPz[50];
  double mTempTreeGenJetsEta[50];
  double mTempTreeGenJetsPhi[50];
  double mTempTreeGenHt;
  double mTempTreeGenMHt;

  int    mTempTreeNelec;
  double mTempTreeElecEt[50];
  double mTempTreeElecPt[50];
  double mTempTreeElecPx[50];
  double mTempTreeElecPy[50];
  double mTempTreeElecPz[50];
  double mTempTreeElecE[50];
  double mTempTreeElecEta[50];
  double mTempTreeElecPhi[50];
  double mTempTreeElecTrkIso[50];
  double mTempTreeElecECalIso[50];
  double mTempTreeElecHCalIso[50];
  double mTempTreeElecAllIso[50];
  double mTempTreeElecTrkChiNorm[50];
  double mTempTreeElecCharge[50];

  double mTempTreeElecIdLoose[50];
  double mTempTreeElecIdTight[50];
  double mTempTreeElecIdRobLoose[50];
  double mTempTreeElecIdRobTight[50];
  double mTempTreeElecChargeMode[50];
  double mTempTreeElecPtTrkMode[50];
  double mTempTreeElecQOverPErrTrkMode[50];
  double mTempTreeGenElecPdgId[50];
  double mTempTreeGenElecMother[50];
  double mTempTreeGenElecPx[50];
  double mTempTreeGenElecPy[50];
  double mTempTreeGenElecPz[50];
  double mTempTreeElecCaloEnergy[50];
  double mTempTreeElecHOverE[50];
  double mTempTreeElecVx[50];
  double mTempTreeElecVy[50];
  double mTempTreeElecVz[50];
  double mTempTreeElecD0[50];
  double mTempTreeElecDz[50];
  double mTempTreeElecPtTrk[50];
  double mTempTreeElecQOverPErrTrk[50];
  double mTempTreeElecLostHits[50];
  double mTempTreeElecValidHits[50];
  double mTempTreeElecNCluster[50];
  double mTempTreeElecEtaTrk[50];
  double mTempTreeElecPhiTrk[50];
  double mTempTreeElecWidthClusterEta[50];
  double mTempTreeElecWidthClusterPhi[50];
  double mTempTreeElecPinTrk[50];
  double mTempTreeElecPoutTrk[50];
  double mTempTreeElecNormChi2[50];
  bool mTempTreeccElecAssoc[50];

  double mTempTreeElecECalIsoDeposit[50];
  double mTempTreeElecHCalIsoDeposit[50];

  int    mTempTreeNmuon;
  double mTempTreeMuonEt[50];
  double mTempTreeMuonPt[50];
  double mTempTreeMuonPx[50];
  double mTempTreeMuonPy[50];
  double mTempTreeMuonPz[50];
  double mTempTreeMuonE[50];
  double mTempTreeMuonEta[50];
  double mTempTreeMuonPhi[50];
  double mTempTreeMuonTrkIso[50];
  double mTempTreeMuonECalIso[50];
  double mTempTreeMuonHCalIso[50];
  double mTempTreeMuonAllIso[50];
  double mTempTreeMuonTrkChiNorm[50];
  double mTempTreeMuonCharge[50];
  bool mTempTreeMuonIsGlobal[50];
  bool mTempTreeMuonIsStandAlone[50];
  bool mTempTreeMuonIsTracker[50]; 
  bool mTempTreeMuonIsGlobalTight[50];
  bool mTempTreeMuonIsTMLastStationLoose[50];
  bool mTempTreeMuonTMLastStationTight[50];
  bool mTempTreeMuonTM2DCompatibilityLoose[50];
  bool mTempTreeMuonTM2DCompatibilityTight[50];
  bool mTempTreeccMuonAssoc[50];

  double mTempTreeMuonECalIsoDeposit[50];
  double mTempTreeMuonHCalIsoDeposit[50];
  
  double mTempTreeMuonCombChi2[50];
  double mTempTreeMuonCombNdof[50];
  double mTempTreeMuonTrkD0[50];
  
  double  mTempTreeMuonId[50];
  double mTempTreeMuonCombVx[50];
  double mTempTreeMuonCombVy[50];
  double mTempTreeMuonCombVz[50];
  double mTempTreeMuonCombD0[50];
  double mTempTreeMuonCombDz[50];

  double mTempTreeMuonStandValidHits[50];
  double mTempTreeMuonStandLostHits[50];
  double mTempTreeMuonStandPt[50];
  double mTempTreeMuonStandPz[50];
  double mTempTreeMuonStandP[50];
  double mTempTreeMuonStandEta[50];
  double mTempTreeMuonStandPhi[50];
  double mTempTreeMuonStandChi[50];
  double mTempTreeMuonStandCharge[50];
  double mTempTreeMuonStandQOverPError[50];

  double mTempTreeMuonTrkValidHits[50];
  double mTempTreeMuonTrkLostHits[50];
  double mTempTreeMuonTrkPt[50];
  double mTempTreeMuonTrkPz[50];
  double mTempTreeMuonTrkP[50];
  double mTempTreeMuonTrkEta[50];
  double mTempTreeMuonTrkPhi[50];
  double mTempTreeMuonTrkChi[50];
  double mTempTreeMuonTrkCharge[50];
  double mTempTreeMuonTrkQOverPError[50];
  double mTempTreeMuonTrkOuterZ[50];
  double mTempTreeMuonTrkOuterR[50];

  double mTempTreeGenMuonPdgId[50];
  double mTempTreeGenMuonMother[50];
  double mTempTreeGenMuonPx[50];
  double mTempTreeGenMuonPy[50];
  double mTempTreeGenMuonPz[50];

  int mTempAlpIdTest;
  double mTempAlpPtScale;

  double mTempMuonPairMass;
  int mTempMuonPairIndex[2];

  //PF jets
  int    mTempTreeNPFjet;
  double mTempTreePFHt;
  double mTempTreePFMHt;
  double mTempTreePFjetEta[50];
  double mTempTreePFjetPhi[50];
  double mTempTreePFjetE[50];
  double mTempTreePFjetPx[50];
  double mTempTreePFjetPy[50];
  double mTempTreePFjetPz[50];
  double mTempTreePFjetPt[50];
  double mTempTreePFjetCharge[50];
  double mTempTreePFjetFem[50];

  // Generated MET
  double mTempTreepfMET_Gen[3];
  double mTempTreetcMET_Gen[3];
  double mTempTreeMET_Gen[3];

  int nFullMET;
 
  // Do the MET save for non cc pfMET
  double mTempTreepfMET_Fullcorr_nocc[3];
  double mTempTreepfMETphi_Fullcorr_nocc;
  double mTempTreepfMET_Fullcorr_nocc_significance;
  double mTempTreepfMET_Nocorr_nocc[2];
  double mTempTreepfMETphi_Nocorr_nocc;
  double mTempTreepfMET_Muoncorr_nocc[2];
  double mTempTreepfMETphi_Muoncorr_nocc;
  double mTempTreepfMET_JECcorr_nocc[2];
  double mTempTreepfMETphi_JECcorr_nocc;

  // Do the MET save for non cc tcMET
  double mTempTreetcMET_Fullcorr_nocc[3];
  double mTempTreetcMETphi_Fullcorr_nocc;
  double mTempTreetcMET_Fullcorr_nocc_significance;
  double mTempTreetcMET_Nocorr_nocc[2];
  double mTempTreetcMETphi_Nocorr_nocc;
  double mTempTreetcMET_Muoncorr_nocc[2];
  double mTempTreetcMETphi_Muoncorr_nocc;
  double mTempTreetcMET_JECcorr_nocc[2];
  double mTempTreetcMETphi_JECcorr_nocc;

  // Do the MET save for non cc calo MET
  double mTempTreeMET_Fullcorr_nocc[3];
  double mTempTreeMETphi_Fullcorr_nocc;
  double mTempTreeMET_Fullcorr_nocc_significance;
  double mTempTreeMET_Nocorr_nocc[2];
  double mTempTreeMETphi_Nocorr_nocc;
  double mTempTreeMET_Muoncorr_nocc[2];
  double mTempTreeMETphi_Muoncorr_nocc;
  double mTempTreeMET_JECcorr_nocc[2];
  double mTempTreeMETphi_JECcorr_nocc;

  // Do the MET save for cc MET
  int nUncorrMET;
  double mTempTreeMET_Fullcorr_cc[3];
  double mTempTreeMETphi_Fullcorr_cc;
  double mTempTreeMET_Nocorr_cc[2];
  double mTempTreeMET_Muoncorr_cc[2];
  double mTempTreeMET_JECcorr_cc[2];

  int mTempTreeNhemispheres;
  double mTempTreeHemispheresEt[2];
  double mTempTreeHemispheresPt[2];
  double mTempTreeHemispheresPx[2];
  double mTempTreeHemispheresPy[2];
  double mTempTreeHemispheresPz[2];
  double mTempTreeHemispheresE[2];
  double mTempTreeHemispheresEta[2];
  double mTempTreeHemispheresPhi[2]; 
  double mTempTreeHemispheresdPhi;

  double mTempTreeMPTPhi;
  double mTempTreeMPTPx;
  double mTempTreeMPTPy;
  double mTempTreeMPTPz;

  int length;
  int ids[1000];
  int refs[1000];
  float genE[1000];
  float genPx[1000];
  float genPy[1000];
  float genPz[1000];
  int genStatus[1000];

  int genLepLength;
  int genLepIds[100];
  int genLepRefs[100];
  float genLepE[100];
  float genLepPx[100];
  float genLepPy[100];
  float genLepPz[100];
  int genLepStatus[100];

  double mTempTreeEventWeight;
  int    mTempTreeProcID;
  double mTempTreePthat;
  int    mGlobalDecision;

  // Data tags
  edm::InputTag triggerResults_; 
  std::vector<std::string> pathNames_;

  edm::TriggerNames triggerNames_;     // TriggerNames class

  unsigned int  nEvents_;              // number of events processed

  unsigned int  nWasRun_;              // # where at least one HLT was run
  unsigned int  nAccept_;              // # of accepted events
  unsigned int  nErrors_;              // # where at least one HLT had error
  std::vector<unsigned int> hlWasRun_; // # where HLT[i] was run
  std::vector<unsigned int> hlAccept_; // # of events accepted by HLT[i]
  std::vector<unsigned int> hlErrors_; // # of events with error in HLT[i]
  bool init_;                          // vectors initialised or not

  edm::InputTag vtxTag_;

  //met tags
  edm::InputTag tcmetTag_;
  edm::InputTag pfmetTag_;
  edm::InputTag metTag_;
  //edm::InputTag mhtTag_;

  edm::InputTag elecTag_;
  edm::InputTag muonTag_;
  edm::InputTag pfelecTag_;
  edm::InputTag pfmuonTag_;
  edm::InputTag genTag_;

  //jet tags
  bool usePfjets_;
  edm::InputTag pfjetTag_;
  edm::InputTag jetTag_;
  edm::InputTag jptTag_;

  std::string outputFileName_;

  //input from .cfg
  bool theSoup;
  double fileWeight;

  double localPi;
  unsigned int *mSelectorResults;

};























#endif















