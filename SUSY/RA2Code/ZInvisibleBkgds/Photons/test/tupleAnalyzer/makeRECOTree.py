import sys,os
import ROOT as r
from array import array
import math
import string
import re
#import argparse
import optparse

import itertools

import analysisCuts as cutF
def main() :

    parser = optparse.OptionParser(description="Switch for data/MC running")
    parser.add_option('-m', action="store_true", default=False, dest="isMC")
    parser.add_option('-d', action="store_true", default=False, dest="debug")
    parser.add_option('-s', action="store",      default="gjets400v1",dest="sample", type="string")
    parser.add_option('-f', action="store",      default="0",       dest="subsec", type="int")
    parser.add_option('-t', action="store",      default="analysisIDPFIso",dest="treeName", type="string")
    parser.add_option('-p', action="store",      default=100.0,dest="minpt",    type="float")
    parser.add_option('-r', action="store",      default=0.0,dest="cutDR",      type="float")
    options, args = parser.parse_args()

    r.gROOT.SetBatch(True)
    
    subfiles = [
        "*",#0
        "*_?_?_???",#1
        "*_1?_?_???","*_2?_?_???","*_3?_?_???","*_4?_?_???","*_5?_?_???","*_6?_?_???","*_7?_?_???","*_8?_?_???","*_9?_?_???",#10
        "*_11?_?_???","*_12?_?_???","*_13?_?_???","*_14?_?_???","*_15?_?_???","*_16?_?_???","*_17?_?_???","*_18?_?_???","*_19?_?_???",#19
        "*_21?_?_???","*_22?_?_???","*_23?_?_???","*_24?_?_???","*_25?_?_???","*_26?_?_???","*_27?_?_???","*_28?_?_???","*_29?_?_???",#28
        "*_31?_?_???","*_32?_?_???","*_33?_?_???","*_34?_?_???","*_35?_?_???","*_36?_?_???","*_37?_?_???","*_38?_?_???","*_39?_?_???",#37
        "*_41?_?_???","*_42?_?_???","*_43?_?_???","*_44?_?_???","*_45?_?_???","*_46?_?_???","*_47?_?_???","*_48?_?_???","*_49?_?_???",#46
        "*_51?_?_???","*_52?_?_???","*_53?_?_???","*_54?_?_???","*_55?_?_???","*_56?_?_???","*_57?_?_???","*_58?_?_???","*_59?_?_???",#55
        "*_61?_?_???","*_62?_?_???","*_63?_?_???","*_64?_?_???","*_65?_?_???","*_66?_?_???","*_67?_?_???","*_68?_?_???","*_69?_?_???",#64
        "*_71?_?_???","*_72?_?_???","*_73?_?_???","*_74?_?_???","*_75?_?_???","*_76?_?_???","*_77?_?_???","*_78?_?_???","*_79?_?_???",#73
        "*_81?_?_???","*_82?_?_???","*_83?_?_???","*_84?_?_???","*_85?_?_???","*_86?_?_???","*_87?_?_???","*_88?_?_???","*_89?_?_???",#82
        "*_91?_?_???","*_92?_?_???","*_93?_?_???","*_94?_?_???","*_95?_?_???","*_96?_?_???","*_97?_?_???","*_98?_?_???","*_99?_?_???",#91
        "*_101?_?_???","*_102?_?_???","*_103?_?_???","*_104?_?_???","*_105?_?_???","*_106?_?_???","*_107?_?_???","*_108?_?_???","*_109?_?_???",#100
        ]

    outFileName = "recoTreeDR%2.1f_%s_%s_%d.root"%(options.cutDR,options.sample,options.treeName,options.subsec)
    if options.debug:
        outFileName = "recoTreeDR%2.1f_%s_%s_test.root"%(options.cutDR,options.treeName,options.sample)

    print outFileName
    sys.stdout.flush()
    outFile  = r.TFile(outFileName,"RECREATE")
    tree     = r.TTree( 'reco', 'tree for reco info ' )
    
    ra2_run     = array( 'i', [ 0 ] )
    ra2_event   = array( 'i', [ 0 ] )
    ra2_lumi    = array( 'i', [ 0 ] )
    nVtx        = array( 'i', [ 0 ] )
    nJetsHT     = array( 'i', [ 0 ] )
    nJetsHTMInv = array( 'i', [ 0 ] )
    nJetsMHT    = array( 'i', [ 0 ] )

    nPhotonsID      = array( 'i', [ 0 ] )
    #nPhotonsLoose   = array( 'i', [ 0 ] )
    nPhotonsTight   = array( 'i', [ 0 ] )
    #nPhotonsCombR03 = array( 'i', [ 0 ] )
    #nPhotonsCombR04 = array( 'i', [ 0 ] )

    htVal       = array( 'd', [ 0. ] )
    htMInvVal   = array( 'd', [ 0. ] )
    mhtVal      = array( 'd', [ 0. ] )
    dphi1       = array( 'd', [ 0. ] )
    dphi2       = array( 'd', [ 0. ] )
    dphi3       = array( 'd', [ 0. ] )
    dphiMin     = array( 'd', [ 0. ] )
    #bosonPt     = array( 'd', [ 0. ] )
    
#    photonIsTightID   = array( 'b', [ 0 ] )
    #photonIsLooseIso  = array( 'b', [ 0 ] )
    photonIsTightIso  = array( 'b', [ 0 ] )
    #photonIsCombIso   = array( 'b', [ 0 ] )
    photonPt  = array( 'd', [ 0. ] )
    photonEta = array( 'd', [ 0. ] )
    photonPhi = array( 'd', [ 0. ] )
    photonMinDR = array( 'd', [ 0. ] )
    photonJet1DR= array( 'd', [ 0. ] )
    photonpfCH = array( 'd', [ 0. ] )
    photonpfNU = array( 'd', [ 0. ] )
    photonpfGA = array( 'd', [ 0. ] )
    photonSieie  = array( 'd', [ 0. ] )
    photonHoverE = array( 'd', [ 0. ] )
    photonEConvVeto    = array( 'b', [ 0 ] )
    photonPixelVeto    = array( 'b', [ 0 ] )

    muon1Pt     = array( 'd', [ 0. ] )
    muon1Eta    = array( 'd', [ 0. ] )
    muon1MinDR  = array( 'd', [ 0. ] )
    muon1Jet1DR = array( 'd', [ 0. ] )
    muon2Pt     = array( 'd', [ 0. ] )
    muon2Eta    = array( 'd', [ 0. ] )
    muon2MinDR  = array( 'd', [ 0. ] )
    muon2Jet1DR = array( 'd', [ 0. ] )
    dimuonPt    = array( 'd', [ 0. ] )
    dimuonEta   = array( 'd', [ 0. ] )
    dimuonMinDR = array( 'd', [ 0. ] )
    dimuonJet1DR= array( 'd', [ 0. ] )
    dimuonM     = array( 'd', [ 0. ] )

    jet1Pt    = array( 'd', [ 0. ] )
    jet1Eta   = array( 'd', [ 0. ] )
    jet2Pt    = array( 'd', [ 0. ] )
    jet2Eta   = array( 'd', [ 0. ] )
    jet3Pt    = array( 'd', [ 0. ] )
    jet3Eta   = array( 'd', [ 0. ] )
    jet4Pt    = array( 'd', [ 0. ] )
    jet4Eta   = array( 'd', [ 0. ] )

    eventWt   = array( 'd', [ 0. ] )
    puWt      = array( 'd', [ 0. ] )

    passRA2ElVeto    = array( 'b', [ 0 ] )
    passRA2MuVeto    = array( 'b', [ 0 ] )
    passDirIsoElVeto = array( 'b', [ 0 ] )
    passDirIsoMuVeto = array( 'b', [ 0 ] )
    passIsoTrkVeto   = array( 'b', [ 0 ] )
    ##################


    tree.Branch( 'ra2_run',    ra2_run,     'ra2_run/I' )
    tree.Branch( 'ra2_event',  ra2_event,   'ra2_event/I' )
    tree.Branch( 'ra2_lumi',   ra2_lumi,    'ra2_lumi/I' )
    tree.Branch( 'nVtx',       nVtx,        'nVtx/I' )
    tree.Branch( 'nJetsHT',    nJetsHT,     'nJetsHT/I' )
    tree.Branch( 'nJetsHTMInv',nJetsHTMInv, 'nJetsHTMInv/I' )
    tree.Branch( 'nJetsMHT',   nJetsMHT,    'nJetsMHT/I' )

    tree.Branch( 'nPhotonsID',      nPhotonsID,      'nPhotonsID/I' )
#    tree.Branch( 'nPhotonsLoose',   nPhotonsLoose,   'nPhotonsLoose/I' )
    tree.Branch( 'nPhotonsTight',   nPhotonsTight,   'nPhotonsTight/I' )
#    tree.Branch( 'nPhotonsCombR03', nPhotonsCombR03, 'nPhotonsCombR03/I' )
#    tree.Branch( 'nPhotonsCombR04', nPhotonsCombR04, 'nPhotonsCombR04/I' )

    tree.Branch( 'htVal',      htVal,       'htVal/D' )
    tree.Branch( 'htMInvVal',  htMInvVal,   'htMInvVal/D' )
    tree.Branch( 'mhtVal',     mhtVal,      'mhtVal/D' )
    tree.Branch( 'dphi1',      dphi1,       'dphi1/D' )
    tree.Branch( 'dphi2',      dphi2,       'dphi2/D' )
    tree.Branch( 'dphi3',      dphi3,       'dphi3/D' )
    tree.Branch( 'dphiMin',    dphiMin,     'dphiMin/D' )
    #tree.Branch( 'bosonPt',    bosonPt,   'bosonPt/D' )

#    tree.Branch( 'photonIsTightID',   photonIsTightID,    'photonIsTightID/O' )
#    tree.Branch( 'photonIsLooseIso',  photonIsLooseIso,   'photonIsLooseIso/O' )
    tree.Branch( 'photonIsTightIso',  photonIsTightIso,   'photonIsTightIso/O' )
#    tree.Branch( 'photonIsCombIso',   photonIsCombIso,    'photonIsCombIso/O' )
    tree.Branch( 'photonPt',   photonPt,    'photonPt/D' )
    tree.Branch( 'photonEta',  photonEta,   'photonEta/D' )
    tree.Branch( 'photonPhi',  photonPhi,   'photonPhi/D' )
    tree.Branch( 'photonMinDR',photonMinDR, 'photonMinDR/D' )
    tree.Branch( 'photonJet1DR',photonJet1DR, 'photonJet1DR/D' )
    tree.Branch( 'photonpfCH', photonpfCH,  'photonpfCH/D' )
    tree.Branch( 'photonpfNU', photonpfNU,  'photonpfNU/D' )
    tree.Branch( 'photonpfGA', photonpfGA,  'photonpfGA/D' )
    tree.Branch( 'photonSieie', photonSieie,  'photonSieie/D' )
    tree.Branch( 'photonHoverE', photonHoverE,  'photonHoverE/D' )
    tree.Branch( 'photonEConvVeto', photonEConvVeto,  'photonEConvVeto/O' )
    tree.Branch( 'photonPixelVeto', photonPixelVeto,  'photonPixelVeto/O' )
    
    tree.Branch( 'muon1Pt',    muon1Pt,    'muon1Pt/D' )
    tree.Branch( 'muon1Eta',   muon1Eta,   'muon1Eta/D' )
    tree.Branch( 'muon1MinDR', muon1MinDR, 'muon1MinDR/D' )
    tree.Branch( 'muon1Jet1DR', muon1Jet1DR, 'muon1Jet1DR/D' )
    tree.Branch( 'muon2Pt',    muon2Pt,    'muon2Pt/D' )
    tree.Branch( 'muon2Eta',   muon2Eta,   'muon2Eta/D' )
    tree.Branch( 'muon2MinDR', muon2MinDR, 'muon2MinDR/D' )
    tree.Branch( 'muon2Jet1DR', muon2Jet1DR, 'muon2Jet1DR/D' )
    tree.Branch( 'dimuonPt',   dimuonPt,   'dimuonPt/D' )
    tree.Branch( 'dimuonEta',  dimuonEta,  'dimuonEta/D' )
    tree.Branch( 'dimuonMinDR',dimuonMinDR,'dimuonMinDR/D' )
    tree.Branch( 'dimuonJet1DR',dimuonJet1DR,'dimuonJet1DR/D' )
    tree.Branch( 'dimuonM',    dimuonM,    'dimuonM/D' )

    tree.Branch( 'jet1Pt',   jet1Pt,    'jet1Pt/D' )
    tree.Branch( 'jet1Eta',  jet1Eta,   'jet1Eta/D' )
    tree.Branch( 'jet2Pt',   jet2Pt,    'jet2Pt/D' )
    tree.Branch( 'jet2Eta',  jet2Eta,   'jet2Eta/D' )
    tree.Branch( 'jet3Pt',   jet3Pt,    'jet3Pt/D' )
    tree.Branch( 'jet3Eta',  jet3Eta,   'jet3Eta/D' )
    tree.Branch( 'jet4Pt',   jet4Pt,    'jet4Pt/D' )
    tree.Branch( 'jet4Eta',  jet4Eta,   'jet4Eta/D' )

    tree.Branch( 'eventWt',  eventWt,   'eventWt/D' )
    tree.Branch( 'puWt',     puWt,      'puWt/D' )
    
    tree.Branch( 'passRA2ElVeto'   ,passRA2ElVeto   ,'passRA2ElVeto/O' )
    tree.Branch( 'passRA2MuVeto'   ,passRA2MuVeto   ,'passRA2MuVeto/O' )
    tree.Branch( 'passDirIsoElVeto',passDirIsoElVeto,'passDirIsoElVeto/O' )
    tree.Branch( 'passDirIsoMuVeto',passDirIsoMuVeto,'passDirIsoMuVeto/O' )
    tree.Branch( 'passIsoTrkVeto'  ,passIsoTrkVeto  ,'passIsoTrkVeto/O' )

    ##################
    myChain = r.TChain('%s/RA2Values'%(options.treeName))
    ##################
    sfCorr = 1.

    if options.debug:
        if options.sample=="photon_ttbar":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_ttjets_semilep_trees_mar7/*_?_?_???.root")
        elif options.sample=="wjets400":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_wjets400_trees_mar7/*_?_?_???.root")
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_wjets400_v2_trees_mar7/*_?_?_???.root")
        elif options.sample=="wjets400v1":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_wjets400_trees_mar7/*_?_?_???.root")
        elif options.sample=="wjets400v2":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_wjets400_v2_trees_mar7/*_?_?_???.root")
        elif options.sample=="gjets200":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_gjets200_trees_mar7/*_?_?_???.root")
        elif options.sample=="gjets400":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_gjets400_trees_mar7/*_?_?_???.root")
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_gjets400_v1_trees_mar7/*_?_?_???.root")
        elif options.sample=="gjets400v1":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_gjets400_trees_mar7/*_?_?_???.root")
        elif options.sample=="gjets400v2":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_gjets400_v1_trees_mar7/*_?_?_???.root")
        elif options.sample=="qcd":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_qcd250_trees_mar7/*_?_?_???.root")
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_qcd500_trees_mar7/*_?_?_???.root")
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_qcd1000_trees_mar7/*_?_?_???.root")
        elif options.sample=="photon2012a":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_run2012a_v1_trees_mar7/*_?_?_???.root")
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_run2012a_recoverv1_trees_mar7/*_?_?_???.root")
        elif options.sample=="photon2012b":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_run2012b_trees_mar7/*_?_?_???.root")
        elif options.sample=="photon2012c":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_run2012c_v1_trees_mar7/*_?_?_???.root")
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_run2012c_v2_trees_mar7/*_?_?_???.root")
        elif options.sample=="photon2012d":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_run2012d_v1_trees_mar7/*_?_?_???.root")
        elif options.sample=="zinv200":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zinvjets200_trees_mar7/*_?_?_???.root")
        elif options.sample=="zinv400":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zinvjets400_trees_mar7/*_?_?_???.root")
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zinvjets400_ext_trees_mar7/*_?_?_???.root")
        elif options.sample=="zinv400v1":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zinvjets400_trees_mar7/*_?_?_???.root")
        elif options.sample=="zinv400v2":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zinvjets400_ext_trees_mar7/*_?_?_???.root")
        elif options.sample=="dimuon_ttbar":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/dimuon_ttjets_semilep_trees_mar7/*_?_?_???.root")
        elif options.sample=="zmumu200":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zmumu200_trees_mar7/*_job?.root")
        elif options.sample=="zmumu400":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/dimuon_zmumu400_trees_mar7/*_?_?_???.root")
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/dimuon_zmumu400_ext_trees_mar7/*_?_?_???.root")
        elif options.sample=="zmumu400v1":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/dimuon_zmumu400_trees_mar7/*_?_?_???.root")
        elif options.sample=="zmumu400v2":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/dimuon_zmumu400_ext_trees_mar7/*_?_?_???.root")
        elif options.sample=="muon2012a":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/doublemu2012Arecoverv1_trees_mar7/*_job?.root")
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/doublemu2012Av1_trees_mar7/*_job?.root")
        elif options.sample=="muon2012b":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/doublemu2012Bv4_trees_mar7/*_job?.root")
        elif options.sample=="muon2012c":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/doublemu2012Cv1_trees_mar7/*_job?.root")
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/doublemu2012Cv2_trees_mar7/*_job?.root")
        elif options.sample=="muon2012d":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/doublemu2012Dv1_trees_mar7/*_job?.root")

    else:
        if options.sample=="photon_ttbar":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_ttjets_semilep_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="wjets400":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_wjets400_trees_mar7/%s.root"%(subfiles[options.subsec]))
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_wjets400_v2_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="wjets400v1":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_wjets400_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="wjets400v2":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_wjets400_v2_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="gjets200":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_gjets200_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="gjets400":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_gjets400_trees_mar7/%s.root"%(subfiles[options.subsec]))
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_gjets400_v2_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="gjets400v1":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_gjets400_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="gjets400v2":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_gjets400_v2_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="qcd":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_qcd250_trees_mar7/%s.root"%(subfiles[options.subsec]))
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_qcd500_trees_mar7/%s.root"%(subfiles[options.subsec]))
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_qcd1000_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="photon2012a":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_run2012a_v1_trees_mar7/%s.root"%(subfiles[options.subsec]))
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_run2012a_recoverv1_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="photon2012b":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_run2012b_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="photon2012c":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_run2012c_v1_trees_mar7/%s.root"%(subfiles[options.subsec]))
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_run2012c_v2_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="photon2012d":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_run2012d_v1_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zinv200":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zinvjets200_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zinv400":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zinvjets400_trees_mar7/%s.root"%(subfiles[options.subsec]))
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zinvjets400_ext_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zinv400v1":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zinvjets400_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zinv400v2":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zinvjets400_ext_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="dimuon_ttbar":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/dimuon_ttjets_semilep_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zmumu200":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zmumu200_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zmumu400":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/dimuon_zmumu400_trees_mar7/%s.root"%(subfiles[options.subsec]))
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/dimuon_zmumu400_ext_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zmumu400v1":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/dimuon_zmumu400_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zmumu400v2":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/dimuon_zmumu400_ext_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="muon2012a":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/doublemu2012Arecoverv1_trees_mar7/%s.root"%(subfiles[options.subsec]))
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/doublemu2012Av1_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="muon2012b":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/doublemu2012Bv4_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="muon2012c":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/doublemu2012Cv1_trees_mar7/%s.root"%(subfiles[options.subsec]))
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/doublemu2012Cv2_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="muon2012d":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/doublemu2012Dv1_trees_mar7/%s.root"%(subfiles[options.subsec]))

    fChain = myChain
##    fChain.SetBranchStatus("*",0)
##    fChain.SetBranchStatus("ra2_Event",1)
##    fChain.SetBranchStatus("ra2_Run",1)
##    fChain.SetBranchStatus("ra2_Lumi",1)
##    fChain.SetBranchStatus("ra2_Vertices",1)
##    fChain.SetBranchStatus("ra2_nJetsPt50Eta25",1)
##    fChain.SetBranchStatus("ra2_nJetsPt50Eta25MInv",1)
##    ###fChain.SetBranchStatus("ra2_nJetsCSVM",1)
##    ###fChain.SetBranchStatus("ra2_nJetsCSVT",1)
##    fChain.SetBranchStatus("ra2_nJetsPt30Eta50",1)
##    fChain.SetBranchStatus("ra2_HT",1)
##    fChain.SetBranchStatus("ra2_HTMInv",1)
##    fChain.SetBranchStatus("ra2_MHT",1)
##    ###fChain.SetBranchStatus("ra2_MET",1)
##    fChain.SetBranchStatus("ra2_dPhiMHT1",1)
##    fChain.SetBranchStatus("ra2_dPhiMHT2",1)
##    fChain.SetBranchStatus("ra2_dPhiMHT3",1)
##    ###fChain.SetBranchStatus("ra2_dPhiMHTMin",1)
##    ###fChain.SetBranchStatus("ra2_dPhiMHTMinBCSVM",1)
##    ###fChain.SetBranchStatus("ra2_dPhiMHTMinBCSVT",1)
##    fChain.SetBranchStatus("ra2_Jet1Pt",1)
##    fChain.SetBranchStatus("ra2_Jet1Eta",1)
##    fChain.SetBranchStatus("ra2_Jet2Pt",1)
##    fChain.SetBranchStatus("ra2_Jet2Eta",1)
##    fChain.SetBranchStatus("ra2_Jet3Pt",1)
##    fChain.SetBranchStatus("ra2_Jet3Eta",1)
##    fChain.SetBranchStatus("ra2_Jet4Pt",1)
##    fChain.SetBranchStatus("ra2_Jet4Eta",1)
##    fChain.SetBranchStatus("ra2_PUWt",1)
##    fChain.SetBranchStatus("ra2_EventWt",1)
##    ###
##    ###fChain.SetBranchStatus("ra2_passLooseTopTagger",1)
##    ###fChain.SetBranchStatus("ra2_loose_bestTopJetMass",1)
##    ###fChain.SetBranchStatus("ra2_loose_MTbJet",1)
##    ###fChain.SetBranchStatus("ra2_loose_MTbestTopJet",1)
##    ###fChain.SetBranchStatus("ra2_loose_MT2",1)
##    ###fChain.SetBranchStatus("ra2_loose_MTbestWJet",1)
##    ###fChain.SetBranchStatus("ra2_loose_MTbestbJet",1)
##    ###fChain.SetBranchStatus("ra2_loose_MTremainingTopJet",1)
##    ###fChain.SetBranchStatus("ra2_loose_linearCombMTbJetPlusMTbestTopJet",1)
##    ###
##    ###fChain.SetBranchStatus("ra2_passNominalTopTagger",1)
##    ###fChain.SetBranchStatus("ra2_nominal_bestTopJetMass",1)
##    ###fChain.SetBranchStatus("ra2_nominal_MTbJet",1)
##    ###fChain.SetBranchStatus("ra2_nominal_MTbestTopJet",1)
##    ###fChain.SetBranchStatus("ra2_nominal_MT2",1)
##    ###fChain.SetBranchStatus("ra2_nominal_MTbestWJet",1)
##    ###fChain.SetBranchStatus("ra2_nominal_MTbestbJet",1)
##    ###fChain.SetBranchStatus("ra2_nominal_MTremainingTopJet",1)
##    ###fChain.SetBranchStatus("ra2_nominal_linearCombMTbJetPlusMTbestTopJet",1)
##    ###
##    ###
##    if options.sample == "photon_ttbar" or options.sample == "wjets200" or options.sample == "wjets400" or options.sample == "wjets400v1" or options.sample == "wjets400v2" or options.sample == "gjets200" or options.sample == "gjets400"  or options.sample == "gjets400v1"  or options.sample == "gjets400v2"  or options.sample == "qcd" or options.sample == "photon2012a" or options.sample == "photon2012b" or options.sample == "photon2012c":
##        print "setting photon branches active"
##        fChain.SetBranchStatus("ra2_nPhotonsIso",1)
##        fChain.SetBranchStatus("ra2_nPhotonsTightIso",1)
##        fChain.SetBranchStatus("ra2_Photon1Eta",1)
##        fChain.SetBranchStatus("ra2_Photon1MinDR",1)
##        fChain.SetBranchStatus("ra2_Photon1DRJet1",1)
##        fChain.SetBranchStatus("ra2_Photon1pfCH",1)
##        fChain.SetBranchStatus("ra2_Photon1pfNU",1)
##        fChain.SetBranchStatus("ra2_Photon1pfGA",1)
##        fChain.SetBranchStatus("ra2_Photon1SigmaIetaIeta",1)
##        fChain.SetBranchStatus("ra2_Photon1HadTowOverEm",1)
##        if not options.isMC:
##            fChain.SetBranchStatus("ra2_Photon70PFHT400",1)
##            fChain.SetBranchStatus("ra2_Photon70PFNoPUHT400",1)
##            
##    if options.sample == "dimuon_ttbar" or options.sample == "zmumu200" or options.sample == "zmumu400" or options.sample == "zmumu400v1" or options.sample == "zmumu400v2" or options.sample == "muon2012a" or options.sample == "muon2012b" or options.sample == "muon2012c" or options.sample == "muon2012d":
##        print "setting dimuon branches active"
##        fChain.SetBranchStatus("ra2_Muon1Pt",1)
##        fChain.SetBranchStatus("ra2_Muon1Eta",1)
##        fChain.SetBranchStatus("ra2_Muon1MinDR",1)
##        fChain.SetBranchStatus("ra2_Muon1DRJet1",1)
##        fChain.SetBranchStatus("ra2_Muon2Pt",1)
##        fChain.SetBranchStatus("ra2_Muon2Eta",1)
##        fChain.SetBranchStatus("ra2_Muon2MinDR",1)
##        fChain.SetBranchStatus("ra2_Muon2DRJet1",1)
##        fChain.SetBranchStatus("ra2_DiMuonPt",1)
##        fChain.SetBranchStatus("ra2_DiMuonEta",1)
##        fChain.SetBranchStatus("ra2_DiMuonMinDR",1)
##        fChain.SetBranchStatus("ra2_DiMuonDRJet1",1)
##        fChain.SetBranchStatus("ra2_DiMuonInvM",1)
##        if not options.isMC:
##            fChain.SetBranchStatus("ra2_Mu13_Mu8",1)
##            fChain.SetBranchStatus("ra2_Mu17_Mu8",1)
#####
    ###Timing information
    decade  = 0
    century = 0
    tsw = r.TStopwatch()
    tenpcount = 1
    onepcount = 1


    nentries = fChain.GetEntries()
    print "nentries %d"%(nentries)
    sys.stdout.flush()
    i = 0
    for event in fChain:
        # ==============print number of events done == == == == == == == =
        if ( i==0):
            tsw.Start()
            #print('.', end='')
            sys.stdout.write('.')
            sys.stdout.flush()
        if ((i*10)/nentries == tenpcount ) :
            tsw.Stop() 
            time = tsw.RealTime() 
            tsw.Start(r.kFALSE) 
            finTime = 0.
            frac = (i*1.0)/(nentries*1.0) 
            if (frac>0):
                finTime = time / frac - time 
                finMin = finTime / 60. 
                sys.stdout.write("%d%% done.  "%(tenpcount*10))
                # sys.stdout.write("t=7.2f"%(time))
                sys.stdout.write("t="+str(time))
                sys.stdout.write(" projected finish=%7d s("%(finTime))
                sys.stdout.write("%2.2f min).   "%(finMin))
                sys.stdout.write("\n")
                sys.stdout.flush()
                tenpcount = tenpcount + 1
        
        elif ( (i*100)/nentries == onepcount ) :
            #print('.', end='')
            sys.stdout.write('.')
            sys.stdout.flush()
            onepcount = onepcount + 1

        ra2_run[0]     = event.ra2_Run
        ra2_event[0]   = event.ra2_Event
        ra2_lumi[0]    = event.ra2_Lumi
        nVtx[0]        = event.ra2_Vertices
        nJetsHT[0]     = event.ra2_nJetsPt50Eta25
        nJetsHTMInv[0] = event.ra2_nJetsPt50Eta25MInv
        nJetsMHT[0]    = event.ra2_nJetsPt30Eta50

        htVal[0]       = event.ra2_HT
        htMInvVal[0]   = event.ra2_HTMInv
        mhtVal[0]      = event.ra2_MHT

        dphi1[0]       = event.ra2_dPhiMHT1
        dphi2[0]       = event.ra2_dPhiMHT2
        dphi3[0]       = event.ra2_dPhiMHT3
        dphiMin[0]     = event.ra2_dPhiMHTMin

        jet1Pt[0]      = event.ra2_Jet1Pt
        jet1Eta[0]     = event.ra2_Jet1Eta
        jet2Pt[0]      = event.ra2_Jet2Pt
        jet2Eta[0]     = event.ra2_Jet2Eta
        jet3Pt[0]      = event.ra2_Jet3Pt
        jet3Eta[0]     = event.ra2_Jet3Eta
        jet4Pt[0]      = event.ra2_Jet4Pt
        jet4Eta[0]     = event.ra2_Jet4Eta

        puWt[0]        = event.ra2_PUWt
        eventWt[0]     = event.ra2_EventWt*sfCorr

        passRA2ElVeto[0]     = event.ra2_passRA2ElVeto
        passRA2MuVeto[0]     = event.ra2_passRA2MuVeto
        passDirIsoElVeto[0]  = event.ra2_passDirIsoElVeto
        passDirIsoMuVeto[0]  = event.ra2_passDirIsoMuVeto
        passIsoTrkVeto[0]    = event.ra2_passIsoTrkVeto

        triggers = True

        nPhotonsID[0]    = 0
#        nPhotonsLoose[0] = 0
        nPhotonsTight[0] = 0
#        nPhotonsCombR03[0]  = 0
#        nPhotonsCombR04[0]  = 0
#        photonIsTightID[0]  = 0
#        photonIsLooseIso[0] = 0
        photonIsTightIso[0] = 0
#        photonIsCombIso[0]  = 0
        photonPt[0]    = -10.
        photonEta[0]   = -10.
        photonPhi[0]   = -10.
        photonMinDR[0] = 10.
        photonJet1DR[0]= 10.
        photonpfCH[0]  = -10.
        photonpfNU[0]  = -10.
        photonpfGA[0]  = -10.
        photonSieie[0] = -10.
        photonHoverE[0]= -10.
        photonEConvVeto[0]= 0
        photonPixelVeto[0]= 0

        muon1Pt[0]     = -10.
        muon1Eta[0]    = 10.
        muon1MinDR[0]  = 10.
        muon1Jet1DR[0] = 10.
        muon2Pt[0]     = -10.
        muon2Eta[0]    = 10.
        muon2MinDR[0]  = 10.
        muon2Jet1DR[0] = 10.
        dimuonPt[0]    = -10.
        dimuonEta[0]   = 10.
        dimuonMinDR[0] = 10.
        dimuonJet1DR[0]= 10.
        dimuonM[0]     = -10.

        if options.sample in ["photon_ttbar", "wjets200", "wjets400", "wjets400v1", "wjets400v2",
                              "gjets200", "gjets400", "gjets400v1", "gjets400v2", "qcd",
                              "photon2012a", "photon2012b", "photon2012c", "photon2012d"]:
            nPhotonsID[0]      = event.ra2_nPhotonsID
#            nPhotonsLoose[0]   = event.ra2_nPhotonsLooseIso
            nPhotonsTight[0]   = event.ra2_nPhotonsTightIso
#            nPhotonsCombR03[0] = event.ra2_nPhotonsCombIsoR03
#            nPhotonsCombR04[0] = event.ra2_nPhotonsCombIsoR04
#            photonIsTightID[0]  = event.ra2_Photon1IsTightID
#            photonIsLooseIso[0] = event.ra2_Photon1IsLoosePFIso
            photonIsTightIso[0] = event.ra2_Photon1IsTightPFIso
#            photonIsCombIso[0]  = event.ra2_Photon1IsCombIsoR03
            photonPt[0]    = event.ra2_Photon1Pt
            photonPhi[0]   = event.ra2_Photon1Phi
            photonEta[0]   = event.ra2_Photon1Eta
            photonMinDR[0] = event.ra2_Photon1MinDR
            photonJet1DR[0]= event.ra2_Photon1DRJet1
            photonpfCH[0]  = event.ra2_Photon1pfCH
            photonpfNU[0]  = event.ra2_Photon1pfNU
            photonpfGA[0]  = event.ra2_Photon1pfGA
            photonSieie[0]  = event.ra2_Photon1SigmaIetaIeta
            photonHoverE[0] = event.ra2_Photon1HadTowOverEm
            photonEConvVeto[0] = event.ra2_Photon1EConvVeto
            photonPixelVeto[0] = event.ra2_Photon1PixelVeto
            if not options.isMC:
                triggers = (event.ra2_Photon70PFHT400 or event.ra2_Photon70PFNoPUHT400)
                
        elif options.sample in ["dimuon_ttbar", "zmumu200", "zmumu400", "zmumu400v1", "zmumu400v2",
                                "muon2012a", "muon2012b", "muon2012c", "muon2012d"]:
            muon1Pt[0]     = event.ra2_Muon1Pt
            muon1Eta[0]    = event.ra2_Muon1Eta
            muon1MinDR[0]  = event.ra2_Muon1MinDR
            muon1Jet1DR[0] = event.ra2_Muon1DRJet1
            muon2Pt[0]     = event.ra2_Muon2Pt
            muon2Eta[0]    = event.ra2_Muon2Eta
            muon2MinDR[0]  = event.ra2_Muon2MinDR
            muon2Jet1DR[0] = event.ra2_Muon2DRJet1
            dimuonPt[0]    = event.ra2_DiMuonPt
            dimuonEta[0]   = event.ra2_DiMuonEta
            dimuonMinDR[0] = event.ra2_DiMuonMinDR
            dimuonJet1DR[0]= event.ra2_DiMuonDRJet1
            dimuonM[0]     = event.ra2_DiMuonInvM
            if not options.isMC:
                triggers = (event.ra2_Mu13_Mu8 or event.ra2_Mu17_Mu8)

        if options.isMC:
            triggers = True

        extra = True
        
        if options.sample in ["photon_ttbar", "wjets200", "wjets400", "wjets400v1", "wjets400v2",
                              "gjets200", "gjets400", "gjets400v1", "gjets400v2", "qcd",
                              "photon2012a", "photon2012b", "photon2012c", "photon2012d"]:
            extra = photonPt[0] > 100

        elif options.sample in ["dimuon_ttbar", "zmumu200", "zmumu400", "zmumu400v1", "zmumu400v2",
                                "muon2012a", "muon2012b", "muon2012c", "muon2012d"]:
            extra = dimuonPt[0] > 100

        if (triggers) and extra and mhtVal[0] > 100 and htVal[0] > 300 and nJetsHT[0] > 1:
            tree.Fill()

        #########
        i = i + 1
    #tree.Write()
    outFile.Write()
    outFile.Close()

    ##########################
if __name__ == '__main__':
    main()
    print "done with __main__!"
    sys.stdout.flush()
    

