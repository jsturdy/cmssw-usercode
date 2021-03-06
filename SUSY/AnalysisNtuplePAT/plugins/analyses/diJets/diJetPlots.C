#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <time.h>
#include "TFile.h"
#include "TROOT.h" 
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TTree.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TImage.h"
#include "TRandom.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TCut.h"
#include "TPostScript.h"
#include "TF1.h"
#include "TLeaf.h"
#include <TLegend.h>
#include <TKey.h>

using namespace std;

void diJetPlots(std::string mettype="CaloMET", std::string jettype="CaloJets", std::string energy="7TeV")
{
  static const int NUMHISTS  = 24;//23;
  static const int NUMHISTS2 = 4;
  static const int NUMVARS   = 4;
  static const int NUMFILES  = 7;
  static const int NUMFILES2 = 5;
  static const double localpi  = acos(-1);

  //gROOT->SetStyle("Plain");

  TFile* file[NUMFILES];
  TFile* file2[NUMFILES2];
  TString filenames[NUMFILES];
  TString filenames2[NUMFILES2];

  Color_t linecolor[NUMFILES];
  Color_t linecolor2[NUMFILES];
  Style_t linestyle[NUMFILES];
  UInt_t  linewidth;

  Color_t fillcolor[NUMFILES];
  Color_t fillcolor2[NUMFILES];
  Style_t fillstyle[NUMFILES];

  Style_t markerstyle[NUMFILES];
  UInt_t  markersize;

  TCanvas *mycanvas[NUMHISTS][NUMVARS];
  TCanvas *mycanvas2[NUMHISTS2][2];
 
  Double_t max[NUMHISTS][NUMVARS] = {{0.0}};
  Double_t max2[NUMHISTS2][2]     = {{0.0}};

  TString histnames[NUMHISTS][NUMVARS];
  TString histnames2[NUMHISTS2][2];
  TLegend *leg;
  TString histtitle[NUMHISTS];
  TString histtitle2[NUMHISTS2+1];

  TH1F *hist[NUMFILES][NUMHISTS][NUMVARS];
  TH1F *hist2[NUMFILES][NUMHISTS2][2];

  TH1F *histbac[NUMFILES][NUMHISTS][NUMVARS];
  TH1F *histbac2[NUMFILES][NUMHISTS2][2];

  THStack* thestack[NUMHISTS][NUMVARS];
  THStack* thestack2[NUMHISTS2][NUMVARS];
  
  //TString suffix_ = "_met_jCalomCaloTypeIl";
  TString suffix_ = "_mht_jCalomCaloTypeIl";
  //TString suffix_ = "_full_jCalomCaloTypeIl";
  filenames[0]  = "LM0";
  filenames[1]  = "LM1";
  filenames[2]  = "LM5";
  filenames[3]  = "LM6";
  filenames[4]  = "LM12";
  filenames[5]  = "LM13";
  filenames[6]  = "7TeV_Data";

  filenames2[0]  = "ZJets-madgraph";
  filenames2[1]  = "TTbarJets-madgraph";  
  filenames2[2]  = "ZInvisibleJets";
  filenames2[3]  = "WJets-madgraph";
  filenames2[4]  = "QCD_MadGraph_Pt50toInf";
  //filenames2[5]  = "VectorBosons";
  //filenames2[6]  = "SM_Background";


  //filenames[1]  = "Zmumu";
  //filenames[1]  = "Wmunu";

  histnames[0][0]  = "h_pre_cuts_01_jet1et";       histnames[0][1]  = "h_individual_cuts_01_jet1et";       histnames[0][2]  = "h_N1_cuts_01_jet1et";            histnames[0][3]  = "h_post_cuts_01_jet1et";
  histnames[1][0]  = "h_pre_cuts_02_jet2et";       histnames[1][1]  = "h_individual_cuts_02_jet2et";       histnames[1][2]  = "h_N1_cuts_02_jet2et";            histnames[1][3]  = "h_post_cuts_02_jet2et";
  histnames[2][0]  = "h_pre_cuts_03_jetallet";     histnames[2][1]  = "h_individual_cuts_03_jetallet";     histnames[2][2]  = "h_N1_cuts_03_jetallet";          histnames[2][3]  = "h_post_cuts_03_jetallet";
  histnames[3][0]  = "h_pre_cuts_11_MET";          histnames[3][1]  = "h_individual_cuts_11_MET";	   histnames[3][2]  = "h_N1_cuts_11_MET";	        histnames[3][3]  = "h_post_cuts_11_MET";
  histnames[4][0]  = "h_pre_cuts_12_HT";           histnames[4][1]  = "h_individual_cuts_12_HT";	   histnames[4][2]  = "h_N1_cuts_12_HT";	        histnames[4][3]  = "h_post_cuts_12_HT";
  histnames[5][0]  = "h_pre_cuts_13_MHT";          histnames[5][1]  = "h_individual_cuts_13_MHT";	   histnames[5][2]  = "h_N1_cuts_13_MHT";	        histnames[5][3]  = "h_post_cuts_13_MHT";
  histnames[6][0]  = "h_pre_cuts_14_Meff";         histnames[6][1]  = "h_individual_cuts_14_Meff";  	   histnames[6][2]  = "h_N1_cuts_14_Meff";	        histnames[6][3]  = "h_post_cuts_14_Meff";
  histnames[7][0]  = "h_pre_cuts_21_jet1metdphi";  histnames[7][1]  = "h_individual_cuts_21_jet1metdphi";  histnames[7][2]  = "h_N1_cuts_21_jet1metdphi";       histnames[7][3]  = "h_post_cuts_21_jet1metdphi";
  histnames[8][0]  = "h_pre_cuts_22_jet2metdphi";  histnames[8][1]  = "h_individual_cuts_22_jet2metdphi";  histnames[8][2]  = "h_N1_cuts_22_jet2metdphi";       histnames[8][3]  = "h_post_cuts_22_jet2metdphi";
  histnames[9][0]  = "h_pre_cuts_23_jet12dphi";    histnames[9][1]  = "h_individual_cuts_23_jet12dphi";    histnames[9][2]  = "h_N1_cuts_23_jet12dphi";         histnames[9][3]  = "h_post_cuts_23_jet12dphi";
  histnames[10][0] = "h_pre_cuts_24_dphistar";     histnames[10][1] = "h_individual_cuts_24_dphistar";     histnames[10][2] = "h_N1_cuts_24_dphistar";          histnames[10][3]  = "h_post_cuts_24_dphistar";
  histnames[11][0] = "h_pre_cuts_30_Njets";        histnames[11][1] = "h_individual_cuts_30_Njets";	   histnames[11][2] = "h_N1_cuts_30_Njets";	        histnames[11][3] = "h_post_cuts_30_Njets";
  histnames[12][0] = "h_pre_cuts_31_Ngoodjets";    histnames[12][1] = "h_individual_cuts_31_Ngoodjets";    histnames[12][2] = "h_N1_cuts_31_Ngoodjets";         histnames[12][3] = "h_post_cuts_31_Ngoodjets";
  histnames[13][0] = "h_pre_cuts_32_jet1eta";      histnames[13][1] = "h_individual_cuts_32_jet1eta";      histnames[13][2] = "h_N1_cuts_32_jet1eta";           histnames[13][3] = "h_post_cuts_32_jet1eta";
  histnames[14][0] = "h_pre_cuts_33_jet2eta";      histnames[14][1] = "h_individual_cuts_33_jet2eta";      histnames[14][2] = "h_N1_cuts_33_jet2eta";           histnames[14][3] = "h_post_cuts_33_jet2eta";
  histnames[15][0] = "h_pre_cuts_41_jetFem";       histnames[15][1] = "h_individual_cuts_41_jetFem";       histnames[15][2] = "h_N1_cuts_41_jetFem";            histnames[15][3] = "h_post_cuts_41_jetFem";
  histnames[16][0] = "h_pre_cuts_42_jet1emfrac";   histnames[16][1] = "h_individual_cuts_42_jet1emfrac";   histnames[16][2] = "h_N1_cuts_42_jet1emfrac";        histnames[16][3] = "h_post_cuts_42_jet1emfrac";
  histnames[17][0] = "h_pre_cuts_43_jet2emfrac";   histnames[17][1] = "h_individual_cuts_43_jet2emfrac";   histnames[17][2] = "h_N1_cuts_43_jet2emfrac";        histnames[17][3] = "h_post_cuts_43_jet2emfrac";
  histnames[18][0] = "h_pre_cuts_51a_Nelecs";      histnames[18][1] = "h_individual_cuts_51a_Nelecs"; 	   histnames[18][2] = "h_N1_cuts_51a_Nelecs";           histnames[18][3] = "h_post_cuts_51a_Nelecs";
  histnames[19][0] = "h_pre_cuts_52a_Nmuons";      histnames[19][1] = "h_individual_cuts_52a_Nmuons"; 	   histnames[19][2] = "h_N1_cuts_52a_Nmuons";           histnames[19][3] = "h_post_cuts_52a_Nmuons";
  histnames[20][0] = "h_pre_cuts_51b_Ngoodelecs";  histnames[20][1] = "h_individual_cuts_51b_Ngoodelecs";  histnames[20][2] = "h_N1_cuts_51b_Ngoodelecs";        histnames[20][3] = "h_post_cuts_51b_Ngoodelecs";
  histnames[21][0] = "h_pre_cuts_52b_Ngoodmuons";  histnames[21][1] = "h_individual_cuts_52b_Ngoodmuons";  histnames[21][2] = "h_N1_cuts_52b_Ngoodmuons";        histnames[21][3] = "h_post_cuts_52b_Ngoodmuons";
  //histnames[19][0] = "h_pre_cuts_53_eleceta";      histnames[19][1] = "h_individual_cuts_53_eleceta";	   histnames[19][2] = "h_N1_cuts_53_eleceta";	   histnames[19][3] = "h_post_cuts_53_eleceta";
  //histnames[20][0] = "h_pre_cuts_54_muoneta";      histnames[20][1] = "h_individual_cuts_54_muoneta";	   histnames[20][2] = "h_N1_cuts_54_muoneta";	   histnames[20][3] = "h_post_cuts_54_muoneta";
  histnames[22][0] = "h_pre_cuts_55_elecet";       histnames[22][1] = "h_individual_cuts_55_elecet";	   histnames[22][2] = "h_N1_cuts_55_elecet";	        histnames[22][3] = "h_post_cuts_55_elecet";
  histnames[23][0] = "h_pre_cuts_56_muonet";       histnames[23][1] = "h_individual_cuts_56_muonet";	   histnames[23][2] = "h_N1_cuts_56_muonet";	        histnames[23][3] = "h_post_cuts_56_muonet";

  //histnames2[0][0] = "h_pre_cuts_9_METphi";   histnames2[0][1] = "h_post_cuts_9_METphi";
  //histnames2[1][0] = "h_pre_cuts_9_jet1phi";  histnames2[1][1] = "h_post_cuts_9_jet1phi";
  //histnames2[2][0] = "h_pre_cuts_9_jet2phi";  histnames2[2][1] = "h_post_cuts_9_jet2phi";
  histnames2[0][0] = "h_pre_cuts_9_MT";       histnames2[0][1] = "h_post_cuts_9_MT";
  histnames2[1][0] = "h_pre_cuts_9_Minv";     histnames2[1][1] = "h_post_cuts_9_Minv";
  histnames2[2][0] = "h_pre_cuts_9_SumEt";    histnames2[2][1] = "h_post_cuts_9_SumEt";
  histnames2[3][0] = "h_selections";          histnames2[3][1] = "h_N1_selections";

  histtitle[0]  = "E_{T}^{J_{1}}";
  histtitle[1]  = "E_{T}^{J_{2}}";
  histtitle[2]  = "E_{T}^{J_{3+}}";
  histtitle[3]  = "#slash E_{T}";
  histtitle[4]  = "H_{T}";
  histtitle[5]  = "#slash H_{T}";
  histtitle[6]  = "M_{eff}";
  histtitle[7]  = "#Delta#phi(J_{1};#slash E_{T})";
  histtitle[8]  = "#Delta#phi(J_{2};#slash E_{T})";
  histtitle[9]  = "#Delta#phi(J_{1}; J_{2})";
  histtitle[10] = "#Delta#phi*";
  histtitle[11] = "N_{jets}";
  histtitle[12] = "N^{good}_{jets}";
  histtitle[13] = "#eta^{J_{1}}";
  histtitle[14] = "#eta^{J_{2}} ";
  histtitle[15] = "E^{J_{3+}}_{EM}/E^{J_{3+}}_{tot}";
  histtitle[16] = "E^{J_{1}}_{EM}/E^{J_{1}}_{tot}";
  histtitle[17] = "E^{J_{2}}_{EM}/E^{J_{2}}_{tot}";
  histtitle[18] = "N_{e}";
  histtitle[19] = "N_{#mu}";
  histtitle[20] = "N^{good}_{e}";
  histtitle[21] = "N^{good}_{#mu}";
  //histtitle[19] = "#eta_{e}";
  //histtitle[20] = "#eta_{#mu}";
  histtitle[22] = "E^{e}_{T}";
  histtitle[23] = "E^{#mu}_{T}";
  //histtitle[25] = "E^{#mu}_{T}";			       
  //histtitle[26] = "E^{#mu}_{T}";

  //histtitle2[0] = "#slashE_{T}#phi";
  //histtitle2[1] = "#phi_{J_{1}}";
  //histtitle2[2] = "#phi_{J_{2}}";
  histtitle2[0] = "M_{T}";
  histtitle2[1] = "M_{Inv}";
  histtitle2[2] = "#SigmaE_{T}";
  histtitle2[3] = "selections";
  histtitle2[4] = "N-1 selections";

 
  linecolor[0] = kBlue;//LM0
  linecolor[1] = kRed+1;//LM1
  linecolor[2] = kGreen+1;//LM5
  linecolor[3] = kOrange+7;//SM
  linecolor[4] = kBlue+3;//LQ to CMu M300
  linecolor[5] = kYellow+3;//LQ to CMu M400
  linecolor[6] = kBlack+3;//LQ to CMu M400

  linecolor2[0] = 8;//kGreen-2;//Z+Jets
  linecolor2[1] = 9;//kBlue-2;//TTbar+Jets
  linecolor2[2] = 26;//kTeal-7;//Z Invi +Jets
  linecolor2[3] = 30;//kRed-3;//W+Jets
  linecolor2[4] = 46;//kOrange+3;//QCD

  linestyle[0] = 4;
  linestyle[1] = 5;
  linestyle[2] = 6;
  linestyle[3] = 7;
  linestyle[4] = 8;
  linestyle[5] = 9;
  linestyle[6] = 10;

  fillcolor[0] = kBlue-4;//LM0
  fillcolor[1] = kRed;//LM1
  fillcolor[2] = kGreen;//LM5
  fillcolor[3] = kOrange;//SM
  fillcolor[4] = kBlue-1;//LQ to CMu M300
  fillcolor[5] = kYellow+1;//LQ to CMu M400
  fillcolor[6] = kBlack+1;//LQ to CMu M400

  fillcolor2[0] = 8;//kGreen-2;//Z+Jets
  fillcolor2[1] = 9;//kBlue-2;//TTbar+Jets
  fillcolor2[2] = 28;//kTeal-7;//Z Invi +Jets
  fillcolor2[3] = 30;//kRed-3;//W+Jets
  fillcolor2[4] = 46;//kOrange+3;//QCD

  fillstyle[0] = 4;
  fillstyle[1] = 5;
  fillstyle[2] = 6;
  fillstyle[3] = 7;
  fillstyle[4] = 8;
  fillstyle[5] = 9;
  fillstyle[6] = 10;

  markerstyle[0] = 20;//kFullCircle;
  markerstyle[1] = 21;//kFullSquare;
  markerstyle[2] = 22;//kFullTriangelUp;
  markerstyle[3] = 23;//kFullTriangleDown;
  markerstyle[4] = 27;//kFullDiamond;
  markerstyle[5] = 29;//kFullStar;
  markerstyle[6] = 28;//kFullCross;

  markersize = 1.25;

  linewidth = 1;

  for (int tt = 0; tt < NUMVARS; tt++) {
    for (int z = 0; z < NUMHISTS; z++) {
      mycanvas[z][tt]  = new TCanvas(histnames[z][tt], "", 800,600);}}
  
  for (int ss = 0; ss < 2; ss++) {
    for (int y = 0; y < NUMHISTS2; y++) {
      mycanvas2[y][ss] = new TCanvas(histnames2[y][ss], "", 800,600);}}

  TString metTag    = mettype;
  TString jetTag    = jettype;
  TString energyTag = energy;

  for (int j = 0; j < NUMFILES; j++) {
    //TString filepath = "./MHT_Analysis/CaloJets/"+filenames[j]+suffix_+".root";
    TString filepath = "./MHT_Analysis/"+filenames[j]+suffix_+".root";
    file[j] = new TFile(filepath);
    
    for (int tt = 0; tt < NUMVARS; tt++) {
      for (int z = 0; z < NUMHISTS; z++) {
	hist[j][z][tt] = (TH1F*)gDirectory->Get(histnames[z][tt]);
	int nbins =  hist[j][z][tt]->GetNbinsX();
	double overflows = hist[j][z][tt]->GetBinContent(nbins+1);
	double lastbin   = hist[j][z][tt]->GetBinContent(nbins);
	hist[j][z][tt]->SetBinContent(nbins,overflows+lastbin);
	double histmax = hist[j][z][tt]->GetMaximum();
	max[z][tt] = (max[z][tt]>histmax)?max[z][tt]:histmax;}}
    
    for (int ss = 0; ss < 2; ss++) {
      for (int y = 0; y < NUMHISTS2; y++) {
	hist2[j][y][ss] = (TH1F*)gDirectory->Get(histnames2[y][ss]);
	int nbins =  hist2[j][y][ss]->GetNbinsX();
	double overflows = hist2[j][y][ss]->GetBinContent(nbins+1);
	double lastbin   = hist2[j][y][ss]->GetBinContent(nbins);
	hist2[j][y][ss]->SetBinContent(nbins,overflows+lastbin);
	Double_t hist2max = hist2[j][y][ss]->GetMaximum();
	max2[y][ss] = (max2[y][ss]>hist2max)?max2[y][ss]:hist2max;}}}

  double tmpmax[NUMHISTS][NUMVARS] = {{0.}};
  double tmpmax2[NUMHISTS2][2] = {{0.}};
    
  for (int tt = 0; tt < NUMVARS; tt++) {
    for (int z = 0; z < NUMHISTS; z++) {
      thestack[z][tt] = new THStack();
      for (int j = 0; j < NUMFILES2; j++) {
	//TString filepath = "./MHT_Analysis/CaloJets/"+filenames[j]+suffix_+".root";
	TString filepath = "./MHT_Analysis/"+filenames2[j]+suffix_+".root";
	file2[j] = new TFile(filepath);
	histbac[j][z][tt] = (TH1F*)gDirectory->Get(histnames[z][tt]);
	int nbins =  histbac[j][z][tt]->GetNbinsX();
	double overflows = histbac[j][z][tt]->GetBinContent(nbins+1);
	double lastbin   = histbac[j][z][tt]->GetBinContent(nbins);
	histbac[j][z][tt]->SetBinContent(nbins,overflows+lastbin);
	tmpmax[z][tt] += histbac[j][z][tt]->GetMaximum();}
      max[z][tt] = (max[z][tt]>tmpmax[z][tt])?max[z][tt]:tmpmax[z][tt];
    }
  }
  
  for (int ss = 0; ss < 2; ss++) {
    for (int y = 0; y < NUMHISTS2; y++) {
      thestack2[y][ss] = new THStack();
      for (int j = 0; j < NUMFILES2; j++) {
	histbac2[j][y][ss] = (TH1F*)file2[j]->Get(histnames2[y][ss]);
	int nbins =  histbac2[j][y][ss]->GetNbinsX();
	double overflows = histbac2[j][y][ss]->GetBinContent(nbins+1);
	double lastbin   = histbac2[j][y][ss]->GetBinContent(nbins);
	histbac2[j][y][ss]->SetBinContent(nbins,overflows+lastbin);
	tmpmax2[y][ss] += histbac2[j][y][ss]->GetMaximum();}
      max2[y][ss] = (max2[y][ss]>tmpmax2[y][ss])?max2[y][ss]:tmpmax2[y][ss];
    }
  }
  
  for (int tt = 0; tt < NUMVARS; tt++) {
    for (int z = 0; z < NUMHISTS; z++ ) max[z][tt] *= 1.25;}
  
  for (int ss = 0; ss < 2; ss++) {
    for (int y = 0; y < NUMHISTS2; y++ ) max2[y][ss] *= 1.25;}

  for (int j = 0; j < NUMFILES; j++) {
    for (int tt = 0; tt < NUMVARS; tt++) {
      for (int z = 0; z < NUMHISTS; z++) {
	//hist[j][z][tt]->SetFillStyle(3001);
	hist[j][z][tt]->SetLineColor(linecolor[j]);
	hist[j][z][tt]->SetLineStyle(linestyle[j]);
	hist[j][z][tt]->SetLineWidth(linewidth);
	hist[j][z][tt]->SetMarkerColor(linecolor[j]);
	hist[j][z][tt]->SetMarkerStyle(markerstyle[j]);
	hist[j][z][tt]->SetMarkerSize(markersize);
	//hist[j][z][tt]->SetLineStyle();
	//hist[j][z][tt]->SetFillColor(fillcolor[j]);
	hist[j][z][tt]->SetStats(kFALSE);}}
    
    for (int ss = 0; ss < 2; ss++) {
      for (int y = 0; y < NUMHISTS2; y++) {
	//hist2[j][y][ss]->SetFillStyle(3001);
	hist2[j][y][ss]->SetLineColor(linecolor[j]);
	hist2[j][y][ss]->SetLineStyle(linestyle[j]);
	hist2[j][y][ss]->SetLineWidth(linewidth);
	hist2[j][y][ss]->SetMarkerColor(linecolor[j]);
	hist2[j][y][ss]->SetMarkerStyle(markerstyle[j]);
	hist2[j][y][ss]->SetMarkerSize(markersize);
	//hist2[j][y][ss]->SetLineStyle();
	//hist2[j][y][ss]->SetFillColor(fillcolor[j]);
	hist2[j][y][ss]->SetStats(kFALSE);}}}

  
  //backgrounds in the stack
  for (int tt = 0; tt < NUMVARS; tt++) {
    for (int z = 0; z < NUMHISTS; z++) {
      for (int j = 0; j < NUMFILES2; j++) {
	histbac[j][z][tt]->SetFillStyle(1001);
	histbac[j][z][tt]->SetFillColor(fillcolor2[j]);
	histbac[j][z][tt]->SetLineColor(linecolor2[j]);
	//histbac[j][z][tt]->SetLineStyle(linestyle[j]);
	//histbac[j][z][tt]->SetLineWidth(linewidth);
	//histbac[j][z][tt]->SetMarkerColor(linecolor2[j]);
	//histbac[j][z][tt]->SetMarkerStyle(markerstyle[j]);
	//histbac[j][z][tt]->SetMarkerSize(markersize);
	histbac[j][z][tt]->SetStats(kFALSE);
	thestack[z][tt]->Add(histbac[j][z][tt],"bar");}
      thestack[z][tt]->BuildStack();}}
  
  for (int ss = 0; ss < 2; ss++) {
    for (int y = 0; y < NUMHISTS2; y++) {
      for (int j = 0; j < NUMFILES2; j++) {
	histbac2[j][y][ss]->SetFillStyle(1001);
	histbac2[j][y][ss]->SetFillColor(fillcolor2[j]);
	histbac2[j][y][ss]->SetLineColor(linecolor2[j]);
	//histbac2[j][y][ss]->SetLineStyle(linestyle[j]);
	//histbac2[j][y][ss]->SetLineWidth(linewidth);
	//histbac2[j][y][ss]->SetMarkerColor(linecolor2[j]);
	//histbac2[j][y][ss]->SetMarkerStyle(markerstyle[j]);
	//histbac2[j][y][ss]->SetMarkerSize(markersize);
	histbac2[j][y][ss]->SetStats(kFALSE);
	thestack2[y][ss]->Add(histbac2[j][y][ss],"bar");}
      thestack2[y][ss]->BuildStack();}}

  //leg = new TLegend(.4, .6, 1, 1);
  leg = new TLegend(0., 0., 1, 0.8);
  //leg->SetDrawOption("Plain");
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);
  for (int j = 0; j < NUMFILES; j++) {
    TString extra = "";
    if (j == 6)
      extra = " #intLdt = 71.7 nb^{-1}";
    leg->AddEntry(hist[j][0][0], filenames[j]+extra,"lep");
    //leg->SetTextSize(0.05);
    //leg->GetEntry()->SetTextSize(2);
    //leg->GetEntry()->SetEntrySeparation(0.075);
  }
  for (int j = 0; j < NUMFILES2; j++) {
    leg->AddEntry(histbac[j][0][0], filenames2[j],"fel");
    //leg->SetTextSize(0.05);
    //leg->SetTextSize(0.2);
    //leg->GetEntry()->SetTextSize(2);
    //leg->GetEntry()->SetEntrySeparation(0.075);
  }

  for(int tt = 0; tt < NUMVARS; tt++) {
    for (int z = 0; z < NUMHISTS; z++) {
      mycanvas[z][tt]->cd();

      TPad *subpad1e = new TPad("subpad1e","Sub-pad number 1e",0.00,0.00,0.8,1.00);
      TPad *subpad1f = new TPad("subpad1f","Sub-pad number 1f",0.725,0.5,1,1.00);
      subpad1e->Draw();
      subpad1f->Draw();
      subpad1e->cd();      
      gStyle->SetOptTitle(kFALSE);
      hist[0][z][tt]->SetMaximum(max[z][tt]);
      //hist[0][z][tt]->SetMinimum(0.0001);
      hist[0][z][tt]->GetXaxis()->SetTitle(histtitle[z]);
      hist[0][z][tt]->GetXaxis()->SetLabelSize(0.03);
      hist[0][z][tt]->GetYaxis()->SetTitleOffset(1.2);
      //hist[0][z][tt]->GetYaxis()->SetLabelOffset(-0.05);
      //if (z==3) {
      //  char ytitle[128];
      //  sprintf(ytitle,"Events / 25 GeV / %2.0f pb^{-1}",100);
      //  hist[0][z][tt]->GetYaxis()->SetTitle(ytitle);}
      //if (z==19||z==20) {
      //  if (tt==1||tt==3) {
      //    char ytitle[128];
      //    sprintf(ytitle,"Events / 1 GeV / %2.0f pb^{-1}",100);
      //    hist[0][z][tt]->GetYaxis()->SetTitle(ytitle);}}

      hist[0][z][tt]->GetYaxis()->SetLabelSize(0.03);
      gPad->SetLogy(kTRUE);
      hist[0][z][tt]->Draw("e0p0");
      thestack[z][tt]->Draw("same");
      hist[0][z][tt]->Draw("e0p0same");
      for (int j = 1; j < NUMFILES; j++) hist[j][z][tt]->Draw("e0p0same");

      subpad1f->cd();
      leg->Draw();
      mycanvas[z][tt]->Update();
      char outputimage[128];
      std::string outhistname = (std::string)histnames[z][tt];
      //sprintf(outputimage,"./MHT_Analysis/CaloJets/%s.png",outhistname.c_str());
      sprintf(outputimage,"./MHT_Analysis/%s.png",outhistname.c_str());
      mycanvas[z][tt]->SaveAs(outputimage);}}
  //mycanvas[z][tt]->SaveAs();}}

  for(int ss = 0; ss <2; ss++) {
    for (int y = 0; y < NUMHISTS2; y++) {
      mycanvas2[y][ss]->cd();
      TPad *subpad1e = new TPad("subpad1e","Sub-pad number 1e",0.00,0.00,0.8,1.00);
      TPad *subpad1f = new TPad("subpad1f","Sub-pad number 1f",0.725,0.5,1,1.00);
      subpad1e->Draw();
      subpad1f->Draw();
      subpad1e->cd();
      gStyle->SetOptTitle(kFALSE);

      hist2[0][y][ss]->SetMaximum(max2[y][ss]);
      //hist2[0][y][ss]->SetMinimum(0.0001);
      hist2[0][y][ss]->GetXaxis()->SetTitle(histtitle2[y]);
      hist2[0][y][ss]->GetXaxis()->SetLabelSize(0.03);
      hist2[0][y][ss]->GetYaxis()->SetTitleOffset(1.2);
      //hist2[0][y][ss]->GetYaxis()->SetLabelOffset(-0.05);
      hist2[0][y][ss]->GetYaxis()->SetLabelSize(0.03);
      gPad->SetLogy(kTRUE);
      hist2[0][y][ss]->Draw("e0p0");
      thestack2[y][ss]->Draw("same");
      hist2[0][y][ss]->Draw("e0p0same");
      for (int j = 1; j < NUMFILES; j++) hist2[j][y][ss]->Draw("e0p0same");

      subpad1f->cd();
      leg->Draw();
      mycanvas[y][ss]->Update();
      char outputimage2[128];
      std::string outhistname2 = (std::string)histnames2[y][ss];
      //sprintf(outputimage2,"./MHT_Analysis/CaloJets/%s.png",outhistname2.c_str());
      sprintf(outputimage2,"./MHT_Analysis/%s.png",outhistname2.c_str());
      mycanvas2[y][ss]->SaveAs(outputimage2);}}  
  //mycanvas2[y][ss]->SaveAs();}}  
}
