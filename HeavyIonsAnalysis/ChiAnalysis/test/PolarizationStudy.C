/////////////////////////////////////////////////
/// Code to study different polarization scenarios 
///////////////////////////////////////////////

/// Based on AcceptanceEfficiency, but removing constraints on muons - requiring only J/psi acceptance in the denominators

//// dealing with the prefilter that is present in the official MC: a new private MC with the same setting as the official, used to correct muon acceptance


#include <iostream>
#include <string.h>

#include "TSystem.h"
#include "TTree.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "TPave.h"
#include "TText.h"

#include "TROOT.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"

#include "TStyle.h"
#include "TLatex.h"
#include "TDirectory.h"
#include "TCollection.h"
#include "TPostScript.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "../../Macros/tdrstyle.C"
#include "../../Macros/CMS_lumi.C"
#include "../../Macros/ChiTreeInit.C"
#include "../../Macros/ChiFitterInit.h"


using namespace std;

//const int PythCode_chic0 = 10441; //Pythia codes
//const int PythCode_chic1 = 20443;
//const int PythCode_chic2 = 445;

double binsChiEffpT[] = { 0.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 15.0, 17.0, 20.0, 23.0, 26.0, 30.0 };
int  nbinsChiEffpT = sizeof(binsChiEffpT) / sizeof(double) - 1;

double binsConvEffpT[] = { 0.0, 0.3, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.7, 3.0, 3.5, 4.0, 4.5, 5.0, 7.0, 10.0 };
int  nbinsConvEffpT = sizeof(binsConvEffpT) / sizeof(double) - 1;
double binsConvEffy[] = { -5.0, -4.0, -3.0, -2.6, -2.5, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.5, 2.6, 3.0, 4.0, 5.0 };
int  nbinsConvEffy = sizeof(binsConvEffy) / sizeof(double) - 1;

int weird_decay_counter = 0;

double bins_Q_pT[] = { 6.5, 9, 12, 18, 30 };
int  nbins_Q_pT = sizeof(bins_Q_pT) / sizeof(double) - 1;

double bins_Q_y[] = { -2.4, -1.6, -1.0, 0, 1.0, 1.6, 2.4 };
int  nbins_Q_y = sizeof(bins_Q_y) / sizeof(double) - 1;

double binsWeightChi_pT[] = { 6.5, 9, 12, 18, 30 };
int  nbinsWeightChi_pT = sizeof(binsWeightChi_pT) / sizeof(double) - 1;
double binsWeightChi_absy[] = { 0.0, 1.0, 1.6, 2.4 };
int  nbinsWeightChi_absy = sizeof(binsWeightChi_absy) / sizeof(double) - 1;

double binsWeightChi_nTrk[] = { 0, 50, 100, 150, 250, 400 };
int  nbinsWeightChi_nTrk = sizeof(binsWeightChi_nTrk) / sizeof(double) - 1;


int PolarizationStudy(double lambdaTheta1 = 0.50, double lambdaTheta2 = -0.39, const char* fileOut = "Chic_PolarizationStudy_vTest-bothDir.root", int PhotSystIdx = 0, const char* fileInMC = "/eos/cms/store/group/phys_heavyions/okukral/Chi_c/Chi_c_pPb8TeV_MC-Official_v3-bothDir.root", const char* fileInMCNoFilter = "/eos/cms/store/group/phys_heavyions/okukral/Chi_c/Chi_c_pPb8TeV_MC_noFilter.root", const char* fileMCWeight = "MCWeight_v2.root")
{
	//int PhotSystIdx = 0;
	

	gStyle->SetOptStat(1111);
	gROOT->ProcessLine("#include <vector>");

	// load the MC weights
	TH1D* h_weightPVtrk;
	TFile* fMCWeight = new TFile(fileMCWeight, "READ");
	h_weightPVtrk= (TH1D*)fMCWeight->Get("h_weightPVtrk");
	h_weightPVtrk->SetDirectory(0);

	fMCWeight->Close();

	// Ota : Ton of histograms, but probably easier at the moment than adding any general holding structure

	/////////////////////////////////////
	// MUON ACCEPTANCE CORRECTIONS  //////
	////////////////////////////////////

	TH1D* h_chiCorrectionMuonAcceptance1D_pT_all_den = new TH1D("h_chiCorrectionMuonAcceptance1D_pT_all_den", "h_chiCorrectionMuonAcceptance1D_pT_all_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrectionMuonAcceptance1D_pT_all_num = new TH1D("h_chiCorrectionMuonAcceptance1D_pT_all_num", "h_chiCorrectionMuonAcceptance1D_pT_all_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrectionMuonAcceptance1D_pT_all_rat = new TH1D("h_chiCorrectionMuonAcceptance1D_pT_all_rat", "h_chiCorrectionMuonAcceptance1D_pT_all_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiCorrectionMuonAcceptance1D_y_den = new TH1D("h_chiCorrectionMuonAcceptance1D_y_den", "h_chiCorrectionMuonAcceptance1D_y_den; y", nbins_y, bins_y);
	TH1D* h_chiCorrectionMuonAcceptance1D_y_num = new TH1D("h_chiCorrectionMuonAcceptance1D_y_num", "h_chiCorrectionMuonAcceptance1D_y_num; y", nbins_y, bins_y);
	TH1D* h_chiCorrectionMuonAcceptance1D_y_rat = new TH1D("h_chiCorrectionMuonAcceptance1D_y_rat", "h_chiCorrectionMuonAcceptance1D_y_rat; y", nbins_y, bins_y);

	TH1D* h_chiCorrectionMuonAcceptance1D_nTrk_all_den = new TH1D("h_chiCorrectionMuonAcceptance1D_nTrk_all_den", "h_chiCorrectionMuonAcceptance1D_nTrk_all_den; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiCorrectionMuonAcceptance1D_nTrk_all_num = new TH1D("h_chiCorrectionMuonAcceptance1D_nTrk_all_num", "h_chiCorrectionMuonAcceptance1D_nTrk_all_num; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiCorrectionMuonAcceptance1D_nTrk_all_rat = new TH1D("h_chiCorrectionMuonAcceptance1D_nTrk_all_rat", "h_chiCorrectionMuonAcceptance1D_nTrk_all_rat; nTrk", nbins_nTrk, bins_nTrk);

	TH1D* h_chiCorrectionMuonAcceptance1D_pT_midCMS_den = new TH1D("h_chiCorrectionMuonAcceptance1D_pT_midCMS_den", "h_chiCorrectionMuonAcceptance1D_pT_midCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrectionMuonAcceptance1D_pT_midCMS_num = new TH1D("h_chiCorrectionMuonAcceptance1D_pT_midCMS_num", "h_chiCorrectionMuonAcceptance1D_pT_midCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrectionMuonAcceptance1D_pT_midCMS_rat = new TH1D("h_chiCorrectionMuonAcceptance1D_pT_midCMS_rat", "h_chiCorrectionMuonAcceptance1D_pT_midCMS_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiCorrectionMuonAcceptance1D_pT_fwdCMS_den = new TH1D("h_chiCorrectionMuonAcceptance1D_pT_fwdCMS_den", "h_chiCorrectionMuonAcceptance1D_pT_fwdCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrectionMuonAcceptance1D_pT_fwdCMS_num = new TH1D("h_chiCorrectionMuonAcceptance1D_pT_fwdCMS_num", "h_chiCorrectionMuonAcceptance1D_pT_fwdCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrectionMuonAcceptance1D_pT_fwdCMS_rat = new TH1D("h_chiCorrectionMuonAcceptance1D_pT_fwdCMS_rat", "h_chiCorrectionMuonAcceptance1D_pT_fwdCMS_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiCorrectionMuonAcceptance1D_pT_bkwCMS_den = new TH1D("h_chiCorrectionMuonAcceptance1D_pT_bkwCMS_den", "h_chiCorrectionMuonAcceptance1D_pT_bkwCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrectionMuonAcceptance1D_pT_bkwCMS_num = new TH1D("h_chiCorrectionMuonAcceptance1D_pT_bkwCMS_num", "h_chiCorrectionMuonAcceptance1D_pT_bkwCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrectionMuonAcceptance1D_pT_bkwCMS_rat = new TH1D("h_chiCorrectionMuonAcceptance1D_pT_bkwCMS_rat", "h_chiCorrectionMuonAcceptance1D_pT_bkwCMS_rat; p_{T}", nbins_pT, bins_pT);

	h_chiCorrectionMuonAcceptance1D_pT_all_den->Sumw2();
	h_chiCorrectionMuonAcceptance1D_y_den->Sumw2();
	h_chiCorrectionMuonAcceptance1D_nTrk_all_den->Sumw2();
	h_chiCorrectionMuonAcceptance1D_pT_midCMS_den->Sumw2();
	h_chiCorrectionMuonAcceptance1D_pT_fwdCMS_den->Sumw2();
	h_chiCorrectionMuonAcceptance1D_pT_bkwCMS_den->Sumw2();

	h_chiCorrectionMuonAcceptance1D_pT_all_num->Sumw2();
	h_chiCorrectionMuonAcceptance1D_y_num->Sumw2();
	h_chiCorrectionMuonAcceptance1D_nTrk_all_num->Sumw2();
	h_chiCorrectionMuonAcceptance1D_pT_midCMS_num->Sumw2();
	h_chiCorrectionMuonAcceptance1D_pT_fwdCMS_num->Sumw2();
	h_chiCorrectionMuonAcceptance1D_pT_bkwCMS_num->Sumw2();


	//  MUON ACCEPTANCE FOR RATIO (CHIC1/CHIC2) //////
	

	/////////////////HISTOGRAM OF CHIC1////////////////////////

	TH1D* h_chiCorrMuonAcceptanceChic1_1D_pT_all_den = new TH1D("h_chiCorrMuonAcceptanceChic1_1D_pT_all_den", "h_chiCorrMuonAcceptanceChic1_1D_pT_all_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMuonAcceptanceChic1_1D_pT_all_num = new TH1D("h_chiCorrMuonAcceptanceChic1_1D_pT_all_num", "h_chiCorrMuonAcceptanceChic1_1D_pT_all_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMuonAcceptanceChic1_1D_pT_all_rat = new TH1D("h_chiCorrMuonAcceptanceChic1_1D_pT_all_rat", "h_chiCorrMuonAcceptanceChic1_1D_pT_all_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiCorrMuonAcceptanceChic1_1D_y_den = new TH1D("h_chiCorrMuonAcceptanceChic1_1D_y_den", "h_chiCorrMuonAcceptanceChic1_1D_y_den; y", nbins_y, bins_y);
	TH1D* h_chiCorrMuonAcceptanceChic1_1D_y_num = new TH1D("h_chiCorrMuonAcceptanceChic1_1D_y_num", "h_chiCorrMuonAcceptanceChic1_1D_y_num; y", nbins_y, bins_y);
	TH1D* h_chiCorrMuonAcceptanceChic1_1D_y_rat = new TH1D("h_chiCorrMuonAcceptanceChic1_1D_y_rat", "h_chiCorrMuonAcceptanceChic1_1D_y_rat; y", nbins_y, bins_y);

	TH1D* h_chiCorrMuonAcceptanceChic1_1D_nTrk_all_den = new TH1D("h_chiCorrMuonAcceptanceChic1_1D_nTrk_all_den", "h_chiCorrMuonAcceptanceChic1_1D_nTrk_all_den; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiCorrMuonAcceptanceChic1_1D_nTrk_all_num = new TH1D("h_chiCorrMuonAcceptanceChic1_1D_nTrk_all_num", "h_chiCorrMuonAcceptanceChic1_1D_nTrk_all_num; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiCorrMuonAcceptanceChic1_1D_nTrk_all_rat = new TH1D("h_chiCorrMuonAcceptanceChic1_1D_nTrk_all_rat", "h_chiCorrMuonAcceptanceChic1_1D_nTrk_all_rat; nTrk", nbins_nTrk, bins_nTrk);

	TH1D* h_chiCorrMuonAcceptanceChic1_1D_pT_midCMS_den = new TH1D("h_chiCorrMuonAcceptanceChic1_1D_pT_midCMS_den", "h_chiCorrMuonAcceptanceChic1_1D_pT_midCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMuonAcceptanceChic1_1D_pT_midCMS_num = new TH1D("h_chiCorrMuonAcceptanceChic1_1D_pT_midCMS_num", "h_chiCorrMuonAcceptanceChic1_1D_pT_midCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMuonAcceptanceChic1_1D_pT_midCMS_rat = new TH1D("h_chiCorrMuonAcceptanceChic1_1D_pT_midCMS_rat", "h_chiCorrMuonAcceptanceChic1_1D_pT_midCMS_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiCorrMuonAcceptanceChic1_1D_pT_fwdCMS_den = new TH1D("h_chiCorrMuonAcceptanceChic1_1D_pT_fwdCMS_den", "h_chiCorrMuonAcceptanceChic1_1D_pT_fwdCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMuonAcceptanceChic1_1D_pT_fwdCMS_num = new TH1D("h_chiCorrMuonAcceptanceChic1_1D_pT_fwdCMS_num", "h_chiCorrMuonAcceptanceChic1_1D_pT_fwdCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMuonAcceptanceChic1_1D_pT_fwdCMS_rat = new TH1D("h_chiCorrMuonAcceptanceChic1_1D_pT_fwdCMS_rat", "h_chiCorrMuonAcceptanceChic1_1D_pT_fwdCMS_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiCorrMuonAcceptanceChic1_1D_pT_bkwCMS_den = new TH1D("h_chiCorrMuonAcceptanceChic1_1D_pT_bkwCMS_den", "h_chiCorrMuonAcceptanceChic1_1D_pT_bkwCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMuonAcceptanceChic1_1D_pT_bkwCMS_num = new TH1D("h_chiCorrMuonAcceptanceChic1_1D_pT_bkwCMS_num", "h_chiCorrMuonAcceptanceChic1_1D_pT_bkwCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMuonAcceptanceChic1_1D_pT_bkwCMS_rat = new TH1D("h_chiCorrMuonAcceptanceChic1_1D_pT_bkwCMS_rat", "h_chiCorrMuonAcceptanceChic1_1D_pT_bkwCMS_rat; p_{T}", nbins_pT, bins_pT);

	h_chiCorrMuonAcceptanceChic1_1D_pT_all_den->Sumw2();
	h_chiCorrMuonAcceptanceChic1_1D_y_den->Sumw2();
	h_chiCorrMuonAcceptanceChic1_1D_nTrk_all_den->Sumw2();
	h_chiCorrMuonAcceptanceChic1_1D_pT_midCMS_den->Sumw2();
	h_chiCorrMuonAcceptanceChic1_1D_pT_fwdCMS_den->Sumw2();
	h_chiCorrMuonAcceptanceChic1_1D_pT_bkwCMS_den->Sumw2();

	h_chiCorrMuonAcceptanceChic1_1D_pT_all_num->Sumw2();
	h_chiCorrMuonAcceptanceChic1_1D_y_num->Sumw2();
	h_chiCorrMuonAcceptanceChic1_1D_nTrk_all_num->Sumw2();
	h_chiCorrMuonAcceptanceChic1_1D_pT_midCMS_num->Sumw2();
	h_chiCorrMuonAcceptanceChic1_1D_pT_fwdCMS_num->Sumw2();
	h_chiCorrMuonAcceptanceChic1_1D_pT_bkwCMS_num->Sumw2();


	/////////////////HISTOGRAM OF CHIC2//////////////////////// 

	TH1D* h_chiCorrMuonAcceptanceChic2_1D_pT_all_den = new TH1D("h_chiCorrMuonAcceptanceChic2_1D_pT_all_den", "h_chiCorrMuonAcceptanceChic2_1D_pT_all_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMuonAcceptanceChic2_1D_pT_all_num = new TH1D("h_chiCorrMuonAcceptanceChic2_1D_pT_all_num", "h_chiCorrMuonAcceptanceChic2_1D_pT_all_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMuonAcceptanceChic2_1D_pT_all_rat = new TH1D("h_chiCorrMuonAcceptanceChic2_1D_pT_all_rat", "h_chiCorrMuonAcceptanceChic2_1D_pT_all_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiCorrMuonAcceptanceChic2_1D_y_den = new TH1D("h_chiCorrMuonAcceptanceChic2_1D_y_den", "h_chiCorrMuonAcceptanceChic2_1D_y_den; y", nbins_y, bins_y);
	TH1D* h_chiCorrMuonAcceptanceChic2_1D_y_num = new TH1D("h_chiCorrMuonAcceptanceChic2_1D_y_num", "h_chiCorrMuonAcceptanceChic2_1D_y_num; y", nbins_y, bins_y);
	TH1D* h_chiCorrMuonAcceptanceChic2_1D_y_rat = new TH1D("h_chiCorrMuonAcceptanceChic2_1D_y_rat", "h_chiCorrMuonAcceptanceChic2_1D_y_rat; y", nbins_y, bins_y);

	TH1D* h_chiCorrMuonAcceptanceChic2_1D_nTrk_all_den = new TH1D("h_chiCorrMuonAcceptanceChic2_1D_nTrk_all_den", "h_chiCorrMuonAcceptanceChic2_1D_nTrk_all_den; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiCorrMuonAcceptanceChic2_1D_nTrk_all_num = new TH1D("h_chiCorrMuonAcceptanceChic2_1D_nTrk_all_num", "h_chiCorrMuonAcceptanceChic2_1D_nTrk_all_num; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiCorrMuonAcceptanceChic2_1D_nTrk_all_rat = new TH1D("h_chiCorrMuonAcceptanceChic2_1D_nTrk_all_rat", "h_chiCorrMuonAcceptanceChic2_1D_nTrk_all_rat; nTrk", nbins_nTrk, bins_nTrk);

	TH1D* h_chiCorrMuonAcceptanceChic2_1D_pT_midCMS_den = new TH1D("h_chiCorrMuonAcceptanceChic2_1D_pT_midCMS_den", "h_chiCorrMuonAcceptanceChic2_1D_pT_midCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMuonAcceptanceChic2_1D_pT_midCMS_num = new TH1D("h_chiCorrMuonAcceptanceChic2_1D_pT_midCMS_num", "h_chiCorrMuonAcceptanceChic2_1D_pT_midCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMuonAcceptanceChic2_1D_pT_midCMS_rat = new TH1D("h_chiCorrMuonAcceptanceChic2_1D_pT_midCMS_rat", "h_chiCorrMuonAcceptanceChic2_1D_pT_midCMS_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiCorrMuonAcceptanceChic2_1D_pT_fwdCMS_den = new TH1D("h_chiCorrMuonAcceptanceChic2_1D_pT_fwdCMS_den", "h_chiCorrMuonAcceptanceChic2_1D_pT_fwdCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMuonAcceptanceChic2_1D_pT_fwdCMS_num = new TH1D("h_chiCorrMuonAcceptanceChic2_1D_pT_fwdCMS_num", "h_chiCorrMuonAcceptanceChic2_1D_pT_fwdCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMuonAcceptanceChic2_1D_pT_fwdCMS_rat = new TH1D("h_chiCorrMuonAcceptanceChic2_1D_pT_fwdCMS_rat", "h_chiCorrMuonAcceptanceChic2_1D_pT_fwdCMS_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiCorrMuonAcceptanceChic2_1D_pT_bkwCMS_den = new TH1D("h_chiCorrMuonAcceptanceChic2_1D_pT_bkwCMS_den", "h_chiCorrMuonAcceptanceChic2_1D_pT_bkwCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMuonAcceptanceChic2_1D_pT_bkwCMS_num = new TH1D("h_chiCorrMuonAcceptanceChic2_1D_pT_bkwCMS_num", "h_chiCorrMuonAcceptanceChic2_1D_pT_bkwCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMuonAcceptanceChic2_1D_pT_bkwCMS_rat = new TH1D("h_chiCorrMuonAcceptanceChic2_1D_pT_bkwCMS_rat", "h_chiCorrMuonAcceptanceChic2_1D_pT_bkwCMS_rat; p_{T}", nbins_pT, bins_pT);

	h_chiCorrMuonAcceptanceChic2_1D_pT_all_den->Sumw2();
	h_chiCorrMuonAcceptanceChic2_1D_y_den->Sumw2();
	h_chiCorrMuonAcceptanceChic2_1D_nTrk_all_den->Sumw2();
	h_chiCorrMuonAcceptanceChic2_1D_pT_midCMS_den->Sumw2();
	h_chiCorrMuonAcceptanceChic2_1D_pT_fwdCMS_den->Sumw2();
	h_chiCorrMuonAcceptanceChic2_1D_pT_bkwCMS_den->Sumw2();

	h_chiCorrMuonAcceptanceChic2_1D_pT_all_num->Sumw2();
	h_chiCorrMuonAcceptanceChic2_1D_y_num->Sumw2();
	h_chiCorrMuonAcceptanceChic2_1D_nTrk_all_num->Sumw2();
	h_chiCorrMuonAcceptanceChic2_1D_pT_midCMS_num->Sumw2();
	h_chiCorrMuonAcceptanceChic2_1D_pT_fwdCMS_num->Sumw2();
	h_chiCorrMuonAcceptanceChic2_1D_pT_bkwCMS_num->Sumw2();


	/////////////////HISTOGRAM OF CHIC1/CHIC2//////////////////////// 

	TH1D* h_chiCorrMuonAcceptanceChic1toChic2_1D_pT_all_rat = new TH1D("h_chiCorrMuonAcceptanceChic1toChic2_1D_pT_all_rat", "h_chiCorrMuonAcceptanceChic1toChic2_1D_pT_all_rat; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMuonAcceptanceChic1toChic2_1D_y_rat = new TH1D("h_chiCorrMuonAcceptanceChic1toChic2_1D_y_rat", "h_chiCorrMuonAcceptanceChic1toChic2_1D_y_rat; y", nbins_y, bins_y);
	TH1D* h_chiCorrMuonAcceptanceChic1toChic2_1D_nTrk_all_rat = new TH1D("h_chiCorrMuonAcceptanceChic1toChic2_1D_nTrk_all_rat", "h_chiCorrMuonAcceptanceChic1toChic2_1D_nTrk_all_rat; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiCorrMuonAcceptanceChic1toChic2_1D_pT_midCMS_rat = new TH1D("h_chiCorrMuonAcceptanceChic1toChic2_1D_pT_midCMS_rat", "h_chiCorrMuonAcceptanceChic1toChic2_1D_pT_midCMS_rat; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMuonAcceptanceChic1toChic2_1D_pT_fwdCMS_rat = new TH1D("h_chiCorrMuonAcceptanceChic1toChic2_1D_pT_fwdCMS_rat", "h_chiCorrMuonAcceptanceChic1toChic2_1D_pT_fwdCMS_rat; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMuonAcceptanceChic1toChic2_1D_pT_bkwCMS_rat = new TH1D("h_chiCorrMuonAcceptanceChic1toChic2_1D_pT_bkwCMS_rat", "h_chiCorrMuonAcceptanceChic1toChic2_1D_pT_bkwCMS_rat; p_{T}", nbins_pT, bins_pT);





	/////////////////////////////
	//  Major Part of  CORRECTION  //////
	///////////////////////////

	TH1D* h_chiCorrectionMajorPart1D_pT_all_den = new TH1D("h_chiCorrectionMajorPart1D_pT_all_den", "h_chiCorrectionMajorPart1D_pT_all_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrectionMajorPart1D_pT_all_num = new TH1D("h_chiCorrectionMajorPart1D_pT_all_num", "h_chiCorrectionMajorPart1D_pT_all_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrectionMajorPart1D_pT_all_rat = new TH1D("h_chiCorrectionMajorPart1D_pT_all_rat", "h_chiCorrectionMajorPart1D_pT_all_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiCorrectionMajorPart1D_y_den = new TH1D("h_chiCorrectionMajorPart1D_y_den", "h_chiCorrectionMajorPart1D_y_den; y", nbins_y, bins_y);
	TH1D* h_chiCorrectionMajorPart1D_y_num = new TH1D("h_chiCorrectionMajorPart1D_y_num", "h_chiCorrectionMajorPart1D_y_num; y", nbins_y, bins_y);
	TH1D* h_chiCorrectionMajorPart1D_y_rat = new TH1D("h_chiCorrectionMajorPart1D_y_rat", "h_chiCorrectionMajorPart1D_y_rat; y", nbins_y, bins_y);

	TH1D* h_chiCorrectionMajorPart1D_nTrk_all_den = new TH1D("h_chiCorrectionMajorPart1D_nTrk_all_den", "h_chiCorrectionMajorPart1D_nTrk_all_den; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiCorrectionMajorPart1D_nTrk_all_num = new TH1D("h_chiCorrectionMajorPart1D_nTrk_all_num", "h_chiCorrectionMajorPart1D_nTrk_all_num; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiCorrectionMajorPart1D_nTrk_all_rat = new TH1D("h_chiCorrectionMajorPart1D_nTrk_all_rat", "h_chiCorrectionMajorPart1D_nTrk_all_rat; nTrk", nbins_nTrk, bins_nTrk);

	TH1D* h_chiCorrectionMajorPart1D_pT_midCMS_den = new TH1D("h_chiCorrectionMajorPart1D_pT_midCMS_den", "h_chiCorrectionMajorPart1D_pT_midCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrectionMajorPart1D_pT_midCMS_num = new TH1D("h_chiCorrectionMajorPart1D_pT_midCMS_num", "h_chiCorrectionMajorPart1D_pT_midCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrectionMajorPart1D_pT_midCMS_rat = new TH1D("h_chiCorrectionMajorPart1D_pT_midCMS_rat", "h_chiCorrectionMajorPart1D_pT_midCMS_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiCorrectionMajorPart1D_pT_fwdCMS_den = new TH1D("h_chiCorrectionMajorPart1D_pT_fwdCMS_den", "h_chiCorrectionMajorPart1D_pT_fwdCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrectionMajorPart1D_pT_fwdCMS_num = new TH1D("h_chiCorrectionMajorPart1D_pT_fwdCMS_num", "h_chiCorrectionMajorPart1D_pT_fwdCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrectionMajorPart1D_pT_fwdCMS_rat = new TH1D("h_chiCorrectionMajorPart1D_pT_fwdCMS_rat", "h_chiCorrectionMajorPart1D_pT_fwdCMS_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiCorrectionMajorPart1D_pT_bkwCMS_den = new TH1D("h_chiCorrectionMajorPart1D_pT_bkwCMS_den", "h_chiCorrectionMajorPart1D_pT_bkwCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrectionMajorPart1D_pT_bkwCMS_num = new TH1D("h_chiCorrectionMajorPart1D_pT_bkwCMS_num", "h_chiCorrectionMajorPart1D_pT_bkwCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrectionMajorPart1D_pT_bkwCMS_rat = new TH1D("h_chiCorrectionMajorPart1D_pT_bkwCMS_rat", "h_chiCorrectionMajorPart1D_pT_bkwCMS_rat; p_{T}", nbins_pT, bins_pT);

	h_chiCorrectionMajorPart1D_pT_all_den->Sumw2();
	h_chiCorrectionMajorPart1D_y_den->Sumw2();
	h_chiCorrectionMajorPart1D_nTrk_all_den->Sumw2();
	h_chiCorrectionMajorPart1D_pT_midCMS_den->Sumw2();
	h_chiCorrectionMajorPart1D_pT_fwdCMS_den->Sumw2();
	h_chiCorrectionMajorPart1D_pT_bkwCMS_den->Sumw2();

	h_chiCorrectionMajorPart1D_pT_all_num->Sumw2();
	h_chiCorrectionMajorPart1D_y_num->Sumw2();
	h_chiCorrectionMajorPart1D_nTrk_all_num->Sumw2();
	h_chiCorrectionMajorPart1D_pT_midCMS_num->Sumw2();
	h_chiCorrectionMajorPart1D_pT_fwdCMS_num->Sumw2();
	h_chiCorrectionMajorPart1D_pT_bkwCMS_num->Sumw2();


	///////////////////////////////////////////////
	//  MAJOR PART CORRECTION RATIO (CHIC1/CHIC2) //////
	//////////////////////////////////////////////

	/////////////////HISTOGRAM OF CHIC1////////////////////////

	TH1D* h_chiCorrMajorPartChic1_1D_pT_all_den = new TH1D("h_chiCorrMajorPartChic1_1D_pT_all_den", "h_chiCorrMajorPartChic1_1D_pT_all_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMajorPartChic1_1D_pT_all_num = new TH1D("h_chiCorrMajorPartChic1_1D_pT_all_num", "h_chiCorrMajorPartChic1_1D_pT_all_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMajorPartChic1_1D_pT_all_rat = new TH1D("h_chiCorrMajorPartChic1_1D_pT_all_rat", "h_chiCorrMajorPartChic1_1D_pT_all_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiCorrMajorPartChic1_1D_y_den = new TH1D("h_chiCorrMajorPartChic1_1D_y_den", "h_chiCorrMajorPartChic1_1D_y_den; y", nbins_y, bins_y);
	TH1D* h_chiCorrMajorPartChic1_1D_y_num = new TH1D("h_chiCorrMajorPartChic1_1D_y_num", "h_chiCorrMajorPartChic1_1D_y_num; y", nbins_y, bins_y);
	TH1D* h_chiCorrMajorPartChic1_1D_y_rat = new TH1D("h_chiCorrMajorPartChic1_1D_y_rat", "h_chiCorrMajorPartChic1_1D_y_rat; y", nbins_y, bins_y);

	TH1D* h_chiCorrMajorPartChic1_1D_nTrk_all_den = new TH1D("h_chiCorrMajorPartChic1_1D_nTrk_all_den", "h_chiCorrMajorPartChic1_1D_nTrk_all_den; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiCorrMajorPartChic1_1D_nTrk_all_num = new TH1D("h_chiCorrMajorPartChic1_1D_nTrk_all_num", "h_chiCorrMajorPartChic1_1D_nTrk_all_num; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiCorrMajorPartChic1_1D_nTrk_all_rat = new TH1D("h_chiCorrMajorPartChic1_1D_nTrk_all_rat", "h_chiCorrMajorPartChic1_1D_nTrk_all_rat; nTrk", nbins_nTrk, bins_nTrk);

	TH1D* h_chiCorrMajorPartChic1_1D_pT_midCMS_den = new TH1D("h_chiCorrMajorPartChic1_1D_pT_midCMS_den", "h_chiCorrMajorPartChic1_1D_pT_midCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMajorPartChic1_1D_pT_midCMS_num = new TH1D("h_chiCorrMajorPartChic1_1D_pT_midCMS_num", "h_chiCorrMajorPartChic1_1D_pT_midCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMajorPartChic1_1D_pT_midCMS_rat = new TH1D("h_chiCorrMajorPartChic1_1D_pT_midCMS_rat", "h_chiCorrMajorPartChic1_1D_pT_midCMS_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiCorrMajorPartChic1_1D_pT_fwdCMS_den = new TH1D("h_chiCorrMajorPartChic1_1D_pT_fwdCMS_den", "h_chiCorrMajorPartChic1_1D_pT_fwdCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMajorPartChic1_1D_pT_fwdCMS_num = new TH1D("h_chiCorrMajorPartChic1_1D_pT_fwdCMS_num", "h_chiCorrMajorPartChic1_1D_pT_fwdCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMajorPartChic1_1D_pT_fwdCMS_rat = new TH1D("h_chiCorrMajorPartChic1_1D_pT_fwdCMS_rat", "h_chiCorrMajorPartChic1_1D_pT_fwdCMS_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiCorrMajorPartChic1_1D_pT_bkwCMS_den = new TH1D("h_chiCorrMajorPartChic1_1D_pT_bkwCMS_den", "h_chiCorrMajorPartChic1_1D_pT_bkwCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMajorPartChic1_1D_pT_bkwCMS_num = new TH1D("h_chiCorrMajorPartChic1_1D_pT_bkwCMS_num", "h_chiCorrMajorPartChic1_1D_pT_bkwCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMajorPartChic1_1D_pT_bkwCMS_rat = new TH1D("h_chiCorrMajorPartChic1_1D_pT_bkwCMS_rat", "h_chiCorrMajorPartChic1_1D_pT_bkwCMS_rat; p_{T}", nbins_pT, bins_pT);

	h_chiCorrMajorPartChic1_1D_pT_all_den->Sumw2();
	h_chiCorrMajorPartChic1_1D_y_den->Sumw2();
	h_chiCorrMajorPartChic1_1D_nTrk_all_den->Sumw2();
	h_chiCorrMajorPartChic1_1D_pT_midCMS_den->Sumw2();
	h_chiCorrMajorPartChic1_1D_pT_fwdCMS_den->Sumw2();
	h_chiCorrMajorPartChic1_1D_pT_bkwCMS_den->Sumw2();

	h_chiCorrMajorPartChic1_1D_pT_all_num->Sumw2();
	h_chiCorrMajorPartChic1_1D_y_num->Sumw2();
	h_chiCorrMajorPartChic1_1D_nTrk_all_num->Sumw2();
	h_chiCorrMajorPartChic1_1D_pT_midCMS_num->Sumw2();
	h_chiCorrMajorPartChic1_1D_pT_fwdCMS_num->Sumw2();
	h_chiCorrMajorPartChic1_1D_pT_bkwCMS_num->Sumw2();


	/////////////////HISTOGRAM OF CHIC2//////////////////////// 

	TH1D* h_chiCorrMajorPartChic2_1D_pT_all_den = new TH1D("h_chiCorrMajorPartChic2_1D_pT_all_den", "h_chiCorrMajorPartChic2_1D_pT_all_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMajorPartChic2_1D_pT_all_num = new TH1D("h_chiCorrMajorPartChic2_1D_pT_all_num", "h_chiCorrMajorPartChic2_1D_pT_all_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMajorPartChic2_1D_pT_all_rat = new TH1D("h_chiCorrMajorPartChic2_1D_pT_all_rat", "h_chiCorrMajorPartChic2_1D_pT_all_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiCorrMajorPartChic2_1D_y_den = new TH1D("h_chiCorrMajorPartChic2_1D_y_den", "h_chiCorrMajorPartChic2_1D_y_den; y", nbins_y, bins_y);
	TH1D* h_chiCorrMajorPartChic2_1D_y_num = new TH1D("h_chiCorrMajorPartChic2_1D_y_num", "h_chiCorrMajorPartChic2_1D_y_num; y", nbins_y, bins_y);
	TH1D* h_chiCorrMajorPartChic2_1D_y_rat = new TH1D("h_chiCorrMajorPartChic2_1D_y_rat", "h_chiCorrMajorPartChic2_1D_y_rat; y", nbins_y, bins_y);

	TH1D* h_chiCorrMajorPartChic2_1D_nTrk_all_den = new TH1D("h_chiCorrMajorPartChic2_1D_nTrk_all_den", "h_chiCorrMajorPartChic2_1D_nTrk_all_den; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiCorrMajorPartChic2_1D_nTrk_all_num = new TH1D("h_chiCorrMajorPartChic2_1D_nTrk_all_num", "h_chiCorrMajorPartChic2_1D_nTrk_all_num; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiCorrMajorPartChic2_1D_nTrk_all_rat = new TH1D("h_chiCorrMajorPartChic2_1D_nTrk_all_rat", "h_chiCorrMajorPartChic2_1D_nTrk_all_rat; nTrk", nbins_nTrk, bins_nTrk);

	TH1D* h_chiCorrMajorPartChic2_1D_pT_midCMS_den = new TH1D("h_chiCorrMajorPartChic2_1D_pT_midCMS_den", "h_chiCorrMajorPartChic2_1D_pT_midCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMajorPartChic2_1D_pT_midCMS_num = new TH1D("h_chiCorrMajorPartChic2_1D_pT_midCMS_num", "h_chiCorrMajorPartChic2_1D_pT_midCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMajorPartChic2_1D_pT_midCMS_rat = new TH1D("h_chiCorrMajorPartChic2_1D_pT_midCMS_rat", "h_chiCorrMajorPartChic2_1D_pT_midCMS_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiCorrMajorPartChic2_1D_pT_fwdCMS_den = new TH1D("h_chiCorrMajorPartChic2_1D_pT_fwdCMS_den", "h_chiCorrMajorPartChic2_1D_pT_fwdCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMajorPartChic2_1D_pT_fwdCMS_num = new TH1D("h_chiCorrMajorPartChic2_1D_pT_fwdCMS_num", "h_chiCorrMajorPartChic2_1D_pT_fwdCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMajorPartChic2_1D_pT_fwdCMS_rat = new TH1D("h_chiCorrMajorPartChic2_1D_pT_fwdCMS_rat", "h_chiCorrMajorPartChic2_1D_pT_fwdCMS_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiCorrMajorPartChic2_1D_pT_bkwCMS_den = new TH1D("h_chiCorrMajorPartChic2_1D_pT_bkwCMS_den", "h_chiCorrMajorPartChic2_1D_pT_bkwCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMajorPartChic2_1D_pT_bkwCMS_num = new TH1D("h_chiCorrMajorPartChic2_1D_pT_bkwCMS_num", "h_chiCorrMajorPartChic2_1D_pT_bkwCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMajorPartChic2_1D_pT_bkwCMS_rat = new TH1D("h_chiCorrMajorPartChic2_1D_pT_bkwCMS_rat", "h_chiCorrMajorPartChic2_1D_pT_bkwCMS_rat; p_{T}", nbins_pT, bins_pT);

	h_chiCorrMajorPartChic2_1D_pT_all_den->Sumw2();
	h_chiCorrMajorPartChic2_1D_y_den->Sumw2();
	h_chiCorrMajorPartChic2_1D_nTrk_all_den->Sumw2();
	h_chiCorrMajorPartChic2_1D_pT_midCMS_den->Sumw2();
	h_chiCorrMajorPartChic2_1D_pT_fwdCMS_den->Sumw2();
	h_chiCorrMajorPartChic2_1D_pT_bkwCMS_den->Sumw2();

	h_chiCorrMajorPartChic2_1D_pT_all_num->Sumw2();
	h_chiCorrMajorPartChic2_1D_y_num->Sumw2();
	h_chiCorrMajorPartChic2_1D_nTrk_all_num->Sumw2();
	h_chiCorrMajorPartChic2_1D_pT_midCMS_num->Sumw2();
	h_chiCorrMajorPartChic2_1D_pT_fwdCMS_num->Sumw2();
	h_chiCorrMajorPartChic2_1D_pT_bkwCMS_num->Sumw2();


	/////////////////HISTOGRAM OF CHIC1/CHIC2//////////////////////// 

	TH1D* h_chiCorrMajorPartChic1toChic2_1D_pT_all_rat = new TH1D("h_chiCorrMajorPartChic1toChic2_1D_pT_all_rat", "h_chiCorrMajorPartChic1toChic2_1D_pT_all_rat; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMajorPartChic1toChic2_1D_y_rat = new TH1D("h_chiCorrMajorPartChic1toChic2_1D_y_rat", "h_chiCorrMajorPartChic1toChic2_1D_y_rat; y", nbins_y, bins_y);
	TH1D* h_chiCorrMajorPartChic1toChic2_1D_nTrk_all_rat = new TH1D("h_chiCorrMajorPartChic1toChic2_1D_nTrk_all_rat", "h_chiCorrMajorPartChic1toChic2_1D_nTrk_all_rat; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiCorrMajorPartChic1toChic2_1D_pT_midCMS_rat = new TH1D("h_chiCorrMajorPartChic1toChic2_1D_pT_midCMS_rat", "h_chiCorrMajorPartChic1toChic2_1D_pT_midCMS_rat; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMajorPartChic1toChic2_1D_pT_fwdCMS_rat = new TH1D("h_chiCorrMajorPartChic1toChic2_1D_pT_fwdCMS_rat", "h_chiCorrMajorPartChic1toChic2_1D_pT_fwdCMS_rat; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiCorrMajorPartChic1toChic2_1D_pT_bkwCMS_rat = new TH1D("h_chiCorrMajorPartChic1toChic2_1D_pT_bkwCMS_rat", "h_chiCorrMajorPartChic1toChic2_1D_pT_bkwCMS_rat; p_{T}", nbins_pT, bins_pT);






















	/////////////////////////////
	//  TOTAL CORRECTION  //////
	///////////////////////////

	TH1D* h_chiTotalCorrection1D_pT_all_den = new TH1D("h_chiTotalCorrection1D_pT_all_den", "h_chiTotalCorrection1D_pT_all_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotalCorrection1D_pT_all_num = new TH1D("h_chiTotalCorrection1D_pT_all_num", "h_chiTotalCorrection1D_pT_all_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotalCorrection1D_pT_all_rat = new TH1D("h_chiTotalCorrection1D_pT_all_rat", "h_chiTotalCorrection1D_pT_all_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiTotalCorrection1D_y_den = new TH1D("h_chiTotalCorrection1D_y_den", "h_chiTotalCorrection1D_y_den; y", nbins_y, bins_y);
	TH1D* h_chiTotalCorrection1D_y_num = new TH1D("h_chiTotalCorrection1D_y_num", "h_chiTotalCorrection1D_y_num; y", nbins_y, bins_y);
	TH1D* h_chiTotalCorrection1D_y_rat = new TH1D("h_chiTotalCorrection1D_y_rat", "h_chiTotalCorrection1D_y_rat; y", nbins_y, bins_y);

	TH1D* h_chiTotalCorrection1D_nTrk_den = new TH1D("h_chiTotalCorrection1D_nTrk_den", "h_chiTotalCorrection1D_nTrk_den; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiTotalCorrection1D_nTrk_num = new TH1D("h_chiTotalCorrection1D_nTrk_num", "h_chiTotalCorrection1D_nTrk_num; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiTotalCorrection1D_nTrk_rat = new TH1D("h_chiTotalCorrection1D_nTrk_rat", "h_chiTotalCorrection1D_nTrk_rat; nTrk", nbins_nTrk, bins_nTrk);

	TH1D* h_chiTotalCorrection1D_nTrk_all_den = new TH1D("h_chiTotalCorrection1D_nTrk_all_den", "h_chiTotalCorrection1D_nTrk_all_den; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiTotalCorrection1D_nTrk_all_num = new TH1D("h_chiTotalCorrection1D_nTrk_all_num", "h_chiTotalCorrection1D_nTrk_all_num; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiTotalCorrection1D_nTrk_all_rat = new TH1D("h_chiTotalCorrection1D_nTrk_all_rat", "h_chiTotalCorrection1D_nTrk_all_rat; nTrk", nbins_nTrk, bins_nTrk);

	TH1D* h_chiTotalCorrection1D_pT_mid_den = new TH1D("h_chiTotalCorrection1D_pT_mid_den", "h_chiTotalCorrection1D_pT_mid_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotalCorrection1D_pT_mid_num = new TH1D("h_chiTotalCorrection1D_pT_mid_num", "h_chiTotalCorrection1D_pT_mid_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotalCorrection1D_pT_mid_rat = new TH1D("h_chiTotalCorrection1D_pT_mid_rat", "h_chiTotalCorrection1D_pT_mid_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiTotalCorrection1D_pT_fwd_den = new TH1D("h_chiTotalCorrection1D_pT_fwd_den", "h_chiTotalCorrection1D_pT_fwd_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotalCorrection1D_pT_fwd_num = new TH1D("h_chiTotalCorrection1D_pT_fwd_num", "h_chiTotalCorrection1D_pT_fwd_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotalCorrection1D_pT_fwd_rat = new TH1D("h_chiTotalCorrection1D_pT_fwd_rat", "h_chiTotalCorrection1D_pT_fwd_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiTotalCorrection1D_pT_fwdOnly_den = new TH1D("h_chiTotalCorrection1D_pT_fwdOnly_den", "h_chiTotalCorrection1D_pT_fwdOnly_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotalCorrection1D_pT_fwdOnly_num = new TH1D("h_chiTotalCorrection1D_pT_fwdOnly_num", "h_chiTotalCorrection1D_pT_fwdOnly_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotalCorrection1D_pT_fwdOnly_rat = new TH1D("h_chiTotalCorrection1D_pT_fwdOnly_rat", "h_chiTotalCorrection1D_pT_fwdOnly_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiTotalCorrection1D_pT_bkwOnly_den = new TH1D("h_chiTotalCorrection1D_pT_bkwOnly_den", "h_chiTotalCorrection1D_pT_bkwOnly_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotalCorrection1D_pT_bkwOnly_num = new TH1D("h_chiTotalCorrection1D_pT_bkwOnly_num", "h_chiTotalCorrection1D_pT_bkwOnly_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotalCorrection1D_pT_bkwOnly_rat = new TH1D("h_chiTotalCorrection1D_pT_bkwOnly_rat", "h_chiTotalCorrection1D_pT_bkwOnly_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiTotalCorrection1D_pT_fwdOnlyWide_den = new TH1D("h_chiTotalCorrection1D_pT_fwdOnlyWide_den", "h_chiTotalCorrection1D_pT_fwdOnlyWide_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotalCorrection1D_pT_fwdOnlyWide_num = new TH1D("h_chiTotalCorrection1D_pT_fwdOnlyWide_num", "h_chiTotalCorrection1D_pT_fwdOnlyWide_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotalCorrection1D_pT_fwdOnlyWide_rat = new TH1D("h_chiTotalCorrection1D_pT_fwdOnlyWide_rat", "h_chiTotalCorrection1D_pT_fwdOnlyWide_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiTotalCorrection1D_pT_bkwOnlyWide_den = new TH1D("h_chiTotalCorrection1D_pT_bkwOnlyWide_den", "h_chiTotalCorrection1D_pT_bkwOnlyWide_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotalCorrection1D_pT_bkwOnlyWide_num = new TH1D("h_chiTotalCorrection1D_pT_bkwOnlyWide_num", "h_chiTotalCorrection1D_pT_bkwOnlyWide_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotalCorrection1D_pT_bkwOnlyWide_rat = new TH1D("h_chiTotalCorrection1D_pT_bkwOnlyWide_rat", "h_chiTotalCorrection1D_pT_bkwOnlyWide_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiTotalCorrection1D_pT_midCMS_den = new TH1D("h_chiTotalCorrection1D_pT_midCMS_den", "h_chiTotalCorrection1D_pT_midCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotalCorrection1D_pT_midCMS_num = new TH1D("h_chiTotalCorrection1D_pT_midCMS_num", "h_chiTotalCorrection1D_pT_midCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotalCorrection1D_pT_midCMS_rat = new TH1D("h_chiTotalCorrection1D_pT_midCMS_rat", "h_chiTotalCorrection1D_pT_midCMS_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiTotalCorrection1D_pT_fwdCMS_den = new TH1D("h_chiTotalCorrection1D_pT_fwdCMS_den", "h_chiTotalCorrection1D_pT_fwdCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotalCorrection1D_pT_fwdCMS_num = new TH1D("h_chiTotalCorrection1D_pT_fwdCMS_num", "h_chiTotalCorrection1D_pT_fwdCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotalCorrection1D_pT_fwdCMS_rat = new TH1D("h_chiTotalCorrection1D_pT_fwdCMS_rat", "h_chiTotalCorrection1D_pT_fwdCMS_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiTotalCorrection1D_pT_bkwCMS_den = new TH1D("h_chiTotalCorrection1D_pT_bkwCMS_den", "h_chiTotalCorrection1D_pT_bkwCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotalCorrection1D_pT_bkwCMS_num = new TH1D("h_chiTotalCorrection1D_pT_bkwCMS_num", "h_chiTotalCorrection1D_pT_bkwCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotalCorrection1D_pT_bkwCMS_rat = new TH1D("h_chiTotalCorrection1D_pT_bkwCMS_rat", "h_chiTotalCorrection1D_pT_bkwCMS_rat; p_{T}", nbins_pT, bins_pT);

	h_chiTotalCorrection1D_pT_all_den->Sumw2();
	h_chiTotalCorrection1D_y_den->Sumw2();
	h_chiTotalCorrection1D_nTrk_den->Sumw2();
	h_chiTotalCorrection1D_nTrk_all_den->Sumw2();
	h_chiTotalCorrection1D_pT_mid_den->Sumw2();
	h_chiTotalCorrection1D_pT_fwd_den->Sumw2();
	h_chiTotalCorrection1D_pT_fwdOnly_den->Sumw2();
	h_chiTotalCorrection1D_pT_bkwOnly_den->Sumw2();
	h_chiTotalCorrection1D_pT_fwdOnlyWide_den->Sumw2();
	h_chiTotalCorrection1D_pT_bkwOnlyWide_den->Sumw2();
	h_chiTotalCorrection1D_pT_midCMS_den->Sumw2();
	h_chiTotalCorrection1D_pT_fwdCMS_den->Sumw2();
	h_chiTotalCorrection1D_pT_bkwCMS_den->Sumw2();

	h_chiTotalCorrection1D_pT_all_num->Sumw2();
	h_chiTotalCorrection1D_y_num->Sumw2();
	h_chiTotalCorrection1D_nTrk_num->Sumw2();
	h_chiTotalCorrection1D_nTrk_all_num->Sumw2();
	h_chiTotalCorrection1D_pT_mid_num->Sumw2();
	h_chiTotalCorrection1D_pT_fwd_num->Sumw2();
	h_chiTotalCorrection1D_pT_fwdOnly_num->Sumw2();
	h_chiTotalCorrection1D_pT_bkwOnly_num->Sumw2();
	h_chiTotalCorrection1D_pT_fwdOnlyWide_num->Sumw2();
	h_chiTotalCorrection1D_pT_bkwOnlyWide_num->Sumw2();
	h_chiTotalCorrection1D_pT_midCMS_num->Sumw2();
	h_chiTotalCorrection1D_pT_fwdCMS_num->Sumw2();
	h_chiTotalCorrection1D_pT_bkwCMS_num->Sumw2();


	///////////////////////////////////////////////
	//  TOTAL CORRECTION RATIO (CHIC1/CHIC2) //////
	//////////////////////////////////////////////

	/////////////////HISTOGRAM OF CHIC1////////////////////////

	TH1D* h_chiTotCorrChic1_1D_pT_all_den = new TH1D("h_chiTotCorrChic1_1D_pT_all_den", "h_chiTotCorrChic1_1D_pT_all_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotCorrChic1_1D_pT_all_num = new TH1D("h_chiTotCorrChic1_1D_pT_all_num", "h_chiTotCorrChic1_1D_pT_all_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotCorrChic1_1D_pT_all_rat = new TH1D("h_chiTotCorrChic1_1D_pT_all_rat", "h_chiTotCorrChic1_1D_pT_all_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiTotCorrChic1_1D_y_den = new TH1D("h_chiTotCorrChic1_1D_y_den", "h_chiTotCorrChic1_1D_y_den; y", nbins_y, bins_y);
	TH1D* h_chiTotCorrChic1_1D_y_num = new TH1D("h_chiTotCorrChic1_1D_y_num", "h_chiTotCorrChic1_1D_y_num; y", nbins_y, bins_y);
	TH1D* h_chiTotCorrChic1_1D_y_rat = new TH1D("h_chiTotCorrChic1_1D_y_rat", "h_chiTotCorrChic1_1D_y_rat; y", nbins_y, bins_y);

	TH1D* h_chiTotCorrChic1_1D_nTrk_all_den = new TH1D("h_chiTotCorrChic1_1D_nTrk_all_den", "h_chiTotCorrChic1_1D_nTrk_all_den; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiTotCorrChic1_1D_nTrk_all_num = new TH1D("h_chiTotCorrChic1_1D_nTrk_all_num", "h_chiTotCorrChic1_1D_nTrk_all_num; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiTotCorrChic1_1D_nTrk_all_rat = new TH1D("h_chiTotCorrChic1_1D_nTrk_all_rat", "h_chiTotCorrChic1_1D_nTrk_all_rat; nTrk", nbins_nTrk, bins_nTrk);

	TH1D* h_chiTotCorrChic1_1D_pT_midCMS_den = new TH1D("h_chiTotCorrChic1_1D_pT_midCMS_den", "h_chiTotCorrChic1_1D_pT_midCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotCorrChic1_1D_pT_midCMS_num = new TH1D("h_chiTotCorrChic1_1D_pT_midCMS_num", "h_chiTotCorrChic1_1D_pT_midCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotCorrChic1_1D_pT_midCMS_rat = new TH1D("h_chiTotCorrChic1_1D_pT_midCMS_rat", "h_chiTotCorrChic1_1D_pT_midCMS_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiTotCorrChic1_1D_pT_fwdCMS_den = new TH1D("h_chiTotCorrChic1_1D_pT_fwdCMS_den", "h_chiTotCorrChic1_1D_pT_fwdCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotCorrChic1_1D_pT_fwdCMS_num = new TH1D("h_chiTotCorrChic1_1D_pT_fwdCMS_num", "h_chiTotCorrChic1_1D_pT_fwdCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotCorrChic1_1D_pT_fwdCMS_rat = new TH1D("h_chiTotCorrChic1_1D_pT_fwdCMS_rat", "h_chiTotCorrChic1_1D_pT_fwdCMS_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiTotCorrChic1_1D_pT_bkwCMS_den = new TH1D("h_chiTotCorrChic1_1D_pT_bkwCMS_den", "h_chiTotCorrChic1_1D_pT_bkwCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotCorrChic1_1D_pT_bkwCMS_num = new TH1D("h_chiTotCorrChic1_1D_pT_bkwCMS_num", "h_chiTotCorrChic1_1D_pT_bkwCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotCorrChic1_1D_pT_bkwCMS_rat = new TH1D("h_chiTotCorrChic1_1D_pT_bkwCMS_rat", "h_chiTotCorrChic1_1D_pT_bkwCMS_rat; p_{T}", nbins_pT, bins_pT);

	h_chiTotCorrChic1_1D_pT_all_den->Sumw2();
	h_chiTotCorrChic1_1D_y_den->Sumw2();
	h_chiTotCorrChic1_1D_nTrk_all_den->Sumw2();
	h_chiTotCorrChic1_1D_pT_midCMS_den->Sumw2();
	h_chiTotCorrChic1_1D_pT_fwdCMS_den->Sumw2();
	h_chiTotCorrChic1_1D_pT_bkwCMS_den->Sumw2();

	h_chiTotCorrChic1_1D_pT_all_num->Sumw2();
	h_chiTotCorrChic1_1D_y_num->Sumw2();
	h_chiTotCorrChic1_1D_nTrk_all_num->Sumw2();
	h_chiTotCorrChic1_1D_pT_midCMS_num->Sumw2();
	h_chiTotCorrChic1_1D_pT_fwdCMS_num->Sumw2();
	h_chiTotCorrChic1_1D_pT_bkwCMS_num->Sumw2();


	/////////////////HISTOGRAM OF CHIC2//////////////////////// 

	TH1D* h_chiTotCorrChic2_1D_pT_all_den = new TH1D("h_chiTotCorrChic2_1D_pT_all_den", "h_chiTotCorrChic2_1D_pT_all_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotCorrChic2_1D_pT_all_num = new TH1D("h_chiTotCorrChic2_1D_pT_all_num", "h_chiTotCorrChic2_1D_pT_all_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotCorrChic2_1D_pT_all_rat = new TH1D("h_chiTotCorrChic2_1D_pT_all_rat", "h_chiTotCorrChic2_1D_pT_all_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiTotCorrChic2_1D_y_den = new TH1D("h_chiTotCorrChic2_1D_y_den", "h_chiTotCorrChic2_1D_y_den; y", nbins_y, bins_y);
	TH1D* h_chiTotCorrChic2_1D_y_num = new TH1D("h_chiTotCorrChic2_1D_y_num", "h_chiTotCorrChic2_1D_y_num; y", nbins_y, bins_y);
	TH1D* h_chiTotCorrChic2_1D_y_rat = new TH1D("h_chiTotCorrChic2_1D_y_rat", "h_chiTotCorrChic2_1D_y_rat; y", nbins_y, bins_y);

	TH1D* h_chiTotCorrChic2_1D_nTrk_all_den = new TH1D("h_chiTotCorrChic2_1D_nTrk_all_den", "h_chiTotCorrChic2_1D_nTrk_all_den; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiTotCorrChic2_1D_nTrk_all_num = new TH1D("h_chiTotCorrChic2_1D_nTrk_all_num", "h_chiTotCorrChic2_1D_nTrk_all_num; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiTotCorrChic2_1D_nTrk_all_rat = new TH1D("h_chiTotCorrChic2_1D_nTrk_all_rat", "h_chiTotCorrChic2_1D_nTrk_all_rat; nTrk", nbins_nTrk, bins_nTrk);

	TH1D* h_chiTotCorrChic2_1D_pT_midCMS_den = new TH1D("h_chiTotCorrChic2_1D_pT_midCMS_den", "h_chiTotCorrChic2_1D_pT_midCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotCorrChic2_1D_pT_midCMS_num = new TH1D("h_chiTotCorrChic2_1D_pT_midCMS_num", "h_chiTotCorrChic2_1D_pT_midCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotCorrChic2_1D_pT_midCMS_rat = new TH1D("h_chiTotCorrChic2_1D_pT_midCMS_rat", "h_chiTotCorrChic2_1D_pT_midCMS_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiTotCorrChic2_1D_pT_fwdCMS_den = new TH1D("h_chiTotCorrChic2_1D_pT_fwdCMS_den", "h_chiTotCorrChic2_1D_pT_fwdCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotCorrChic2_1D_pT_fwdCMS_num = new TH1D("h_chiTotCorrChic2_1D_pT_fwdCMS_num", "h_chiTotCorrChic2_1D_pT_fwdCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotCorrChic2_1D_pT_fwdCMS_rat = new TH1D("h_chiTotCorrChic2_1D_pT_fwdCMS_rat", "h_chiTotCorrChic2_1D_pT_fwdCMS_rat; p_{T}", nbins_pT, bins_pT);

	TH1D* h_chiTotCorrChic2_1D_pT_bkwCMS_den = new TH1D("h_chiTotCorrChic2_1D_pT_bkwCMS_den", "h_chiTotCorrChic2_1D_pT_bkwCMS_den; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotCorrChic2_1D_pT_bkwCMS_num = new TH1D("h_chiTotCorrChic2_1D_pT_bkwCMS_num", "h_chiTotCorrChic2_1D_pT_bkwCMS_num; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotCorrChic2_1D_pT_bkwCMS_rat = new TH1D("h_chiTotCorrChic2_1D_pT_bkwCMS_rat", "h_chiTotCorrChic2_1D_pT_bkwCMS_rat; p_{T}", nbins_pT, bins_pT);

	h_chiTotCorrChic2_1D_pT_all_den->Sumw2();
	h_chiTotCorrChic2_1D_y_den->Sumw2();
	h_chiTotCorrChic2_1D_nTrk_all_den->Sumw2();
	h_chiTotCorrChic2_1D_pT_midCMS_den->Sumw2();
	h_chiTotCorrChic2_1D_pT_fwdCMS_den->Sumw2();
	h_chiTotCorrChic2_1D_pT_bkwCMS_den->Sumw2();

	h_chiTotCorrChic2_1D_pT_all_num->Sumw2();
	h_chiTotCorrChic2_1D_y_num->Sumw2();
	h_chiTotCorrChic2_1D_nTrk_all_num->Sumw2();
	h_chiTotCorrChic2_1D_pT_midCMS_num->Sumw2();
	h_chiTotCorrChic2_1D_pT_fwdCMS_num->Sumw2();
	h_chiTotCorrChic2_1D_pT_bkwCMS_num->Sumw2();


	/////////////////HISTOGRAM OF CHIC1/CHIC2//////////////////////// 

	TH1D* h_chiTotCorrChic1toChic2_1D_pT_all_rat = new TH1D("h_chiTotCorrChic1toChic2_1D_pT_all_rat", "h_chiTotCorrChic1toChic2_1D_pT_all_rat; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotCorrChic1toChic2_1D_y_rat = new TH1D("h_chiTotCorrChic1toChic2_1D_y_rat", "h_chiTotCorrChic1toChic2_1D_y_rat; y", nbins_y, bins_y);
	TH1D* h_chiTotCorrChic1toChic2_1D_nTrk_all_rat = new TH1D("h_chiTotCorrChic1toChic2_1D_nTrk_all_rat", "h_chiTotCorrChic1toChic2_1D_nTrk_all_rat; nTrk", nbins_nTrk, bins_nTrk);
	TH1D* h_chiTotCorrChic1toChic2_1D_pT_midCMS_rat = new TH1D("h_chiTotCorrChic1toChic2_1D_pT_midCMS_rat", "h_chiTotCorrChic1toChic2_1D_pT_midCMS_rat; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotCorrChic1toChic2_1D_pT_fwdCMS_rat = new TH1D("h_chiTotCorrChic1toChic2_1D_pT_fwdCMS_rat", "h_chiTotCorrChic1toChic2_1D_pT_fwdCMS_rat; p_{T}", nbins_pT, bins_pT);
	TH1D* h_chiTotCorrChic1toChic2_1D_pT_bkwCMS_rat = new TH1D("h_chiTotCorrChic1toChic2_1D_pT_bkwCMS_rat", "h_chiTotCorrChic1toChic2_1D_pT_bkwCMS_rat; p_{T}", nbins_pT, bins_pT);




	////////////////////////////////////////////
	/////// STARTING RUNNING    ///////////////
	////////////////////////////////////////////


	//////////// First, we run over the small MC, to fill in the effects of the acceptance of the muon


	TFile* fMCNoFilter = new TFile(fileInMCNoFilter, "READ");


	TTree* event_treeNoFilter = (TTree*)fMCNoFilter->Get("ChiRootuple/event_tree");
	if (!event_treeNoFilter) {
		cout << "Problem with event Tree";
		return 1;
	}

	LoadChiBranches(event_treeNoFilter, true);

	Long64_t nentriesNoFilter = event_treeNoFilter->GetEntries();
	cout << "n entries in the tree without pre-filter: " << nentriesNoFilter << endl;
	//if (nentriesNoFilter > 50000) { nentriesNoFilter = 10000; }


	for (Long64_t i = 0; i < nentriesNoFilter; i++) {

		event_treeNoFilter->GetEntry(i);

		if (i % 10000 == 0) { cout << "event: " << i << " done: " << 100 * i / nentriesNoFilter << "%" << endl; }

		if (gen_Jpsi_pt->size() != 1 || gen_muon_pt->size() != 2 || gen_phot_pt->size() != 1 || gen_isGoodChicDecay->at(0) == false)
		{
			weird_decay_counter++;
			continue;
		}


		double eta1 = gen_muon_eta->at(0);
		double pt1 = gen_muon_pt->at(0);
		double eta2 = gen_muon_eta->at(1);
		double pt2 = gen_muon_pt->at(1);
		double eta_phot = gen_phot_eta->at(0);
		double pt_phot = gen_phot_pt->at(0);
		double eta_Jpsi = gen_Jpsi_eta->at(0);
		double pt_Jpsi = gen_Jpsi_pt->at(0);
		TLorentzVector* LVJpsi = (TLorentzVector*)gen_Jpsi_p4->At(0);
		double rap_Jpsi = LVJpsi->Rapidity();
		int muonPos1 = MuonMCMatched(0);
		int muonPos2 = MuonMCMatched(1);
		int JpsiPos = DimuonPassAllCutsMC(0);
		int matchJpsiPos = DimuonMCMatched(0); //includes those that don't pass the cuts, but reco exists
		int matchPositionConv = PhotMCMatched(0);
		int chicPos = ChiPassAllCutsMC(0);
		TLorentzVector* LVmuon1 = (TLorentzVector*)gen_muon_p4->At(0);

		bool JpsiInAcceptance = DimuonAcceptance(rap_Jpsi, pt_Jpsi);
		bool muon1InAcceptance = MuonAcceptance(eta1, pt1);
		bool muon2InAcceptance = MuonAcceptance(eta2, pt2);

		if (JpsiInAcceptance == true) // J/psi has to be in acceptance
		{

			double nTrack_inPV = ntracks_inEvent; //here use total ntrack, since unless there is reconstructed J/psi, we don't know what vertex to chose (total is the greatest PV)
			double MCweight = h_weightPVtrk->GetBinContent(h_weightPVtrk->FindBin(nTrack_inPV));
			double rap_Jpsi_corrections = rap_Jpsi; // rap_Jpsi not flipped, corrections yes
			if (ispPb == true) { //direction dependent stuff
				MCweight = 1.28*MCweight*WeightForMC_pTpart(pt_Jpsi)*WeightPhotonAcceptanceSystematic(pt_phot, PhotSystIdx);
			}
			else
			{
				MCweight = 0.72*MCweight*WeightForMC_pTpart(pt_Jpsi)*WeightPhotonAcceptanceSystematic(pt_phot, PhotSystIdx);
				rap_Jpsi_corrections = -rap_Jpsi_corrections;
			}

			///////////////////////////////////////
			// polarization reweighting based on MC
			///////////////////////////////////////
			double polarizedWeight = 1.0;

			// check pdgId : whether chic_1 or chic_2
			int pdgId = gen_pdgId->at(0);
			// get LorentzVector of positive charged muon from J/Psi decays
			TLorentzVector* LVmuonPositive = new TLorentzVector();
			double charge1 = gen_muon_charge->at(0);
			double charge2 = gen_muon_charge->at(1);
			if (charge1 > 0) {
				LVmuonPositive = (TLorentzVector*)gen_muon_p4->At(0);
			}
			else {
				LVmuonPositive = (TLorentzVector*)gen_muon_p4->At(1);
			}
			// get the polarization reweighting based on lambdaTheta for chic_1 or chic_2
			if (pdgId == PythCode_chic1) {
				polarizedWeight = PolarizationWeight(LVJpsi, LVmuonPositive, lambdaTheta1);
			}
			if (pdgId == PythCode_chic2) {
				polarizedWeight = PolarizationWeight(LVJpsi, LVmuonPositive, lambdaTheta2);
			}
			double MCweightPol = MCweight * polarizedWeight;


			///////////////////////////////////////
			// acceptance correction for chic/Jpsi ratio
			///////////////////////////////////////


			h_chiCorrectionMuonAcceptance1D_pT_all_den->Fill(pt_Jpsi, MCweightPol);
			h_chiCorrectionMuonAcceptance1D_y_den->Fill(rap_Jpsi, MCweightPol);
			h_chiCorrectionMuonAcceptance1D_nTrk_all_den->Fill(nTrack_inPV, MCweightPol);

			if (rap_Jpsi > rapCM_Edge2 && rap_Jpsi < rapCM_Edge3) { h_chiCorrectionMuonAcceptance1D_pT_midCMS_den->Fill(pt_Jpsi, MCweightPol); }
			if (rap_Jpsi > rapCM_Edge3 && rap_Jpsi < rapCM_Edge4) { h_chiCorrectionMuonAcceptance1D_pT_fwdCMS_den->Fill(pt_Jpsi, MCweightPol); }
			if (rap_Jpsi > rapCM_Edge1 && rap_Jpsi < rapCM_Edge2) { h_chiCorrectionMuonAcceptance1D_pT_bkwCMS_den->Fill(pt_Jpsi, MCweightPol); }

			if (muon1InAcceptance==true && muon2InAcceptance == true) //here, both muons in acceptance is the ratio we care about
			{
				h_chiCorrectionMuonAcceptance1D_pT_all_num->Fill(pt_Jpsi, MCweightPol);
				h_chiCorrectionMuonAcceptance1D_y_num->Fill(rap_Jpsi, MCweightPol);
				h_chiCorrectionMuonAcceptance1D_nTrk_all_num->Fill(nTrack_inPV, MCweightPol);

				if (rap_Jpsi > rapCM_Edge2 && rap_Jpsi < rapCM_Edge3) { h_chiCorrectionMuonAcceptance1D_pT_midCMS_num->Fill(pt_Jpsi, MCweightPol); }
				if (rap_Jpsi > rapCM_Edge3 && rap_Jpsi < rapCM_Edge4) { h_chiCorrectionMuonAcceptance1D_pT_fwdCMS_num->Fill(pt_Jpsi, MCweightPol); }
				if (rap_Jpsi > rapCM_Edge1 && rap_Jpsi < rapCM_Edge2) { h_chiCorrectionMuonAcceptance1D_pT_bkwCMS_num->Fill(pt_Jpsi, MCweightPol); }
			}


			//////////////////////////////////
			// Chic1 to Chic2 ratio 
			///////////////////////////////////
			// Chic1 
			if (gen_pdgId->at(0) == PythCode_chic1) {

				h_chiCorrMuonAcceptanceChic1_1D_pT_all_den->Fill(pt_Jpsi, MCweightPol);
				h_chiCorrMuonAcceptanceChic1_1D_y_den->Fill(rap_Jpsi, MCweightPol);
				h_chiCorrMuonAcceptanceChic1_1D_nTrk_all_den->Fill(nTrack_inPV, MCweightPol);
				if (rap_Jpsi > rapCM_Edge2 && rap_Jpsi < rapCM_Edge3) { h_chiCorrMuonAcceptanceChic1_1D_pT_midCMS_den->Fill(pt_Jpsi, MCweightPol); }
				if (rap_Jpsi > rapCM_Edge3 && rap_Jpsi < rapCM_Edge4) { h_chiCorrMuonAcceptanceChic1_1D_pT_fwdCMS_den->Fill(pt_Jpsi, MCweightPol); }
				if (rap_Jpsi > rapCM_Edge1 && rap_Jpsi < rapCM_Edge2) { h_chiCorrMuonAcceptanceChic1_1D_pT_bkwCMS_den->Fill(pt_Jpsi, MCweightPol); }

				if (muon1InAcceptance == true && muon2InAcceptance == true) {

					h_chiCorrMuonAcceptanceChic1_1D_pT_all_num->Fill(pt_Jpsi, MCweightPol);
					h_chiCorrMuonAcceptanceChic1_1D_y_num->Fill(rap_Jpsi, MCweightPol);
					h_chiCorrMuonAcceptanceChic1_1D_nTrk_all_num->Fill(nTrack_inPV, MCweightPol);
					if (rap_Jpsi > rapCM_Edge2 && rap_Jpsi < rapCM_Edge3) { h_chiCorrMuonAcceptanceChic1_1D_pT_midCMS_num->Fill(pt_Jpsi, MCweightPol); }
					if (rap_Jpsi > rapCM_Edge3 && rap_Jpsi < rapCM_Edge4) { h_chiCorrMuonAcceptanceChic1_1D_pT_fwdCMS_num->Fill(pt_Jpsi, MCweightPol); }
					if (rap_Jpsi > rapCM_Edge1 && rap_Jpsi < rapCM_Edge2) { h_chiCorrMuonAcceptanceChic1_1D_pT_bkwCMS_num->Fill(pt_Jpsi, MCweightPol); }

				}

			}
			//
			// Chic2
			if (gen_pdgId->at(0) == PythCode_chic2) {

				h_chiCorrMuonAcceptanceChic2_1D_pT_all_den->Fill(pt_Jpsi, MCweightPol);
				h_chiCorrMuonAcceptanceChic2_1D_y_den->Fill(rap_Jpsi, MCweightPol);
				h_chiCorrMuonAcceptanceChic2_1D_nTrk_all_den->Fill(nTrack_inPV, MCweightPol);
				if (rap_Jpsi > rapCM_Edge2 && rap_Jpsi < rapCM_Edge3) { h_chiCorrMuonAcceptanceChic2_1D_pT_midCMS_den->Fill(pt_Jpsi, MCweightPol); }
				if (rap_Jpsi > rapCM_Edge3 && rap_Jpsi < rapCM_Edge4) { h_chiCorrMuonAcceptanceChic2_1D_pT_fwdCMS_den->Fill(pt_Jpsi, MCweightPol); }
				if (rap_Jpsi > rapCM_Edge1 && rap_Jpsi < rapCM_Edge2) { h_chiCorrMuonAcceptanceChic2_1D_pT_bkwCMS_den->Fill(pt_Jpsi, MCweightPol); }

				if (muon1InAcceptance == true && muon2InAcceptance == true) {

					h_chiCorrMuonAcceptanceChic2_1D_pT_all_num->Fill(pt_Jpsi, MCweightPol);
					h_chiCorrMuonAcceptanceChic2_1D_y_num->Fill(rap_Jpsi, MCweightPol);
					h_chiCorrMuonAcceptanceChic2_1D_nTrk_all_num->Fill(nTrack_inPV, MCweightPol);
					if (rap_Jpsi > rapCM_Edge2 && rap_Jpsi < rapCM_Edge3) { h_chiCorrMuonAcceptanceChic2_1D_pT_midCMS_num->Fill(pt_Jpsi, MCweightPol); }
					if (rap_Jpsi > rapCM_Edge3 && rap_Jpsi < rapCM_Edge4) { h_chiCorrMuonAcceptanceChic2_1D_pT_fwdCMS_num->Fill(pt_Jpsi, MCweightPol); }
					if (rap_Jpsi > rapCM_Edge1 && rap_Jpsi < rapCM_Edge2) { h_chiCorrMuonAcceptanceChic2_1D_pT_bkwCMS_num->Fill(pt_Jpsi, MCweightPol); }

				}


			}

		}


	

	}

	fMCNoFilter->Close();

	// calculate the ratios for muon acceptance

	// quick muon acceptance:
	h_chiCorrectionMuonAcceptance1D_pT_all_rat->Divide(h_chiCorrectionMuonAcceptance1D_pT_all_num, h_chiCorrectionMuonAcceptance1D_pT_all_den, 1, 1, "B");
	h_chiCorrectionMuonAcceptance1D_y_rat->Divide(h_chiCorrectionMuonAcceptance1D_y_num, h_chiCorrectionMuonAcceptance1D_y_den, 1, 1, "B");
	h_chiCorrectionMuonAcceptance1D_nTrk_all_rat->Divide(h_chiCorrectionMuonAcceptance1D_nTrk_all_num, h_chiCorrectionMuonAcceptance1D_nTrk_all_den, 1, 1, "B");
	h_chiCorrectionMuonAcceptance1D_pT_midCMS_rat->Divide(h_chiCorrectionMuonAcceptance1D_pT_midCMS_num, h_chiCorrectionMuonAcceptance1D_pT_midCMS_den, 1, 1, "B");
	h_chiCorrectionMuonAcceptance1D_pT_fwdCMS_rat->Divide(h_chiCorrectionMuonAcceptance1D_pT_fwdCMS_num, h_chiCorrectionMuonAcceptance1D_pT_fwdCMS_den, 1, 1, "B");
	h_chiCorrectionMuonAcceptance1D_pT_bkwCMS_rat->Divide(h_chiCorrectionMuonAcceptance1D_pT_bkwCMS_num, h_chiCorrectionMuonAcceptance1D_pT_bkwCMS_den, 1, 1, "B");

	/////////Efficiency of Chic1///////////////////
	h_chiCorrMuonAcceptanceChic1_1D_pT_all_rat->Divide(h_chiCorrMuonAcceptanceChic1_1D_pT_all_num, h_chiCorrMuonAcceptanceChic1_1D_pT_all_den, 1, 1, "B");
	h_chiCorrMuonAcceptanceChic1_1D_y_rat->Divide(h_chiCorrMuonAcceptanceChic1_1D_y_num, h_chiCorrMuonAcceptanceChic1_1D_y_den, 1, 1, "B");
	h_chiCorrMuonAcceptanceChic1_1D_nTrk_all_rat->Divide(h_chiCorrMuonAcceptanceChic1_1D_nTrk_all_num, h_chiCorrMuonAcceptanceChic1_1D_nTrk_all_den, 1, 1, "B");
	h_chiCorrMuonAcceptanceChic1_1D_pT_midCMS_rat->Divide(h_chiCorrMuonAcceptanceChic1_1D_pT_midCMS_num, h_chiCorrMuonAcceptanceChic1_1D_pT_midCMS_den, 1, 1, "B");
	h_chiCorrMuonAcceptanceChic1_1D_pT_fwdCMS_rat->Divide(h_chiCorrMuonAcceptanceChic1_1D_pT_fwdCMS_num, h_chiCorrMuonAcceptanceChic1_1D_pT_fwdCMS_den, 1, 1, "B");
	h_chiCorrMuonAcceptanceChic1_1D_pT_bkwCMS_rat->Divide(h_chiCorrMuonAcceptanceChic1_1D_pT_bkwCMS_num, h_chiCorrMuonAcceptanceChic1_1D_pT_bkwCMS_den, 1, 1, "B");

	/////////Efficiency of Chic2/////////////////// 
	h_chiCorrMuonAcceptanceChic2_1D_pT_all_rat->Divide(h_chiCorrMuonAcceptanceChic2_1D_pT_all_num, h_chiCorrMuonAcceptanceChic2_1D_pT_all_den, 1, 1, "B");
	h_chiCorrMuonAcceptanceChic2_1D_y_rat->Divide(h_chiCorrMuonAcceptanceChic2_1D_y_num, h_chiCorrMuonAcceptanceChic2_1D_y_den, 1, 1, "B");
	h_chiCorrMuonAcceptanceChic2_1D_nTrk_all_rat->Divide(h_chiCorrMuonAcceptanceChic2_1D_nTrk_all_num, h_chiCorrMuonAcceptanceChic2_1D_nTrk_all_den, 1, 1, "B");
	h_chiCorrMuonAcceptanceChic2_1D_pT_midCMS_rat->Divide(h_chiCorrMuonAcceptanceChic2_1D_pT_midCMS_num, h_chiCorrMuonAcceptanceChic2_1D_pT_midCMS_den, 1, 1, "B");
	h_chiCorrMuonAcceptanceChic2_1D_pT_fwdCMS_rat->Divide(h_chiCorrMuonAcceptanceChic2_1D_pT_fwdCMS_num, h_chiCorrMuonAcceptanceChic2_1D_pT_fwdCMS_den, 1, 1, "B");
	h_chiCorrMuonAcceptanceChic2_1D_pT_bkwCMS_rat->Divide(h_chiCorrMuonAcceptanceChic2_1D_pT_bkwCMS_num, h_chiCorrMuonAcceptanceChic2_1D_pT_bkwCMS_den, 1, 1, "B");

	/////////Efficiency of Chic1/Chic2/////////////////// 
	h_chiCorrMuonAcceptanceChic1toChic2_1D_pT_all_rat->Divide(h_chiCorrMuonAcceptanceChic1_1D_pT_all_rat, h_chiCorrMuonAcceptanceChic2_1D_pT_all_rat, 1, 1);
	h_chiCorrMuonAcceptanceChic1toChic2_1D_y_rat->Divide(h_chiCorrMuonAcceptanceChic1_1D_y_rat, h_chiCorrMuonAcceptanceChic2_1D_y_rat, 1, 1);
	h_chiCorrMuonAcceptanceChic1toChic2_1D_nTrk_all_rat->Divide(h_chiCorrMuonAcceptanceChic1_1D_nTrk_all_rat, h_chiCorrMuonAcceptanceChic2_1D_nTrk_all_rat, 1, 1);
	h_chiCorrMuonAcceptanceChic1toChic2_1D_pT_midCMS_rat->Divide(h_chiCorrMuonAcceptanceChic1_1D_pT_midCMS_rat, h_chiCorrMuonAcceptanceChic2_1D_pT_midCMS_rat, 1, 1);
	h_chiCorrMuonAcceptanceChic1toChic2_1D_pT_fwdCMS_rat->Divide(h_chiCorrMuonAcceptanceChic1_1D_pT_fwdCMS_rat, h_chiCorrMuonAcceptanceChic2_1D_pT_fwdCMS_rat, 1, 1);
	h_chiCorrMuonAcceptanceChic1toChic2_1D_pT_bkwCMS_rat->Divide(h_chiCorrMuonAcceptanceChic1_1D_pT_bkwCMS_rat, h_chiCorrMuonAcceptanceChic2_1D_pT_bkwCMS_rat, 1, 1);
















	///////////////////////////////////////////////////////
	////  M A I N   M C   /////////////////////////////
	/////////////////////////////////////////////////////

	// Here we correct for polarization effects for everything else - "MAJOR PART"


	TFile* fMC = new TFile(fileInMC, "READ");


	TTree* event_tree = (TTree*)fMC->Get("ChiRootuple/event_tree");
	if (!event_tree) { 
		cout << "Problem with event Tree";
		return 1;
	}

	LoadChiBranches(event_tree, true);

	
	Long64_t nentries = event_tree->GetEntries();
	cout << "n entries: "<<nentries << endl;
	if (nentries > 50000) { nentries = 10000; }


	for (Long64_t i = 0; i < nentries; i++) {

		event_tree->GetEntry(i);

		if (i % 10000 == 0) { cout << "event: " << i << " done: " << 100 * i / nentries << "%" << endl; }

		//in the latest MC, a few chic decay weirdly (no phot, or no J/psi) - skip those
		if (gen_Jpsi_pt->size() != 1 || gen_muon_pt->size() != 2 || gen_phot_pt->size() != 1 || gen_isGoodChicDecay->at(0) == false)
		{
			weird_decay_counter++;
			continue;
		}


		double eta1 = gen_muon_eta->at(0);
		double pt1 = gen_muon_pt->at(0);
		double eta2 = gen_muon_eta->at(1);
		double pt2 = gen_muon_pt->at(1);
		double eta_phot = gen_phot_eta->at(0);
		double pt_phot = gen_phot_pt->at(0);
		double eta_Jpsi = gen_Jpsi_eta->at(0);
		double pt_Jpsi = gen_Jpsi_pt->at(0);
		TLorentzVector* LVJpsi = (TLorentzVector*)gen_Jpsi_p4->At(0);
		double rap_Jpsi = LVJpsi->Rapidity();
		int muonPos1 = MuonMCMatched(0);
		int muonPos2 = MuonMCMatched(1);
		int JpsiPos = DimuonPassAllCutsMC(0);
		int matchJpsiPos = DimuonMCMatched(0); //includes those that don't pass the cuts, but reco exists
		int matchPositionConv = PhotMCMatched(0);
		int chicPos = ChiPassAllCutsMC(0);
		TLorentzVector* LVmuon1 = (TLorentzVector*)gen_muon_p4->At(0);



		bool JpsiInAcceptance = DimuonAcceptance(rap_Jpsi, pt_Jpsi);  
		bool muon1InAcceptance = MuonAcceptance(eta1, pt1);
		bool muon2InAcceptance = MuonAcceptance(eta2, pt2);

		if (JpsiInAcceptance==true && muon1InAcceptance==true && muon2InAcceptance == true) // J/psi has to be in acceptance and muon also in acceptance (we fixed that in previous step)
		{
			double nTrack_inPV = ntracks_inEvent; //here use total ntrack, since unless there is reconstructed J/psi, we don't know what vertex to chose (total is the greatest PV)
			double MCweight = h_weightPVtrk->GetBinContent(h_weightPVtrk->FindBin(nTrack_inPV));
			double rap_Jpsi_corrections = rap_Jpsi; // rap_Jpsi not flipped, corrections yes
			if (ispPb == true) { //direction dependent stuff
				MCweight = 1.28*MCweight*WeightForMC_pTpart(pt_Jpsi)*WeightPhotonAcceptanceSystematic(pt_phot, PhotSystIdx);
			}
			else
			{
				MCweight = 0.72*MCweight*WeightForMC_pTpart(pt_Jpsi)*WeightPhotonAcceptanceSystematic(pt_phot, PhotSystIdx);
				rap_Jpsi_corrections = -rap_Jpsi_corrections;
			}

			///////////////////////////////////////
			// polarization reweighting based on MC
			///////////////////////////////////////
			double polarizedWeight = 1.0;

			// check pdgId : whether chic_1 or chic_2
			int pdgId = gen_pdgId->at(0);
			// get LorentzVector of positive charged muon from J/Psi decays
			TLorentzVector* LVmuonPositive = new TLorentzVector();
			double charge1 = gen_muon_charge->at(0);
			double charge2 = gen_muon_charge->at(1);
			if (charge1 > 0) {
				LVmuonPositive = (TLorentzVector*)gen_muon_p4->At(0);
			}
			else {
				LVmuonPositive = (TLorentzVector*)gen_muon_p4->At(1);
			}
			// get the polarization reweighting based on lambdaTheta for chic_1 or chic_2
			if (pdgId == PythCode_chic1) {
				polarizedWeight = PolarizationWeight(LVJpsi, LVmuonPositive, lambdaTheta1);
			}
			if (pdgId == PythCode_chic2) {
				polarizedWeight = PolarizationWeight(LVJpsi, LVmuonPositive, lambdaTheta2);
			}
			double MCweightPol = MCweight * polarizedWeight;

		///////////////////////////////////////
		// most of the correction for chic/Jpsi ratio
		///////////////////////////////////////
		
					   
			h_chiCorrectionMajorPart1D_pT_all_den->Fill(pt_Jpsi, MCweightPol);
			h_chiCorrectionMajorPart1D_y_den->Fill(rap_Jpsi, MCweightPol);
			h_chiCorrectionMajorPart1D_nTrk_all_den->Fill(nTrack_inPV, MCweightPol);

			if (rap_Jpsi > rapCM_Edge2 && rap_Jpsi < rapCM_Edge3) { h_chiCorrectionMajorPart1D_pT_midCMS_den->Fill(pt_Jpsi, MCweightPol); }
			if (rap_Jpsi > rapCM_Edge3 && rap_Jpsi < rapCM_Edge4) { h_chiCorrectionMajorPart1D_pT_fwdCMS_den->Fill(pt_Jpsi, MCweightPol); }
			if (rap_Jpsi > rapCM_Edge1 && rap_Jpsi < rapCM_Edge2) { h_chiCorrectionMajorPart1D_pT_bkwCMS_den->Fill(pt_Jpsi, MCweightPol); }
				
			if (chicPos > -1) //found chic -> numerator
			{
				h_chiCorrectionMajorPart1D_pT_all_num->Fill(pt_Jpsi, MCweightPol);
				h_chiCorrectionMajorPart1D_y_num->Fill(rap_Jpsi, MCweightPol);
				h_chiCorrectionMajorPart1D_nTrk_all_num->Fill(nTrack_inPV, MCweightPol);
	
				if (rap_Jpsi > rapCM_Edge2 && rap_Jpsi < rapCM_Edge3) { h_chiCorrectionMajorPart1D_pT_midCMS_num->Fill(pt_Jpsi, MCweightPol); }
				if (rap_Jpsi > rapCM_Edge3 && rap_Jpsi < rapCM_Edge4) { h_chiCorrectionMajorPart1D_pT_fwdCMS_num->Fill(pt_Jpsi, MCweightPol); }
				if (rap_Jpsi > rapCM_Edge1 && rap_Jpsi < rapCM_Edge2) { h_chiCorrectionMajorPart1D_pT_bkwCMS_num->Fill(pt_Jpsi, MCweightPol); }
			}


			//////////////////////////////////
			// Chic1 to Chic2 ratio 
			///////////////////////////////////
			// Chic1 
			if (gen_pdgId->at(0) == PythCode_chic1) {

				h_chiCorrMajorPartChic1_1D_pT_all_den->Fill(pt_Jpsi, MCweightPol);
				h_chiCorrMajorPartChic1_1D_y_den->Fill(rap_Jpsi, MCweightPol);
				h_chiCorrMajorPartChic1_1D_nTrk_all_den->Fill(nTrack_inPV, MCweightPol);
				if (rap_Jpsi > rapCM_Edge2 && rap_Jpsi < rapCM_Edge3) { h_chiCorrMajorPartChic1_1D_pT_midCMS_den->Fill(pt_Jpsi, MCweightPol); }
				if (rap_Jpsi > rapCM_Edge3 && rap_Jpsi < rapCM_Edge4) { h_chiCorrMajorPartChic1_1D_pT_fwdCMS_den->Fill(pt_Jpsi, MCweightPol); }
				if (rap_Jpsi > rapCM_Edge1 && rap_Jpsi < rapCM_Edge2) { h_chiCorrMajorPartChic1_1D_pT_bkwCMS_den->Fill(pt_Jpsi, MCweightPol); }

				if (chicPos > -1) {

					h_chiCorrMajorPartChic1_1D_pT_all_num->Fill(pt_Jpsi, MCweightPol);
					h_chiCorrMajorPartChic1_1D_y_num->Fill(rap_Jpsi, MCweightPol);
					h_chiCorrMajorPartChic1_1D_nTrk_all_num->Fill(nTrack_inPV, MCweightPol);
					if (rap_Jpsi > rapCM_Edge2 && rap_Jpsi < rapCM_Edge3) { h_chiCorrMajorPartChic1_1D_pT_midCMS_num->Fill(pt_Jpsi, MCweightPol); }
					if (rap_Jpsi > rapCM_Edge3 && rap_Jpsi < rapCM_Edge4) { h_chiCorrMajorPartChic1_1D_pT_fwdCMS_num->Fill(pt_Jpsi, MCweightPol); }
					if (rap_Jpsi > rapCM_Edge1 && rap_Jpsi < rapCM_Edge2) { h_chiCorrMajorPartChic1_1D_pT_bkwCMS_num->Fill(pt_Jpsi, MCweightPol); }

				}

			}
			//
			// Chic2
			if (gen_pdgId->at(0) == PythCode_chic2) {

				h_chiCorrMajorPartChic2_1D_pT_all_den->Fill(pt_Jpsi, MCweightPol);
				h_chiCorrMajorPartChic2_1D_y_den->Fill(rap_Jpsi, MCweightPol);
				h_chiCorrMajorPartChic2_1D_nTrk_all_den->Fill(nTrack_inPV, MCweightPol);
				if (rap_Jpsi > rapCM_Edge2 && rap_Jpsi < rapCM_Edge3) { h_chiCorrMajorPartChic2_1D_pT_midCMS_den->Fill(pt_Jpsi, MCweightPol); }
				if (rap_Jpsi > rapCM_Edge3 && rap_Jpsi < rapCM_Edge4) { h_chiCorrMajorPartChic2_1D_pT_fwdCMS_den->Fill(pt_Jpsi, MCweightPol); }
				if (rap_Jpsi > rapCM_Edge1 && rap_Jpsi < rapCM_Edge2) { h_chiCorrMajorPartChic2_1D_pT_bkwCMS_den->Fill(pt_Jpsi, MCweightPol); }

				if (chicPos > -1) {

					h_chiCorrMajorPartChic2_1D_pT_all_num->Fill(pt_Jpsi, MCweightPol);
					h_chiCorrMajorPartChic2_1D_y_num->Fill(rap_Jpsi, MCweightPol);
					h_chiCorrMajorPartChic2_1D_nTrk_all_num->Fill(nTrack_inPV, MCweightPol);
					if (rap_Jpsi > rapCM_Edge2 && rap_Jpsi < rapCM_Edge3) { h_chiCorrMajorPartChic2_1D_pT_midCMS_num->Fill(pt_Jpsi, MCweightPol); }
					if (rap_Jpsi > rapCM_Edge3 && rap_Jpsi < rapCM_Edge4) { h_chiCorrMajorPartChic2_1D_pT_fwdCMS_num->Fill(pt_Jpsi, MCweightPol); }
					if (rap_Jpsi > rapCM_Edge1 && rap_Jpsi < rapCM_Edge2) { h_chiCorrMajorPartChic2_1D_pT_bkwCMS_num->Fill(pt_Jpsi, MCweightPol); }

				}


			}
		
		}

		

	} //end of event loop

	cout << "weird decays: " << weird_decay_counter << "  out of  " << nentries << endl;


	// CALCULATE THE RATIOS

	// quick  correction:
	h_chiCorrectionMajorPart1D_pT_all_rat->Divide(h_chiCorrectionMajorPart1D_pT_all_num, h_chiCorrectionMajorPart1D_pT_all_den, 1, 1, "B");
	h_chiCorrectionMajorPart1D_y_rat->Divide(h_chiCorrectionMajorPart1D_y_num, h_chiCorrectionMajorPart1D_y_den, 1, 1, "B");
	h_chiCorrectionMajorPart1D_nTrk_all_rat->Divide(h_chiCorrectionMajorPart1D_nTrk_all_num, h_chiCorrectionMajorPart1D_nTrk_all_den, 1, 1, "B");
	h_chiCorrectionMajorPart1D_pT_midCMS_rat->Divide(h_chiCorrectionMajorPart1D_pT_midCMS_num, h_chiCorrectionMajorPart1D_pT_midCMS_den, 1, 1, "B");
	h_chiCorrectionMajorPart1D_pT_fwdCMS_rat->Divide(h_chiCorrectionMajorPart1D_pT_fwdCMS_num, h_chiCorrectionMajorPart1D_pT_fwdCMS_den, 1, 1, "B");
	h_chiCorrectionMajorPart1D_pT_bkwCMS_rat->Divide(h_chiCorrectionMajorPart1D_pT_bkwCMS_num, h_chiCorrectionMajorPart1D_pT_bkwCMS_den, 1, 1, "B");

	/////////Efficiency of Chic1///////////////////
	h_chiCorrMajorPartChic1_1D_pT_all_rat->Divide(h_chiCorrMajorPartChic1_1D_pT_all_num, h_chiCorrMajorPartChic1_1D_pT_all_den, 1, 1, "B");
	h_chiCorrMajorPartChic1_1D_y_rat->Divide(h_chiCorrMajorPartChic1_1D_y_num, h_chiCorrMajorPartChic1_1D_y_den, 1, 1, "B");
	h_chiCorrMajorPartChic1_1D_nTrk_all_rat->Divide(h_chiCorrMajorPartChic1_1D_nTrk_all_num, h_chiCorrMajorPartChic1_1D_nTrk_all_den, 1, 1, "B");
	h_chiCorrMajorPartChic1_1D_pT_midCMS_rat->Divide(h_chiCorrMajorPartChic1_1D_pT_midCMS_num, h_chiCorrMajorPartChic1_1D_pT_midCMS_den, 1, 1, "B");
	h_chiCorrMajorPartChic1_1D_pT_fwdCMS_rat->Divide(h_chiCorrMajorPartChic1_1D_pT_fwdCMS_num, h_chiCorrMajorPartChic1_1D_pT_fwdCMS_den, 1, 1, "B");
	h_chiCorrMajorPartChic1_1D_pT_bkwCMS_rat->Divide(h_chiCorrMajorPartChic1_1D_pT_bkwCMS_num, h_chiCorrMajorPartChic1_1D_pT_bkwCMS_den, 1, 1, "B");

	/////////Efficiency of Chic2/////////////////// 
	h_chiCorrMajorPartChic2_1D_pT_all_rat->Divide(h_chiCorrMajorPartChic2_1D_pT_all_num, h_chiCorrMajorPartChic2_1D_pT_all_den, 1, 1, "B");
	h_chiCorrMajorPartChic2_1D_y_rat->Divide(h_chiCorrMajorPartChic2_1D_y_num, h_chiCorrMajorPartChic2_1D_y_den, 1, 1, "B");
	h_chiCorrMajorPartChic2_1D_nTrk_all_rat->Divide(h_chiCorrMajorPartChic2_1D_nTrk_all_num, h_chiCorrMajorPartChic2_1D_nTrk_all_den, 1, 1, "B");
	h_chiCorrMajorPartChic2_1D_pT_midCMS_rat->Divide(h_chiCorrMajorPartChic2_1D_pT_midCMS_num, h_chiCorrMajorPartChic2_1D_pT_midCMS_den, 1, 1, "B");
	h_chiCorrMajorPartChic2_1D_pT_fwdCMS_rat->Divide(h_chiCorrMajorPartChic2_1D_pT_fwdCMS_num, h_chiCorrMajorPartChic2_1D_pT_fwdCMS_den, 1, 1, "B");
	h_chiCorrMajorPartChic2_1D_pT_bkwCMS_rat->Divide(h_chiCorrMajorPartChic2_1D_pT_bkwCMS_num, h_chiCorrMajorPartChic2_1D_pT_bkwCMS_den, 1, 1, "B");

	/////////Efficiency of Chic1/Chic2/////////////////// 
	h_chiCorrMajorPartChic1toChic2_1D_pT_all_rat->Divide(h_chiCorrMajorPartChic1_1D_pT_all_rat, h_chiCorrMajorPartChic2_1D_pT_all_rat, 1, 1);
	h_chiCorrMajorPartChic1toChic2_1D_y_rat->Divide(h_chiCorrMajorPartChic1_1D_y_rat, h_chiCorrMajorPartChic2_1D_y_rat, 1, 1);
	h_chiCorrMajorPartChic1toChic2_1D_nTrk_all_rat->Divide(h_chiCorrMajorPartChic1_1D_nTrk_all_rat, h_chiCorrMajorPartChic2_1D_nTrk_all_rat, 1, 1);
	h_chiCorrMajorPartChic1toChic2_1D_pT_midCMS_rat->Divide(h_chiCorrMajorPartChic1_1D_pT_midCMS_rat, h_chiCorrMajorPartChic2_1D_pT_midCMS_rat, 1, 1);
	h_chiCorrMajorPartChic1toChic2_1D_pT_fwdCMS_rat->Divide(h_chiCorrMajorPartChic1_1D_pT_fwdCMS_rat, h_chiCorrMajorPartChic2_1D_pT_fwdCMS_rat, 1, 1);
	h_chiCorrMajorPartChic1toChic2_1D_pT_bkwCMS_rat->Divide(h_chiCorrMajorPartChic1_1D_pT_bkwCMS_rat, h_chiCorrMajorPartChic2_1D_pT_bkwCMS_rat, 1, 1);




	/////////////////////////////////////////////////////////
	//  C O M B I N E   T H E   C O R R E C T I O N S   /////
	////////////////////////////////////////////////////////

	//Ota: add some loop to combine the histograms to obtain total correction
	// where the "total polarization correction" = "correction due to muon acceptance" * "correction due to the rest"
	// something like h_chiTotalCorrection1D_pT_all_rat = h_chiCorrectionMuonAcceptance1D_pT_all_rat * h_chiCorrectionMajorPart1D_pT_all_rat

	//Then, delete useless correction histograms (h_chiTotalCorrection1D_pT_all_den, h_chiTotalCorrection1D_pT_all_num)
	
	/// !!! WARNING - the small MC is non-embedded. You need to fix that for the ntrack plot - Just use integrated ratio across the board (since the muon acceptance does not depend on ntrack)


	//////////////////////////////////////////////////////////////
	// WRITE IN THE FILE
	///////////////////////

	TFile* fout = new TFile(fileOut, "RECREATE");

	



	//// Write the muon acceptance



	h_chiCorrectionMuonAcceptance1D_pT_all_den->Write();
	h_chiCorrectionMuonAcceptance1D_y_den->Write();
	h_chiCorrectionMuonAcceptance1D_nTrk_all_den->Write();
	h_chiCorrectionMuonAcceptance1D_pT_midCMS_den->Write();
	h_chiCorrectionMuonAcceptance1D_pT_fwdCMS_den->Write();
	h_chiCorrectionMuonAcceptance1D_pT_bkwCMS_den->Write();

	h_chiCorrectionMuonAcceptance1D_pT_all_num->Write();
	h_chiCorrectionMuonAcceptance1D_y_num->Write();
	h_chiCorrectionMuonAcceptance1D_nTrk_all_num->Write();
	h_chiCorrectionMuonAcceptance1D_pT_midCMS_num->Write();
	h_chiCorrectionMuonAcceptance1D_pT_fwdCMS_num->Write();
	h_chiCorrectionMuonAcceptance1D_pT_bkwCMS_num->Write();

	h_chiCorrectionMuonAcceptance1D_pT_all_rat->Write();
	h_chiCorrectionMuonAcceptance1D_y_rat->Write();
	h_chiCorrectionMuonAcceptance1D_nTrk_all_rat->Write();
	h_chiCorrectionMuonAcceptance1D_pT_midCMS_rat->Write();
	h_chiCorrectionMuonAcceptance1D_pT_fwdCMS_rat->Write();
	h_chiCorrectionMuonAcceptance1D_pT_bkwCMS_rat->Write();


	//Chic1/Chic2 Efficiency

	h_chiCorrMuonAcceptanceChic1_1D_pT_all_den->Write();
	h_chiCorrMuonAcceptanceChic1_1D_y_den->Write();
	h_chiCorrMuonAcceptanceChic1_1D_nTrk_all_den->Write();
	h_chiCorrMuonAcceptanceChic1_1D_pT_midCMS_den->Write();
	h_chiCorrMuonAcceptanceChic1_1D_pT_fwdCMS_den->Write();
	h_chiCorrMuonAcceptanceChic1_1D_pT_bkwCMS_den->Write();

	h_chiCorrMuonAcceptanceChic1_1D_pT_all_num->Write();
	h_chiCorrMuonAcceptanceChic1_1D_y_num->Write();
	h_chiCorrMuonAcceptanceChic1_1D_nTrk_all_num->Write();
	h_chiCorrMuonAcceptanceChic1_1D_pT_midCMS_num->Write();
	h_chiCorrMuonAcceptanceChic1_1D_pT_fwdCMS_num->Write();
	h_chiCorrMuonAcceptanceChic1_1D_pT_bkwCMS_num->Write();


	h_chiCorrMuonAcceptanceChic2_1D_pT_all_den->Write();
	h_chiCorrMuonAcceptanceChic2_1D_y_den->Write();
	h_chiCorrMuonAcceptanceChic2_1D_nTrk_all_den->Write();
	h_chiCorrMuonAcceptanceChic2_1D_pT_midCMS_den->Write();
	h_chiCorrMuonAcceptanceChic2_1D_pT_fwdCMS_den->Write();
	h_chiCorrMuonAcceptanceChic2_1D_pT_bkwCMS_den->Write();

	h_chiCorrMuonAcceptanceChic2_1D_pT_all_num->Write();
	h_chiCorrMuonAcceptanceChic2_1D_y_num->Write();
	h_chiCorrMuonAcceptanceChic2_1D_nTrk_all_num->Write();
	h_chiCorrMuonAcceptanceChic2_1D_pT_midCMS_num->Write();
	h_chiCorrMuonAcceptanceChic2_1D_pT_fwdCMS_num->Write();
	h_chiCorrMuonAcceptanceChic2_1D_pT_bkwCMS_num->Write();

	h_chiCorrMuonAcceptanceChic1_1D_pT_all_rat->Write();
	h_chiCorrMuonAcceptanceChic1_1D_y_rat->Write();
	h_chiCorrMuonAcceptanceChic1_1D_nTrk_all_rat->Write();
	h_chiCorrMuonAcceptanceChic1_1D_pT_midCMS_rat->Write();
	h_chiCorrMuonAcceptanceChic1_1D_pT_fwdCMS_rat->Write();
	h_chiCorrMuonAcceptanceChic1_1D_pT_bkwCMS_rat->Write();

	h_chiCorrMuonAcceptanceChic2_1D_pT_all_rat->Write();
	h_chiCorrMuonAcceptanceChic2_1D_y_rat->Write();
	h_chiCorrMuonAcceptanceChic2_1D_nTrk_all_rat->Write();
	h_chiCorrMuonAcceptanceChic2_1D_pT_midCMS_rat->Write();
	h_chiCorrMuonAcceptanceChic2_1D_pT_fwdCMS_rat->Write();
	h_chiCorrMuonAcceptanceChic2_1D_pT_bkwCMS_rat->Write();

	h_chiCorrMuonAcceptanceChic1toChic2_1D_pT_all_rat->Write();
	h_chiCorrMuonAcceptanceChic1toChic2_1D_y_rat->Write();
	h_chiCorrMuonAcceptanceChic1toChic2_1D_nTrk_all_rat->Write();
	h_chiCorrMuonAcceptanceChic1toChic2_1D_pT_midCMS_rat->Write();
	h_chiCorrMuonAcceptanceChic1toChic2_1D_pT_fwdCMS_rat->Write();
	h_chiCorrMuonAcceptanceChic1toChic2_1D_pT_bkwCMS_rat->Write();






	//// write the major part of the corrections


	h_chiCorrectionMajorPart1D_pT_all_den->Write();
	h_chiCorrectionMajorPart1D_y_den->Write();
	h_chiCorrectionMajorPart1D_nTrk_all_den->Write();
	h_chiCorrectionMajorPart1D_pT_midCMS_den->Write();
	h_chiCorrectionMajorPart1D_pT_fwdCMS_den->Write();
	h_chiCorrectionMajorPart1D_pT_bkwCMS_den->Write();

	h_chiCorrectionMajorPart1D_pT_all_num->Write();
	h_chiCorrectionMajorPart1D_y_num->Write();
	h_chiCorrectionMajorPart1D_nTrk_all_num->Write();
	h_chiCorrectionMajorPart1D_pT_midCMS_num->Write();
	h_chiCorrectionMajorPart1D_pT_fwdCMS_num->Write();
	h_chiCorrectionMajorPart1D_pT_bkwCMS_num->Write();

	h_chiCorrectionMajorPart1D_pT_all_rat->Write();
	h_chiCorrectionMajorPart1D_y_rat->Write();
	h_chiCorrectionMajorPart1D_nTrk_all_rat->Write();
	h_chiCorrectionMajorPart1D_pT_midCMS_rat->Write();
	h_chiCorrectionMajorPart1D_pT_fwdCMS_rat->Write();
	h_chiCorrectionMajorPart1D_pT_bkwCMS_rat->Write();


	h_chiCorrMajorPartChic1_1D_pT_all_den->Write();
	h_chiCorrMajorPartChic1_1D_y_den->Write();
	h_chiCorrMajorPartChic1_1D_nTrk_all_den->Write();
	h_chiCorrMajorPartChic1_1D_pT_midCMS_den->Write();
	h_chiCorrMajorPartChic1_1D_pT_fwdCMS_den->Write();
	h_chiCorrMajorPartChic1_1D_pT_bkwCMS_den->Write();

	h_chiCorrMajorPartChic1_1D_pT_all_num->Write();
	h_chiCorrMajorPartChic1_1D_y_num->Write();
	h_chiCorrMajorPartChic1_1D_nTrk_all_num->Write();
	h_chiCorrMajorPartChic1_1D_pT_midCMS_num->Write();
	h_chiCorrMajorPartChic1_1D_pT_fwdCMS_num->Write();
	h_chiCorrMajorPartChic1_1D_pT_bkwCMS_num->Write();


	h_chiCorrMajorPartChic2_1D_pT_all_den->Write();
	h_chiCorrMajorPartChic2_1D_y_den->Write();
	h_chiCorrMajorPartChic2_1D_nTrk_all_den->Write();
	h_chiCorrMajorPartChic2_1D_pT_midCMS_den->Write();
	h_chiCorrMajorPartChic2_1D_pT_fwdCMS_den->Write();
	h_chiCorrMajorPartChic2_1D_pT_bkwCMS_den->Write();

	h_chiCorrMajorPartChic2_1D_pT_all_num->Write();
	h_chiCorrMajorPartChic2_1D_y_num->Write();
	h_chiCorrMajorPartChic2_1D_nTrk_all_num->Write();
	h_chiCorrMajorPartChic2_1D_pT_midCMS_num->Write();
	h_chiCorrMajorPartChic2_1D_pT_fwdCMS_num->Write();
	h_chiCorrMajorPartChic2_1D_pT_bkwCMS_num->Write();

	h_chiCorrMajorPartChic1_1D_pT_all_rat->Write();
	h_chiCorrMajorPartChic1_1D_y_rat->Write();
	h_chiCorrMajorPartChic1_1D_nTrk_all_rat->Write();
	h_chiCorrMajorPartChic1_1D_pT_midCMS_rat->Write();
	h_chiCorrMajorPartChic1_1D_pT_fwdCMS_rat->Write();
	h_chiCorrMajorPartChic1_1D_pT_bkwCMS_rat->Write();

	h_chiCorrMajorPartChic2_1D_pT_all_rat->Write();
	h_chiCorrMajorPartChic2_1D_y_rat->Write();
	h_chiCorrMajorPartChic2_1D_nTrk_all_rat->Write();
	h_chiCorrMajorPartChic2_1D_pT_midCMS_rat->Write();
	h_chiCorrMajorPartChic2_1D_pT_fwdCMS_rat->Write();
	h_chiCorrMajorPartChic2_1D_pT_bkwCMS_rat->Write();

	h_chiCorrMajorPartChic1toChic2_1D_pT_all_rat->Write();
	h_chiCorrMajorPartChic1toChic2_1D_y_rat->Write();
	h_chiCorrMajorPartChic1toChic2_1D_nTrk_all_rat->Write();
	h_chiCorrMajorPartChic1toChic2_1D_pT_midCMS_rat->Write();
	h_chiCorrMajorPartChic1toChic2_1D_pT_fwdCMS_rat->Write();
	h_chiCorrMajorPartChic1toChic2_1D_pT_bkwCMS_rat->Write();








	/// write the final plots





	h_chiTotalCorrection1D_pT_all_den->Write();
	h_chiTotalCorrection1D_y_den->Write();
	h_chiTotalCorrection1D_nTrk_den->Write();
	h_chiTotalCorrection1D_nTrk_all_den->Write();
	h_chiTotalCorrection1D_pT_mid_den->Write();
	h_chiTotalCorrection1D_pT_fwd_den->Write();
	h_chiTotalCorrection1D_pT_fwdOnly_den->Write();
	h_chiTotalCorrection1D_pT_bkwOnly_den->Write();
	h_chiTotalCorrection1D_pT_fwdOnlyWide_den->Write();
	h_chiTotalCorrection1D_pT_bkwOnlyWide_den->Write();
	h_chiTotalCorrection1D_pT_midCMS_den->Write();
	h_chiTotalCorrection1D_pT_fwdCMS_den->Write();
	h_chiTotalCorrection1D_pT_bkwCMS_den->Write();

	h_chiTotalCorrection1D_pT_all_num->Write();
	h_chiTotalCorrection1D_y_num->Write();
	h_chiTotalCorrection1D_nTrk_num->Write();
	h_chiTotalCorrection1D_nTrk_all_num->Write();
	h_chiTotalCorrection1D_pT_mid_num->Write();
	h_chiTotalCorrection1D_pT_fwd_num->Write();
	h_chiTotalCorrection1D_pT_fwdOnly_num->Write();
	h_chiTotalCorrection1D_pT_bkwOnly_num->Write();
	h_chiTotalCorrection1D_pT_fwdOnlyWide_num->Write();
	h_chiTotalCorrection1D_pT_bkwOnlyWide_num->Write();
	h_chiTotalCorrection1D_pT_midCMS_num->Write();
	h_chiTotalCorrection1D_pT_fwdCMS_num->Write();
	h_chiTotalCorrection1D_pT_bkwCMS_num->Write();

	h_chiTotalCorrection1D_pT_all_rat->Write();
	h_chiTotalCorrection1D_y_rat->Write();
	h_chiTotalCorrection1D_nTrk_rat->Write();
	h_chiTotalCorrection1D_nTrk_all_rat->Write();
	h_chiTotalCorrection1D_pT_mid_rat->Write();
	h_chiTotalCorrection1D_pT_fwd_rat->Write();
	h_chiTotalCorrection1D_pT_fwdOnly_rat->Write();
	h_chiTotalCorrection1D_pT_bkwOnly_rat->Write();
	h_chiTotalCorrection1D_pT_fwdOnlyWide_rat->Write();
	h_chiTotalCorrection1D_pT_bkwOnlyWide_rat->Write();
	h_chiTotalCorrection1D_pT_midCMS_rat->Write();
	h_chiTotalCorrection1D_pT_fwdCMS_rat->Write();
	h_chiTotalCorrection1D_pT_bkwCMS_rat->Write();


	//write histograms for Chic1/Chic2 Efficiency

	h_chiTotCorrChic1_1D_pT_all_den->Write();
	h_chiTotCorrChic1_1D_y_den->Write();
	h_chiTotCorrChic1_1D_nTrk_all_den->Write();
	h_chiTotCorrChic1_1D_pT_midCMS_den->Write();
	h_chiTotCorrChic1_1D_pT_fwdCMS_den->Write();
	h_chiTotCorrChic1_1D_pT_bkwCMS_den->Write();

	h_chiTotCorrChic1_1D_pT_all_num->Write();
	h_chiTotCorrChic1_1D_y_num->Write();
	h_chiTotCorrChic1_1D_nTrk_all_num->Write();
	h_chiTotCorrChic1_1D_pT_midCMS_num->Write();
	h_chiTotCorrChic1_1D_pT_fwdCMS_num->Write();
	h_chiTotCorrChic1_1D_pT_bkwCMS_num->Write();


	h_chiTotCorrChic2_1D_pT_all_den->Write();
	h_chiTotCorrChic2_1D_y_den->Write();
	h_chiTotCorrChic2_1D_nTrk_all_den->Write();
	h_chiTotCorrChic2_1D_pT_midCMS_den->Write();
	h_chiTotCorrChic2_1D_pT_fwdCMS_den->Write();
	h_chiTotCorrChic2_1D_pT_bkwCMS_den->Write();

	h_chiTotCorrChic2_1D_pT_all_num->Write();
	h_chiTotCorrChic2_1D_y_num->Write();
	h_chiTotCorrChic2_1D_nTrk_all_num->Write();
	h_chiTotCorrChic2_1D_pT_midCMS_num->Write();
	h_chiTotCorrChic2_1D_pT_fwdCMS_num->Write();
	h_chiTotCorrChic2_1D_pT_bkwCMS_num->Write();

	h_chiTotCorrChic1_1D_pT_all_rat->Write();
	h_chiTotCorrChic1_1D_y_rat->Write();
	h_chiTotCorrChic1_1D_nTrk_all_rat->Write();
	h_chiTotCorrChic1_1D_pT_midCMS_rat->Write();
	h_chiTotCorrChic1_1D_pT_fwdCMS_rat->Write();
	h_chiTotCorrChic1_1D_pT_bkwCMS_rat->Write();

	h_chiTotCorrChic2_1D_pT_all_rat->Write();
	h_chiTotCorrChic2_1D_y_rat->Write();
	h_chiTotCorrChic2_1D_nTrk_all_rat->Write();
	h_chiTotCorrChic2_1D_pT_midCMS_rat->Write();
	h_chiTotCorrChic2_1D_pT_fwdCMS_rat->Write();
	h_chiTotCorrChic2_1D_pT_bkwCMS_rat->Write();

	h_chiTotCorrChic1toChic2_1D_pT_all_rat->Write();
	h_chiTotCorrChic1toChic2_1D_y_rat->Write();
	h_chiTotCorrChic1toChic2_1D_nTrk_all_rat->Write();
	h_chiTotCorrChic1toChic2_1D_pT_midCMS_rat->Write();
	h_chiTotCorrChic1toChic2_1D_pT_fwdCMS_rat->Write();
	h_chiTotCorrChic1toChic2_1D_pT_bkwCMS_rat->Write();




	fout->Close();
	fMC->Close();
	cout << "END" << endl;
	return 0;
}
