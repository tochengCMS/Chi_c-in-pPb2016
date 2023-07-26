#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include <iostream>
#include "Riostream.h"
#include <math.h>
#include "TStyle.h"
#include "TCanvas.h"
#include <TROOT.h>
#include <cmath>
#include <TMatrixTSym.h>
#include <TGraphAsymmErrors.h>

#include "tdrstyle.C"
#include "CMS_lumi.C"
#include "../ChiFitterInit.h"
#include "plottingHelper.C"

const TString sEffName= "Nominal";
const float cLowX = 0, cHighX = 30;
const TString sNameTag = "_pT_CMS";
double systErrorWidth = 0.5;
double pointXShift = 0.3;
//const float cLowX = -2.4, cHighX = 2.4;
//const TString sNameTag = "_y";
//const float cLowX = 0.0, cHighX = 250;
//const TString sNameTag = "_nTrk";

const float _markerSize = 2.4;


int plotNominal_pTCMS()
{
	
	gStyle->SetOptStat(0);
	setTDRStyle();
	gStyle->SetEndErrorSize(8);
	gStyle->SetPadRightMargin(0.03);

	//TCanvas* cankres1 = new TCanvas("cankres1", "Canvas with results1", 640, 480);
	TCanvas* cankres1 = new TCanvas("cankres1", "Canvas with results1", 800, 600);
	cankres1->SetBottomMargin(0.145);
	TH1F* hframe = new TH1F("hframe", "", 1, cLowX, cHighX);
	hframe->Draw();
	hframe->GetYaxis()->SetRangeUser(0.0, 0.45);
	hframe->GetXaxis()->SetTitle("p_{T} (J/#psi) [GeV/c]");
	//hframe->GetXaxis()->SetTitle("y (J/#psi)");
	//hframe->GetXaxis()->SetTitle("N_{tracks}");
	//hframe->GetXaxis()->SetTitleSize(0.05);
	hframe->GetXaxis()->SetTitleOffset(1.03);
	hframe->GetYaxis()->SetTitle("(#chi_{c1}+#chi_{c2}) / J/#psi");
	hframe->GetYaxis()->SetTitleSize(0.06);
	//hframe->GetYaxis()->SetTitleOffset(1.12);







	//TF1* one = new TF1("one", "1", cLowX, cHighX);
	//one->SetLineStyle(2);
	//one->SetLineWidth(2);
	//one->SetLineColor(kBlack);

	// Load data

	TFile* myFile1 = new TFile("Chi_c_output_Nominal_vDissertation-bothDir_DCB_NewBins.root", "READ");
	//TCanvas* c1 = (TCanvas*)myFile1->Get("tpTreeTrk/Trk_ntracksdep/fit_eff_plots/tag_hiBin_PLOT_TrackCuts_true"); //Get canvas with the plot
	//TGraphAsymmErrors* gAS_PbPb= (TGraphAsymmErrors*)c1->GetListOfPrimitives()->FindObject("hxy_fit_eff");

	//TGraphAsymmErrors* gAS_Result = (TGraphAsymmErrors*)myFile1->Get("gAsRatio_pT");
	//TGraphAsymmErrors* gAS_Result2 = (TGraphAsymmErrors*)myFile1->Get("gAsRatio_pT_mid");
	//TGraphAsymmErrors* gAS_Result3 = (TGraphAsymmErrors*)myFile1->Get("gAsRatio_pT_fwd");

	TGraphAsymmErrors* gAS_Result = (TGraphAsymmErrors*)myFile1->Get("gAsRatio_pT_midCMS");
	TGraphAsymmErrors* gAS_Result2 = (TGraphAsymmErrors*)myFile1->Get("gAsRatio_pT_fwdCMS");
	TGraphAsymmErrors* gAS_Result3 = (TGraphAsymmErrors*)myFile1->Get("gAsRatio_pT_bkwCMS");

	//TGraphAsymmErrors* gAS_Result = (TGraphAsymmErrors*)myFile1->Get("gAsRatio_y");
	//TGraphAsymmErrors* gAS_Result = (TGraphAsymmErrors*)myFile1->Get("gAsRatio_nTrk");
	//TGraphAsymmErrors* gAS_Result2 = (TGraphAsymmErrors*)myFile1->Get("gAsRatio_nTrk_all");



	/////////// SYSTEMATICS

	TFile* myFile2 = new TFile("TotalSyst_vDissertation.root", "READ");

	TGraph* g_Systematics = (TGraph*)myFile2->Get("gSystOutputArrayTotal_pt_midCMS");
	TGraph* g_Systematics2 = (TGraph*)myFile2->Get("gSystOutputArrayTotal_pt_fwdCMS");
	TGraph* g_Systematics3 = (TGraph*)myFile2->Get("gSystOutputArrayTotal_pt_bkwCMS");

	TGraphAsymmErrors* gAS_Result_Syst = (TGraphAsymmErrors*)gAS_Result->Clone();
	TGraphAsymmErrors* gAS_Result_Syst2 = (TGraphAsymmErrors*)gAS_Result2->Clone();
	TGraphAsymmErrors* gAS_Result_Syst3 = (TGraphAsymmErrors*)gAS_Result3->Clone();

	PrepareSystPlotting(gAS_Result_Syst, gAS_Result, g_Systematics, systErrorWidth);
	PrepareSystPlotting(gAS_Result_Syst2, gAS_Result2, g_Systematics2, systErrorWidth);
	PrepareSystPlotting(gAS_Result_Syst3, gAS_Result3, g_Systematics3, systErrorWidth);


	cout << "Loading done" << endl;



	//Setup

	////kill x errors
	//RemoveXError(gAS_pp);
	//RemoveXError(gAS_pPb);
	//RemoveXError(gAS_PbPb);

	ShiftXPosition(gAS_Result2, -pointXShift);
	ShiftXPosition(gAS_Result3, pointXShift);
	ShiftXPosition(gAS_Result_Syst2, -pointXShift);
	ShiftXPosition(gAS_Result_Syst3, pointXShift);

	gAS_Result_Syst->SetLineWidth(2);
	gAS_Result_Syst->SetLineColor(kOrange + 1);
	gAS_Result_Syst->SetMarkerSize(_markerSize);
	gAS_Result_Syst->SetFillColor(kOrange + 1);
	//gAS_Result_Syst->SetFillStyle(3002);//3002 looks good as png, but poor in pdf
	gAS_Result_Syst->SetFillStyle(3013);

	gAS_Result_Syst2->SetLineWidth(2);
	gAS_Result_Syst2->SetLineColor(kBlue);
	gAS_Result_Syst2->SetMarkerSize(_markerSize);
	gAS_Result_Syst2->SetFillColor(kBlue);
	//gAS_Result_Syst2->SetFillStyle(3002);//3002 looks good as png, but poor in pdf
	gAS_Result_Syst2->SetFillStyle(3013);

	gAS_Result_Syst3->SetLineWidth(2);
	gAS_Result_Syst3->SetLineColor(kGreen + 3);
	gAS_Result_Syst3->SetMarkerSize(_markerSize);
	gAS_Result_Syst3->SetFillColor(kGreen + 3);
	//gAS_Result_Syst3->SetFillStyle(3002);
	gAS_Result_Syst3->SetFillStyle(3013);






	gAS_Result->SetMarkerSize(0.9*_markerSize);
	gAS_Result->SetMarkerColor(kOrange + 1);
	gAS_Result->SetMarkerStyle(71);
	gAS_Result->SetLineColor(kOrange + 1);
	gAS_Result->SetLineWidth(2);

	gAS_Result2->SetMarkerSize(0.9*_markerSize);
	gAS_Result2->SetMarkerColor(kBlue);
	gAS_Result2->SetMarkerStyle(72);
	gAS_Result2->SetLineColor(kBlue);
	gAS_Result2->SetLineWidth(2);

	gAS_Result3->SetMarkerSize(1.1*_markerSize);
	gAS_Result3->SetMarkerColor(kGreen+3);
	gAS_Result3->SetMarkerStyle(76);
	gAS_Result3->SetLineColor(kGreen+3);
	gAS_Result3->SetLineWidth(2);






	cankres1->cd();
	//cankres1->SetLogx();
	//one->Draw("Same");

	//gAS_Result_Syst->Draw("[]");
	//gAS_Result_Syst2->Draw("[]");
	//gAS_Result_Syst3->Draw("[]");

	//gAS_Result->Draw("ZPE");
	//gAS_Result2->Draw("ZPE");
	//gAS_Result3->Draw("ZPE");

	gAS_Result_Syst->Draw("2");
	gAS_Result_Syst2->Draw("2");
	gAS_Result_Syst3->Draw("2");

	gAS_Result->Draw("ZPE");
	gAS_Result2->Draw("ZPE");
	gAS_Result3->Draw("ZPE");

	// add CMS text
	//CMS_lumi(cankres1, 0, 22); //mid top
	CMS_lumi(cankres1, 0, 10); //left top

	//TLegend*leg = new TLegend(0.50, 0.20, 0.88, 0.45, sEffName+" efficiency");
	//TLegend*leg = new TLegend(0.58, 0.18, 0.88, 0.4, "");
	TLegend*leg = new TLegend(0.45, 0.19, 0.88, 0.39, "");
	leg->SetFillColor(kWhite);
	leg->SetBorderSize(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.045);

	//leg->AddEntry(gAS_Result, "6.5<p_{T}<30 GeV/c", "p");
	//leg->AddEntry(gAS_Result, "|y|<1.0, 6.5<p_{T}<30 GeV/c", "p");
	//leg->AddEntry(gAS_Result2, "|y|<2.4, 6.5<p_{T}<30 GeV/c", "p");

	//leg->AddEntry(gAS_Result, "Midrapidity: |y|<1.0", "p");
	//leg->AddEntry(gAS_Result2, "Forward: 1.6<y<2.4", "p");
	//leg->AddEntry(gAS_Result3, "Backward: -2.4<y<-1.6", "p");

	leg->AddEntry(gAS_Result3, "Backward:   -2.0<y_{CM}<-1.0", "p");
	leg->AddEntry(gAS_Result,  "Midrapidity: -1.0<y_{CM}<1.0", "p");
	leg->AddEntry(gAS_Result2, "Forward:       1.0<y_{CM}<1.9", "p");

	//leg->AddEntry(gAS_Result, "|y|<2.4", "p");
	//leg->AddEntry(gAS_Result2, "|y|<1.0", "p");
	//leg->AddEntry(gAS_Result3, "1.0<|y|<2.4", "p");
	leg->Draw("same");

	//TPaveText* pText1 = new TPaveText(0.62, 0.4, 0.88, 0.48, "NDC NB");
	//pText1->AddText("6.5<p_{T}(J/#psi)<30.0 GeV/c");
	//pText1->SetTextSize(0.05);
	//pText1->SetFillColor(0);
	//pText1->Draw("");

	PrintOutTable(gAS_Result3, gAS_Result_Syst3);
	PrintOutTable(gAS_Result, gAS_Result_Syst);
	PrintOutTable(gAS_Result2, gAS_Result_Syst2);


	cankres1->SaveAs("NominalResultDCB_" + sEffName + sNameTag + ".root");
	cankres1->SaveAs("NominalResultDCB_" + sEffName + sNameTag + ".pdf");
	cankres1->SaveAs("NominalResultDCB_" + sEffName + sNameTag + ".png");

	return 0;
}

