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

#include "../../tdrstyle.C"
#include "../../CMS_lumi.C"
#include "../../ChiFitterInit.h"


const TString sEffName= "Chic_polarization";
//const float cLowX = 5, cHighX = 30;
//const TString sNameTag = "_pT_all_Acceptance";
//const float cLowX = -2.4, cHighX = 2.4;
//const TString sNameTag = "_y";
//const float cLowX = 0.0, cHighX = 270;
//const TString sNameTag = "_ntrk";

//std::vector<std::string> fittingSets = { "pt_all", "y", "nTrack", "pt_mid", "pt_fwd", "pt_fwdOnly", "pt_bkwOnly", "pt_fwdOnlyWide", "pt_bkwOnlyWide", "pt_midCMS", "pt_fwdCMS", "pt_bkwCMS", "nTrack_all" };
std::vector<std::string> setsToPlot = { "pT_all", "y", "pT_midCMS", "pT_fwdCMS", "pT_bkwCMS", "nTrk_all" }; //the name in the files are different for historic reasons
int nSets = 7;

TGraph* gSystOutputArrayMC[nFittingSets]; // stores the systematic uncertainty - total
const string fileOutCon = "Chi_c_MCSyst_vDissertation.root";


const float _markerSize = 2.4;
void RemoveXError(TGraphAsymmErrors* gAS);

int GetSystematicPolarization()  // produces MC systematic from the different MC trees, makes plots, and stores the actual value in a root file
{

	setTDRStyle();
	gStyle->SetPadRightMargin(0.03);

	TCanvas* cankres1 = new TCanvas("cankres1","Canvas with results1",960,960);


	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.360, 0.98, 1.0);
	pad1->SetBottomMargin(0.02);
	pad1->SetTicks(1, 1);
	pad1->Draw();

	cankres1->cd();
	//pull pad
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.006, 0.98, 0.360);
	pad2->SetTopMargin(0.035); // Upper and lower plot are joined
	pad2->SetBottomMargin(0.28);
	pad2->Draw();
	//pad2->SetTicks(1, 1);




	// Load data
	TFile* myFile1 = new TFile("Chi_c_WeightsMC_Official_vDissertation-bothDir_lambdaTheta1_0.00_lambdaTheta2_0.00.root", "READ");

	TFile* myFile2 = new TFile("Chi_c_WeightsMC_Official_vDissertation-bothDir_lambdaTheta1_-0.33_lambdaTheta2_1.00.root", "READ");

	TFile* myFile3 = new TFile("Chi_c_WeightsMC_Official_vDissertation-bothDir_lambdaTheta1_-0.33_lambdaTheta2_-0.33.root", "READ");

        TFile* myFile4 = new TFile("Chi_c_WeightsMC_Official_vDissertation-bothDir_lambdaTheta1_-0.33_lambdaTheta2_-0.60.root", "READ");

        TFile* myFile5 = new TFile("Chi_c_WeightsMC_Official_vDissertation-bothDir_lambdaTheta1_1.00_lambdaTheta2_1.00.root", "READ");

        TFile* myFile6 = new TFile("Chi_c_WeightsMC_Official_vDissertation-bothDir_lambdaTheta1_1.00_lambdaTheta2_-0.33.root", "READ");

        TFile* myFile7 = new TFile("Chi_c_WeightsMC_Official_vDissertation-bothDir_lambdaTheta1_1.00_lambdaTheta2_-0.60.root", "READ");

	// loop over all the results
	for (int iResult = 0; iResult < nFittingSets; iResult++)
	{
		//std::string hName = "h_photAcceptanceCor1D_Q_"+ setsToPlot.at(iResult)+"_rat";
		std::string hName = "h_chiTotalCorrection1D_" + setsToPlot.at(iResult) + "_rat";

		TH1D* gAS_Result = (TH1D*)myFile1->Get(hName.c_str());
		TH1D* gAS_Result2 = (TH1D*)myFile2->Get(hName.c_str());
		TH1D* gAS_Result3 = (TH1D*)myFile3->Get(hName.c_str());
		TH1D* gAS_Result4 = (TH1D*)myFile4->Get(hName.c_str());
		TH1D* gAS_Result5 = (TH1D*)myFile5->Get(hName.c_str());
                TH1D* gAS_Result6 = (TH1D*)myFile6->Get(hName.c_str());
                TH1D* gAS_Result7 = (TH1D*)myFile7->Get(hName.c_str());

		if (gAS_Result == nullptr || gAS_Result2 == nullptr || gAS_Result3 == nullptr || gAS_Result4 == nullptr || gAS_Result5 == nullptr) { cout << "FAILED TO OPEN, WILL SKIP THE BIN: " << hName<< endl; continue; }



		cout << "Loading maybe done" << endl;


		//Plotting setup

		TString sNameTag = setsToPlot.at(iResult);
		float cLowX = 5, cHighX = 30;
		string xAxisTitle = "p_{T}(J/#psi) [GeV/c]";
		if (setsToPlot.at(iResult) == "y") {
			cLowX = -2.4, cHighX = 2.4;
			xAxisTitle = "y_{Lab, p}(J/#psi)";
		}
		if (setsToPlot.at(iResult) == "nTrk_all" || setsToPlot.at(iResult) == "nTrk") {
			cLowX = 0.0, cHighX = 250;
			xAxisTitle = "N_{tracks}";
		}
		

		pad1->cd();
		TH1F* hframe = new TH1F("hframe", "", 1, cLowX, cHighX);
		hframe->Draw();
		hframe->GetYaxis()->SetRangeUser(0.0, 0.08);
		hframe->GetXaxis()->SetTitle(xAxisTitle.c_str());
		//hframe->GetXaxis()->SetTitleSize(0.05);
		//hframe->GetXaxis()->SetTitleOffset(1.12);
		hframe->GetYaxis()->SetTitle("Efficiency ((#chi_{c1}+#chi_{c2}) / J/#psi)");
		hframe->GetYaxis()->SetTitleSize(0.06);
		hframe->GetXaxis()->SetLabelSize(0); //remove it from the top plot

		gAS_Result->SetMarkerSize(0.9*_markerSize);
		gAS_Result->SetMarkerColor(kBlack);
		gAS_Result->SetMarkerStyle(20);
		gAS_Result->SetLineColor(kBlack);
		gAS_Result->SetLineWidth(2);

		gAS_Result2->SetMarkerSize(0.9*_markerSize);
		gAS_Result2->SetMarkerColor(kGreen);
		gAS_Result2->SetMarkerStyle(72);
		gAS_Result2->SetLineColor(kGreen);
		gAS_Result2->SetLineWidth(2);

		gAS_Result3->SetMarkerSize(1.2*_markerSize);
		gAS_Result3->SetMarkerColor(kGreen + 3);
		gAS_Result3->SetMarkerStyle(74);
		gAS_Result3->SetLineColor(kGreen + 3);
		gAS_Result3->SetLineWidth(2);

		gAS_Result4->SetMarkerSize(1.2*_markerSize);
		gAS_Result4->SetMarkerColor(kOrange);
		gAS_Result4->SetMarkerStyle(76);
		gAS_Result4->SetLineColor(kOrange);
		gAS_Result4->SetLineWidth(2);

		gAS_Result5->SetMarkerSize(1.2*_markerSize);
		gAS_Result5->SetMarkerColor(kCyan);
		gAS_Result5->SetMarkerStyle(83);
		gAS_Result5->SetLineColor(kCyan);
		gAS_Result5->SetLineWidth(2);

                gAS_Result6->SetMarkerSize(1.2*_markerSize);
                gAS_Result6->SetMarkerColor(kBlue);
                gAS_Result6->SetMarkerStyle(83);
                gAS_Result6->SetLineColor(kBlue);
                gAS_Result6->SetLineWidth(2);

                gAS_Result7->SetMarkerSize(1.2*_markerSize);
                gAS_Result7->SetMarkerColor(kRed);
                gAS_Result7->SetMarkerStyle(83);
                gAS_Result7->SetLineColor(kRed);
                gAS_Result7->SetLineWidth(2);

		cout << "Here" << endl;

		gAS_Result->Draw("same");
		gAS_Result2->Draw("same");
		gAS_Result3->Draw("same");
		gAS_Result4->Draw("same");
		gAS_Result5->Draw("same");
                gAS_Result6->Draw("same");
                gAS_Result7->Draw("same");

		TLatex latex;
		latex.SetNDC();
		latex.SetTextAngle(0);
		latex.SetTextColor(kBlack);
		latex.SetTextFont(42);
		latex.SetTextAlign(11);
		latex.SetTextSize(0.095);
		latex.DrawLatex(0.2, 0.83, "CMS");

		latex.SetTextFont(52);
		latex.SetTextSize(0.06);
		latex.DrawLatex(0.2, 0.76, "Preliminary");

		//TLegend*leg = new TLegend(0.58, 0.18, 0.88, 0.4, "");
		//leg->SetFillColor(kWhite);
		//leg->SetBorderSize(0);
		//leg->SetTextFont(42);
		//leg->SetTextSize(0.05);

		//TLegend*leg = new TLegend(0.58, 0.68, 0.88, 0.9, "");
		TLegend*leg = new TLegend(0.6, 0.3, 0.88, 0.75, "");
		leg->SetFillColor(kWhite);
		leg->SetBorderSize(0);
		leg->SetTextFont(42);
		leg->SetTextSize(0.03);

                /*
		if (setsToPlot.at(iResult) == "nTrk_all" || setsToPlot.at(iResult) == "nTrk") { // plotting details
			leg->SetY1(0.58);
			leg->SetY2(0.84);
		}
		if (setsToPlot.at(iResult) == "y") { // plotting details
			leg->SetY1(0.58);
			leg->SetX1(0.5);
			leg->SetY2(0.84);
			leg->SetX2(0.78);
		}
                */

		//leg->AddEntry(gAS_Result, "|y|<1.0, 6.5<p_{T}<30 GeV/c", "p");
		//leg->AddEntry(gAS_Result, "6.5<p_{T}<30 GeV/c", "p");

		//leg->AddEntry(gAS_Result, "Midrapidity: |y|<1.0", "p");
		//leg->AddEntry(gAS_Result2, "Forward: 1.6<y<2.4", "p");
		//leg->AddEntry(gAS_Result3, "Backward: -2.4<y<-1.6", "p");



		leg->AddEntry(gAS_Result, "unpolarized", "p");
		leg->AddEntry(gAS_Result2, "#chi_{c1} #pm 1, #chi_{c2} #pm 2", "l");
                leg->AddEntry(gAS_Result3, "#chi_{c1} #pm 1, #chi_{c2} #pm 1", "l");
                leg->AddEntry(gAS_Result4, "#chi_{c1} #pm 1, #chi_{c2} 0", "l");
                leg->AddEntry(gAS_Result5, "#chi_{c1} 0, #chi_{c2} #pm 2", "l");
                leg->AddEntry(gAS_Result6, "#chi_{c1} 0, #chi_{c2} #pm 1", "l");
                leg->AddEntry(gAS_Result7, "#chi_{c1} 0, #chi_{c2} 0", "l");

		leg->Draw("same");

		TLatex latex_text;
		latex_text.SetNDC();
		latex_text.SetTextColor(kBlack);
		latex_text.SetTextFont(42);
		latex_text.SetTextAlign(11);
		latex_text.SetTextSize(0.07);
		if (setsToPlot.at(iResult) == "pT_midCMS") { // plotting details
			latex_text.DrawLatex(0.6,0.8, "-1<y_{CM}(J/#psi)<1");
		}
		if (setsToPlot.at(iResult) == "pT_bkwCMS") { // plotting details
			latex_text.DrawLatex(0.6, 0.8, "-2<y_{CM}(J/#psi)<-1");
		}
		if (setsToPlot.at(iResult) == "pT_fwdCMS") { // plotting details
			latex_text.DrawLatex(0.6, 0.8, "1<y_{CM}(J/#psi)<1.9");
		}
		pad1->Update();


		// ADD DIFFERENCE PLOT
		pad2->cd();
		TH1D* h_diffToNom2 = (TH1D*)gAS_Result2->Clone(); //needed to initialize hist
		h_diffToNom2->Divide(gAS_Result2, gAS_Result);
		TH1D* h_diffToNom3 = (TH1D*)gAS_Result3->Clone(); //needed to initialize hist
		h_diffToNom3->Divide(gAS_Result3, gAS_Result);
		TH1D* h_diffToNom4 = (TH1D*)gAS_Result4->Clone(); //needed to initialize hist
		h_diffToNom4->Divide(gAS_Result4, gAS_Result);
		TH1D* h_diffToNom5 = (TH1D*)gAS_Result5->Clone(); //needed to initialize hist
		h_diffToNom5->Divide(gAS_Result5, gAS_Result);
                TH1D* h_diffToNom6 = (TH1D*)gAS_Result5->Clone(); //needed to initialize hist
                h_diffToNom6->Divide(gAS_Result5, gAS_Result);
                TH1D* h_diffToNom7 = (TH1D*)gAS_Result5->Clone(); //needed to initialize hist
                h_diffToNom7->Divide(gAS_Result5, gAS_Result);

		TH1F* hframe2 = new TH1F("hframe2", "", 1, cLowX, cHighX);
		//hframe2->Draw();
		hframe2->GetYaxis()->SetRangeUser(0.95, 1.05);
		hframe2->GetXaxis()->SetTitle(xAxisTitle.c_str());
		//hframe2->GetXaxis()->SetTitle("rapidity (p-going)");
		//hframe2->GetXaxis()->SetTitle("N_{tracks}");
		hframe2->GetXaxis()->SetTitleSize(0.11);
		hframe2->GetXaxis()->SetTitleOffset(1.0);
		hframe2->GetXaxis()->SetLabelSize(0.11);

		hframe2->GetYaxis()->SetLabelSize(0.09);
		hframe2->GetYaxis()->SetTitle("Ratio to unpolarized");
		hframe2->GetYaxis()->SetTitleOffset(0.75);
		hframe2->GetYaxis()->SetTitleSize(0.10);

		hframe2->Draw();

		//one
		TLine *l1 = new TLine(cLowX, 1, cHighX, 1);
		l1->SetLineStyle(9);
		l1->Draw("same");

		h_diffToNom2->Draw("same");
		h_diffToNom3->Draw("same");
		h_diffToNom4->Draw("same");
		h_diffToNom5->Draw("same");
                h_diffToNom6->Draw("same");
                h_diffToNom7->Draw("same");

		//// Add text to frame
		//TText* txt = new TText(2, 100, "Signal");
		//txt->SetTextSize(0.04);
		//txt->SetTextColor(kRed);
		//txt->Draw("same");
		////massframeBin->addObject(txt);




		//TPaveText* pText1 = new TPaveText(0.62, 0.4, 0.88, 0.48, "NDC NB");
		//pText1->AddText("6.5<p_{T}(J/#psi)<30.0 GeV/c");
		//pText1->SetTextSize(0.05);
		//pText1->SetFillColor(0);
		//pText1->Draw("");

		//*/
		cankres1->SaveAs("PlotsActualResults/Systematics_" + sEffName + "_" + sNameTag + ".root");
		cankres1->SaveAs("PlotsActualResults/Systematics_" + sEffName + "_" + sNameTag + ".pdf");
		cankres1->SaveAs("PlotsActualResults/Systematics_" + sEffName + "_" + sNameTag + ".png");


		// READOUT TO FILE

		gSystOutputArrayMC[iResult] = new TGraph(gAS_Result->GetNbinsX()); // create the graph with the same number of points

		cout << "gSystOutputArrayMC_" + fittingSets.at(iResult) << endl;
		gSystOutputArrayMC[iResult]->SetNameTitle(("gSystOutputArrayMC_" + fittingSets.at(iResult)).c_str(), ("MC syst uncertainty for " + fittingSets.at(iResult)).c_str());

		for (int iPoint = 0; iPoint < gAS_Result->GetNbinsX(); iPoint++)
		{
			double systValue = std::max({ abs(h_diffToNom2->GetBinContent(iPoint+1)-1), abs(h_diffToNom3->GetBinContent(iPoint+1) - 1), abs(h_diffToNom4->GetBinContent(iPoint+1) - 1), abs(h_diffToNom5->GetBinContent(iPoint+1) - 1), abs(h_diffToNom6->GetBinContent(iPoint+1) - 1), abs(h_diffToNom6->GetBinContent(iPoint+1) - 1), abs(h_diffToNom7->GetBinContent(iPoint+1) - 1) });
			systValue = systValue * 100; //make it percents
                        cout<<"systValue "<<systValue<<endl;
			gSystOutputArrayMC[iResult]->SetPoint(iPoint, gAS_Result->GetBinCenter(iPoint+1), systValue);
		}








		delete hframe;
		delete hframe2;
	}

	TFile* fout = new TFile(fileOutCon.c_str(), "RECREATE");


	for (int i = 0; i < nFittingSets; i++)
	{

		if (gSystOutputArrayMC[i] != NULL)
		{

			gSystOutputArrayMC[i]->Write();
		}

	}


	return 0;
}


void RemoveXError(TGraphAsymmErrors* gAS)
{
	for (int i = 0; i < gAS->GetN(); i++)
	{
		gAS->SetPointEXlow(i, 0);
		gAS->SetPointEXhigh(i, 0);
	}
}
