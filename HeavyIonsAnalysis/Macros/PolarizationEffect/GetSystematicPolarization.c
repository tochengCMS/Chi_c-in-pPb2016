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

TGraph* gSystOutputArrayMC[nFittingSets]; // stores the polarization effect - chic over J/Psi
TGraph* gSystOutputArrayMC_chic1[nFittingSets]; // stores the systematic uncertainty - chic1 over J/Psi
TGraph* gSystOutputArrayMC_chic2[nFittingSets]; // stores the systematic uncertainty - chic2 over J/Psi
TGraph* gSystOutputArrayMC_chic1_over_chic2[nFittingSets]; // stores the systematic uncertainty - chic1 over chic2

const string fileOutCon = "Chi_c_Polarization_vTest.root";


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
	TFile* myFile1 = new TFile("Chi_c_WeightsMC_Official_vTest-bothDir_lambdaTheta1_0.00_lambdaTheta2_0.00.root", "READ");
	TFile* myFile2 = new TFile("Chi_c_WeightsMC_Official_vTest-bothDir_lambdaTheta1_0.50_lambdaTheta2_-0.39.root", "READ");

	// loop over all the results
	for (int iResult = 0; iResult < nFittingSets; iResult++)
	{
		//std::string hName = "h_photAcceptanceCor1D_Q_"+ setsToPlot.at(iResult)+"_rat";
		std::string hName = "h_chiTotalCorrection1D_" + setsToPlot.at(iResult) + "_rat";

		TH1D* gAS_Result = (TH1D*)myFile1->Get(hName.c_str());
		TH1D* gAS_Result2 = (TH1D*)myFile2->Get(hName.c_str());

		if (gAS_Result == nullptr || gAS_Result2 == nullptr) { cout << "FAILED TO OPEN, WILL SKIP THE BIN: " << hName<< endl; continue; }


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

		cout << "Here" << endl;

		gAS_Result->Draw("same");
		gAS_Result2->Draw("same");

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

		TLegend*leg = new TLegend(0.1,0.7,0.48,0.9);
		leg->SetFillColor(kWhite);
	    leg->SetFillStyle(0);	
		leg->SetBorderSize(0);
		leg->SetTextFont(42);
		leg->SetTextSize(0.03);
		leg->AddEntry(gAS_Result, "unpolarized", "p");
		leg->AddEntry(gAS_Result2, "#lambda_{c1} = 0.5, #lambda_{c2} = -0.39", "p");

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


		cankres1->SaveAs("PolarizationEffect/"+TString(type)+ "_" + sEffName + "_" + sNameTag + ".root");
		cankres1->SaveAs("PolarizationEffect/"+TString(type)+ "_" + sEffName + "_" + sNameTag + ".pdf");
		cankres1->SaveAs("PolarizationEffect/"+TString(type)+ "_" + sEffName + "_" + sNameTag + ".png");


		// READOUT TO FILE

		gSystOutputArrayMC[iResult] = new TGraph(gAS_Result->GetNbinsX()); // create the graph with the same number of points

		cout << "gSystOutputArrayMC_" + fittingSets.at(iResult) << endl;
		gSystOutputArrayMC[iResult]->SetNameTitle(("gSystOutputArrayMC_" + fittingSets.at(iResult)).c_str(), ("MC syst uncertainty for " + fittingSets.at(iResult)).c_str());

		for (int iPoint = 0; iPoint < gAS_Result->GetNbinsX(); iPoint++)
		{
			double systValue = abs(h_diffToNom2->GetBinContent(iPoint+1)-1);
			systValue = systValue * 100; //make it percents
            cout<<"systValue "<<systValue<<endl;
			gSystOutputArrayMC[iResult]->SetPoint(iPoint, gAS_Result->GetBinCenter(iPoint+1), systValue);
		}


		delete hframe;
		delete hframe2;
	}


	// loop over all the results
	for (int iResult = 0; iResult < nFittingSets; iResult++)
	{
		std::string hName = "h_chiTotCorrChic1_1D_" + setsToPlot.at(iResult) + "_rat";

		TH1D* gAS_Result = (TH1D*)myFile1->Get(hName.c_str());
		TH1D* gAS_Result2 = (TH1D*)myFile2->Get(hName.c_str());

		if (gAS_Result == nullptr || gAS_Result2 == nullptr) { cout << "FAILED TO OPEN, WILL SKIP THE BIN: " << hName<< endl; continue; }

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
		if (type==)
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

		cout << "Here" << endl;

		gAS_Result->Draw("same");
		gAS_Result2->Draw("same");

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

		TLegend*leg = new TLegend(0.1,0.7,0.48,0.9);
		leg->SetFillColor(kWhite);
	    leg->SetFillStyle(0);	
		leg->SetBorderSize(0);
		leg->SetTextFont(42);
		leg->SetTextSize(0.03);
		leg->AddEntry(gAS_Result, "unpolarized", "p");
		leg->AddEntry(gAS_Result2, "#lambda_{c1} = 0.5, #lambda_{c2} = -0.39", "p");

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

		TH1F* hframe2 = new TH1F("hframe2", "", 1, cLowX, cHighX);
		//hframe2->Draw();
		hframe2->GetYaxis()->SetRangeUser(0.95, 1.05);
		hframe2->GetXaxis()->SetTitle(xAxisTitle.c_str());
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

		cankres1->SaveAs("PolarizationEffect/Systematics_" + sEffName + "_" + sNameTag + ".root");
		cankres1->SaveAs("PolarizationEffect/Systematics_" + sEffName + "_" + sNameTag + ".pdf");
		cankres1->SaveAs("PolarizationEffect/Systematics_" + sEffName + "_" + sNameTag + ".png");


		// READOUT TO FILE
		gSystOutputArrayMC[iResult] = new TGraph(gAS_Result->GetNbinsX()); // create the graph with the same number of points

		cout << "gSystOutputArrayMC_chic1" + fittingSets.at(iResult) << endl;
		gSystOutputArrayMC_chic1[iResult]->SetNameTitle(("gSystOutputArrayMC_chic1_" + fittingSets.at(iResult)).c_str(), ("MC syst uncertainty for " + fittingSets.at(iResult)).c_str());

		for (int iPoint = 0; iPoint < gAS_Result->GetNbinsX(); iPoint++)
		{
			double systValue = abs(h_diffToNom2->GetBinContent(iPoint+1)-1);
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

		if (gSystOutputArrayMC_chic1[i] != NULL)
		{
			gSystOutputArrayMC_chic1[i]->Write();
		}

		if (gSystOutputArrayMC_chic2[i] != NULL)
		{
			gSystOutputArrayMC_chic2[i]->Write();
		}

		if (gSystOutputArrayMC_chic1_over_chic2[i] != NULL)
		{
			gSystOutputArrayMC_chic1_over_chic2[i]->Write();
		}		

	}


	return 0;
}


void plotSystematics(TH1D* myFile1, TH1D* myFile2, TGraph* gSystOutputArrayMC, std:string type)
{

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

	for (int iResult = 0; iResult < nFittingSets; iResult++)
	{
		std::string hName = "h_"+type+"_1D_" + setsToPlot.at(iResult) + "_rat";

		TH1D* gAS_Result = (TH1D*)myFile1->Get(hName.c_str());
		TH1D* gAS_Result2 = (TH1D*)myFile2->Get(hName.c_str());

		if (gAS_Result == nullptr || gAS_Result2 == nullptr) { cout << "FAILED TO OPEN, WILL SKIP THE BIN: " << hName<< endl; continue; }

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

		cout << "Here" << endl;

		gAS_Result->Draw("same");
		gAS_Result2->Draw("same");

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

		TLegend*leg = new TLegend(0.1,0.7,0.48,0.9);
		leg->SetFillColor(kWhite);
	    leg->SetFillStyle(0);	
		leg->SetBorderSize(0);
		leg->SetTextFont(42);
		leg->SetTextSize(0.03);
		leg->AddEntry(gAS_Result, "unpolarized", "p");
		leg->AddEntry(gAS_Result2, "#lambda_{c1} = 0.5, #lambda_{c2} = -0.39", "p");

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

		TH1F* hframe2 = new TH1F("hframe2", "", 1, cLowX, cHighX);
		//hframe2->Draw();
		hframe2->GetYaxis()->SetRangeUser(0.95, 1.05);
		hframe2->GetXaxis()->SetTitle(xAxisTitle.c_str());
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

		cankres1->SaveAs("PolarizationEffect/Systematics_" + sEffName + "_" + sNameTag + ".root");
		cankres1->SaveAs("PolarizationEffect/Systematics_" + sEffName + "_" + sNameTag + ".pdf");
		cankres1->SaveAs("PolarizationEffect/Systematics_" + sEffName + "_" + sNameTag + ".png");


		// READOUT TO FILE
		gSystOutputArrayMC[iResult] = new TGraph(gAS_Result->GetNbinsX()); // create the graph with the same number of points

		cout << "gSystOutputArrayMC_"+ type +"_"+ fittingSets.at(iResult) << endl;
		gSystOutputArrayMC[iResult]->SetNameTitle(("gSystOutputArrayMC_"+ TString(type) + "_" + fittingSets.at(iResult)).c_str(), ("MC syst uncertainty for " + fittingSets.at(iResult)).c_str());

		for (int iPoint = 0; iPoint < gAS_Result->GetNbinsX(); iPoint++)
		{
			double systValue = abs(h_diffToNom2->GetBinContent(iPoint+1)-1);
			systValue = systValue * 100; //make it percents
            cout<<"systValue "<<systValue<<endl;
			gSystOutputArrayMC[iResult]->SetPoint(iPoint, gAS_Result->GetBinCenter(iPoint+1), systValue);
		}


		delete hframe;
		delete hframe2;
	}

}

void RemoveXError(TGraphAsymmErrors* gAS)
{
	for (int i = 0; i < gAS->GetN(); i++)
	{
		gAS->SetPointEXlow(i, 0);
		gAS->SetPointEXhigh(i, 0);
	}
}
