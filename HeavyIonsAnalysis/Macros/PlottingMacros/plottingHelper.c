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

#include "../ChiFitterInit.h"

void RemoveXError(TGraphAsymmErrors* gAS);
void PrepareSystPlotting(TGraphAsymmErrors* gAS_Result_Syst, TGraphAsymmErrors* gAS_Result, TGraph* g_Systematics, double errorWidth); // adds the systematic uncertainty stored in percents as TGraph to our results (gAS_Result), and stores it in the separate gAS for plotting
void ShiftXPosition(TGraphAsymmErrors* gAS, double shiftSize); // move the points to avoid overlap
void ApplyPolarization(TGraphAsymmErrors* gASPolarized, TGraphAsymmErrors* gASNominal, const char* corrName, bool isInverseEfficiency = false, const char* fileCorrection = "Chi_c_Polarization_vTest.root"); // if the efficiency is flipped (e.g. chic1/chic2), use true
void PrintOutTable(TGraphAsymmErrors* gAS, TGraphAsymmErrors* gAS_Syst);

void RemoveXError(TGraphAsymmErrors* gAS)
{
	for (int i = 0; i < gAS->GetN(); i++)
	{
		gAS->SetPointEXlow(i, 0);
		gAS->SetPointEXhigh(i, 0);
	}
}

void PrepareSystPlotting(TGraphAsymmErrors* gAS_Result_Syst, TGraphAsymmErrors* gAS_Result, TGraph* g_Systematics, double errorWidth) {
	if (gAS_Result_Syst->GetN() <= g_Systematics->GetN()) { //smaller, because of nTrk dependence which has an extra bin
		for (int i = 0; i < gAS_Result_Syst->GetN(); i++) {
			if (gAS_Result_Syst->GetPointX(i) != g_Systematics->GetPointX(i)) {
				cout << "SYSTEMATICS AND NOMINAL HAVE DIFFERENT BINNING: " << gAS_Result_Syst->GetPointX(i) << " " << g_Systematics->GetPointX(i) << endl;
			}

			gAS_Result_Syst->SetPointEYhigh(i, g_Systematics->GetPointY(i)*0.01*gAS_Result->GetPointY(i)); //change from percent, and multiply by value
			gAS_Result_Syst->SetPointEYlow(i, g_Systematics->GetPointY(i)*0.01*gAS_Result->GetPointY(i)); //change from percent, and multiply by value
			gAS_Result_Syst->SetPointEXhigh(i, errorWidth);
			gAS_Result_Syst->SetPointEXlow(i, errorWidth);
		}

	}
	else { cout << "DIFFERENT NUMBER OF BINS BETWEEN SYST AND NOMINAL " << gAS_Result_Syst->GetN() << " " << g_Systematics->GetN() << endl; }
}

void ShiftXPosition(TGraphAsymmErrors* gAS, double shiftSize)
{
	for (int i = 0; i < gAS->GetN(); i++) {
		gAS->SetPointX(i, gAS->GetPointX(i) + shiftSize);
	}

}

void ApplyPolarization(TGraphAsymmErrors* gASPolarized, TGraphAsymmErrors* gASNominal, const char* corrName, bool isInverseEfficiency, const char* fileCorrection)
{
	if (gASPolarized->GetN() == gASNominal->GetN()) { //polarized and nominal to be declared above in the main file, since we might want to use in a loop, etc.
		TFile* myFile1 = new TFile(fileCorrection, "READ");
		cout << "Using polarization: " << corrName << endl;
		TH1D* hCor = (TH1D*)myFile1->Get(corrName);
		if (hCor != NULL)
		{
			cout << "Opened" << endl;
		}
		else {
			cout << "Problem opening histogram, exiting" << endl; return;
		}
		for (int i = 0; i < gASNominal->GetN(); i++) {
			double x, y;
			gASNominal->GetPoint(i, x, y); // Doing this the old way, should work for root 5.
			double polarizationShift = 1.0 / hCor->GetBinContent(i + 1);

			if (isInverseEfficiency) {
				polarizationShift = 1.0 / polarizationShift;					
			}
			cout << "Tests: Polarization shift: " << polarizationShift << endl;

			gASPolarized->SetPoint(i, x, y*polarizationShift);
			//gASPolarized->SetPointEYhigh(i, gASNominal->GetErrorYhigh(i)*polarizationShift); //not needed
			//gASPolarized->SetPointEYlow(i, gASNominal->GetErrorYlow(i)*polarizationShift); //not needed

			gASPolarized->SetPointEYhigh(i,0); //set to 0
			gASPolarized->SetPointEYlow(i,0); //set to 0

		}

	}
	else { cout << "DIFFERENT NUMBER OF BINS BETWEEN POLARIZED AND NOMINAL " << gASPolarized->GetN() << " " << gASNominal->GetN() << endl; }

}

void PrintOutTable(TGraphAsymmErrors* gAS, TGraphAsymmErrors* gAS_Syst)
{
	cout << endl << "******************************" << endl;
	cout << " Printing out table: " << endl;
	cout << "Graph name: " << gAS->GetName() << " and title: " << gAS->GetTitle() << endl;
	cout << "******************************" << endl << endl;
	//cout.precision(precision);
	//cout<<std::fixed;

	cout << "Bin  " << " Values" << endl;
	for (int i = 0; i < gAS->GetN(); i++)
	{
		cout << endl << gAS->GetPointX(i); printf(" & $ %.3f \\pm %.3f \\pm %.3f $ ", gAS->GetPointY(i), gAS->GetErrorY(i), gAS_Syst->GetErrorY(i));
	}


	cout << endl << endl;
}