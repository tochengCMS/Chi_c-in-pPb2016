#include "ChiTreeInit.h"
//#include "ChiFitterInit.h"
#include <iostream>
using namespace std;

/////////////////////////
/////// functions ///////
//////////////////////////

int LoadChiBranches(TTree* tree, bool isMC, bool minimalOnly) {
	if (!tree) {
		cout << "Problem with event tree";
		return 1;
	}

	if (minimalOnly == true) //changed as needed, effort to speed up simple loops. Needs to be checked, as of now content not guaranteed
	{
		tree->SetBranchAddress("ntracks_inEvent", &ntracks_inEvent);
		tree->SetBranchAddress("pvtx_nTracks", &pvtx_nTracks);
		tree->SetBranchAddress("pvtx_nTracksUncut", &pvtx_nTracksUncut);
		tree->SetBranchAddress("dimuon_p4", &dimuon_p4); //TLorentzVector
		tree->SetBranchAddress("dimuon_eta", &dimuon_eta);
		tree->SetBranchAddress("dimuon_pt", &dimuon_pt);
		tree->SetBranchAddress("dimuon_charge", &dimuon_charge);
		tree->SetBranchAddress("dimuon_vtxProb", &dimuon_vtxProb);
		tree->SetBranchAddress("dimuon_pvtx_indexFromOniaMuMu", &dimuon_pvtx_indexFromOniaMuMu);
		tree->SetBranchAddress("dimuon_muon1_position", &dimuon_muon1_position); //stores position of first muon in muon collection (probably one with higher pT due to ordering in the collection)
		tree->SetBranchAddress("dimuon_muon2_position", &dimuon_muon2_position); //stores position of second muon in muon collection 
		tree->SetBranchAddress("dimuon_ctpv", &dimuon_ctpv);
		tree->SetBranchAddress("dimuon_ctpvError", &dimuon_ctpvError);

		tree->SetBranchAddress("muonIsHLTDoubleMuOpen", &muonIsHLTDoubleMuOpen);
		tree->SetBranchAddress("muonIsHLTDoubleMuOpenFilter", &muonIsHLTDoubleMuOpenFilter);
		tree->SetBranchAddress("muonIsSoft", &muonIsSoft);
		tree->SetBranchAddress("muon_eta", &muon_eta);
		tree->SetBranchAddress("muon_pt", &muon_pt);

		if (isMC == true)
		{
			tree->SetBranchAddress("gen_Jpsi_pt", &gen_Jpsi_pt);
			tree->SetBranchAddress("gen_muon_pt", &gen_muon_pt);
			tree->SetBranchAddress("gen_phot_pt", &gen_phot_pt);
			tree->SetBranchAddress("gen_isGoodChicDecay", &gen_isGoodChicDecay);
		}
	}
	else
	{
		//general
		tree->SetBranchAddress("ispPb", &ispPb);
		tree->SetBranchAddress("runNumber", &runNumber);
		tree->SetBranchAddress("eventNumber", &eventNumber);
		tree->SetBranchAddress("nPrimVertices", &nPrimVertices);
		tree->SetBranchAddress("muonPerEvent_noCuts", &muonPerEvent);
		tree->SetBranchAddress("convPerEvent_noCuts", &convPerTriggeredEvent);
		tree->SetBranchAddress("dimuonPerEvent", &dimuonPerEvent);
		tree->SetBranchAddress("chiCandPerEvent", &chiCandPerEvent);

		tree->SetBranchAddress("ntracks_inEvent", &ntracks_inEvent);
		tree->SetBranchAddress("hfTowerSum_inEvent", &hfTowerSum_inEvent);
		tree->SetBranchAddress("Trig_Event_HLTDoubleMuOpen", &Trig_Event_HLTDoubleMuOpen);


		// vertex
		tree->SetBranchAddress("pvtx_z", &pvtx_z);
		tree->SetBranchAddress("pvtx_zError", &pvtx_zError);
		tree->SetBranchAddress("pvtx_x", &pvtx_x);
		tree->SetBranchAddress("pvtx_y", &pvtx_y);
		tree->SetBranchAddress("pvtx_nTracksUncut", &pvtx_nTracksUncut);
		tree->SetBranchAddress("pvtx_nTracks", &pvtx_nTracks);
		tree->SetBranchAddress("pvtx_isFake", &pvtx_isFake);

		//muon info
		tree->SetBranchAddress("muonIsHLTDoubleMuOpen", &muonIsHLTDoubleMuOpen);
		tree->SetBranchAddress("muonIsHLTDoubleMuOpenFilter", &muonIsHLTDoubleMuOpenFilter);
		tree->SetBranchAddress("muonIsGlobal", &muonIsGlobal);
		tree->SetBranchAddress("muonIsTracker", &muonIsTracker);
		tree->SetBranchAddress("muonIsPF", &muonIsPF);
		tree->SetBranchAddress("muonIsSoft", &muonIsSoft);
		tree->SetBranchAddress("muonIsTight", &muonIsTight);
		tree->SetBranchAddress("muonIsNotGlobalNorTracker", &muonIsNotGlobalNorTracker);
		tree->SetBranchAddress("muonIDHas_TMOneStationTight", &muonIDHas_TMOneStationTight);
		tree->SetBranchAddress("muon_pvtx_index", &muon_pvtx_index);
		tree->SetBranchAddress("muonInnerTrack_dxy", &muonInnerTrack_dxy);
		tree->SetBranchAddress("muonInnerTrack_dz", &muonInnerTrack_dz);
		tree->SetBranchAddress("muonTrackerLayersWithMeasurement", &muonTrackerLayersWithMeasurement);
		tree->SetBranchAddress("muonPixelLayersWithMeasurement", &muonPixelLayersWithMeasurement);
		tree->SetBranchAddress("muonQuality_isHighPurity", &muonQuality_isHighPurity);
		tree->SetBranchAddress("muon_charge", &muon_charge);
		tree->SetBranchAddress("muon_eta", &muon_eta);
		tree->SetBranchAddress("muon_pt", &muon_pt);
		tree->SetBranchAddress("muon_p4", &muon_p4); //TLorentzVector
		//tree->SetBranchAddress("patMuonStored", &patMuonStored);

		//muon MC
		if (isMC == true && flag_saveExtraThings == true) {
			tree->SetBranchAddress("muon_isMatchedMC", &muon_isMatchedMC);
			tree->SetBranchAddress("muonGen_eta", &muonGen_eta);
			tree->SetBranchAddress("muonGen_pt", &muonGen_pt);
			tree->SetBranchAddress("muonGen_p4", &muonGen_p4); //TLorentzVector
			tree->SetBranchAddress("muonGen_rDelta", &muonGen_rDelta);
			tree->SetBranchAddress("muonGen_ptDelta", &muonGen_ptDelta);
			tree->SetBranchAddress("muonGen_ptDeltaRel", &muonGen_ptDeltaRel);
		}

		//dimuon info

		tree->SetBranchAddress("dimuon_p4", &dimuon_p4); //TLorentzVector
		tree->SetBranchAddress("dimuon_eta", &dimuon_eta);
		tree->SetBranchAddress("dimuon_pt", &dimuon_pt);
		tree->SetBranchAddress("dimuon_charge", &dimuon_charge);
		tree->SetBranchAddress("dimuon_vtx", &dimuon_vtx); //TVector3
		tree->SetBranchAddress("dimuon_pvtx_indexFromOniaMuMu", &dimuon_pvtx_indexFromOniaMuMu);
		tree->SetBranchAddress("dimuon_pvtx_index", &dimuon_pvtx_index);

		tree->SetBranchAddress("dimuon_dz_dimuonvtx_pvtx", &dimuon_dz_dimuonvtx_pvtx);
		tree->SetBranchAddress("dimuon_vtxProb", &dimuon_vtxProb);
		//tree->SetBranchAddress("dimuonStored", &dimuonStored);
		tree->SetBranchAddress("dimuon_muon1_position", &dimuon_muon1_position); //stores position of first muon in muon collection (probably one with higher pT due to ordering in the collection)
		tree->SetBranchAddress("dimuon_muon2_position", &dimuon_muon2_position); //stores position of second muon in muon collection 
		tree->SetBranchAddress("dimuon_ctpv", &dimuon_ctpv);
		tree->SetBranchAddress("dimuon_ctpvError", &dimuon_ctpvError);



		//conversion info

		tree->SetBranchAddress("convRaw_duplicityStatus", &convRaw_duplicityStatus);
		tree->SetBranchAddress("convRaw_duplicityStatus_AV", &convRaw_duplicityStatus_AV);
		tree->SetBranchAddress("convRaw_splitDR", &convRaw_splitDR);
		tree->SetBranchAddress("convRaw_splitDpT", &convRaw_splitDpT);
		tree->SetBranchAddress("conv_positionRaw", &conv_positionRaw);

		tree->SetBranchAddress("conv_tk1ValidHits", &conv_tk1ValidHits);
		tree->SetBranchAddress("conv_tk2ValidHits", &conv_tk2ValidHits);
		tree->SetBranchAddress("convQuality_isHighPurity", &convQuality_isHighPurity);
		tree->SetBranchAddress("convQuality_isGeneralTracksOnly", &convQuality_isGeneralTracksOnly);
		tree->SetBranchAddress("conv_vtx", &conv_vtx); //TVector3
		tree->SetBranchAddress("conv_vertexPositionRho", &conv_vertexPositionRho);
		tree->SetBranchAddress("conv_sigmaTkVtx1", &conv_sigmaTkVtx1);
		tree->SetBranchAddress("conv_sigmaTkVtx2", &conv_sigmaTkVtx2);
		tree->SetBranchAddress("conv_tkVtxCompatibilityOK", &conv_tkVtxCompatibilityOK);
		tree->SetBranchAddress("conv_tkVtxCompatible_bestVertex", &conv_tkVtxCompatible_bestVertex);
		tree->SetBranchAddress("conv_tkVtxCompatible_secondBestVertexA", &conv_tkVtxCompatible_secondBestVertexA);
		tree->SetBranchAddress("conv_tkVtxCompatible_secondBestVertexB", &conv_tkVtxCompatible_secondBestVertexB);
		//tree->SetBranchAddress("conv_tkVtxCompatibilityOK_test", &conv_tkVtxCompatibilityOK_test); //test ones
		//tree->SetBranchAddress("conv_tkVtxCompatible_bestVertex_test", &conv_tkVtxCompatible_bestVertex_test); //test ones
		//tree->SetBranchAddress("conv_tkVtxCompatible_secondBestVertexA_test", &conv_tkVtxCompatible_secondBestVertexA_test); //test
		//tree->SetBranchAddress("conv_tkVtxCompatible_secondBestVertexB_test", &conv_tkVtxCompatible_secondBestVertexB_test); //test

		tree->SetBranchAddress("conv_compatibleInnerHitsOK", &conv_compatibleInnerHitsOK); //-1: less than 2 tracks, 0: not compatible, 1: yes
		//tree->SetBranchAddress("conv_hitPat1", &conv_hitPat1);
		//tree->SetBranchAddress("conv_hitPat2", &conv_hitPat2);
		tree->SetBranchAddress("conv_vertexChi2Prob", &conv_vertexChi2Prob);
		tree->SetBranchAddress("conv_pvtx_index", &conv_pvtx_index);
		tree->SetBranchAddress("conv_zOfPriVtx", &conv_zOfPriVtx); // z of primary vertex that is used in the conversions (could be obtained also from pvtx_z)
		tree->SetBranchAddress("conv_zOfPriVtxFromTracks", &conv_zOfPriVtxFromTracks);
		tree->SetBranchAddress("conv_dzToClosestPriVtx", &conv_dzToClosestPriVtx);
		tree->SetBranchAddress("conv_dxyPriVtx_Tr1", &conv_dxyPriVtx_Tr1);
		tree->SetBranchAddress("conv_dxyPriVtx_Tr2", &conv_dxyPriVtx_Tr2);
		tree->SetBranchAddress("conv_dxyPriVtxTimesCharge_Tr1", &conv_dxyPriVtxTimesCharge_Tr1);
		tree->SetBranchAddress("conv_dxyPriVtxTimesCharge_Tr2", &conv_dxyPriVtxTimesCharge_Tr2);
		tree->SetBranchAddress("conv_dxyError_Tr1", &conv_dxyError_Tr1);
		tree->SetBranchAddress("conv_dxyError_Tr2", &conv_dxyError_Tr2);

		tree->SetBranchAddress("conv_tk1NumOfDOF", &conv_tk1NumOfDOF);
		tree->SetBranchAddress("conv_tk2NumOfDOF", &conv_tk2NumOfDOF);
		tree->SetBranchAddress("conv_track1Chi2", &conv_track1Chi2);
		tree->SetBranchAddress("conv_track2Chi2", &conv_track2Chi2);
		tree->SetBranchAddress("conv_Tr1_pt", &conv_Tr1_pt);
		tree->SetBranchAddress("conv_Tr2_pt", &conv_Tr2_pt);

		tree->SetBranchAddress("conv_minDistanceOfApproach", &conv_minDistanceOfApproach);
		tree->SetBranchAddress("conv_p4", &conv_p4); //TLorentzVector
		tree->SetBranchAddress("conv_eta", &conv_eta);
		tree->SetBranchAddress("conv_pt", &conv_pt);


		//conv MC
		if (isMC == true && flag_saveExtraThings == true) {
			tree->SetBranchAddress("conv_isMatchedMC", &conv_isMatchedMC);
			tree->SetBranchAddress("convGen_eta", &convGen_eta);
			tree->SetBranchAddress("convGen_pt", &convGen_pt);
			tree->SetBranchAddress("convGen_p4", &convGen_p4); //TLorentzVector
			tree->SetBranchAddress("convGen_rDelta", &convGen_rDelta);
			tree->SetBranchAddress("convGen_ptDelta", &convGen_ptDelta);
			tree->SetBranchAddress("convGen_ptDeltaRel", &convGen_ptDeltaRel);
			tree->SetBranchAddress("convGen_motherCode", &convGen_motherCode);
		}

		// MC general
		if (isMC == true) {
			tree->SetBranchAddress("gen_isGoodChicDecay", &gen_isGoodChicDecay);
			tree->SetBranchAddress("gen_pdgId", &gen_pdgId);
			tree->SetBranchAddress("gen_chic_pt", &gen_chic_pt);
			tree->SetBranchAddress("gen_chic_eta", &gen_chic_eta);
			tree->SetBranchAddress("gen_chic_p4", &gen_chic_p4); //TLorentzVector
			tree->SetBranchAddress("gen_chic_matchPosition", &gen_chic_matchPosition);
			tree->SetBranchAddress("gen_chic_nMatches", &gen_chic_nMatches);
			tree->SetBranchAddress("gen_chic_rDelta", &gen_chic_rDelta); //in principle duplicates information
			tree->SetBranchAddress("gen_chic_ptDeltaRel", &gen_chic_ptDeltaRel);//in principle duplicates information

			tree->SetBranchAddress("gen_Jpsi_pt", &gen_Jpsi_pt);
			tree->SetBranchAddress("gen_Jpsi_eta", &gen_Jpsi_eta);
			tree->SetBranchAddress("gen_Jpsi_matchPosition", &gen_Jpsi_matchPosition);
			tree->SetBranchAddress("gen_Jpsi_nMatches", &gen_Jpsi_nMatches);
			tree->SetBranchAddress("gen_Jpsi_rDelta", &gen_Jpsi_rDelta); //in principle duplicates information
			tree->SetBranchAddress("gen_Jpsi_ptDeltaRel", &gen_Jpsi_ptDeltaRel);//in principle duplicates information
			tree->SetBranchAddress("gen_Jpsi_p4", &gen_Jpsi_p4); //TLorentzVector

			tree->SetBranchAddress("gen_Jpsi_photon_n", &gen_Jpsi_photon_n);
			tree->SetBranchAddress("gen_Jpsi_photon_pt", &gen_Jpsi_photon_pt);
			tree->SetBranchAddress("gen_Jpsi_photon_p4", &gen_Jpsi_photon_p4); //TLorentzVector

			tree->SetBranchAddress("gen_muon_charge", &gen_muon_charge);
			tree->SetBranchAddress("gen_muon_pt", &gen_muon_pt);
			tree->SetBranchAddress("gen_muon_eta", &gen_muon_eta);
			tree->SetBranchAddress("gen_muon_matchPosition", &gen_muon_matchPosition);
			tree->SetBranchAddress("gen_muon_nMatches", &gen_muon_nMatches);
			tree->SetBranchAddress("gen_muon_rDelta", &gen_muon_rDelta); //in principle duplicates information
			tree->SetBranchAddress("gen_muon_ptDeltaRel", &gen_muon_ptDeltaRel);//in principle duplicates information
			tree->SetBranchAddress("gen_muon_p4", &gen_muon_p4); //TLorentzVector

			tree->SetBranchAddress("gen_phot_pt", &gen_phot_pt);
			tree->SetBranchAddress("gen_phot_eta", &gen_phot_eta);
			tree->SetBranchAddress("gen_phot_p4", &gen_phot_p4); //TLorentzVector
			tree->SetBranchAddress("gen_conv_matchPosition", &gen_conv_matchPosition);
			tree->SetBranchAddress("gen_conv_nMatches", &gen_conv_nMatches);
			tree->SetBranchAddress("gen_conv_rDelta", &gen_conv_rDelta); //in principle duplicates information
			tree->SetBranchAddress("gen_conv_ptDeltaRel", &gen_conv_ptDeltaRel);//in principle duplicates information
		}

		// chi

		tree->SetBranchAddress("chi_p4", &chi_p4); //TLorentzVector
		tree->SetBranchAddress("chi_eta", &chi_eta);
		tree->SetBranchAddress("chi_pt", &chi_pt);
		tree->SetBranchAddress("chi_daughterJpsi_position", &chi_daughterJpsi_position); //stores position of daughter Jpsi in dimuon collection
		tree->SetBranchAddress("chi_daughterConv_position", &chi_daughterConv_position); //stores position of daughter photon (conversion)
		tree->SetBranchAddress("chi_dzPhotToDimuonVtx", &chi_dzPhotToDimuonVtx); //z distance of photon to dimuon vertex when dxy is minimal
		tree->SetBranchAddress("chi_dxyPhotToDimuonVtx", &chi_dxyPhotToDimuonVtx); //dxy distance of photon to dimuon vertex when dz is 0 - probably not too good for very midrapidity conversions
		//tree->SetBranchAddress("chiStored", &chiStored);
		tree->SetBranchAddress("chi_kinematicRefitFlag", &chi_kinematicRefitFlag); // -1 kinematic refit not done, 1 done: +2 needed extra refit for photon +4 something wrong with photon at the end +8 something wrong with the final fit 
		tree->SetBranchAddress("chi_refit_origChicPosition", &chi_refit_origChicPosition); //stores position of the chic candidate for the refit (there will be gaps if refit fails) 
		//tree->SetBranchAddress("chi_refitStored", &chi_refitStored);
		tree->SetBranchAddress("chi_refit_vprob", &chi_refit_vprob);
		tree->SetBranchAddress("chi_refit_ctauPV", &chi_refit_ctauPV);
		tree->SetBranchAddress("chi_refit_ctauErrPV", &chi_refit_ctauErrPV);
		tree->SetBranchAddress("chi_refit_ctauPV3D", &chi_refit_ctauPV3D);
		tree->SetBranchAddress("chi_refit_pvtxFromPVwithMuons_x", &chi_refit_pvtxFromPVwithMuons_x);
		tree->SetBranchAddress("chi_refit_pvtxFromPVwithMuons_y", &chi_refit_pvtxFromPVwithMuons_y);
		tree->SetBranchAddress("chi_refit_pvtxFromPVwithMuons_z", &chi_refit_pvtxFromPVwithMuons_z);

	}
	return 0;
}

//muon  - acceptance

bool MuonAcceptanceLoose(double eta, double pt)
{
	if (fabs(eta) > 2.4) return false;  //2.4
	if (fabs(eta) < 0.3 && pt < 3.4) return false;
	if (fabs(eta) < 1.1 && pt < 3.3) return false;
	if (fabs(eta) >= 1.1 && fabs(eta) < 2.1 && pt < 5.5 - 2 * fabs(eta)) return false;
	if (fabs(eta) >= 2.1 && pt < 1.3) return false;
	return true;
}
bool MuonAcceptanceTight(double eta, double pt)
{
	if (fabs(eta) > 2.4) return false;  //2.4
	if (fabs(eta) < 0.3 && pt < 3.4) return false;
	if (fabs(eta) < 2.4 && pt < 3.3) return false;
	return true;
}

bool MuonAcceptance(double eta, double pt) {
	return MuonAcceptanceLoose(eta, pt);
}
// muon selection

bool MuonSelectionPass(int muonPos)  //uses variables loaded in main function
{
	if (muonIsSoft->at(muonPos) != 1) return false;
	if (muonIsHLTDoubleMuOpen->at(muonPos) != 1) return false;
	return true;
}

int MuonMCMatched(int muonMCPos)
{
	int matchPosition = gen_muon_matchPosition->at(muonMCPos);
	if (matchPosition < -0.5) return -1; //not matched
	if (gen_muon_rDelta->at(muonMCPos) > muon_maxDeltaR_analysis || gen_muon_ptDeltaRel->at(muonMCPos) > muon_maxDPtRel_analysis) return -1;//not matched
	return matchPosition;
}

bool MuonSelectionPassMC(int muonMCPos)  //uses variables loaded in main function
{
	int matchPosition = gen_muon_matchPosition->at(muonMCPos);
	if (matchPosition < -0.5) return false; //not matched
	if (gen_muon_rDelta->at(muonMCPos)> muon_maxDeltaR_analysis || gen_muon_ptDeltaRel->at(muonMCPos) > muon_maxDPtRel_analysis ) return false;//not matched
	if (muonIsSoft->at(matchPosition) != 1) return false;
	if (muonIsHLTDoubleMuOpen->at(matchPosition) != 1) return false;
	return true;
}


// conversion

bool PhotAcceptance(double eta, double pt)
{
	if (fabs(eta) > 2.4) return false; //2.5
	if (pt < 0.5) return false; // 0.2
	return true;
}



bool PhotSelectionPassTight(int photPos)  //uses variables loaded in main function //nominal by Alberto and Jhovanny (8.2021) (with updated conversion removal)
{
	//if (conv_duplicityStatus->at(photPos) != 0 && conv_duplicityStatus->at(photPos) != 1 && conv_duplicityStatus->at(photPos) != 2) return false;
	if (convQuality_isHighPurity->at(photPos) != 1) return false;
	if (convQuality_isGeneralTracksOnly->at(photPos) != 1) return false;
	if (conv_vertexPositionRho->at(photPos) <= 1.5) return false;
	if (conv_sigmaTkVtx1->at(photPos) > 5) return false;
	if (conv_sigmaTkVtx2->at(photPos) > 5) return false;
	//if (conv_tkVtxCompatibilityOK->at(photPos) != 1) return false;
	if (conv_compatibleInnerHitsOK->at(photPos) != 1) return false;
	if (conv_vertexChi2Prob->at(photPos) <= 0.0005) return false;
	//if (fabs(conv_zOfPriVtx->at(photPos)) >= 20) return false;
	//if (fabs(conv_dzToClosestPriVtx->at(photPos)) >= 10) return false;
	if (conv_tk1NumOfDOF->at(photPos) < 2.5) return false;
	if (conv_tk2NumOfDOF->at(photPos) < 2.5) return false;
	if (conv_track1Chi2->at(photPos) >= 10) return false;
	if (conv_track2Chi2->at(photPos) >= 10) return false;
	if (conv_minDistanceOfApproach->at(photPos) <= -0.25) return false; //removes -1000 values
	if (conv_minDistanceOfApproach->at(photPos) >= 1.00) return false;

	return true;
}

bool PhotSelectionPassMedium(int photPos)  //uses variables loaded in main function
{
	//if (conv_duplicityStatus->at(photPos) != 0 && conv_duplicityStatus->at(photPos) != 1 && conv_duplicityStatus->at(photPos) != 2) return false;
	//if (convQuality_isHighPurity->at(photPos) != 1) return false;
	if (convQuality_isGeneralTracksOnly->at(photPos) != 1) return false;
	//if (conv_vertexPositionRho->at(photPos) <= 1.5) return false;
	//if (conv_sigmaTkVtx1->at(photPos) > 10) return false;
	//if (conv_sigmaTkVtx2->at(photPos) > 10) return false;
	//if (conv_tkVtxCompatibilityOK->at(photPos) != 1) return false;
	if (conv_compatibleInnerHitsOK->at(photPos) != 1) return false;
	if (conv_vertexChi2Prob->at(photPos) <= 0.0005) return false;
	//if (fabs(conv_zOfPriVtx->at(photPos)) >= 20) return false;
	if (fabs(conv_dzToClosestPriVtx->at(photPos)) >= 10) return false;
	if (conv_tk1NumOfDOF->at(photPos) < 2.5) return false;
	if (conv_tk2NumOfDOF->at(photPos) < 2.5) return false;
	//if (conv_track1Chi2->at(photPos) >= 10) return false;
	//if (conv_track2Chi2->at(photPos) >= 10) return false;
	if (conv_minDistanceOfApproach->at(photPos) <= -10) return false; //removes -1000 values
	//if (conv_minDistanceOfApproach->at(photPos) >= 1.00) return false;

	return true;
}

bool PhotSelectionPassLoose(int photPos)  //uses variables loaded in main function
{
	//if (conv_duplicityStatus->at(photPos) != 0 && conv_duplicityStatus->at(photPos) != 1 && conv_duplicityStatus->at(photPos) != 2) return false;
	//if (convQuality_isHighPurity->at(photPos) != 1) return false;
	if (convQuality_isGeneralTracksOnly->at(photPos) != 1) return false;
	//if (conv_vertexPositionRho->at(photPos) <= 1.5) return false;
	//if (conv_sigmaTkVtx1->at(photPos) > 10) return false;
	//if (conv_sigmaTkVtx2->at(photPos) > 10) return false;
	//if (conv_tkVtxCompatibilityOK->at(photPos) != 1) return false;
	if (conv_compatibleInnerHitsOK->at(photPos) != 1) return false;
	//if (conv_vertexChi2Prob->at(photPos) <= 0.001) return false;
	//if (fabs(conv_zOfPriVtx->at(photPos)) >= 20) return false;
	if (fabs(conv_dzToClosestPriVtx->at(photPos)) >= 10) return false;
	//if (conv_tk1NumOfDOF->at(photPos) < 2.5) return false;
	//if (conv_tk2NumOfDOF->at(photPos) < 2.5) return false;
	//if (conv_track1Chi2->at(photPos) >= 10) return false;
	//if (conv_track2Chi2->at(photPos) >= 10) return false;
	if (conv_minDistanceOfApproach->at(photPos) <= -10) return false; //removes -1000 values
	//if (conv_minDistanceOfApproach->at(photPos) >= 1.00) return false;

	return true;
}

bool PhotSelectionPassVeryLoose(int photPos)  //uses variables loaded in main function
{
	if (convQuality_isGeneralTracksOnly->at(photPos) != 1) return false;
	return true;
}

bool PhotSelectionPass(int photPos) // nominal version
{
	return PhotSelectionPassLoose(photPos);
}

int PhotMCMatched(int photMCPos)
{
	int matchPosition = gen_conv_matchPosition->at(photMCPos);
	if (matchPosition < -0.5) return -1; //not matched
	if (gen_conv_rDelta->at(photMCPos) > conv_maxDeltaR_analysis || gen_conv_ptDeltaRel->at(photMCPos) > conv_maxDPtRel_analysis) return -1;//not matched
	return matchPosition;
}

bool PhotSelectionPassMC(int photMCPos)  //uses variables loaded in main function
{
	int matchPosition = PhotMCMatched(photMCPos);
	if (matchPosition < -0.5) return false; //not matched
	return PhotSelectionPass(matchPosition);
}

bool PhotSelectionPassVeryLooseMC(int photMCPos)  //uses variables loaded in main function
{
	int matchPosition = PhotMCMatched(photMCPos);
	if (matchPosition < -0.5) return false; //not matched
	return PhotSelectionPassVeryLoose(matchPosition);
}

bool PhotSelectionPassMediumMC(int photMCPos)  //uses variables loaded in main function
{
	int matchPosition = gen_conv_matchPosition->at(photMCPos);
	if (matchPosition < -0.5) return false; //not matched
	if (gen_conv_rDelta->at(photMCPos) > conv_maxDeltaR_analysis || gen_conv_ptDeltaRel->at(photMCPos) > conv_maxDPtRel_analysis) return false;//not matched
	return PhotSelectionPassMedium(matchPosition);
}

bool PhotSelectionPassTightMC(int photMCPos)  //uses variables loaded in main function
{
	int matchPosition = gen_conv_matchPosition->at(photMCPos);
	if (matchPosition < -0.5) return false; //not matched
	if (gen_conv_rDelta->at(photMCPos) > conv_maxDeltaR_analysis || gen_conv_ptDeltaRel->at(photMCPos) > conv_maxDPtRel_analysis) return false;//not matched
	return PhotSelectionPassTight(matchPosition);
}

/////////////////////
//  D I M U O N /////
/////////////////////

// dimuon - acceptance

bool DimuonAcceptanceLoose(double rap, double pt)
{
	if (fabs(rap) > 2.4) return false;
	if (pt < 6.5) return false;
	if (pt > 30.0) return false;
	return true;
}
bool DimuonAcceptanceTight(double rap, double pt)
{
	if (fabs(rap) > 1.0) return false;  //2.4
	if (pt < 6.5) return false;
	if (pt > 25) return false;
	return true;
}

bool DimuonAcceptance(double rap, double pt) {
	return DimuonAcceptanceLoose(rap, pt);
}



// dimuon selection

bool DimuonSelectionPass(int dimuonPos)  //uses variables loaded in main function // contains acceptance cut for historic reasons
{
	if (dimuon_charge->at(dimuonPos) != 0) return false;
	if (dimuon_vtxProb->at(dimuonPos) < 0.01) return false;
	double ctauSigJpsi = dimuon_ctpv->at(dimuonPos) / dimuon_ctpvError->at(dimuonPos);
	if (ctauSigJpsi > 3 || ctauSigJpsi < -3) return false;
	double rap = ((TLorentzVector*)dimuon_p4->At(dimuonPos))->Rapidity();
	return DimuonAcceptance(rap, dimuon_pt->at(dimuonPos));
}

bool DimuonSelectionPassTight(int dimuonPos)  //uses variables loaded in main function // contains acceptance cut for historic reasons
{
	if (dimuon_charge->at(dimuonPos) != 0) return false;
	if (dimuon_vtxProb->at(dimuonPos) < 0.01) return false;
	double ctauSigJpsi = dimuon_ctpv->at(dimuonPos) / dimuon_ctpvError->at(dimuonPos);
	if (ctauSigJpsi > 3 || ctauSigJpsi < -3) return false;
	double rap = ((TLorentzVector*)dimuon_p4->At(dimuonPos))->Rapidity();
	return DimuonAcceptanceTight(rap, dimuon_pt->at(dimuonPos));
}

bool DimuonSelectionPassMC(int dimuonMCPos)
{
	int matchPosition = DimuonMCMatched(dimuonMCPos);
	if (matchPosition < -0.5) return false; //not matched
	return DimuonSelectionPass(matchPosition);
}

bool DimuonSelectionPassNoCharge(int dimuonPos)  //to allow selecting same sign bkg (but no sign cut directly)
{
	//if (dimuon_charge->at(dimuonPos) != 0) return false;
	//if (dimuon_vtxProb->at(dimuonPos) < 0.01) return false;
	double ctauSigJpsi = dimuon_ctpv->at(dimuonPos) / dimuon_ctpvError->at(dimuonPos);
	if (ctauSigJpsi > 3 || ctauSigJpsi < -3) return false;
	double rap = ((TLorentzVector*)dimuon_p4->At(dimuonPos))->Rapidity();
	return DimuonAcceptance(rap, dimuon_pt->at(dimuonPos));
}

int DimuonMCMatched(int dimuonMCPos)// check if the gen jpsi was matched to reco. Usually one Jpsi per event (thus index 0). New version - matching to the muons only: 
//-1 no muon match, -2 muons matched, but the dimuon doesn't exist (probably removed by dimuon preselection - confirmed for most)
{
	int dimuonPos = -2; //set to default (no match)

	int iMuon1MC = gen_muon_matchPosition->at(2 * dimuonMCPos);
	int iMuon2MC = gen_muon_matchPosition->at(2 * dimuonMCPos + 1);

	if (iMuon1MC < 0 || gen_muon_rDelta->at(2 * dimuonMCPos) > muon_maxDeltaR_analysis || gen_muon_ptDeltaRel->at(2 * dimuonMCPos) > muon_maxDPtRel_analysis) return -1; //muon not matched
	if (iMuon2MC < 0 || gen_muon_rDelta->at(2 * dimuonMCPos + 1) > muon_maxDeltaR_analysis || gen_muon_ptDeltaRel->at(2 * dimuonMCPos + 1) > muon_maxDPtRel_analysis) return -1; //muon not matched

	// all the final objects were matched, let's find the dimuon that is the proper match 


	for (int iJpsi = 0; iJpsi < dimuon_p4->GetEntriesFast(); iJpsi++) // Jpsi loop
	{
		if (dimuon_muon1_position->at(iJpsi) != iMuon1MC && (dimuon_muon1_position->at(iJpsi) != iMuon2MC)) { continue; } //first muon don't agree with either gen
		if (dimuon_muon2_position->at(iJpsi) != iMuon1MC && (dimuon_muon2_position->at(iJpsi) != iMuon2MC)) { continue; } //second muon don't agree with either gen
		if (dimuonPos < 0) { dimuonPos = iJpsi; break; } //we found the one matched
		else { cout << "Something wrong, two Jpsi shouldn't be both matched" << endl; } //we already found one, this is a crosscheck (once tested, add break above)
	}

	return dimuonPos; //either matched, or -2 if such reco doesn't exist
}

// overall dimuon pass

bool DimuonPassAllCuts(int dimuonPos)
{
	if (DimuonSelectionPass(dimuonPos) == false) return false;
	int muon1Pos = dimuon_muon1_position->at(dimuonPos);
	int muon2Pos = dimuon_muon2_position->at(dimuonPos);
	if (MuonSelectionPass(muon1Pos) == false) return false;
	if (MuonSelectionPass(muon2Pos) == false) return false;
	if (MuonAcceptance(muon_eta->at(muon1Pos), muon_pt->at(muon1Pos)) == false) return false;
	if (MuonAcceptance(muon_eta->at(muon2Pos), muon_pt->at(muon2Pos)) == false) return false;
	return true;
}


int DimuonPassAllCutsMC(int dimuonMCPos)  // -1 if failed, else returns the position of good reco
{
	int dimuonPos = DimuonMCMatched(dimuonMCPos);
	if (dimuonPos < 0) { return -1; }
	else if (DimuonPassAllCuts(dimuonPos) == false) { return -1; }
	else return dimuonPos;
}




/////////////////////
////  C H I C   /////
/////////////////////


bool ChiPassAllCuts(int chiPos)
{
	int dimuonPos = chi_daughterJpsi_position->at(chiPos);
	if (DimuonPassAllCuts(dimuonPos) == false) return false;
	double dimuonM = ((TLorentzVector*)dimuon_p4->At(dimuonPos))->M();
	if (dimuonM< mass_cutoffJpsi_l || dimuonM > mass_cutoffJpsi_h) return false;

	// photon
	int convPos = chi_daughterConv_position->at(chiPos);
	if (PhotAcceptance(conv_eta->at(convPos), conv_pt->at(convPos)) == false) return false;
	if (PhotSelectionPass(convPos) == false) return false;
	return true;

//	//cuts or MVA for photons   //NO TMVA FOR NOW
//#ifdef UsesTMVA
//	convQuality_isHighPurityValue = convQuality_isHighPurity->at(convPos);
//	convQuality_isGeneralTracksOnlyValue = convQuality_isGeneralTracksOnly->at(convPos);
//	conv_vertexPositionRhoValue = conv_vertexPositionRho->at(convPos);
//	conv_sigmaTkVtx1Value = conv_sigmaTkVtx1->at(convPos);
//	conv_sigmaTkVtx2Value = conv_sigmaTkVtx2->at(convPos);
//	conv_tkVtxCompatibilityOKValue = conv_tkVtxCompatibilityOK->at(convPos);
//	conv_compatibleInnerHitsOKValue = conv_compatibleInnerHitsOK->at(convPos);
//	conv_vertexChi2ProbValue = conv_vertexChi2Prob->at(convPos);
//	conv_dzToClosestPriVtxValue = conv_dzToClosestPriVtx->at(convPos);
//	conv_dxyPriVtxTimesCharge_Tr1Value = conv_dxyPriVtxTimesCharge_Tr1->at(convPos);
//	conv_dxyPriVtxTimesCharge_Tr2Value = conv_dxyPriVtxTimesCharge_Tr2->at(convPos);
//	conv_tk1NumOfDOFValue = conv_tk1NumOfDOF->at(convPos);
//	conv_tk2NumOfDOFValue = conv_tk2NumOfDOF->at(convPos);
//	conv_minDistanceOfApproachValue = conv_minDistanceOfApproach->at(convPos);
//	conv_etaValue = conv_eta->at(convPos);
//	conv_ptValue = conv_pt->at(convPos);
//
//
//	double photMVA = TMWAreader->EvaluateMVA("BDT");
//	hphoton_MVA_response->Fill(photMVA);
//	//cout << "Values MVA: HP: " << convQuality_isHighPurityValue << "   TracksOnly: " << convQuality_isGeneralTracksOnlyValue << "    Vprob: " << conv_vertexChi2ProbValue << "    and the response: "<< photMVA << endl;
//	//if (photMVA < -0.2) continue;
//#endif
//


}

bool ChiPassAllCutsVeryLooseConversion(int chiPos)
{
	int dimuonPos = chi_daughterJpsi_position->at(chiPos);
	if (DimuonPassAllCuts(dimuonPos) == false) return false;
	double dimuonM = ((TLorentzVector*)dimuon_p4->At(dimuonPos))->M();
	if (dimuonM< mass_cutoffJpsi_l || dimuonM > mass_cutoffJpsi_h) return false;

	// photon
	int convPos = chi_daughterConv_position->at(chiPos);
	if (PhotAcceptance(conv_eta->at(convPos), conv_pt->at(convPos)) == false) return false;
	if (PhotSelectionPassVeryLoose(convPos) == false) return false;
	return true;
}

bool ChiPassAllCutsMediumConversion(int chiPos)
{
	int dimuonPos = chi_daughterJpsi_position->at(chiPos);
	if (DimuonPassAllCuts(dimuonPos) == false) return false;
	double dimuonM = ((TLorentzVector*)dimuon_p4->At(dimuonPos))->M();
	if (dimuonM< mass_cutoffJpsi_l || dimuonM > mass_cutoffJpsi_h) return false;

	// photon
	int convPos = chi_daughterConv_position->at(chiPos);
	if (PhotAcceptance(conv_eta->at(convPos), conv_pt->at(convPos)) == false) return false;
	if (PhotSelectionPassMedium(convPos) == false) return false;
	return true;
}

bool ChiPassAllCutsTightConversion(int chiPos)
{
	int dimuonPos = chi_daughterJpsi_position->at(chiPos);
	if (DimuonPassAllCuts(dimuonPos) == false) return false;
	double dimuonM = ((TLorentzVector*)dimuon_p4->At(dimuonPos))->M();
	if (dimuonM< mass_cutoffJpsi_l || dimuonM > mass_cutoffJpsi_h) return false;

	// photon
	int convPos = chi_daughterConv_position->at(chiPos);
	if (PhotAcceptance(conv_eta->at(convPos), conv_pt->at(convPos)) == false) return false;
	if (PhotSelectionPassTight(convPos) == false) return false;
	return true;
}

int ChiMCMatched(int chiMCPos)// check if the gen chic was matched to reco. Usually one chic per event (thus index 0). New version - matching to the muons and conversions only
//-1 no matches, -2 conversions and muons matched, but the chic doesn't exist (probably removed by dimuon preselection - confirmed for most, but could maybe rarely be caused by chic selection)
{
	int chiPos = -2; //set to default (no match)
	
	// check if the matching to the gen objects happened (this we can get directly from the stored information, without checking the actual matched reco object)
	if (gen_conv_matchPosition->at(chiMCPos) < 0 || gen_conv_rDelta->at(chiMCPos) > conv_maxDeltaR_analysis || gen_conv_ptDeltaRel->at(chiMCPos) > conv_maxDPtRel_analysis) return -1; //conv was not matched
	if (gen_muon_matchPosition->at(2 * chiMCPos) < 0 || gen_muon_rDelta->at(2 * chiMCPos) > muon_maxDeltaR_analysis || gen_muon_ptDeltaRel->at(2 * chiMCPos) > muon_maxDPtRel_analysis) return -1; //muon not matched
	if (gen_muon_matchPosition->at(2 * chiMCPos + 1) < 0 || gen_muon_rDelta->at(2 * chiMCPos + 1) > muon_maxDeltaR_analysis || gen_muon_ptDeltaRel->at(2 * chiMCPos + 1) > muon_maxDPtRel_analysis) return -1; //muon not matched

	// all the final objects were matched, let's find the chic that is the proper match 


	for (int iChi = 0; iChi < chi_p4->GetEntriesFast(); iChi++) //loop over chic to find the one matched
	{
		if (chi_daughterConv_position->at(iChi) != gen_conv_matchPosition->at(chiMCPos)) { continue; }
		int dimuonPos = chi_daughterJpsi_position->at(iChi);
		int iMuon1MC = gen_muon_matchPosition->at(2 * chiMCPos);
		int iMuon2MC = gen_muon_matchPosition->at(2 * chiMCPos + 1);

		if (dimuon_muon1_position->at(dimuonPos) != iMuon1MC && (dimuon_muon1_position->at(dimuonPos) != iMuon2MC)) { continue; } //first muon don't agree with either gen
		if (dimuon_muon2_position->at(dimuonPos) != iMuon1MC && (dimuon_muon2_position->at(dimuonPos) != iMuon2MC)) { continue; } //second muon don't agree with either gen
		if (chiPos < 0) { chiPos = iChi; } //we found the one matched
		else { cout << "Something wrong, two chic shouldn't be both matched" << endl; } //we already found one, this is a crosscheck (once tested, add break above)
	}

	return chiPos; //either matched, or -2 if such reco doesn't exist
}



int ChiPassAllCutsMC(int chiMCPos)
{
	int chiPos = ChiMCMatched(chiMCPos); // we have 1 chic, check as which one it was potentially reconstructed
	if (chiPos < -0.5) { return -1; }// -1 not matched
	else if (ChiPassAllCuts(chiPos) == false) { return -1; }
	else return chiPos;
}

int ChiPassAllCutsVeryLooseConversionMC(int chiMCPos)
{
	int chiPos = ChiMCMatched(chiMCPos); // we have 1 chic, check as which one it was potentially reconstructed
	if (chiPos < -0.5) { return -1; }// -1 not matched
	else if (ChiPassAllCutsVeryLooseConversion(chiPos) == false) { return -1; }
	else return chiPos;
}

int ChiPassAllCutsMediumConversionMC(int chiMCPos)
{
	int chiPos = ChiMCMatched(chiMCPos); // we have 1 chic, check as which one it was potentially reconstructed
	if (chiPos < -0.5) { return -1; }// -1 not matched
	else if (ChiPassAllCutsMediumConversion(chiPos) == false) { return -1; }
	else return chiPos;
}

int ChiPassAllCutsTightConversionMC(int chiMCPos)
{
	int chiPos = ChiMCMatched(chiMCPos); // we have 1 chic, check as which one it was potentially reconstructed
	if (chiPos < -0.5) { return -1; }// -1 not matched
	else if (ChiPassAllCutsTightConversion(chiPos) == false) { return -1; }
	else return chiPos;
}


bool ChiIsMatchedAllDaughters(int chiPos, int chiMCPos)// check if the gen chic was matched to reco, as well as J/psi and photon. Usually one chic per event (thus index 0)
// mostly obsolete
{
	int iChiMC = gen_chic_matchPosition->at(chiMCPos); // we have 1 chic, check as which one it was potentially reconstructed
	if (iChiMC < -0.5) { return false; }// -1 not matched
	if (chiPos != iChiMC) { return false; }
	if (chi_daughterConv_position->at(chiPos) != gen_conv_matchPosition->at(chiMCPos)) { return false; } //conversions don't agree
	if (gen_conv_matchPosition->at(chiMCPos)<0 || gen_conv_rDelta->at(chiMCPos) > conv_maxDeltaR_analysis || gen_conv_ptDeltaRel->at(chiMCPos) > conv_maxDPtRel_analysis) return false;//not matched
	int dimuonPos = chi_daughterJpsi_position->at(chiPos);
	if (dimuonPos != gen_Jpsi_matchPosition->at(chiMCPos)) { return false; } //dimuons don't agree
	if (gen_Jpsi_matchPosition->at(chiMCPos) < 0 || gen_Jpsi_rDelta->at(chiMCPos) > jpsi_maxDeltaR_analysis || gen_Jpsi_ptDeltaRel->at(chiMCPos) > jpsi_maxDPtRel_analysis) return false;//not matched
	int iMuon1MC = gen_muon_matchPosition->at(2 * chiMCPos);
	int iMuon2MC = gen_muon_matchPosition->at(2 * chiMCPos + 1);
	if (dimuon_muon1_position->at(dimuonPos) != iMuon1MC && (dimuon_muon1_position->at(dimuonPos) != iMuon2MC)) { return false; } //first muon don't agree with either gen
	if (dimuon_muon2_position->at(dimuonPos) != iMuon1MC && (dimuon_muon2_position->at(dimuonPos) != iMuon2MC)) { return false; } //second muon don't agree with either gen
	if (gen_muon_matchPosition->at(2 * chiMCPos) < 0 || gen_muon_rDelta->at(2 * chiMCPos) > muon_maxDeltaR_analysis || gen_muon_ptDeltaRel->at(2 * chiMCPos) > muon_maxDPtRel_analysis) return false;//not matched
	if (gen_muon_matchPosition->at(2 * chiMCPos +1) < 0 || gen_muon_rDelta->at(2 * chiMCPos +1) > muon_maxDeltaR_analysis || gen_muon_ptDeltaRel->at(2 * chiMCPos +1) > muon_maxDPtRel_analysis) return false;//not matched

	return true;
}



double WeightForMC_pTpart(double pt) //weights official MC to match data better, to be used for nominal results (and pretty much everywhere). Input - pT of J/psi from chic decay in GeV; output - weight (RD/MC) with which to multiply MC so it matches data 
{
	if (pt < 6.5 || pt>30) {
		cout << "WARNING: pT is smaller or larger than our acceptance cut-off, no weight calculated for those cases. 1 (no additional weight) returned" << endl;
		return 1;
	}
	else {
		return (0.104*pt - 0.07); //parametrization from the study comparing the RD MC spectra
	}
}

double WeightPhotonAcceptanceSystematic(double photon_pt, int idx) // weights official MC with variations obtained from alternative MC settings. To be used as systematics. idx: 1->pthat=3; 2->pthat=6; 3->cmass=1.43; 4->renormalization+factorization variation 
// currently idx 1 and 2 are largest variations, others need not to be considered
{
	if (idx == 0 || photon_pt > 5.0) return 1; //no additional weight for 0 (nominal) and above 5 GeV (the statistics there is close to none anyway, better no weight than randomly extrapolated function (moreover, all functions are close to 1 at 5GeV))
	// for following, I copied the result of the fit with all the digits, however the actual precision is lower, and documented in the AN (last digits mean nothing)
	if (idx == 1) { return (0.025861 * photon_pt*photon_pt - 0.144868 * photon_pt + 1.1045); } //pthat 3
	if (idx == 2) { return (-0.0391252 * photon_pt*photon_pt + 0.239757 * photon_pt + 0.809756); } //pthat 6
	if (idx == 3) { return (-0.0105452 * photon_pt*photon_pt + 0.0164734 * photon_pt + 1.0015); } //cmass =1.43 (nominal seems to be 1.5)
	if (idx == 4) { return (-0.0308825 * photon_pt*photon_pt + 0.186884 * photon_pt + 0.852952); } //renorm and factorization scale change
	cout << "Warning, in function WeightPhotonAcceptanceSystematic idx doesn't correspond to any selection, returning no weight." << endl;
	return 1;

}



double PolarizationCosTheta(TLorentzVector* LVdimuon, TLorentzVector* LVmuon) // calculates cos theta for the J/psi, based on Jeongho
{
	TLorentzVector* LVdimuonCopy = new TLorentzVector(*LVdimuon); // this is here as a quick fix, so we do not change the original LVs
	TLorentzVector* LVmuonCopy = new TLorentzVector(*LVmuon); // this is here as a quick fix, so we do not change the original LVs
	TVector3 TVdimuon = LVdimuonCopy->Vect();
	TVector3 TVUdimuon = TVdimuon.Unit(); //Unit vector of J/psi Z axis
	TVector3 boost = -LVdimuonCopy->BoostVector();
	TLorentzVector LVmuon_prime = LVmuonCopy->Transform(boost); //Transverse to J/psi frame
	TVector3 TVmuon = LVmuon_prime.Vect();
	TVector3 TVUmuon = TVmuon.Unit(); //muon Unit Vector
	double costheta = TVUdimuon.Dot(TVUmuon);
	return costheta;
}


double PolarizationWeight(TLorentzVector* LVdimuon, TLorentzVector* LVmuon, double lambdaTheta) // assigns weight to chic, assuming the polarization axes are the same as for J/psi. Following Arxiv 1103.4882 
{
	double cosTheta = PolarizationCosTheta(LVdimuon, LVmuon);
	double weight = 3/(3-lambdaTheta)*(1 + lambdaTheta * cosTheta*cosTheta);
	return weight;
}

double PolarizationWeight_ChicStateWeighted(TLorentzVector* LVdimuon, TLorentzVector* LVmuon, int gen_pdgId, double lambdaTheta) // created as a new function that can be called instead of "PolarizationWeight", if chic2/chic1 ratio matters
{
	double weight = PolarizationWeight(LVdimuon, LVmuon, lambdaTheta); 
	if (gen_pdgId == 20443) { //chic1
		weight *= 1.12;	//weights based on a study from 12.7.2023
	}
	else if (gen_pdgId == 445) { //chic2
		weight *= 0.8;	//weights based on a study from 12.7.2023
	}
	else {
		cout << "Warning, the chic state Pythia code not recognized. Weight returned but no chic state weighting done" << endl;
	}
	return weight;
}


double Polarizationlambda_pTdependence(double jpsi_pt, TLorentzVector* jpsi_p4, int gen_pdgId)  // copied from Jeongho
{
	double jpsi_mass = jpsi_p4->M();
	double jpsiptM = jpsi_pt / jpsi_mass;
	double lambdaTheta = 1;
	if (jpsiptM < 5 && gen_pdgId == 20443) {
		lambdaTheta = -0.06*jpsiptM + 0.8;
	}
	if (jpsiptM > 5 && gen_pdgId == 20443) {
		lambdaTheta = 0.004*jpsiptM + 0.48;
	}
	if (gen_pdgId == 445) {
		lambdaTheta = 0.075*jpsiptM - 0.95;
	}
	return lambdaTheta;
}


////////////////////////////
//// Chic2/chic1 ratio /////
////////////////////////////

double CalculateChicRatioValue(double x) {
	return x / (1 - x);
}
double CalculateChicRatioError(double x, double delta_x) {
	double function_value = x / (1-x);
	double derivative = 1 / ((1-x)*(1-x));
	double propagated_error = std::abs(derivative) * delta_x;
	return propagated_error;
}
TGraphAsymmErrors* CalculateChicRatioFromC2Ratio(TGraphAsymmErrors* graph) {
	int nPoints = graph->GetN();

	// Create arrays to store the calculated values
	double* xValues = graph->GetX();
	double* yValues = graph->GetY();
	double* xErrorsLow = graph->GetEXlow();
	double* xErrorsHigh = graph->GetEXhigh();
	double* yErrorsLow = graph->GetEYlow();
	double* yErrorsHigh = graph->GetEYhigh();

	// Create arrays to store the calculated values for the new graph
	double* newXValues = new double[nPoints];
	double* newYValues = new double[nPoints];
	double* newXErrorsLow = new double[nPoints];
	double* newXErrorsHigh = new double[nPoints];
	double* newYErrorsLow = new double[nPoints];
	double* newYErrorsHigh = new double[nPoints];

	for (int i = 0; i < nPoints; ++i) {

		newXValues[i] = xValues[i];
		newYValues[i] = CalculateChicRatioValue(yValues[i]);
		newXErrorsLow[i] = CalculateChicRatioError(xValues[i], xErrorsLow[i]);
		newXErrorsHigh[i] = CalculateChicRatioError(xValues[i], xErrorsHigh[i]);
		newYErrorsLow[i] = CalculateChicRatioError(yValues[i], yErrorsLow[i]);
		newYErrorsHigh[i] = CalculateChicRatioError(yValues[i], yErrorsHigh[i]);

	}

	TGraphAsymmErrors* outputGraph = new TGraphAsymmErrors(nPoints, newXValues, newYValues, newXErrorsLow, newXErrorsHigh, newYErrorsLow, newYErrorsHigh);

	outputGraph->SetMarkerColor(kBlue);
	outputGraph->SetLineColor(kBlue);
	outputGraph->SetMarkerStyle(25);
	outputGraph->SetMaximum(1.001);
	outputGraph->SetMinimum(-0.001);
	outputGraph->GetYaxis()->SetTitleSize(0.05);
	outputGraph->GetXaxis()->SetTitleSize(0.05);
	outputGraph->GetYaxis()->SetTitleOffset(1.00);
	outputGraph->GetXaxis()->SetTitleOffset(1.00);
	outputGraph->GetYaxis()->SetTitle("#chi_{c2}/#chi_{c1}");
	return outputGraph;
}

